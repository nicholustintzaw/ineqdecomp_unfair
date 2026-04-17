capture program drop ineqdecomp_unfair
program define ineqdecomp_unfair, rclass
    version 17.0

    /******************************************************************
    ineqdecomp_unfair
    -------------------------------------------------------------------
    Purpose
      Decompose inequality in a binary outcome using unfair factors.

    Core design
      1. Detect omitted / non-estimable unfair variables
      2. Build one fixed common analytic sample
      3. Build unfairness rank from survey-weighted predicted risk
      4. Estimate AMEs from survey-weighted probit + margins
      5. Compute elasticity using observed mean outcome (WB-style)
      6. Compute factor CI using selected CI type:
           - relative   (default)
           - wagstaff
           - erreygers
      7. Compute signed contributions
      8. Report:
           - signed contribution
           - absolute percentage contribution
           - normalized absolute percentage contribution
      9. Report omitted variables separately

    Methodological rationale
      This program implements a hybrid workflow based on the decisions
      made during method development:

      (a) Ranking variable:
          The inequality ranking variable is created from the predicted
          probability of the outcome using unfair factors only.
          This makes the rank a multivariate unfairness score rather
          than a simple SES or wealth rank.

      (b) Elasticity:
          Elasticities are computed using average marginal effects (AMEs)
          from survey-weighted probit + margins, and then scaled using
          the observed weighted mean of the outcome:
              elasticity_x = (AME_x * mean(x)) / mean(y)
          This follows the World Bank-style denominator logic more
          closely than using a standardized predicted outcome mean.

      (c) Contributions:
          Signed contributions are retained for the decomposition
          identity:
              CI_y = sum(contribution_x) + residual
          For reporting only, absolute contributions are converted into
          absolute percentage shares and normalized to sum to 100,
          following a VERSE-style presentation logic.

      (d) Common analytic sample:
          A single fixed analytic sample is created and then used for
          every decomposition component:
              - rank model
              - probit + margins model
              - outcome mean
              - unfair-factor means
              - outcome CI
              - unfair-factor CI
          This avoids mixing a full-sample outcome CI with decomposition
          parameters estimated on a smaller model sample.

    Important assumptions
      - The data must already be svyset before running this program.
      - unfair() should contain analysis-ready numeric variables.
      - Do not use factor-variable notation such as i.var in unfair().
      - For binary outcomes, relative CI is the default and most
        defensible decomposition scale.
      - Wagstaff and Erreygers options are provided as bounded-index
        extensions for sensitivity analysis.
      - The common analytic sample is defined after omission checks
        using outcome, retained unfair variables, and the supplied
        weight variables.

    Output
      The resulting dataset is left in memory only if clear is specified.
      Otherwise, the original dataset is restored after returning scalars
      and macros in r().

    ******************************************************************/

	syntax , ///
		OUTcome(varname numeric) ///
		UNFair(varlist numeric) ///
		WVAR(varname numeric) ///
		[ CIType(string) ///
		  RANKModel(string) ///
		  CLEAR ]

    *---------------------------------------------------------------*
    * 0. Default options and basic checks
    *---------------------------------------------------------------*
    if "`citype'"     == "" local citype relative
    if "`rankmodel'"  == "" local rankmodel logit

    local citype    = lower("`citype'")
    local rankmodel = lower("`rankmodel'")

    if !inlist("`citype'", "relative", "wagstaff", "erreygers") {
        di as error "citype() must be relative, wagstaff, or erreygers"
        exit 198
    }

    if !inlist("`rankmodel'", "logit", "probit") {
        di as error "rankmodel() must be logit or probit"
        exit 198
    }

    marksample touse, strok

    tempvar analytic rank esample_rank esample_probit
    tempfile resultsfile finalout

    *---------------------------------------------------------------*
    * 1. Initial omission / collinearity check
    *    Use survey-weighted probit because final AMEs come from probit
    *---------------------------------------------------------------*
    quietly svy if `touse': probit `outcome' `unfair'

    matrix b0 = e(b)
    local terms : colnames b0

    local X_keep
    foreach term of local terms {
        if "`term'" == "_cons" continue

        capture noisily _ms_parse_parts `term'
        if !_rc {
            if r(omit) == 1 continue
            if r(base) == 1 continue
        }

        local X_keep `X_keep' `term'
    }

    local X_drop
    foreach x of local unfair {
        local found = 0
        foreach k of local X_keep {
            if "`x'" == "`k'" local found = 1
        }
        if `found' == 0 local X_drop `X_drop' `x'
    }

    if "`X_keep'" == "" {
        di as error "No estimable unfair variables remained after omission checks."
        exit 498
    }

    *---------------------------------------------------------------*
    * 2. Build a single common analytic sample
    *    This sample will be used consistently for all later steps:
    *      - rank model
    *      - probit + margins model
    *      - outcome mean
    *      - unfair-factor means
    *      - outcome CI
    *      - unfair-factor CI
    *
    *    This avoids mixing decomposition parameters across different
    *    effective samples.
    *---------------------------------------------------------------*
    gen byte `analytic' = `touse'
    markout `analytic' `outcome' `X_keep' `wvar'

    quietly count if `analytic'
    if r(N) == 0 {
        di as error "No observations remain in the common analytic sample."
        exit 498
    }

    *---------------------------------------------------------------*
    * 3. Build unfairness rank from survey-weighted predicted risk
    *    The rank model is estimated on the common analytic sample.
    *---------------------------------------------------------------*
    if "`rankmodel'" == "logit" {
        quietly svy if `analytic': logit `outcome' `X_keep'
        gen byte `esample_rank' = e(sample)
        quietly predict double `rank' if e(sample), pr
    }
    else {
        quietly svy if `analytic': probit `outcome' `X_keep'
        gen byte `esample_rank' = e(sample)
        quietly predict double `rank' if e(sample), pr
    }

    * Defensive check: rank-model estimation sample should match analytic
    quietly count if `analytic' & `esample_rank' != 1
    if r(N) > 0 {
        di as error "Rank model estimation sample does not fully match the common analytic sample."
        exit 498
    }

    *---------------------------------------------------------------*
    * 4. Final survey-weighted probit for AMEs via margins
    *    The probit model is also estimated on the same common analytic
    *    sample, so the AMEs are aligned with the rest of the
    *    decomposition.
    *---------------------------------------------------------------*
    quietly svy if `analytic': probit `outcome' `X_keep'
    gen byte `esample_probit' = e(sample)

    * Defensive check: probit estimation sample should match analytic
    quietly count if `analytic' & `esample_probit' != 1
    if r(N) > 0 {
        di as error "Final probit estimation sample does not fully match the common analytic sample."
        exit 498
    }

    quietly margins, dydx(`X_keep')
    matrix dfdx = r(b)

    *---------------------------------------------------------------*
    * 5. Outcome mean and outcome CI
    *    All calculations are now based on the same common analytic
    *    sample.
    *
    *    CI type options:
    *      relative  : 2*cov(rank,y)/mean(y)
    *      wagstaff  : CI/(1-mean(y))
    *      erreygers : CI * [4*mean(y)/(max(y)-min(y))]
    *---------------------------------------------------------------*
    quietly summarize `outcome' if `analytic' [aw=`wvar'], meanonly
    scalar mu_y   = r(mean)
    scalar min_y  = r(min)
    scalar max_y  = r(max)

    if abs(mu_y) < 1e-12 {
        di as error "Mean of outcome is zero or too close to zero."
        exit 498
    }

    quietly corr `rank' `outcome' if `analytic' [aw=`wvar'], c
    scalar cov_y = r(cov_12)
    scalar CI_rel_y = 2 * cov_y / mu_y

    if "`citype'" == "relative" {
        scalar CI_y = CI_rel_y
    }
    else if "`citype'" == "wagstaff" {
        if abs(1 - mu_y) < 1e-12 {
            di as error "Wagstaff correction undefined because mean(outcome)=1."
            exit 498
        }
        scalar CI_y = CI_rel_y / (1 - mu_y)
    }
    else if "`citype'" == "erreygers" {
        if abs(max_y - min_y) < 1e-12 {
            di as error "Erreygers correction undefined because outcome has zero range."
            exit 498
        }
        scalar CI_y = CI_rel_y * (4 * mu_y / (max_y - min_y))
    }

    if abs(CI_y) < 1e-12 {
        di as error "Outcome CI is zero or too close to zero; percentage contributions undefined."
        exit 498
    }

    *---------------------------------------------------------------*
    * 6. Initialize result posting structure
    *---------------------------------------------------------------*
    tempname posth
    postfile `posth' ///
        int order_id ///
        str40 var ///
        str18 status ///
        double elasticity ///
        double var_ci ///
        double contribution ///
        double contribution_abs ///
        double contribution_pct_signed ///
        double contribution_pct_abs ///
        double outcome_ci ///
        using `resultsfile', replace

    scalar unfair_sum = 0

    *---------------------------------------------------------------*
    * 7. Loop over retained unfair variables
    *    Each unfair factor uses:
    *      - AME from common-sample probit + margins
    *      - weighted mean from common analytic sample
    *      - CI from common analytic sample
    *---------------------------------------------------------------*
    local i = 1
    foreach x of local X_keep {

        capture scalar ame_x = dfdx[1,"`x'"]
        if _rc {
            di as error "Could not extract AME from margins result for variable `x'."
            di as error "Check variable naming and confirm unfair() contains analysis-ready numeric variables."
            exit 498
        }

        quietly summarize `x' if `analytic' [aw=`wvar'], meanonly
        scalar mu_x  = r(mean)
        scalar min_x = r(min)
        scalar max_x = r(max)

        if abs(mu_x) < 1e-12 {
            scalar elas_x = .
            scalar CI_x   = .
            scalar con_x  = .
            scalar cona_x = .
            scalar pcts_x = .
            scalar pcta_x = .
        }
        else {
            * Elasticity: WB-style denominator uses observed mean outcome
            scalar elas_x = (ame_x * mu_x) / mu_y

            * Relative CI of unfair factor
            quietly corr `rank' `x' if `analytic' [aw=`wvar'], c
            scalar cov_x = r(cov_12)
            scalar CI_rel_x = 2 * cov_x / mu_x

            * Apply selected CI correction to unfair factor
            if "`citype'" == "relative" {
                scalar CI_x = CI_rel_x
            }
            else if "`citype'" == "wagstaff" {
                if abs(1 - mu_x) < 1e-12 {
                    scalar CI_x = .
                }
                else {
                    scalar CI_x = CI_rel_x / (1 - mu_x)
                }
            }
            else if "`citype'" == "erreygers" {
                if abs(max_x - min_x) < 1e-12 {
                    scalar CI_x = .
                }
                else {
                    scalar CI_x = CI_rel_x * (4 * mu_x / (max_x - min_x))
                }
            }

            if missing(CI_x) {
                scalar con_x  = .
                scalar cona_x = .
                scalar pcts_x = .
                scalar pcta_x = .
            }
            else {
                * Signed contribution for decomposition identity
                scalar con_x = elas_x * CI_x

                * Absolute contribution only for reporting
                scalar cona_x = abs(con_x)

                * Signed and absolute percentages
                scalar pcts_x = 100 * (con_x / CI_y)
                scalar pcta_x = abs(pcts_x)

                scalar unfair_sum = unfair_sum + con_x
            }
        }

        post `posth' ///
            (`i') ///
            ("`x'") ///
            ("estimated") ///
            (elas_x) ///
            (CI_x) ///
            (con_x) ///
            (cona_x) ///
            (pcts_x) ///
            (pcta_x) ///
            (CI_y)

        local ++i
    }

    *---------------------------------------------------------------*
    * 8. Add omitted / dropped unfair variables to output
    *---------------------------------------------------------------*
    foreach x of local X_drop {
        post `posth' ///
            (`i') ///
            ("`x'") ///
            ("omitted") ///
            (.) ///
            (.) ///
            (.) ///
            (.) ///
            (.) ///
            (.) ///
            (CI_y)
        local ++i
    }

    *---------------------------------------------------------------*
    * 9. Residual from signed decomposition
    *---------------------------------------------------------------*
    scalar residual_signed = CI_y - unfair_sum
    scalar residual_abs    = abs(residual_signed)
    scalar residual_pct_signed = 100 * (residual_signed / CI_y)
    scalar residual_pct_abs    = abs(residual_pct_signed)

    post `posth' ///
        (`i') ///
        ("Residual") ///
        ("residual") ///
        (.) ///
        (.) ///
        (residual_signed) ///
        (residual_abs) ///
        (residual_pct_signed) ///
        (residual_pct_abs) ///
        (CI_y)

    postclose `posth'

    *---------------------------------------------------------------*
    * 10. Finalize result dataset
    *---------------------------------------------------------------*
	quietly count if `analytic'
	local n_analytic = r(N)

    preserve
        use `resultsfile', clear
        sort order_id

		gen n_analytic = `n_analytic' 
		
        egen contribution_pct_tot_abs = total(contribution_pct_abs)

        gen contribution_pct_norml = 100 * contribution_pct_abs / contribution_pct_tot_abs ///
            if !missing(contribution_pct_abs)

        gen note = ""
        replace note = "Estimated unfair factor" if status == "estimated"
        replace note = "Dropped from model due to omission/base/collinearity" if status == "omitted"
        replace note = "Residual from signed decomposition" if status == "residual"

        label var var                       "Unfair factor variable name"
        label var status                    "Estimation status"
		label var n_analytic 				"Number of observations in common analytic sample"
        label var elasticity                "Elasticity"
        label var var_ci                    "Concentration index of unfair factor"
        label var contribution              "Signed contribution"
        label var contribution_abs          "Absolute contribution"
        label var contribution_pct_signed   "Signed percentage contribution"
        label var contribution_pct_abs      "Absolute percentage contribution"
        label var contribution_pct_tot_abs  "Total absolute percentage contribution"
        label var contribution_pct_norml    "Normalized absolute percentage contribution"
        label var outcome_ci                "Concentration index of outcome"
        label var note                      "Interpretation note"

		order order_id var status n_analytic elasticity var_ci contribution contribution_abs ///
			  contribution_pct_signed contribution_pct_abs contribution_pct_tot_abs ///
			  contribution_pct_norml outcome_ci note

        count if status == "estimated"
        return scalar n_estimated = r(N)

        count if status == "omitted"
        return scalar n_omitted = r(N)

        //quietly count if `analytic'
        //return scalar n_analytic = r(N)
		return scalar n_analytic = `n_analytic'

        return scalar outcome_ci    = CI_y
        return scalar outcome_mean  = mu_y
        return scalar explained_sum = unfair_sum
        return scalar residual      = residual_signed
        return scalar residual_abs  = residual_abs
        return local unfair_kept    "`X_keep'"
        return local unfair_dropped "`X_drop'"
        return local citype         "`citype'"
        return local rankmodel      "`rankmodel'"

        save `finalout', replace
    restore

    if "`clear'" != "" {
        use `finalout', clear
    }

end