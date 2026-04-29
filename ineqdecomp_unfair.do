capture program drop ineqdecomp_unfair
program define ineqdecomp_unfair, rclass
    version 17.0

	/******************************************************************
	ineqdecomp_unfair
	-------------------------------------------------------------------
	Purpose
	  Decompose inequality in a binary outcome using unfair factors
	  within a concentration index (CI) framework.

	Core design
	  1. Detect omitted / non-estimable unfair variables
	  2. Build one fixed common analytic sample
	  3. Construct unfairness rank from predicted risk and convert to
		 weighted fractional rank (conindex-consistent)
	  4. Estimate AMEs from survey-weighted probit + margins
	  5. Compute elasticity using observed mean outcome (WB-style)
	  6. Compute factor CI using weighted summation formula:
		   - relative CI (default, used for decomposition)
		   - wagstaff (sensitivity only)
		   - erreygers (sensitivity only)
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
		  probability of the outcome using unfair factors, and then
		  converted into a weighted fractional rank:

			  R_i = (cum_weight_i − 0.5 × w_i) / total_weight

		  This ensures consistency with standard CI estimation and
		  comparability with conindex results. The rank therefore
		  represents a multivariate unfairness ordering.

	  (b) Elasticity:
		  Elasticities are computed using average marginal effects (AMEs)
		  from survey-weighted probit + margins, and then scaled using
		  the observed weighted mean of the outcome:

			  elasticity_x = (AME_x * mean(x)) / mean(y)

		  This follows the World Bank decomposition logic adapted to
		  nonlinear models.

	  (c) Contributions:
		  Signed contributions are retained for the decomposition identity:

			  CI_y = sum(contribution_x) + residual

		  where:

			  contribution_x = elasticity_x × CI_x

		  For reporting only, absolute contributions are converted into
		  percentage shares and normalized to sum to 100, following a
		  VERSE-style presentation.

	  (d) CI calculation:
		  The concentration index is computed using a weighted summation
		  formula consistent with conindex:

			  CI = (2 / mean(y)) × Σ [w_i × y_i × (R_i − 0.5)]

		  This avoids covariance-based shortcuts and ensures consistency
		  under weighted data and fractional ranking.

	  (e) Common analytic sample:
		  A single fixed analytic sample is created and then used for
		  every decomposition component:
			  - rank model
			  - probit + margins model
			  - outcome mean
			  - unfair-factor means
			  - outcome CI
			  - unfair-factor CI

		  This avoids mixing different estimation samples and ensures
		  internal consistency.

	  (f) CI type and interpretation:
		  The standard (relative) CI is used for decomposition because
		  the additive identity holds under this formulation.

		  Wagstaff and Erreygers corrections are included as optional
		  sensitivity analyses:

			  CI_w = CI / (1 − μ)
			  CI_e = 4μ × CI

		  However, these transformations rescale the index and break the
		  additive decomposition identity, often leading to inflated
		  residual terms. They should therefore not be used for the main
		  decomposition results.

	Important assumptions
	  - The data must already be svyset before running the program.
	  - unfair() should contain analysis-ready numeric variables.
	  - Do not use factor-variable notation such as i.var in unfair().
	  - Relative CI is the appropriate scale for decomposition.
	  - Wagstaff and Erreygers are provided for sensitivity checks only.
	  - The analytic sample is defined after omission checks using
		outcome, retained unfair variables, and the weight variable.

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

    tempvar analytic rankscore rank esample_rank esample_probit ///
			wnorm cumw wtmp obsid rcenter prod
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
	*    The model predicts the unfairness score. This score is then
	*    converted into a weighted fractional rank, following the
	*    concentration-index logic used by conindex.
	*
	*    Important:
	*      - rankscore = predicted probability from unfair factors
	*      - rank      = weighted fractional rank of rankscore
	*---------------------------------------------------------------*
	gen long `obsid' = _n

	if "`rankmodel'" == "logit" {
		quietly svy if `analytic': logit `outcome' `X_keep'
		gen byte `esample_rank' = e(sample)
		quietly predict double `rankscore' if e(sample), pr
	}
	else {
		quietly svy if `analytic': probit `outcome' `X_keep'
		gen byte `esample_rank' = e(sample)
		quietly predict double `rankscore' if e(sample), pr
	}

	* Defensive check: rank-model estimation sample should match analytic
	quietly count if `analytic' & `esample_rank' != 1
	if r(N) > 0 {
		di as error "Rank model estimation sample does not fully match the common analytic sample."
		exit 498
	}

	* Weighted fractional rank, conindex-style
	quietly summarize `wvar' if `analytic', meanonly
	scalar total_w = r(sum)

	if total_w <= 0 | missing(total_w) {
		di as error "Total analytic sample weight is zero or missing."
		exit 498
	}

	gen double `wtmp' = cond(`analytic', `wvar', 0)

	sort `rankscore' `obsid'
	gen double `cumw' = sum(`wtmp')

	gen double `rank' = (`cumw' - 0.5*`wvar') / total_w if `analytic'

	sort `obsid'

	gen double `wnorm' = `wvar' / total_w if `analytic'
	gen double `rcenter' = `rank' - 0.5 if `analytic'

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

	gen double `prod' = `wnorm' * `outcome' * `rcenter' if `analytic'
	quietly summarize `prod' if `analytic', meanonly
	scalar CI_rel_y = 2 * r(sum) / mu_y
	drop `prod'

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
            tempvar prodx
			gen double `prodx' = `wnorm' * `x' * `rcenter' if `analytic'
			quietly summarize `prodx' if `analytic', meanonly
			scalar CI_rel_x = 2 * r(sum) / mu_x
			drop `prodx'

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