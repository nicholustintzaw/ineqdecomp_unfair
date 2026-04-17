# ineqdecomp_unfair


---

# 📊 Inequality Decomposition with Unfair Factors (Stata)

## Overview

This repository provides a Stata program:

```stata
ineqdecomp_unfair
```

for decomposing inequality in a binary outcome using **unfair factors**, based on a hybrid methodological approach combining:

* **World Bank concentration index decomposition framework**
* **Model-based decomposition logic inspired by the VERSE project**
* **Robust CI calculation principles consistent with the `rineq` package**

---

## 🎯 Purpose

The function is designed to:

* Quantify **how much each unfair factor contributes to inequality**
* Provide both:

  * **Signed contributions** (for decomposition identity)
  * **Normalized absolute contributions** (for reporting and interpretation)
* Ensure **internal consistency** by using a **single analytic sample**

---

## 🧠 Methodological Framework

### 1. Concentration Index (CI)

The decomposition is based on the standard concentration index:

[
CI = \frac{2}{\mu} \cdot \text{Cov}(y, R)
]

where:

* ( y ) = outcome
* ( R ) = rank variable (here: predicted probability from unfair factors)
* ( \mu ) = mean of outcome

---

### 2. Decomposition Structure

Each unfair factor contributes:

[
Contribution_k = Elasticity_k \times CI_k
]

and:

[
CI_y = \sum_k Contribution_k + Residual
]

---

### 3. Elasticity (Hybrid WB Approach)

Elasticity is computed as:

[
Elasticity_k = \frac{AME_k \cdot \bar{x}_k}{\bar{y}}
]

where:

* ( AME_k ) = average marginal effect from `svy: probit + margins`
* ( \bar{x}_k ) = weighted mean of factor
* ( \bar{y} ) = weighted mean of outcome

This follows the **World Bank-style denominator** while using nonlinear marginal effects.

---

### 4. Ranking Variable (Unfairness Score)

Instead of SES rank, the function uses:

* predicted probability from:

  ```stata
  svy: logit outcome unfair_factors
  predict rank, pr
  ```

This creates a **multivariate unfairness ranking**, consistent with modern equity analysis.

---

### 5. Common Analytic Sample

All calculations are performed on a **single analytic sample**:

* same observations used for:

  * rank model
  * probit model
  * means
  * CI calculations
  * decomposition

This avoids internal inconsistency.

---

### 6. Reporting Strategy (Hybrid WB + VERSE)

Two types of outputs:

#### ✔ Signed (theoretical decomposition)

* preserves identity
* allows interpretation of direction

#### ✔ Absolute normalized (%)

* easier interpretation
* sums to 100%
* aligns with **VERSE-style reporting**

---

## ⚙️ Features

* ✔ Survey-weighted estimation (`svy:`)
* ✔ Supports:

  * Relative CI (default)
  * Wagstaff correction
  * Erreygers correction
* ✔ Handles omitted / collinear variables
* ✔ Returns:

  * decomposition table
  * analytic sample size
  * CI values
* ✔ Fully reproducible workflow

---

## 🚀 Usage

### Example

```stata
ineqdecomp_unfair , ///
    outcome(skilled_battend) ///
    unfair($X_raw) ///
    wvar(weight_var) ///
    citype(relative) ///
    rankmodel(logit) ///
    clear
```

---

## 📦 Output Variables

| Variable                  | Description                    |
| ------------------------- | ------------------------------ |
| `elasticity`              | Elasticity of unfair factor    |
| `var_ci`                  | CI of unfair factor            |
| `contribution`            | Signed contribution            |
| `contribution_abs`        | Absolute contribution          |
| `contribution_pct_signed` | Signed % contribution          |
| `contribution_pct_abs`    | Absolute %                     |
| `contribution_pct_norml`  | Normalized % (sum = 100)       |
| `outcome_ci`              | CI of outcome                  |
| `n_analytic`              | Number of observations used    |
| `status`                  | estimated / omitted / residual |

---

## ⚠️ Important Notes

* Use **numeric variables only** (no `i.` factor notation)
* Data must be `svyset`
* For binary outcomes:

  * `citype(relative)` is recommended for main analysis
* Wagstaff / Erreygers:

  * recommended as sensitivity analysis

---

## 📚 References

### World Bank (Core Decomposition Framework)

O'Donnell et al. (2008)
**Analyzing Health Equity Using Household Survey Data**
Chapter 13: Decomposition of inequality

🔗 [https://documents.worldbank.org/en/publication/documents-reports/documentdetail/633931468139502235](https://documents.worldbank.org/en/publication/documents-reports/documentdetail/633931468139502235)

---

### VERSE Project (Equity & Decomposition Toolkit)

GitHub:
🔗 [https://github.com/VERSE-Equity/Toolkit-DHS/tree/main](https://github.com/VERSE-Equity/Toolkit-DHS/tree/main)

Paper:

> Verguet et al. (2022).
> *Measuring equity in health service use and health outcomes in LMICs*
> 🔗 [https://www.sciencedirect.com/science/article/pii/S0277953622002854](https://www.sciencedirect.com/science/article/pii/S0277953622002854)

---

### R Package for CI Calculation

`rineq` package (CRAN):

🔗 [https://cran.r-project.org/web/packages/rineq/index.html](https://cran.r-project.org/web/packages/rineq/index.html)

Used as reference for:

* CI formula implementation
* consistent weighting approach
* bounded CI corrections

---

## 🧩 Methodological Positioning

This function represents a **hybrid approach**:

| Component                   | Source              |
| --------------------------- | ------------------- |
| Decomposition identity      | World Bank          |
| Elasticity (AME-based)      | Nonlinear extension |
| Reporting (absolute shares) | VERSE               |
| CI calculation logic        | rineq               |

---

## 📌 When to Use This Tool

Use this when you want to:

* Understand **drivers of inequality**
* Decompose inequality using **nonlinear models**
* Apply **equity-focused analysis beyond SES ranking**
* Produce **policy-relevant decomposition results**

---
