# When Weight Doesn’t Weigh the Same on Everyone

Applying The Lancet’s Clinical Obesity Framework to Map Phenotypic Diagnostic Gaps in NHANES 2021–2023

## Executive Summary

| Section | Content |
|:---|:---|
| **Context** | A 3-notebook analytical series applying The Lancet Commission (2025) Clinical Obesity Framework to the NHANES 2021–2023 data cycle (**N = 1,588** fasting subsample cohort). |
| **Challenge** | Clinicians almost universally rely on Body Mass Index (**BMI ≥ 30 kg/m²**) as the default diagnostic gatekeeper for obesity assessment. This crude anthropometric metric creates a critical diagnostic gap: over-pathologizing metabolically healthy individuals with high BMI while failing to identify normal-BMI patients suffering from undiagnosed end-organ damage. |
| **Strategy** | Synthesize 20 CDC survey modules across 3 standalone analytical notebooks: execute deterministic structural imputation for gatekeeper skip logic, align complex survey architecture using Masked Variance Units (MVUs) and fasting weights, and construct multi-system clinical risk indices (e.g. CKD-EPI eGFR, Fatty Liver Index, METS-IR) to stratify population survey data into obesity phenotypes. |
| **Key Metric** | Robust Prevalence Ratio (PR) for hepatic end-organ damage derived via GEE Poisson Regression (log link) clustered on Masked Variance Units (MVUs) and weighted by Fasting Subsample weights (`WTSAF2YR`). | 
| **Outcome** | Highlighted that elevated METS-IR serves as a important independent risk factor of hepatic end-organ damage, driving distinct prevalence ratios between peripheral and classic phenotypes. |

## 📑 Table of Contents

* [Outcome: METS-IR as a surrogate for hepatic end-organ damage](#-outcome-mets-ir-as-a-surrogate-for-hepatic-end-organ-damage)
* [Technical structure](#-technical-structure)
* [Data Context & Population Flow](#-data-context--population-flow)
* [Complex Survey design](#-complex-survey-design)
* [Risk indices feature engineering](#-risk-indices-feature-engineering)
* [Descriptive Population Stratification by phenotype](#-descriptive-population-stratification-by-phenotype)
* [Regression Framework: Pooled and Stratified Models](#️-regression-framework-pooled-and-stratified-models)
* [Survey-Weighted Trajectory Analysis](#-survey-weighted-trajectory-analysis)

## 📊 Outcome: METS-IR as a driver for hepatic end-organ damage

### Key Analytical Takeaways

1. **The Diagnostic Misclassification of BMI:** Individuals classified as Overweight ($\text{BMI } 25.0\text{--}29.9\text{ kg/m}^2$) bifurcate into the Peripheral Obesity Phenotype, while a metabolically vulnerable subset with normal BMI also transitions into this stratum.

2. **Surrogate Index Performance (METS-IR):** Non-invasive metabolic index (METS-IR) demonstrated superior parameter stability and diagnostic concordance for FibroScan-confirmed hepatic fibrosis compared to other lipid and glycemic ratios.

3. **Phenotypic Risk Tipping Points:** Within the Peripheral Phenotype ($N = 585$), crossing the METS-IR threshold triggers a 3.60-fold increase ($PR = 3.60, p = 0.004$) in hepatic damage risk relative to metabolically normal peripheral controls. In the Classic Phenotype ($N = 731$), high baseline risk is already established, with elevated METS-IR compounding risk by 2.76-fold ($PR = 2.76, p < 0.001$).

## 🔎 Data Context & Population Flow

This project synthesizes data from the NHANES 2021–2023 Cycle released by the CDC National Center for Health Statistics (NCHS). To capture silent end-organ pathology, 20 survey modules spanning Demographics, Examination, Laboratory, and Questionnaire categories were integrated.

![Cohort flowchart](/assets/final_flowchart.svg)

The final target cohort ($n = 1,588$ working-age adults) is established through a multi-stage epidemiological filtering strategy. Rather than executing a destructive flat subsetting of the dataframe, the pipeline isolates this subpopulation to minimize non-response bias, eliminate structural confounding, and guarantee full data density across downstream index calculations:

- **Fasting Subsample Weighting (`WTSAF2YR`):** Restricted strictly to participants dynamically assigned to the morning fasting subsample who met physiological fasting durations and provided valid blood specimens. Utilizing dedicated WTSAF2YR weights is epidemiologically mandatory to counteract selection bias and maintain true population representativeness during GEE variance estimation.

- **Physiological Stabilization Filters:**

    - **Pregnancy Exclusion (RIDEXPRG == 1):** Eliminates acute, non-pathological gestational shifts in body composition, lipid metabolism, and insulin sensitivity that would otherwise confound baseline obesity phenotypes.
    - **Working-Age Restraints ($18 \le \text{Age} < 65$):** Isolates a mature metabolic profile, removing pediatric development variables on the lower bound and minimizing geriatric survivorship bias or competing mortality risks on the upper bound.

- **Complete-Case Feature Integrity:** Ensures zero missingness across mandatory clinical parameters in order to calculate risk indices such as:
    - **AHA PREVENT Equations:** Retaining complete records for cardiovascular-kidney-metabolic (CKM) variables (SBP, anti-hypertensive med history, smoking status, lipid fractions, eGFR).
    - **Metabolic, hepatic and renal indices:** Guaranteeing zero missingness among baseline anthropometric and biochemical components required for composite metabolic markers (e.g., METS-IR, TyG index).
    
### Phenotype Operationalization & Unclassified Exclusion Rule

The obesity_criteria engine operationalizes demographic-adjusted thresholds (WHO Waist Circumference: Men $> 102\text{ cm}$, Women $> 88\text{ cm}$; WHR: Men $> 0.90$, Women $> 0.85$; WHtR $\ge 0.50$) and sums positive central markers (anthro_criteria_count) to assign subjects into discrete clinical tiers:

- **Control Phenotype:** Clean baseline controls with normal BMI ($18.5 \le \text{BMI} < 25$) and zero central adiposity risk criteria ($0$).

- **Peripheral Phenotype:** Normal BMI ($18.5 \le \text{BMI} < 25$) paired with high central visceral deposition ($\ge 2$ criteria).

- **Classic Phenotype:** High BMI ($\text{BMI} \ge 30$) paired with clinical central adiposity ($\ge 1$ criterion).

- **Unclassified Tiers (Excluded from Primary Analysis):** Any observation falling outside these explicit phenotypic boundaries—such as overweight subjects ($\text{BMI } 25.0\text{--}29.9\text{ kg/m}^2$) presenting with 0 or 1 central adiposity criteria—is assigned as "unclassified". These unclassified records are ruled out from primary comparative regression models to prevent boundary blurring and maintain clean contrast against the control baseline.

> **Outliers**  
> Explain why these outliers may exist and their interpretation


## 💻 Technical structure
### Project Schema
```text
├── sources/                    # Raw NHANES 2021-2022 .XPT files
├── notebooks/                  
│   ├── 01_prepare.ipynb        # Batch-loading, data audit and imputation
│   ├── 02_process.ipynb        # Phenotypic Stratification, risk index Construction
│   └── 03_analysis.ipynb       # Survey-Weighted Statistics and GEE Poisson Regressions
├── README.md                   # Project overview
├── requirements.txt            # Python environment dependencies
└── LICENSE                     # Open-source license (e.g., MIT)
```

## 💽 Complex Survey design

To produce statistically valid US population-level inferences, this pipeline embeds probability sampling parameters throughout the data processing and modeling stages.

### Methodological Core: Poisson GEE + Sandwich Variance

Variance estimation in complex survey data requires handling both overdispersion and intra-cluster correlation across Masked Variance Units (MVU). Rather than relying on model-based standard errors, our regression architecture combines a Poisson Mean Link with an Empirical Sandwich Covariance Estimator

### Key Survey Rules Applied

- **Least Common Denominator Weighting:** Because some risk indices like metabolic insulin resistance (METS-IR) and fasting plasma glucose/triglycerides require morning fasting, the 2-year fasting subsample weight (`WTSAF2YR`) was designated as the primary analysis weight for the cohort ($N = 1,588$).

- **Exchangeable Correlation Structure:** Primary sampling unit dependency is modeled using an exchangeable working correlation structure, preventing variance inflation and artificially narrow $p$-values.

## 🧮 Risk indices and feature engineering

To move beyond crude anthropometrics, key variables were transformed into validated clinical risk indices. Mathematical definitions implemented across the pipeline include:

**Metabolic Score for Insulin Resistance (METS-IR)**

$$\text{METS-IR} = \frac{\ln\left((2 \times G_0) + \text{TG}_0\right) \times \text{BMI}}{\ln(\text{HDL}_c)}$$

**Chronic Kidney Disease Epidemiology Collaboration (CKD-EPI 2021)**

$$\text{eGFR} = \mu \times \min\left(\frac{S_{cr}}{\kappa}, 1\right)^{\alpha_1} \times \max\left(\frac{S_{cr}}{\kappa}, 1\right)^{\alpha_2} \times c^{\text{Age}} \times d$$

**Product of Fasting Glucose and Triglycerides (TyG Index)**

$$\text{TyG} = \ln \left( \frac{\text{TG} \times \text{FBG}}{2} \right)$$

### Gatekeeper Skip Logic & Deterministic Imputation

In population survey databases, questionnaire and examination modules utilize Gatekeeper Skip Logic (e.g., if a participant answers "No" to having ever been told they have kidney disease, subsequent questions regarding dialysis treatment are skipped by design).

## 📊 Descriptive Population Stratification by phenotype

To evaluate how health status transitions across body composition phenotypes, population data was stratified across three clinical tiers: Control, Peripheral phenotype, and Classic phenotype.

| Variable | Overall | Control | Peripheral | Classic |
| :- | -: | -: | -: | -: |
| Age (years) | 40.92 | 32.05 | 45.55 | 42.94 |
| Systolic BP (mmHg) | 117.21 | 112.75 | 118.29 | 119.31 |
| Diastolic BP (mmHg) | 74.59 | 68.59 | 74.11 | 78.64 |
| Pulse rate | 70.88 | 69.30 | 70.41 | 72.06 |
| Fasting Glucose (mg/dL) | 105.06 | 98.72 | 103.94 | 111.31 |
| Triglycerides (mg/dL) | 108.60 | 74.17 | 117.46 | 124.36 |
| HDL Cholesterol (mg/dL) | 53.62 | 60.70 | 53.50 |	49.59 |
| Liver Stiffness (kPa) | 5.76 | 4.86 | 5.42 | 6.79 |
| eGFR (mL/min/1.73m²) | 102.51 | 108.14 | 99.57 | 101.41 |
| Hepatic damage (%) | 13.96 | 4.08 | 9.55 | 25.13 |
| Metabolic damage (%) | 6.42 | 0.09 | 7.90 | 9.82 |
| Renal damage (%) | 0.99 | 0.0 | 0.93 | 1.54 |

## 🗃️ Biomarker Distribution & Phenotypic Separation

To evaluate how risk indices separate across the three clinical tiers, survey-weighted mean distributions were compared using design-adjusted Wald F-tests.

![Indices ridgelines](/assets/plot_indices_ridgeline.png)

| Index | Control (Mean ± SD) | Peripheral (Mean ± SD) | Classic (Mean ± SD) | Wald F-Test (p-value) |
| :- | -: | -: | -: | -: |
| Estimated GFR	| 108.14 (±15.76) | 99.57 (±16.33) | 101.41 (±18.15) | < 0.001
| Triglycerides/Glucose | 8.09 (±0.48) | 8.57 (±0.58) | 8.70 (±0.57) | < 0.001 |
| Atherogenic Plasma Index | 0.06 (±0.22) | 0.30 (±0.29) | 0.37 (±0.26) | < 0.001 |
| METS-IR Score | 30.45 (±3.41) | 38.67 (±5.23) | 55.37 (±10.86) | < 0.001 |
| FIB4 | 0.73 (±0.42) | 0.96 (±0.69) | 0.76 (±0.42) | < 0.001 |
| Fatty Liver Index | 9.29 (±7.39) | 41.25 (±22.31)	| 82.72 (±16.98) | < 0.001 |
| Fibrotic NASH Index | 0.06 (±0.09) | 0.10 (±0.12) | 0.13 (±0.16) | < 0.001 |
| PREVENT total risk | 0.96 (±2.02) | 3.38 (±3.94) | 3.37 (±4.14) | < 0.001 |

### Ridgeline density insights

- **Fatty Liver Index (FLI):** Demonstrates near-complete multimodal separation between phenotypes. While Control sits in a tight lower boundary ($< 20$) and Peripheral spans a broad intermediate range, Classic exhibits a dense shift toward the extreme upper boundary ($\ge 80$).

- **METS-IR Score:** Highlights a progressive rightward density migration, shifting from a narrow distribution centered near $30$ (Control) to a wide, right-skewed distribution centered above $55$ (Classic).

*All indices demonstrate statistically significant phenotypic separation ($p < 0.001$ via survey-adjusted Wald F-tests).*

## 🛠️ Regression Framework: Pooled and Stratified Models

### Methodological Justification for Stratification

When evaluating interaction terms between physical phenotypes and metabolic surrogates in pooled regression architectures ($\text{Outcome} \sim \text{Phenotype} \times \text{METS-IR}$), standard models encounter Numerical Non-Convergence:

- **The Control Reference Cell Failure:** In the Control (non-obese) phenotype ($N = 272$), elevated METS-IR ($> 50.39$) exhibits near-zero positivity ($0.4\%$, $n \approx 1$).

- **Variance Explosion:** Because the dummy-coded reference group main effect (METS-IR cutoff) evaluates risk strictly within Controls, the absence of exposed events causes the Hessian matrix to become singular. Standard errors explode to $\text{SE} \approx 10^6$, rendering $p = 1.000$ and invalidating multiplicative interaction parameters.

- **The Stratified Solution:** Stratifying models by clinical phenotype isolates parameter estimation within each target population ($N_{\text{Peripheral}} = 585$, $N_{\text{Classic}} = 731$). This bypasses reference-cell sparsity, allows covariate slopes to vary freely, and produces robust Relative Risk ($\text{RR}$) estimates.

## 🏹 Survey-Weighted Trajectory Analysis

Survey-weighted Sankey diagrams map the population flow from standard anthropometrics through intermediate metabolic drivers to multi-organ damage endpoints.

### Framework 1: Standard BMI $\rightarrow$ Clinical Phenotype $\rightarrow$ End-Organ Targets

![BMI trajectory](/assets/plot_framework_sankey_bmi.png)

- **Phenotypic Bifurcation:** Overweight participants ($\text{BMI } 25.0\text{--}29.9\text{ kg/m}^2$) feed entirely into the Peripheral Obesity Phenotype, accompanied by a metabolically vulnerable sub-stream from the Normal Weight tier.

- **Target Concentration:** While the Control cohort flows exclusively to No Organ Damage, the Classic Phenotype acts as the primary source for Multi-Organ Dysfunction ($\ge 2$ domains), alongside isolated Hepatic Damage (+).

### Framework 2: BMI Category $\rightarrow$ Visual Phenotype $\rightarrow$ METS-IR Risk $\rightarrow$ Hepatic Damage

![METS-IR Trajectory](/assets/plot_framework_sankey_mets-ir.png)

- **Metabolic Gateway Effect:** High METS-IR ($> 50.39$) acts as an important funnel. The Classic Phenotype accounts for the dominant flow into the High METS-IR node, which subsequently captures the majority of downstream Hepatic Damage (+) cases.

- **Peripheral Protection Limit:** The majority of Peripheral subjects remain within the Normal METS-IR band ($\le 50.39$), escaping overt hepatic injury.

### Framework 3: Phenotype $\rightarrow$ Domain Clustering $\rightarrow$ End-Organ Targets

![Domain Trajectory](/assets/plot_framework_sankey.png)

- **Pathology Accumulation:** The Peripheral Phenotype primarily fuels the Single-Domain Risk intermediate tier, terminating in isolated hepatic or metabolic strain. Conversely, the Classic Phenotype feeds the Multi-Domain High Risk intermediate tier, driving systemic end-organ damage.

These visual trajectories reinforce that obesity phenotypes coupled with continuous metabolic insulin resistance escalates risk from isolated tissue strain to multi-system organ damage.