# ABCD-Bio-and-Brain-Age  
Sample Code from Biological and Brain Age Project using the Adolescent Brain and Cognitive Development Study

**NB:** Only includes selected files and no data because of data privacy.

---

## Overview

This repository contains sample code and analysis scripts for the Brain and Biological Age project using data from the Adolescent Brain and Cognitive Development (ABCD) study. The project focuses on cleaning, processing, and modeling biomarkers to generate Brain Development Index (BRDI) and Biological Development Index (BIDI) estimates, and to explore their associations with physical, mental, and academic outcomes in adolescents.

---

## Repository Structure

The analysis pipeline is organized into several stages—from initial data cleaning to advanced modeling and outcome evaluation. The repository is split into two major workflow versions:

- **Initial Workflow:** Prepares the data by cleaning and standardizing biomarkers.
- **Workflow Version 5.0:** Contains the modeling pipelines for both brain and biological age, including training models, feature selection, and outcome predictions.

---

## Workflow Details

### Initial Processing

- **Initial.rmd**  
  - **Purpose:** Reads in biomarkers, selects, cleans, and standardizes them.  
  - **Inputs:** All biomarker files, site ID file, family ID file.  
  - **Output:** `Biomarkers.csv` (list of biomarkers ready for modeling).

- **ParentQuestionnaire.rmd**  
  - **Purpose:** Processes parent CBCL data to split GMV and cortical thickness into healthy and unhealthy groups.  
  - **Inputs:**  
    - `abcd_cbcls01.txt`  
    - `biomarkers.csv`  
    - `smrip_gmv_y2.csv` (GMV)  
    - `smrip_cortthick_y2.csv` (Cortical Thickness)  
  - **Outputs:**  
    - `bm_brain_cbcl_67.csv` (combined GMV and CBCL data)  
    - `data/healthy_cbcl_brain_67.csv` (healthy: below threshold)  
    - `unhealthy_cbcl_brain_67.csv` (unhealthy: above threshold)

- **PhysHealth.rmd**  
  - **Purpose:** Splits GMV and cortical thickness measures into physically healthy and unhealthy groups.  
  - **Inputs:**  
    - `abcd_lpmh01.csv` (physical diagnoses)  
    - `biomarkers.csv`  
    - `smrip_gmv_y2.csv` (GMV)  
  - **Output:** Four files corresponding to healthy and unhealthy classifications for both GMV and CBCL data.

- **PhysMentHealth.rmd**  
  - **Purpose:** Processes GMV and cortical thickness data based on both mental (via CBCL) and physical health criteria.  
  - **Input:** GMV and cortical thickness data categorized by health status.  
  - **Output:** Healthy and unhealthy brain and bio age files.

- **KDM_cbcl.rmd**  
  - **Purpose:** Computes KDM Biological Age for healthy and unhealthy groups separately.  
  - **Inputs:**  
    - `healthy_cbcl_phys.csv`  
    - `unhealthy_cbcl_phys.csv`  
  - **Outputs:**  
    - R model objects with fit statistics  
    - `year2_cs.csv` (includes test set participants for correlation analyses)

---

### Workflow Version 5.0

#### Training Model (Brain Age)
- **BrainAge_GMV_V5.0.Rmd**
- **CBCL_5.0.Rmd**
- **PhysHealth_5.0.Rmd**
- **BrainAge_HealthyUnhealthy_Merge_5.0.Rmd**
- **BrainAge_FeatureSelection_5.0.Rmd**  
  *Uses backwards stepwise selection.*
- **KDM_GMV_5.0_V1.Rmd**
- **KDM_THK_5.0_V1.qmd**
- **KDM_GMV_THK_5.0_V1.qmd**
- **KDM_BrainAge_Projections_5.0.Rmd**  
  *Includes 2- and 4-year projections for the brain age estimates.*

#### Training Model (Biological Age)
- **Initial.Rmd**
- **CBCL_5.0.Rmd**
- **PhysHealth_5.0.Rmd**
- **BioAge_HealthyUnhealthy_Merge_5.0.Rmd**
- **BioAge_FeatureSelection_5.0_V2.Rmd**  
  *Uses backwards stepwise selection and outputs separate files for organ age.*
- **KDM_BBB_5.0_V1.Rmd**
- **KDM_Organ_5.0_V1_Temp.qmd**  
  *Template for running KDM organ age; outputs merged organ age dataframes.*

#### Bioage–Brainage Comparisons
- **Bioage_Brainage_Comparisons_V5.Rmd**  
  *Generates graphs comparing biological and brain ages; separates ages into tertiles.*  
  - **Output:** `gmv_bbb_all.csv`
- **Heat_Scatter_Plots.qmd**  
  *Creates scatterplots for biological, brain, and organ ages.*
- **MAE_BBB_Organ_Comp_Graphs.qmd**  
  *Creates barplots for MAE and correlation metrics for biological, brain, and organ ages.*
- **SEM_Organ_KDM.qmd**  
  *Performs SEM analyses on the calculated KDM methods; outputs merged Bio, Brain, and Organ KDM files.*
- **BNA_Organ_KDM.qmd**

---

### Outcome Analyses

#### KDM Projections
- **KDM_GMV_Projections_5.0.Rmd**
- **KDM_BBB_Projections_5.0.Rmd**
- **KSADS_5.0.Rmd**  
  *Note: Currently only implemented for Brain Age (insufficient data at Year 4; pending ABCD update).*

#### Additional Outcome Files
- **CBCL_Outcome_5.0.Rmd** – Analyzes all CBCL outcomes.
- **Mania_5.0.Rmd** – Analyzes mania outcomes.
- **Barkley_NIHT_Pearson_Outcome_5.0.Rmd**
- **SchoolGrades_Outcome_5.0.qmd**
- **Alcohol_5.0.qmd**
- **MedHistory_Final.Rmd**
- **ScreenTime_Final.Rmd**

#### Outcome Summary Files
- **Neurocog_Summary_5.0.Rmd**  
  *Summarizes outcomes for CBCL, Mania, and KSADS.*
- **Social_Summary_5.0.Rmd**  
  *Summarizes outcomes for school grades, media use, alcohol/substance use, and physical activity.*

---

## How to Use

1. **Setup:**  
   Ensure that all necessary input files (e.g., biomarker data, CBCL files, physical diagnoses) are placed in their respective directories.

2. **Execution:**  
   Run the R Markdown (`.Rmd`) and Quarto (`.qmd`) documents in the order indicated in the workflow. These scripts process the data, build the models, and generate various outputs.

3. **Outputs:**  
   Outputs include CSV files, model objects, and plots, which will be saved in their designated output directories.

---

## Contributing

Contributions are welcome! Please fork the repository and submit pull requests with your improvements or bug fixes. For any major changes, open an issue first to discuss your ideas.

---

## License

[Insert license information here.]

---

This repository is structured to provide a comprehensive pipeline for assessing adolescent brain and biological development using state-of-the-art modeling techniques. For any questions or further assistance, please contact the repository maintainer.
