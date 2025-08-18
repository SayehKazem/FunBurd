# FunBurd 
FunBurd is designed to test the association between variants aggregated across a gene set and a given trait. This repository contains all code and data required to reproduce the figures and statistical analyses presented in the project.

<p align="center">
  <img src="/FunBurd_Logo.png" alt="Project Logo" width="150"/>
</p>



## **Project Title**  
**Determinants of Pleiotropy and Monotonic Gene Dosage Responses Across Human Traits**

## Project Workflow

<p align="center">
  <img src="/Project_Workflow.jpeg" alt="Project Workflow Diagram" width="900"/>
</p>


## Repository Contents

### Codes and Figures

The following scripts and data generate the main figures for the study:

- **Figure 2:** Heatmap of effect sizes for whole-body tissues and cell-type gene sets across traits.
- **Figure 3:** Dissecting pleiotropy, gene function, and genetic constraint.
- **Figure 4:** Rare and common variant architectures across complex traits.
- **Figure 5:** Gene dosage responses across traits.

### Additional Developed Tools
- **FunBurd (Python notebook):** FunBurd is designed to test the association between variants aggregated across a gene set and a given trait. The traits of interest were considered as a function of the number of genes within the gene set disrupted by CNVs. To avoid effect size inflation, due to multigenic CNVs, we adjusted for the number of genes (not members of the gene set) disrupted by the same CNV. We also adjusted for age, sex, and ancestry.
<p align="center">
  <img src="/FunBurd_RegressionModel.png" alt="Regression Model" width="600"/>
</p>
The coefficient B1 represents the effect size of a given gene set on the scaled trait(IRNT) of interest. In this analysis, we examined 14,792 associations, derived from 172 gene sets, 43 traits, and 2 types of variations (deletions and duplications). Our group has recently explored similar burden analysis models for the analysis of a single trait. 


- **Jaccard-Based p-value (Python notebook):** A robust method for assessing gene set associations. The notebook is included for reproducing these calculations.
- **S-LDSC:** Scripts to recompute GWAS enrichments for your geneset of interest using stratified LD score regression.

### File Types

- `*.Rmd` and `*.md` files: Scripts for generating figures and statistical analyses.
- `*.ipynb` files: Python notebooks for Jaccard-based p-value analysis.
- `*.RData`: Binary data files in R format. These files contain all datasets required to reproduce the figures and analyses presented in this project.
  
## Usage

1. Clone the repository:

```bash
git clone https://github.com/SayehKazem/FunBurd.git
