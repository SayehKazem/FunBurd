# FunBurd 
**Project:** *Determinants of pleiotropy and monotonic gene dosage responses across human traits*

<p align="center">
  <img src="/FunBurd_Logo.png" alt="Project Logo" width="150"/>
</p>

FunBurd is designed to test the association between variants aggregated across a gene set and a given trait. This repository contains all code and data required to reproduce the figures and statistical analyses presented in the project.

## Workflow

<p align="center">
  <img src="/Project_Workflow.png" alt="Project Logo" width="500"/>
</p>


## Repository Contents

### Codes and Figures

The following scripts and data generate the main figures for the study:

- **Figure 2:** Heatmap of effect sizes for whole-body tissues and cell-type gene sets across traits.
- **Figure 3:** Dissecting pleiotropy, gene function, and genetic constraint.
- **Figure 4:** Rare and common variant architectures across complex traits.
- **Figure 5:** Gene dosage responses across traits.

### Additional Developed Tools

- **Jaccard-Based p-value (Python notebook):** A robust method for assessing gene set associations. The notebook is included for reproducing these calculations.
- **S-LDSC:** Scripts to recompute GWAS enrichments for your geneset of interest using stratified LD score regression.

### File Types

- `*.Rmd` and `*.md` files: Scripts for generating figures and statistical analyses.
- `*.ipynb` files: Python notebooks for Jaccard-based p-value analysis.
- `data/`: Contains all datasets necessary for reproducing the figures and analyses. All data files are in **RData format**.
  
## Usage

1. Clone the repository:

```bash
git clone https://github.com/SayehKazem/FunBurd.git
