# FunBurd 
FunBurd is designed to test the association between variants aggregated across a gene set and a given trait. This repository contains all code and data required to reproduce the figures and statistical analyses presented in the project.

<p align="center">
  <img src="/FunBurd_Logo.png" alt="Project Logo" width="350"/>
</p>

https://github.com/user-attachments/assets/31558116-877f-43d1-87e5-a253ffa3e561

## **Project Title**  
**Determinants of Pleiotropy and Monotonic Gene Dosage Responses Across Human Traits**

## Project Workflow

<p align="center">
  <img src="/Project_Workflow.jpeg" alt="Project Workflow Diagram" width="900"/>
</p>

## Citation
If you use this project or its code in your research or pipelines, please cite the corresponding paper. Your citation helps others discover the project and acknowledges our work.

- **Paper:** Sayeh Kazem, Kuldeep Kumar, et al. "Gene dosage architecture across complex traits." medRxiv 2025.02.25.25322833; 2025.
- **DOI Badge:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17038354.svg)](https://doi.org/10.5281/zenodo.17038354)

## Repository Contents

### Upstream Analysis Pipelines
- **FunBurd (Python notebook):** FunBurd is designed to test the association between variants aggregated across a gene set and a given trait. The traits of interest were considered as a function of the number of genes within the gene set disrupted by CNVs. To avoid effect size inflation, due to multigenic CNVs, we adjusted for the number of genes (not members of the gene set) disrupted by the same CNV. We also adjusted for age, sex, and ancestry:
<p align="center">
  <img src="/FunBurd_RegressionModel.png" alt="Regression Model" width="600"/>
</p>

- **S-LDSC:** Scripts to recompute GWAS enrichments for the geneset of interest using stratified LD score regression.

### Downstream Analysis & Figures

The following scripts and data generate the main figures for the study:

- **Figure 2:** Heatmap of effect sizes for whole-body tissues and cell-type gene sets across traits.
- **Figure 3:** Dissecting pleiotropy, gene function, and genetic constraint.
- **Figure 4:** Rare and common variant architectures across complex traits.
- **Figure 5:** Gene dosage responses across traits.

### Other Codes
- **Jaccard-Based p-value (Python notebook):** The robust method for assessing gene set associations. The approach conditions on the degree of overlap between gene sets using the Jaccard index. This method helps to avoid inflated significance that can arise from gene sets with high degrees of shared genes.</p>  

### File Types

- Python Notebooks (`*.ipynb`): These notebooks contain the code for the FunBurd pipeline and the Jaccard-based p-value computations.

- R Scripts (`*.Rmd`, `*.md`): These files are used for all downstream analysis, including generating the study's key figures and performing statistical tests.

- R Data Files (`*.RData`): This directory contains all the binary data files needed to reproduce the figures and analyses presented in the project.



## Usage

1. Clone the repository:

```bash
git clone https://github.com/SayehKazem/FunBurd.git




