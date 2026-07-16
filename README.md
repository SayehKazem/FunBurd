# FunBurd 
FunBurd is designed to test the association between variants aggregated across a gene set and a given trait. This repository contains all code and data required to reproduce the figures and statistical analyses presented in the project.

<p align="center">
  <img src="/FunBurd_Logo.png" alt="Project Logo" width="300"/>
</p>

https://github.com/user-attachments/assets/31558116-877f-43d1-87e5-a253ffa3e561

## **Project Title**  
**Determinants of functional burden pleiotropy and gene dosage responses across human traits**

## Project Workflow

<p align="center">
  <img src="/FunBurd_Workflow.jpeg" alt="Project Workflow Diagram" width="1000"/>
</p>

## Citation
If you use this project or its code in your research or pipelines, please cite the corresponding paper. Your citation helps others discover the project and acknowledges our work.

- **Paper:** Sayeh Kazem, Kuldeep Kumar, et al. "Determinants of functional burden pleiotropy and gene dosage responses across human traits" [medRxiv 2025.02.25.25322833; 2025.](https://www.medrxiv.org/content/10.1101/2025.02.25.25322833v3) 
- **DOI Badge:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17369429.svg)](https://doi.org/10.5281/zenodo.17369429)

## Repository Contents

### Upstream Analysis Pipelines
- **[`Functional_Burden_Association_Analysis/`](Functional_Burden_Association_Analysis) — FunBurd:** FunBurd is designed to test the association between variants aggregated across a gene set and a given trait. The traits of interest were considered as a function of the number of genes within the gene set disrupted by CNVs. To avoid effect size inflation, due to multigenic CNVs, we adjusted for the number of genes (not members of the gene set) disrupted by the same CNV. We also adjusted for age, sex, and ancestry:
<p align="center">
  <img src="/FunBurd_RegressionModel.png" alt="Regression Model" width="600"/>
</p>

Provided as both a notebook (`FunBurd_Multitrait_EffectSizes_Computation_Script_On_ToyDatasets.ipynb`) and an equivalent command-line script (`.py`), so the same code can be run interactively or launched as a batch job. Since UK Biobank / All of Us genotypes can't be redistributed, the folder ships with a synthetic **toy dataset** (`Genotypic_info_Toy/`, `Phenotypic_info_Toy/`, `Outcome_Toy/`) so the full pipeline can be run end-to-end, plus the three **functional gene set collections** used throughout the paper (`Gene_Sets/`): 60 whole-body tissue sets, 81 whole-body cell-type sets, and 31 whole-brain cell-type (Siletti) sets. `Slurm_FunBurd.sh` and `Args_File_FunBurd.txt` show how to launch the script as a SLURM/`glost` array job across many trait × gene-set combinations.

### Downstream Analysis & Figures
The following scripts and data ([`Figure_Generation/`](Figure_Generation)) generate the main figures for the study:

- **Figure 2:** Tissue-specific associations of CNVs with complex traits.
- **Figure 3:** Whole-body cell type-specific associations of CNVs with complex traits.
- **Figure 4:** Dissecting pleiotropy, gene function, and genetic constraint.
- **Figure 5:** Rare and common variant architectures across complex traits.
- **Figure 6:** Gene dosage responses across traits and functional gene sets.

Each `FigureN.Rmd` / `FigureN.md` pair is a self-contained, rendered analysis (data loading → statistics → plotting), with output plots saved alongside in the matching `FigureN_files/` folder.

### Other Codes
- **[`Permutation_PreservingJaccardDistance_P_Jaccard/`](Permutation_PreservingJaccardDistance_P_Jaccard) — P-Jaccard:** A permutation-based method for assessing gene set associations that conditions on the degree of overlap between gene sets using the Jaccard index, avoiding inflated significance for gene sets that share many genes. For each trait, the method generates 1,000 Jaccard-preserving surrogate gene sets (via `BrainSMASH`) to build a null distribution for the deletion–duplication burden correlation, producing two p-values per trait — one fixing deletions and surrogating duplications, one fixing duplications and surrogating deletions — which are then FDR-corrected (`PJaccard/`, `JaccardDist/`, `surrogate_maps_DelDup/`).

### File Types

- **Python** (`*.ipynb`, `*.py`): FunBurd effect-size computation and P-Jaccard computation.
- **R Markdown / Markdown** (`*.Rmd`, `*.md`): Downstream statistical analyses and figure generation.
- **R Data Files** (`*.RData`): Processed data underlying each figure.
- **Tabular / text data** (`*.tsv`, `*.csv`, `*.txt`): Toy genotype/phenotype inputs, gene set definitions, and effect-size / P-Jaccard outputs.
- **NumPy arrays** (`*.npy`): Jaccard-preserving surrogate maps used to build null distributions.
- **Cluster scripts** (`*.sh` + args files): SLURM/`glost` templates for running FunBurd across many traits and gene sets in parallel.

## Usage

1. Clone the repository:

```bash
git clone https://github.com/SayehKazem/FunBurd.git
```

2. Try FunBurd on the provided toy dataset:

```bash
cd FunBurd/Functional_Burden_Association_Analysis
python3 FunBurd_Multitrait_EffectSizes_Computation_Script_On_ToyDatasets.py \
    Phenotypic_info_Toy/phenotype1_toy.tsv phenotype1 phenotype1 \
    Gene_Sets/GeneSet_HPA_tissue_fantom_nTPM HPA_tissue_fantom_nTPM
```

   The script takes five positional arguments — phenotype file, phenotype column, phenotype name, gene-set directory, and gene-set name — and writes effect sizes to `Outcome_Toy/`. See `Args_File_FunBurd.txt` for the full list of trait × gene-set combinations, and `Slurm_FunBurd.sh` for running them all as a parallel batch job.

### Requirements
- **Python (FunBurd):** `numpy`, `pandas`, `statsmodels`, `patsy`, `firthlogist`
- **Python (P-Jaccard):** `numpy`, `pandas`, `scipy`, `statsmodels`, `matplotlib`, `nibabel`, `brainsmash`
- **R (Figure Generation):** `ggplot2`, `ggpubr`, `ggprism`, `ggcorrplot`, `ggraph`, `ggrepel`, `ComplexHeatmap`, `circlize`, `igraph`, `vegan`, `dplyr`, `tidyr`, `tibble`, `forcats`, `stringr`, `reshape2`, `patchwork`, `gridExtra`, `openxlsx`, `svglite`, `knitr`, `MASS`, `Rmpfr`

## License
This project is released under the [MIT License](LICENSE).
