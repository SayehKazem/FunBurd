#!/usr/bin/env python
# coding: utf-8

# --- 1. Library Imports ---
import os
import sys
import time
import numpy as np
import pandas as pd
from multiprocessing import Pool
import statsmodels.formula.api as sm
from patsy import dmatrices
from firthlogist import FirthLogisticRegression

os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

# --- 2. Configuration Parameters ---
# PHENOTYPE_FILE = "Phenotypic_info_Toy/phenotype1_toy.tsv"
# PHENOTYPE_COLUMN = "phenotype1"
# PHENOTYPE_NAME = "phenotype1"

PHENOTYPE_FILE = sys.argv[1]
PHENOTYPE_COLUMN = sys.argv[2]
PHENOTYPE_NAME = sys.argv[3]
GENE_SET_DIRECTORY = sys.argv[4]
GENE_SET_NAME = sys.argv[5]

COVARIATES = ' + age + Sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10'

# File paths for input and output data.
CURRENT_DIR = os.getcwd()

PATH_GENE_FILES = os.path.join(CURRENT_DIR, GENE_SET_DIRECTORY)
PATH_GENE_INFO = os.path.join(CURRENT_DIR,"Genotypic_info_Toy", "cnv_gene_df_toy.tsv")
PATH_CNV_INDIVIDUAL = os.path.join(CURRENT_DIR,"Genotypic_info_Toy", "cnv_individual_df_toy.tsv")
PATH_CLEAN_DATA = os.path.join(CURRENT_DIR, "Phenotypic_info_Toy","clean_individual_df_toy.tsv")

OUTPUT_PATH = os.path.join(CURRENT_DIR,'Outcome_Toy', f"Effect_sizes_{GENE_SET_NAME}_{PHENOTYPE_NAME}.tsv")

# Number of CPU cores for parallel processing.
CPU_COUNT = int(os.environ.get("SLURM_CPUS_PER_TASK", 4))

print("Notebook setup complete. Proceed to the next cell to load data.")




print("Loading all required data files...")

cnv_individual_df = pd.read_csv(PATH_CNV_INDIVIDUAL, sep='\t')
cnv_individual_df['individual'] = cnv_individual_df['individual'].str[:7].astype(str)

phenotype_individual_df = pd.read_csv(PHENOTYPE_FILE, sep='\t')
phenotype_individual_df['individual'] = phenotype_individual_df['individual'].astype(str)

cnv_gene_df = pd.read_csv(
    PATH_GENE_INFO,
    sep='\t', usecols=['CHR', 'START', 'STOP', 'TYPE', 'proportion_gene_overlap', 'gene_id']
)

clean_individual_df = pd.read_csv(PATH_CLEAN_DATA, usecols=["individual", "Ancestry"],sep='\t')
clean_individual_df['individual'] = clean_individual_df['individual'].str[:7].astype(str)

# Filter for individuals with 'EUR' ancestry and merge with phenotype data
phenotype_individual_df = pd.merge(phenotype_individual_df, clean_individual_df[clean_individual_df['Ancestry'] == 'EUR'], on='individual', how='inner')

print("Data loading complete.")




# ==================================================================================================
## Data Processing and Aggregation
# ==================================================================================================

def annotate_gene_set_status(data, gene_column_name, gene_list):
    """
    Annotates genes as being either 'inside' or 'outside' a specified gene list.

    Args:
        data (pd.DataFrame): DataFrame containing gene information.
        gene_column_name (str): The name of the column containing gene identifiers.
        gene_list (list): A list of gene identifiers to check against.

    Returns:
        pd.DataFrame: The original DataFrame with two new columns, 'inside' and 'outside',
                      indicating gene set membership (1 for true, 0 for false).
    """
    data["inside"] = data[gene_column_name].isin(gene_list).astype(int)
    data["outside"] = (~data[gene_column_name].isin(gene_list)).astype(int)
    return data

def aggregate_cnv_counts(cnv_individual_df, cnv_gene_df, phenotype_individual_df):
    """
    Counts CNV events for each individual based on genes inside and outside the set.

    Args:
        cnv_individual_df (pd.DataFrame): DataFrame with individual-level CNV data.
        cnv_gene_df (pd.DataFrame): DataFrame with gene-level CNV data and annotations.
        phenotype_individual_df (pd.DataFrame): DataFrame with individual phenotype data.

    Returns:
        pd.DataFrame: A DataFrame with CNV counts (inside/outside for DEL/DUP) aggregated
                      by individual, merged with phenotype data.
    """
    cnv_gene_del_df = cnv_gene_df[cnv_gene_df["TYPE"] == "DEL"]
    cnv_gene_dup_df = cnv_gene_df[cnv_gene_df["TYPE"] == "DUP"]
    
    cnv_individual_del_df = cnv_individual_df[cnv_individual_df["TYPE"] == "DEL"]
    cnv_individual_dup_df = cnv_individual_df[cnv_individual_df["TYPE"] == "DUP"]

    def process_cnv_type(cnv_subset, gene_subset):
        """Helper function to aggregate counts for a single CNV type (DEL or DUP)."""
        gene_subset = gene_subset.copy()
        gene_subset['position'] = gene_subset.apply(lambda row: f"{row['CHR']}--{row['START']}--{row['STOP']}", axis=1)
        gene_subset.set_index('position', inplace=True)
        gene_counts = gene_subset[['inside', 'outside']].groupby(level="position").sum()

        cnv_subset = cnv_subset.copy()
        cnv_subset['position'] = cnv_subset.apply(lambda row: f"{row['CHR']}--{row['START']}--{row['STOP']}", axis=1)
        cnv_positions = cnv_subset[['position', 'individual']]

        merged_counts = pd.merge(cnv_positions, gene_counts, on="position", how="outer")
        individual_counts = merged_counts[['inside', 'outside']]
        individual_counts.index = merged_counts["individual"]
        return individual_counts.groupby(level="individual").sum()

    count_individual_del_cnvs = process_cnv_type(cnv_individual_del_df, cnv_gene_del_df)
    count_individual_dup_cnvs = process_cnv_type(cnv_individual_dup_df, cnv_gene_dup_df)

    final_counts = pd.merge(count_individual_dup_cnvs, count_individual_del_cnvs, on="individual", how="outer", suffixes=('_DUP', '_DEL'))
    
    final_counts.columns = [
        "inside_DUP", "outside_DUP",
        "inside_DEL", "outside_DEL"
    ]

    final_counts = pd.merge(phenotype_individual_df, final_counts, on="individual", how="inner")
    final_counts.fillna(0, inplace=True)
    
    return final_counts


# ==================================================================================================
## Statistical Modeling (Firth Logistic Regression for Binary traits and Linear regression for Continues Traits)
# ==================================================================================================

class FirthLogisticRegressionResults:
    """A class to hold Firth Logistic Regression results for consistent output formatting."""
    def __init__(self, coefficients_series, std_errors_series, p_values_series):
        self.params = coefficients_series
        self.bse = std_errors_series
        self.pvalues = p_values_series

def fit_firth_logistic_regression(formula, data):
    """
    Performs Firth logistic regression with robust data handling and result formatting.

    Args:
        formula (str): A patsy-compatible formula string for the model.
        data (pd.DataFrame): The DataFrame containing the model variables.

    Returns:
        FirthLogisticRegressionResults: An object containing the model's coefficients,
                                        standard errors, and p-values.
    """
    y, X0 = dmatrices(formula, data, return_type='dataframe')
    X = X0.drop('Intercept', axis=1)
    feature_names = X.columns.tolist()
    
    y = y.values.ravel()
    X = X.values
    
    # Firth logistic regression model
    fl0 = FirthLogisticRegression()
    fl = fl0.fit(X, y)
    
    variable_names =  ['Intercept'] + feature_names 
    variable_def =  feature_names +['Intercept']
    
    # Extract coefficients and intercept
    intercept = fl.intercept_
    coefficients = fl.coef_
    coefficients_series = pd.Series(np.concatenate(( [intercept],coefficients)), index=variable_names) 
    # Extract std
    std_errors = fl.bse_
    std_errors_series = pd.Series(std_errors, index=variable_def)
    # Extract pvals
    p_values = fl.pvals_
    p_values_series = pd.Series(p_values, index=variable_def)
    # Reorder the std and pvalues
    std_errors_series = std_errors_series.reindex(variable_names)
    p_values_series = p_values_series.reindex(variable_names)
    
    return FirthLogisticRegressionResults(coefficients_series,std_errors_series,p_values_series)


def fit_model_and_extract_results(data, cnv_type, pheno_col, covariates):
    """
    Fits the appropriate regression model and returns a DataFrame of key results.

    Args:
        data (pd.DataFrame): The DataFrame containing all data for the model.
        cnv_type (str): The CNV type ("DEL" or "DUP").
        pheno_col (str): The name of the phenotype column.
        covariates (str): A string of additional covariates for the formula.

    Returns:
        pd.DataFrame: A DataFrame with the model's effect size, standard error,
                      and p-value. Returns an empty DataFrame on error.
    """
    formula = f"{pheno_col} ~ inside_{cnv_type} + outside_{cnv_type}{covariates}"
    
    try:
        # if pheno_col in ['Final_phenotype_value_info_GP', 'Final_phenotype_info_Depression', 'Final_phenotype_value_Hypertention', 'Final_phenotype_value_MoodSwings', 'Final_phenotype_value_loneliness', 'Final_phenotype_value_RiskTaking', 'Final_phenotype_value_Miserableness', 'Final_phenotype_value_Irritability', 'Final_phenotype_value_GuiltyFeelings', 'Final_phenotype_value_SubjectiveWellBeing']:
        #     regression_results = fit_firth_logistic_regression(formula, data)
        # else:
        #     regression_results = sm.gls(formula, data=data).fit()
        is_binary_trait = np.array_equal(data[pheno_col].dropna().unique(), [0, 1])
        if is_binary_trait:
            # Run Firth logistic regression
            regression_results = fit_firth_logistic_regression(formula, data)
        else:
    # Run GLM for continuous traits
            regression_results = sm.gls(formula, data=data).fit()
        
    except Exception as e:
        print(f"Error running model for {cnv_type}: {e}")
        return pd.DataFrame()

    results_df = pd.DataFrame({
        'Effectsize': regression_results.params.values,
        'se': regression_results.bse.values,
        'pvalue': regression_results.pvalues.values
    })
    
    return results_df

# ==================================================================================================
## Main Analysis Workflow
# ==================================================================================================

def analyze_gene_set_data(list_name):
    """
    Primary worker function that orchestrates the analysis for a single gene list.

    It reads the gene list, processes the data, aggregates CNV counts, and fits
    the appropriate regression models for both DEL and DUP CNV types.

    Args:
        list_name (str): The name of the gene list file to be processed.

    Returns:
        pd.DataFrame: A summary DataFrame of the final regression results for the
                      processed gene list.
    """
    print(f"Processing gene list: {list_name}")
    start_time = time.time()
    
    try:
        gene_list = pd.read_csv(f"{PATH_GENE_FILES}/{list_name}", sep='\t')
        gene_list.columns = ["gene_id"]
        
        gene_info_subset = cnv_gene_df[cnv_gene_df["proportion_gene_overlap"] >= 1].copy()
        gene_info_with_indicators = annotate_gene_set_status(
            gene_info_subset, "gene_id", gene_list.iloc[:, 0].tolist()
        )
        
        combined_counts_df = aggregate_cnv_counts(
            cnv_individual_df,
            gene_info_with_indicators,
            phenotype_individual_df
        )
        
        cleaned_data = combined_counts_df.dropna(subset=[PHENOTYPE_COLUMN])
        
        final_results = []
        
        for cnv_type in ["DEL", "DUP"]:
            results_df = fit_model_and_extract_results(
                data=cleaned_data,
                cnv_type=cnv_type,
                pheno_col=PHENOTYPE_COLUMN,
                covariates=COVARIATES
            )
            
            # Extracts the row corresponding to the 'inside' effect, as per original logic.
            results_df = results_df.iloc[1:2].copy()
            results_df['TYPE'] = cnv_type
            results_df['gene_list_name'] = list_name
            final_results.append(results_df)

        final_list_summary = pd.concat(final_results, ignore_index=True)
        print(f"Finished processing {list_name} in {(time.time() - start_time):.2f} seconds.")
        return final_list_summary
    
    except Exception as e:
        print(f"An error occurred while processing {list_name}: {e}")
        return pd.DataFrame()



if __name__ == '__main__':
    
    print("Beginning the gene set association analysis...")

    list_of_gene_files = [f for f in os.listdir(PATH_GENE_FILES) if f.endswith(".tsv")]
    
    start_all = time.time()
    
    with Pool(CPU_COUNT) as pool:
        results = pool.map(analyze_gene_set_data, list_of_gene_files)
    
    final_results_df = pd.concat(results, ignore_index=True)
    
    if not final_results_df.empty:
        final_results_df.to_csv(OUTPUT_PATH, sep="\t", index=False)
        print(f"Analysis successfully completed. Results saved to: {OUTPUT_PATH}")
    else:
        print("No results were generated. The process may have failed.")
        
    print(f"Total execution time: {(time.time() - start_all):.2f} seconds.")


