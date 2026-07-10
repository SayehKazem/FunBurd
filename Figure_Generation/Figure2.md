Fig2
================
Sayeh Kazem

## Fig. 2: Tissue-specific associations of CNVs with complex traits

#### – Figure legend –

**Legend:** (A) The heatmap displays association effect sizes between
five categories of traits (x-axis) and tissue-specific gene sets
(y-axis) for deletions (left) and duplications (right). Categories are
shown along the x-axis, with gene sets on the y-axis, and their
annotations indicated on the right. Color intensity reflects the
direction and magnitude of the association (blue = negative; red =
positive). Asterisks (*) indicate FDR-significant associations
(two-sided t-test from multiple linear regression; 172 gene sets X 43
traits X 2 CNV types = 14,792 tests). Only a representative subset of
traits is shown (complete heatmaps in Figs. S2-S3). (B-C) Summary of
variant association differences. (B) Whole-Body Comparison: Differences
in aggregate association/enrichment levels across variant types
(Deletions, Duplications, SNPs) for all tissue gene sets. (C) Higher
pleiotropy in brain tissue gene sets: Proportion of significant trait
associations for brain and non-brain functional gene sets. In both
panels, data are presented as the aggregate statistical proportion of
significant enrichments ± standard error. Statistics are derived from
gene sets rather than individual sample observations. For these panels,
n = 60 whole-body tissue gene sets (n = 27 brain, n = 33 non-brain). No
technical or biological replicates are utilized. Asterisks (*) denote
FDR-adjusted statistical significance (two-sided Pearson’s Chi-squared
test). (D) Pearson correlation of gene-set effect sizes between AoU and
UKBB; shaded bands represent 95% confidence intervals, asterisks denote
significant P-Jaccard. (E) Discovery power and cross-cohort replication.
Per-trait concordance of 60 gene-set effect sizes (BMI, Heart Rate,
Platelet Count, Standing Height) vs. number of UKBB FDR-significant gene
sets. Note: Townsend Deprivation Index (TDI) is an area-level measure of
material deprivation, grouped with brain-related traits based on its
established genetic and phenotypic correlations with cognitive and
mental health outcomes. Abbreviations: HDL, high-density lipoprotein;
HbA1c, glycated haemoglobin; BMD, bone mineral density; BMI, body mass
index; EA, educational attainment; FI, fluid intelligence; TDI, Townsend
deprivation index; VSWM, visuospatial working memory; FVC, forced vital
capacity; Del, deletion; Dup, duplication; SNP, single-nucleotide
polymorphism; AoU, All of Us.

#### Libraries

``` r
#### Libraries for Heatmap, Bar-Plots, and Exporting
library(ggplot2)
library(ggprism)
library(grid)
library(ggpubr)
library(gridExtra)
library(tidyr)
library(tibble)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(svglite) # Required for saving ComplexHeatmap as SVG
library(knitr)   # Required for clean markdown table rendering
library(patchwork)
library(stringr)
library(ggrepel)
```

#### Panel A: Tissue-specific associations of CNVs with complex traits (Heatmap)

##### Description: The heatmap displays association effect sizes between five categories of traits (x-axis) and tissue-specific gene sets (y-axis) for deletions (left) and duplications (right).

##### Panel A Data

``` r
# Load the FunBurd associations output for UKBB dataset 
load('df_funburd_ukbb.RData')

# Render a clean markdown table preview for GitHub
kable(head(df_funburd_ukbb), format = "markdown")
```

| Trait       | Trait_Cat | Geneset    | Geneset_Type | Geneset_Cat | Geneset_SubCat | Geneset_Cat_Detailed | CNV_Type | Effectsize |        se |   p_value | FDR_p_value | Effectsize_outside | se_outside | p_value_outside | FDR_p_value_outside |
|:------------|:----------|:-----------|:-------------|:------------|:---------------|:---------------------|:---------|-----------:|----------:|----------:|------------:|-------------------:|-----------:|----------------:|--------------------:|
| Albumin     | NonBrain  | adipocytes | SC_HPA81     | NonBrain    | NA             | Adipocytes           | DUP      | -0.0026390 | 0.0110931 | 0.7730077 |   0.8946567 |         -0.0026465 |  0.0008490 |       0.0018256 |           0.0028934 |
| Albumin     | NonBrain  | adipocytes | SC_HPA81     | NonBrain    | NA             | Adipocytes           | DEL      |  0.0024743 | 0.0149456 | 0.8416826 |   0.9298805 |         -0.0058628 |  0.0015788 |       0.0002045 |           0.0003783 |
| BirthWeight | NonBrain  | adipocytes | SC_HPA81     | NonBrain    | NA             | Adipocytes           | DEL      | -0.0066052 | 0.0149456 | 0.6728936 |   0.8401685 |         -0.0062473 |  0.0020550 |       0.0023652 |           0.0036855 |
| BirthWeight | NonBrain  | adipocytes | SC_HPA81     | NonBrain    | NA             | Adipocytes           | DUP      | -0.0301690 | 0.0110931 | 0.0084059 |   0.0548476 |          0.0001433 |  0.0010699 |       0.8934737 |           0.9049136 |
| BMD         | NonBrain  | adipocytes | SC_HPA81     | NonBrain    | NA             | Adipocytes           | DUP      | -0.0135559 | 0.0110931 | 0.2217055 |   0.4706469 |         -0.0022584 |  0.0010175 |       0.0264508 |           0.0349589 |
| BMD         | NonBrain  | adipocytes | SC_HPA81     | NonBrain    | NA             | Adipocytes           | DEL      |  0.0224277 | 0.0149456 | 0.1334548 |   0.3468137 |         -0.0023796 |  0.0018936 |       0.2088704 |           0.2415268 |

##### Panel A Figure (For Deletion)

``` r
# Load required libraries (ensure these are loaded in your Rmd session)
# library(dplyr)
# library(tidyr)
# library(ComplexHeatmap)
# library(circlize)
# library(grid)

# ==============================================================================
# 1. DEFINE VECTORS FOR CUSTOM ORDERING (35 Tissues x 20 Traits)
# ==============================================================================
tissue_order <- c(
  # Brain Tissues (18)
  "parietal_lobe", "postcentral_gyrus", "occipital_cortex",
  "frontal_lobe", "medial_frontal_gyrus", "temporal_cortex",
  "caudate", "putamen", "globus_pallidus", "thalamus",
  "hippocampus", "amygdala", "substantia_nigra",
  "locus_coeruleus", "pons", "corpus_callosum",
  "cerebellum", "spinal_cord",
  # Non-Brain Tissues (17)
  "tonsil", "spleen", "lymph_node", "heart_muscle",
  "skeletal_muscle", "liver", "kidney", "esophagus",
  "colon", "appendix", "thyroid_gland", "ovary",
  "testis", "endometrium", "placenta", "prostate",
  "seminal_vesicle"
)

trait_order <- c(
  # Blood Assays (5)
  'HDL', 'Platelet', 'Redblood', 'Gamma', 'Triglycerides',
  # Reproductive Factors (2)
  'Menarche', 'Menopause',
  # Physical Measures (4)
  'FVC', 'BMI', 'Hypertention', 'HeartRate',
  # Cognitive Metrics (4)
  'EA', 'FI', 'VSWM1', 'TDI',
  # Mental Health (5)
  'Depression', 'Neuroticism', 'GuiltyFeelings', 'Loneliness', 'MoodAnxiety'
)

# ==============================================================================
# 2. SUBSETTING THE MASTER DATAFRAME
# ==============================================================================
df_FantomTissue = subset(df_funburd_ukbb, Geneset_Type == 'Tissue_Fantom60')
df_FantomTissue_Del = subset(df_FantomTissue, CNV_Type == 'DEL')

df_FantomTissue_Del_EffectSize = subset(df_FantomTissue_Del, select = c('Trait', 'Geneset', 'Effectsize'))
df_FantomTissue_Del_fdrpval = subset(df_FantomTissue_Del, select = c('Trait', 'Geneset', 'FDR_p_value'))

# ==============================================================================
# 3. DATA TRANSFORMATION (LONG TO WIDE FORMAT)
# ==============================================================================
Matrix_FantomTissue_Del_EffectSize <- df_FantomTissue_Del_EffectSize %>%
  pivot_wider(names_from = Trait, values_from = Effectsize) %>%
  column_to_rownames("Geneset")

Matrix_FantomTissue_Del_fdrpval <- df_FantomTissue_Del_fdrpval %>%
  pivot_wider(names_from = Trait, values_from = FDR_p_value) %>%
  column_to_rownames("Geneset")

# ==============================================================================
# 4. ORDERING AND SAFE SUBSETTING (FIXES ROW_SPLIT ERROR)
# ==============================================================================
# Force the matrices to exact dimensions (35x20) using the target vectors.
# Any tissues missing from the raw data will safely become NAs, ensuring row counts match.
Matrix_FantomTissue_Del_EffectSize_ordered <- Matrix_FantomTissue_Del_EffectSize[tissue_order, trait_order]
Matrix_FantomTissue_Del_fdrpval_ordered <- Matrix_FantomTissue_Del_fdrpval[tissue_order, trait_order]

# Convert numerical p-values into string asterisks, safely handling NAs.
Matrix_FantomTissue_Del_fdrpval_ordered <- ifelse(
  !is.na(Matrix_FantomTissue_Del_fdrpval_ordered) & Matrix_FantomTissue_Del_fdrpval_ordered < 0.05, 
  "*", 
  ""
)

# ==============================================================================
# 5. HEATMAP COLOR PALETTE SETUP
# ==============================================================================
colors <- colorRampPalette(c("#386cb0", "white", "#a50f15"))(100)
min_val <- -0.2
breaks <- seq(min_val, -min_val, length.out = length(colors))
color_function <- colorRamp2(breaks, colors)

# ==============================================================================
# 6. ANNOTATIONS AND GROUPING
# ==============================================================================
## ROW ANNOTATIONS (Must be exactly 35 elements to match the forced matrix rows)
row_annotation <- data.frame(
  Tissue = c(rep("Brain", 18), rep("Non-Brain", 17)),
  Cat_Detailed = c(
    rep("Sensory/Visual", 3), rep("Cognitive/Motor", 3), rep("Basal Ganglia", 2), 
    rep("Thalamus/Pallidum", 2), rep("Limbic System", 2), rep("Brainstem", 3), 
    rep("Callosum/Cerebellum/Spinal", 3),#, rep("Spinal Cord", 1),
    rep("Immune", 3), rep("Muscle", 2), rep("Metabolic/Excretory", 2), 
    rep("Digestive", 3), rep("Endocrine/Gonads", 3), rep("Female Repro", 2), 
    rep("Male Repro", 2)
  )
)
row_annotation$Tissue <- factor(row_annotation$Tissue, levels = c("Brain", "Non-Brain"))
row_annotation$Cat_Detailed <- factor(row_annotation$Cat_Detailed, levels = unique(row_annotation$Cat_Detailed))

## COLUMN ANNOTATIONS
column_annotation <- data.frame(
  Trait = c(
    rep("Blood Assays", 5), 
    rep("Reproductive and\n Activity Factors", 2), 
    rep("Physical\n Measure", 4), 
    rep("Cognitive\n Metrics", 4), 
    rep("Mental\n Health", 5)
  )
)
column_annotation$Trait <- factor(column_annotation$Trait, levels = unique(column_annotation$Trait))

column_ha <- HeatmapAnnotation(
  Trait = column_annotation$Trait,
  col = list(Trait = c(
    "Blood Assays" = "#e78ac3", 
    "Reproductive and\n Activity Factors" = "#8da0cb", 
    "Physical\n Measure" = "#a6d854", 
    "Cognitive\n Metrics" = "#66c2a5", 
    "Mental\n Health" = "#fc8d62"
  )),
  annotation_name = NULL, 
  show_legend = FALSE, 
  show_annotation_name = FALSE
)

# ==============================================================================
# 7. HEATMAP CONSTRUCTION
# ==============================================================================

# png("Heatmap_DEL_Fantom_FDR.png",  width = 20000, height = 25000, res = 1600)
# There are 8 detailed categories in Brain and 7 in Non-Brain
num_detailed_splits <- length(unique(row_annotation$Cat_Detailed)) # 15
num_brain_categories <- 7

# Define the gap sizes (adjust these unit values as needed)
small_gap <- unit(0.2, "cm") # Separation between detailed categories
large_gap <- unit(1.2, "cm") # Larger separation between Brain and Non-Brain

# Create the vector of gaps (14 gaps for 15 groups)
row_gaps_list_raw <- list(
  # Gaps within Brain
  rep(small_gap, num_brain_categories - 1), 
  # Large gap between Brain/Non-Brain
  large_gap,                                
  # Gaps within Non-Brain
  rep(small_gap, num_detailed_splits - num_brain_categories - 1) 
)

# Flatten the list and combine all unit objects
row_gap_units <- do.call(unit.c, row_gaps_list_raw)


ht1 <- Heatmap(
  as.matrix(Matrix_FantomTissue_Del_EffectSize_ordered), 
  name = "Effect size",
  na_col = "white",  # Renders any missing tissues cleanly as white boxes
  
  cluster_rows = FALSE,  
  cluster_columns = FALSE,  
  show_row_dend = FALSE,  
  show_column_dend = FALSE,  
  
  top_annotation = column_ha,
  
  row_split = row_annotation$Cat_Detailed, 
  column_split = column_annotation$Trait,  
  
  row_title_rot = 0, 
  row_title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "sans"),
  
  row_names_gp = gpar(col = "black", fontsize = 12, fontfamily = "sans", fontface = "bold"),
  column_title_gp = gpar(col = c('#e78ac3', '#8da0cb', '#a6d854', '#66c2a5', '#fc8d62'), fontsize = 14, fontface = "bold"),
  column_names_gp = gpar(col = "black", fontsize = 12, rot = 45, fontfamily = "sans", fontface = "bold"),
  column_names_rot = 45,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    # Check if the value is not NA before trying to draw an asterisk
    if (!is.na(Matrix_FantomTissue_Del_fdrpval_ordered[i, j])) {
      grid.text(
        Matrix_FantomTissue_Del_fdrpval_ordered[i, j], 
        x = x, 
        y = y - convertHeight(grobHeight(textGrob('*')), "mm"), 
        gp = gpar(fontsize = 36, col = 'grey40', fontfamily = "sans", just = "center", fontface = "bold", hjust=30)
      )
    }
  },
  
  col = color_function, 
  row_gap = row_gap_units,   
  column_gap = unit(0.5, "cm"), 
  show_heatmap_legend = FALSE, 
  rect_gp = gpar(col = "lightgrey", lwd = 0.5) 
)


draw(ht1, padding = unit(c(0.5, 2, 0.5, 0.5), "cm"))
```

![](Fig2_June2026_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

#### Panel B: Proportion of Associations/Enrichmenets

##### Description: Differences in aggregate association/enrichment levels across variant types (Deletions, Duplications, SNPs) for all tissue gene sets.

##### Panel B Data

``` r
# Load Del/Dup/SNP comparison stats 
load('df_compare_del_dup_snp.RData')

# Render a clean markdown table preview for GitHub
kable(head(df_compare_del_dup_snp), format = "markdown")
```

| Geneset_Type | Type1 | Type2 |   p_value |  statistic | estimate_prop_Type1 | estimate_prop_Type2 | p_value_fdr |
|:-------------|:------|:------|----------:|-----------:|--------------------:|--------------------:|------------:|
| Fantom       | SNP   | Del   | 0.0000000 | 185.602972 |           0.0785088 |           0.2201550 |   0.0000000 |
| Fantom       | SNP   | Dup   | 0.0000000 |  41.136684 |           0.0785088 |           0.1364341 |   0.0000000 |
| Fantom       | Del   | Dup   | 0.0000000 |  61.146688 |           0.2201550 |           0.1364341 |   0.0000000 |
| HPA          | SNP   | Del   | 0.0000000 | 197.840003 |           0.0399610 |           0.1418318 |   0.0000000 |
| HPA          | SNP   | Dup   | 0.0000000 | 139.619410 |           0.0399610 |           0.1208728 |   0.0000000 |
| HPA          | Del   | Dup   | 0.0106531 |   6.522292 |           0.1418318 |           0.1208728 |   0.0119847 |

##### Panel B Figure

``` r
# ------------------------------------------------------------------------------
# 1. DEFINE SAMPLE SIZE (N)
# ------------------------------------------------------------------------------
# Based on the Fantom matrix dimensions (43 traits * 60 tissues)
total_pvalues_Fantom <- 43 * 60

# ------------------------------------------------------------------------------
# 2. EXTRACT DATA & COMPUTE STANDARD ERROR
# ------------------------------------------------------------------------------
# Extract the unique proportions from the pairwise comparison table.
df_fantom_props <- bind_rows(
  df_compare_del_dup_snp %>% 
    filter(Geneset_Type == "Fantom") %>% 
    select(Type = Type1, prop = estimate_prop_Type1),
  df_compare_del_dup_snp %>% 
    filter(Geneset_Type == "Fantom") %>% 
    select(Type = Type2, prop = estimate_prop_Type2)
) %>% 
  distinct() %>% 
  mutate(
    # Calculate Standard Error using the provided logic: sqrt(p * (1-p) / N)
    se = sqrt((prop * (1 - prop)) / total_pvalues_Fantom),
    
    # Force the strict X-axis order matching the reference figure
    Type = factor(Type, levels = c("Del", "Dup", "SNP"))
  ) %>%
  arrange(Type)

# ------------------------------------------------------------------------------
# 3. AESTHETICS & COLOR PALETTE
# ------------------------------------------------------------------------------
custom_colors <- c(
  "Del" = "#bd3122",   # Brick red
  "Dup" = "royalblue", # Royal blue
  "SNP" = "#2e6320"    # Dark forest green
)

# ------------------------------------------------------------------------------
# 4. GENERATE VISUALIZATION
# ------------------------------------------------------------------------------
p_panel_b <- ggplot(df_fantom_props, aes(x = Type, y = prop, fill = Type)) +
  
  # Draw the primary bars (no borders)
  geom_bar(stat = "identity", width = 0.5) +
  
  # Overlay Standard Error (SE) bars
  geom_errorbar(
    aes(ymin = prop - se, ymax = prop + se), 
    width = 0.15,      
    linewidth = 0.8,        
    color = "black"    
  ) +
  
  # Apply the custom variant color palette
  scale_fill_manual(values = custom_colors) +
  
  # Define axis titles (X-axis title removed to match reference image)
  labs(
    y = "Proportion of\nAssociations / Enrichments", 
    x = NULL 
  ) + 
  
  # Configure Y-axis scaling (Flush to 0, up to 0.25)
  scale_y_continuous(
    limits = c(0, 0.25), 
    breaks = seq(0, 0.25, by = 0.05), 
    expand = c(0, 0)
  ) +
  
  # Apply a classic theme mimicking GraphPad Prism
  theme_classic(base_size = 16) +
  theme(
    legend.position   = "none",                                  
    axis.line         = element_line(linewidth = 1, color = "black"), 
    axis.ticks        = element_line(linewidth = 1, color = "black"), 
    axis.ticks.length = unit(0.2, "cm"),                         
    axis.text.y       = element_text(size = 16, color = "black", face = "bold"),
    axis.text.x       = element_text(size = 16, color = "black", face = "bold"),
    axis.title.y      = element_text(size = 16, face = "bold", margin = margin(r = 10))
  )

# Render the final plot
print(p_panel_b)
```

![](Fig2_June2026_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

##### Panel B Stats

``` r
# ------------------------------------------------------------------------------
# 1. PREPARE THE DATA
# ------------------------------------------------------------------------------
# Filter for Fantom and keep the EXACT column names from your source data
formatted_stats_table <- df_compare_del_dup_snp %>%
  filter(Geneset_Type == "Fantom") %>%
  select(
    Geneset_Type,
    Type1,
    Type2,
    p_value,
    statistic,
    estimate_prop_Type1,
    estimate_prop_Type2,
    p_value_fdr
  ) %>%
  mutate(
    # Format numerical values to match the clean look of your screenshot
    statistic = round(statistic, 6),
    estimate_prop_Type1 = round(estimate_prop_Type1, 7),
    estimate_prop_Type2 = round(estimate_prop_Type2, 7),
    
    # Format p-values to prevent scientific notation and keep the zeros consistent
    p_value = format(p_value, scientific = FALSE, digits = 7),
    p_value_fdr = format(p_value_fdr, scientific = FALSE, digits = 7)
  )

# ------------------------------------------------------------------------------
# 2. RENDER THE TABLE
# ------------------------------------------------------------------------------
# Generate a clean, Markdown-compatible table that will knit anywhere
kable(
  formatted_stats_table, 
  align = "c",
  caption = "Pairwise Proportion Comparisons (Fantom60)"
)
```

| Geneset_Type | Type1 | Type2 |                    p_value                    | statistic | estimate_prop_Type1 | estimate_prop_Type2 |                 p_value_fdr                  |
|:------------:|:-----:|:-----:|:---------------------------------------------:|:---------:|:-------------------:|:-------------------:|:--------------------------------------------:|
|    Fantom    |  SNP  |  Del  | 0.0000000000000000000000000000000000000000029 | 185.60297 |      0.0785088      |      0.2201550      | 0.000000000000000000000000000000000000000013 |
|    Fantom    |  SNP  |  Dup  | 0.0000000001420000000000000007789750146127198 | 41.13668  |      0.0785088      |      0.1364341      | 0.000000000255999999999999994850981793089438 |
|    Fantom    |  Del  |  Dup  | 0.0000000000000053000000000000001357373025189 | 61.14669  |      0.2201550      |      0.1364341      | 0.000000000000011900000000000000096390308050 |

Pairwise Proportion Comparisons (Fantom60)

#### Panel C: Comparison of Brain vs. Non-Brain Associations

##### Description This section generates a faceted analysis of the Fantom60 geneset, comparing the proportion of significant phenotypic associations between Brain and Non-Brain tissues. The analysis is separated by structural variant type (Deletions, Duplications) and SNPs.

##### Panel C Data

``` r
# Load Del/Dup/SNP comparison stats 
load('df_compare_brain_nonbrain.RData')

# Render a clean markdown table preview for GitHub
kable(head(df_compare_brain_nonbrain), format = "markdown")
```

| Type | Geneset_Type | geneset_cat1 | geneset_cat2 |   p_value |   statistic | estimate_prop_cat1 | estimate_prop_cat2 | p_value_fdr |
|:-----|:-------------|:-------------|:-------------|----------:|------------:|-------------------:|-------------------:|------------:|
| Del  | Fantom       | Brain        | NonBrain     | 0.0000000 | 285.4459309 |          0.3729543 |          0.0951374 |   0.0000000 |
| Del  | HPA          | Brain        | NonBrain     | 0.0036681 |   8.4412379 |          0.1993355 |          0.1363922 |   0.0055021 |
| Del  | Siletti      | Neuron       | NonNeuronal  | 0.7067248 |   0.1415704 |          0.1439646 |          0.1534884 |   0.7067248 |
| Dup  | Fantom       | Brain        | NonBrain     | 0.0000000 |  41.8323609 |          0.1851852 |          0.0965469 |   0.0000000 |
| Dup  | HPA          | Brain        | NonBrain     | 0.2577819 |   1.2806302 |          0.1428571 |          0.1187932 |   0.3314339 |
| Dup  | Siletti      | Neuron       | NonNeuronal  | 0.3472231 |   0.8835796 |          0.1151716 |          0.1348837 |   0.3906260 |

##### Panel C Figure

``` r
# ------------------------------------------------------------------------------
# 1. DEFINE SAMPLE SIZES (N)
# ------------------------------------------------------------------------------
# Fantom Total: 60 Tissues * 43 Traits
# Brain Tissues: 27 Tissues * 43 Traits = 1161
# Non-Brain Tissues: 33 Tissues * 43 Traits = 1419
n_brain <- 27 * 43
n_nonbrain <- 33 * 43

# ------------------------------------------------------------------------------
# 2. EXTRACT DATA & COMPUTE STANDARD ERROR FOR PLOTTING
# ------------------------------------------------------------------------------
# Bind geneset_cat1 (Brain) and geneset_cat2 (NonBrain) into a unified format
df_panel_c_props <- bind_rows(
  
  # Extract Brain proportions 
  df_compare_brain_nonbrain %>% 
    filter(Geneset_Type == "Fantom") %>% 
    select(Type, Tissue_Class = geneset_cat1, prop = estimate_prop_cat1),
  
  # Extract Non-Brain proportions 
  df_compare_brain_nonbrain %>% 
    filter(Geneset_Type == "Fantom") %>% 
    select(Type, Tissue_Class = geneset_cat2, prop = estimate_prop_cat2)
    
) %>% 
  distinct() %>% 
  mutate(
    # Map the correct N based on the Tissue Class
    N = case_when(
      Tissue_Class == "Brain"    ~ n_brain,
      Tissue_Class == "NonBrain" ~ n_nonbrain
    ),
    
    # Calculate Standard Error: sqrt(p * (1-p) / N)
    se = sqrt((prop * (1 - prop)) / N),
    
    # Enforce strict plotting order matching the reference figure
    Type = factor(Type, levels = c("Del", "Dup", "SNP")),
    Tissue_Class = factor(Tissue_Class, levels = c("Brain", "NonBrain"))
  ) %>%
  arrange(Type, Tissue_Class)

# ------------------------------------------------------------------------------
# 3. AESTHETICS & COLOR PALETTE
# ------------------------------------------------------------------------------
# Custom colors carefully matched to the provided Panel C figure
custom_colors_c <- c(
  "Brain"    = "#cf6323", # Rust / Burnt Orange
  "NonBrain" = "#659e38"  # Olive Green
)

# ------------------------------------------------------------------------------
# 4. GENERATE VISUALIZATION
# ------------------------------------------------------------------------------
p_panel_c <- ggplot(df_panel_c_props, aes(x = Tissue_Class, y = prop, fill = Tissue_Class)) +
  
  # Draw the bars (stat="identity" uses our pre-calculated proportions)
  geom_bar(stat = "identity", width = 0.4) +
  
  # Overlay Standard Error (SE) whiskers
  geom_errorbar(
    aes(ymin = prop - se, ymax = prop + se), 
    width = 0.1,      
    linewidth = 0.7,        
    color = "black"    
  ) +
  
  # Apply the custom color palette
  scale_fill_manual(values = custom_colors_c) +
  
  # Split the plot into 3 separate panels based on Variant Type (Del, Dup, SNP)
  # using the exact 'Type' column from your dataframe
  facet_wrap(~ Type) +
  
  # Define axis titles
  labs(
    y = "Proportion of\nAssociations / Enrichments", 
    x = "Genesets" 
  ) + 
  
  # Set Y-Axis limits (Flush to 0, up to 0.45)
  scale_y_continuous(
    limits = c(0, 0.45),      
    breaks = seq(0, 0.4, by = 0.1), 
    expand = c(0, 0)
  ) +
  
  # ----------------------------------------------------------------------------
  # 5. PRISM-STYLE THEME CUSTOMIZATION
  # ----------------------------------------------------------------------------
  theme_classic(base_size = 16) +
  theme(
    legend.position   = "none",                                  
    axis.line         = element_line(linewidth = 1, color = "black"), 
    axis.ticks        = element_line(linewidth = 1, color = "black"), 
    axis.ticks.length = unit(0.2, "cm"),                         
    axis.text.y       = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x       = element_text(size = 14, color = "black", face = "bold"),
    axis.title.y      = element_text(size = 16, face = "bold", margin = margin(r = 10)),
    axis.title.x      = element_text(size = 16, face = "bold", margin = margin(t = 10)),
    
    # Facet Header Styling
    strip.text        = element_text(size = 18, face = "bold"),  
    strip.background  = element_blank()                          
  )

# Render the final plot
print(p_panel_c)
```

![](Fig2_June2026_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

##### Panel C Stats

``` r
# ------------------------------------------------------------------------------
# 1. PREPARE THE DATA
# ------------------------------------------------------------------------------
formatted_stats_c <- df_compare_brain_nonbrain %>%
  
  # Filter strictly for the Fantom geneset
  filter(Geneset_Type == "Fantom") %>%
  
  # Select the exact columns (Removed backticks and spaces to fix parsing error!)
  select(
    Type,                 
    Geneset_Type,         
    Tissue_1 = geneset_cat1, 
    Tissue_2 = geneset_cat2, 
    p_value,
    statistic,
    estimate_prop_cat1,
    estimate_prop_cat2,
    p_value_fdr
  ) %>%
  
  mutate(
    # Clean up numerical formatting for professional presentation
    statistic = round(statistic, 6),
    estimate_prop_cat1 = round(estimate_prop_cat1, 7),
    estimate_prop_cat2 = round(estimate_prop_cat2, 7),
    
    # Format p-values to strictly prevent scientific notation 
    p_value = format(p_value, scientific = FALSE, digits = 7),
    p_value_fdr = format(p_value_fdr, scientific = FALSE, digits = 7)
  )

# ------------------------------------------------------------------------------
# 2. RENDER THE TABLE
# ------------------------------------------------------------------------------
kable(
  formatted_stats_c, 
  align = "c",
  caption = "Brain vs Non-Brain Proportion Comparisons (Fantom60)"
)
```

| Type | Geneset_Type | Tissue_1 | Tissue_2 |                               p_value                                | statistic | estimate_prop_cat1 | estimate_prop_cat2 |                            p_value_fdr                             |
|:----:|:------------:|:--------:|:--------:|:--------------------------------------------------------------------:|:---------:|:------------------:|:------------------:|:------------------------------------------------------------------:|
| Del  |    Fantom    |  Brain   | NonBrain | 0.000000000000000000000000000000000000000000000000000000000000000489 | 285.44593 |     0.3729543      |     0.0951374      | 0.0000000000000000000000000000000000000000000000000000000000000044 |
| Dup  |    Fantom    |  Brain   | NonBrain | 0.000000000099400000000000000545282510228903847041004304685429815436 | 41.83236  |     0.1851852      |     0.0965469      | 0.0000000002240000000000000019569576045237865521975173521695978707 |
| SNP  |    Fantom    |  Brain   | NonBrain | 0.000000000000000000000000000000089099999999999993149546970826837385 | 137.60171 |     0.1520468      |     0.0183413      | 0.0000000000000000000000000000004009999999999999738627199298372546 |

Brain vs Non-Brain Proportion Comparisons (Fantom60)

#### Panel D: FunBurd on All of Us Cohort

##### Description: BMI Concordance Scatter Plots (UKBB vs AoU)

##### Panel D Data (AoU)

``` r
# Load AoU output
load('df_funburd_aou.RData')

kable(head(df_funburd_aou), format = "markdown")
```

| Trait | Geneset         | Geneset_Type    | CNV_Type | Effectsize_AoU |    se_AoU | p_value_AoU | FDR_p_value_AoU |
|:------|:----------------|:----------------|:---------|---------------:|----------:|------------:|----------------:|
| BMI   | caudate         | Tissue_Fantom60 | DEL      |      0.0115241 | 0.0166375 |   0.4885254 |       0.7097988 |
| BMI   | thymus          | Tissue_Fantom60 | DEL      |      0.0711385 | 0.0122504 |   0.0000000 |       0.0000002 |
| BMI   | lymph_node      | Tissue_Fantom60 | DEL      |      0.0307653 | 0.0162883 |   0.0589214 |       0.2062210 |
| BMI   | locus_coeruleus | Tissue_Fantom60 | DEL      |      0.0685991 | 0.0121597 |   0.0000000 |       0.0000006 |
| BMI   | occipital_pole  | Tissue_Fantom60 | DEL      |      0.1095969 | 0.0136947 |   0.0000000 |       0.0000000 |
| BMI   | frontal_lobe    | Tissue_Fantom60 | DEL      |      0.1143842 | 0.0152108 |   0.0000000 |       0.0000000 |

##### Panel D Data (UKBB)

``` r
# Load AoU output
load('df_funburd_ukbb.RData')

kable(head(df_funburd_ukbb), format = "markdown")
```

| Trait       | Trait_Cat | Geneset    | Geneset_Type | Geneset_Cat | Geneset_SubCat | Geneset_Cat_Detailed | CNV_Type | Effectsize |        se |   p_value | FDR_p_value | Effectsize_outside | se_outside | p_value_outside | FDR_p_value_outside |
|:------------|:----------|:-----------|:-------------|:------------|:---------------|:---------------------|:---------|-----------:|----------:|----------:|------------:|-------------------:|-----------:|----------------:|--------------------:|
| Albumin     | NonBrain  | adipocytes | SC_HPA81     | NonBrain    | NA             | Adipocytes           | DUP      | -0.0026390 | 0.0110931 | 0.7730077 |   0.8946567 |         -0.0026465 |  0.0008490 |       0.0018256 |           0.0028934 |
| Albumin     | NonBrain  | adipocytes | SC_HPA81     | NonBrain    | NA             | Adipocytes           | DEL      |  0.0024743 | 0.0149456 | 0.8416826 |   0.9298805 |         -0.0058628 |  0.0015788 |       0.0002045 |           0.0003783 |
| BirthWeight | NonBrain  | adipocytes | SC_HPA81     | NonBrain    | NA             | Adipocytes           | DEL      | -0.0066052 | 0.0149456 | 0.6728936 |   0.8401685 |         -0.0062473 |  0.0020550 |       0.0023652 |           0.0036855 |
| BirthWeight | NonBrain  | adipocytes | SC_HPA81     | NonBrain    | NA             | Adipocytes           | DUP      | -0.0301690 | 0.0110931 | 0.0084059 |   0.0548476 |          0.0001433 |  0.0010699 |       0.8934737 |           0.9049136 |
| BMD         | NonBrain  | adipocytes | SC_HPA81     | NonBrain    | NA             | Adipocytes           | DUP      | -0.0135559 | 0.0110931 | 0.2217055 |   0.4706469 |         -0.0022584 |  0.0010175 |       0.0264508 |           0.0349589 |
| BMD         | NonBrain  | adipocytes | SC_HPA81     | NonBrain    | NA             | Adipocytes           | DEL      |  0.0224277 | 0.0149456 | 0.1334548 |   0.3468137 |         -0.0023796 |  0.0018936 |       0.2088704 |           0.2415268 |

##### Panel D Figure

``` r
# ------------------------------------------------------------------------------
# 1. DEFINE TARGET TRAITS
# ------------------------------------------------------------------------------
target_traits <- c("BMI", "StandingHeight", "HeartRate", "Platelet")

# ------------------------------------------------------------------------------
# 2. MERGE COHORTS & COMPUTE PLOTTING METRICS
# ------------------------------------------------------------------------------
df_panel_d_data <- inner_join(
  df_funburd_ukbb, 
  df_funburd_aou, 
  by = c("Trait", "Geneset", "Geneset_Type", "CNV_Type")
) %>%
  # Rename UKBB effect size so it visually balances with Effectsize_AoU
  rename(Effectsize_UKBB = Effectsize) %>%
  
  # Calculate -log10 p-values
  mutate(
    UKBBLogp = -log10(p_value),
    AoULogp  = -log10(p_value_AoU),
    
    # Enforce strict plotting order for variants
    CNV_Type = factor(CNV_Type, levels = c("DEL", "DUP", "SNP"))
  ) %>%
  
  # Filter for our specific traits
  filter(Trait %in% target_traits) %>%
  mutate(Trait = factor(Trait, levels = target_traits))

# Filter specifically for the plot
df_bmi_concordance <- df_panel_d_data %>% 
  filter(Trait == "BMI", Geneset_Type == "Tissue_Fantom60")

# ------------------------------------------------------------------------------
# 3. PLOTTING FUNCTION
# ------------------------------------------------------------------------------
# Create a reusable function to guarantee identical formatting for Del and Dup
plot_concordance <- function(data, variant_type, title_text) {
  
  # Filter for the specific variant and remove any NAs to prevent cor.test crashes
  sub_df <- data %>% 
    filter(CNV_Type == variant_type) %>%
    filter(!is.na(Effectsize_UKBB) & !is.na(Effectsize_AoU))
  
  # Calculate Pearson r and exactly compute its p-value
  cor_res <- cor.test(sub_df$Effectsize_UKBB, sub_df$Effectsize_AoU)
  r_val <- cor_res$estimate
  p_val <- cor_res$p.value
  
  # Format the correlation label: Append '*' ONLY if p-value < 0.05
  r_label <- paste0("r = ", round(r_val, 2), if(p_val < 0.05) "*" else "")
  
  # Generate Plot
  ggplot(sub_df, aes(x = Effectsize_UKBB, y = Effectsize_AoU)) +
    
    # Background alignment lines
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.8) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60", linewidth = 1) +
    
    # Scatter points
    geom_point(color = "grey70", size = 2.5, alpha = 0.6) +
    
    # Linear Regression Line (Orange)
    geom_smooth(method = "lm", formula = y ~ x,color = "#E69F00", fill = "#E69F00", alpha = 0.2, linewidth = 1.5) +
    
    # Add Correlation text
    annotate("text", x = -Inf, y = Inf, label = r_label, 
             hjust = -0.1, vjust = 1.5, size = 6, fontface = "bold") +
    
    labs(title = title_text, x = "UKBB", y = "AoU") +
    
    # Prism-Style minimal theme
    theme_bw(base_size = 14) +
    theme(
      aspect.ratio = 1,
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(face = "bold", size = 14)
    )
}

# ------------------------------------------------------------------------------
# 4. GENERATE & ASSEMBLE PANELS
# ------------------------------------------------------------------------------
p_del <- plot_concordance(df_bmi_concordance, "DEL", "BMI\nDel EffectSize")
p_dup <- plot_concordance(df_bmi_concordance, "DUP", "BMI\nDup EffectSize")

# Stitch them together using patchwork
p_panel_d_final <- p_del + p_dup

# Render the final figure
print(p_panel_d_final)
```

![](Fig2_June2026_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

##### Panel D Stats

``` r
# ------------------------------------------------------------------------------
# 1. COMPUTE REGRESSION & CORRELATION STATS
# ------------------------------------------------------------------------------
library(dplyr)
library(knitr)

formatted_stats_d <- df_panel_d_data %>%
  filter(Trait == "BMI", Geneset_Type == "Tissue_Fantom60") %>%
  # 1. Clean the data: Keep only rows where both cohorts have an effect size
  filter(!is.na(Effectsize_UKBB) & !is.na(Effectsize_AoU)) %>%
  
  # 2. Group by our primary categories
  group_by(Geneset_Type, Trait, CNV_Type) %>%
  
  # 3. SAFETY NET: cor.test and lm require at least 3 points to calculate degrees of freedom
  filter(n() >= 3) %>%
  
  # 4. Calculate metrics on the fly
  summarize(
    
    # Correlation Metrics
    `Pearson r` = round(cor.test(Effectsize_UKBB, Effectsize_AoU)$estimate, 3),
    
    # Extract EXACT p-value (using scientific notation for readability on tiny values)
    `Corr P-Value` = format(cor.test(Effectsize_UKBB, Effectsize_AoU)$p.value, scientific = TRUE, digits = 4),
    
    # Linear Regression Metrics (AoU on Y-axis, UKBB on X-axis)
    `Slope (β)` = round(coef(lm(Effectsize_AoU ~ Effectsize_UKBB))[2], 3),
    `Intercept` = round(coef(lm(Effectsize_AoU ~ Effectsize_UKBB))[1], 4),
    `R-Squared` = round(summary(lm(Effectsize_AoU ~ Effectsize_UKBB))$r.squared, 3),
    
    # Add an N column so reviewers know how many points went into the model
    `N Variants` = n(),
    
    .groups = "drop"
  ) %>%
  
  # Sort for clean presentation
  arrange(Geneset_Type, Trait, CNV_Type)

# ------------------------------------------------------------------------------
# 2. RENDER THE TABLE
# ------------------------------------------------------------------------------
kable(
  formatted_stats_d, 
  align = "c",
  caption = "Linear Regression & Correlation Stats: UKBB vs AoU Effect Sizes"
)
```

|  Geneset_Type   | Trait | CNV_Type | Pearson r | Corr P-Value | Slope (β) | Intercept | R-Squared | N Variants |
|:---------------:|:-----:|:--------:|:---------:|:------------:|:---------:|:---------:|:---------:|:----------:|
| Tissue_Fantom60 |  BMI  |   DEL    |   0.915   |  1.601e-24   |   1.203   |  -0.0198  |   0.837   |     60     |
| Tissue_Fantom60 |  BMI  |   DUP    |   0.422   |  7.733e-04   |   0.469   |  0.0025   |   0.178   |     60     |

Linear Regression & Correlation Stats: UKBB vs AoU Effect Sizes

#### Panel E: Power vs Reproducibility (UKBB & AoU)

##### Description: Evaluates the relationship between discovery power (number of significant genesets in UKBB) and cross-cohort reproducibility (Pearson correlation of effect sizes between UKBB and AoU)

##### Panel E Figure

``` r
# ------------------------------------------------------------------------------
# 1. CALCULATE CONCORDANCE (Y-Axis)
# ------------------------------------------------------------------------------
# We use the merged dataset from Panel D to calculate the Pearson r for each group
df_concordance <- df_panel_d_data %>%
  # Target Fantom genesets (Using grepl to catch "Fantom" or "Tissue_Fantom60")
  filter(grepl("Fantom", Geneset_Type, ignore.case = TRUE)) %>%
  # Keep only DEL and DUP
  filter(CNV_Type %in% c("DEL", "DUP")) %>%
  group_by(Trait, CNV_Type) %>%
  summarise(
    Concordance = cor(Effectsize_UKBB, Effectsize_AoU, use = "complete.obs"), 
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 2. CALCULATE SIGNIFICANT COUNTS (X-Axis)
# ------------------------------------------------------------------------------
# We pull this directly from the raw UKBB dataset to ensure we count all 
# significant findings, even if they didn't overlap with AoU
df_sig_counts <- df_funburd_ukbb %>%
  filter(grepl("Fantom", Geneset_Type, ignore.case = TRUE)) %>%
  filter(CNV_Type %in% c("DEL", "DUP")) %>%
  # Filter for FDR significance threshold
  filter(FDR_p_value < 0.05) %>%
  group_by(Trait, CNV_Type) %>%
  summarise(N_Sig = n(), .groups = "drop")

# ------------------------------------------------------------------------------
# 3. MERGE & COMPUTE GLOBAL STATS
# ------------------------------------------------------------------------------
df_plot <- df_concordance %>%
  left_join(df_sig_counts, by = c("Trait", "CNV_Type")) %>%
  mutate(
    # If a group had 0 significant genesets, the join puts NA. We replace it with 0.
    N_Sig = ifelse(is.na(N_Sig), 0, N_Sig),
    # Create the clean "Del-BMI" labels
    Label = paste0(str_to_title(CNV_Type), "-", Trait)
  )

# Calculate the global correlation to display on the plot (Concordance vs N_Sig)
cor_test_res <- cor.test(df_plot$Concordance, df_plot$N_Sig)
global_r <- round(cor_test_res$estimate, 2)
global_p <- cor_test_res$p.value

# Format the label dynamically: add an asterisk if p-value < 0.05
global_label <- paste0("r = ", global_r, ifelse(global_p < 0.05, "*", ""))

# Print stats to console for your reference
#message("Global Correlation: r = ", global_r, " | p-value = ", global_p)

# ------------------------------------------------------------------------------
# 4. GENERATE THE HARMONIZED PLOT
# ------------------------------------------------------------------------------
p_summary <- ggplot(df_plot, aes(x = N_Sig, y = Concordance)) +
  
  # Solid grey regression line
  geom_smooth(method = "lm",formula = y ~ x, color = "darkgrey", linetype = "solid", se = FALSE, linewidth = 1.2) +
  
  # Scatter points colored by CNV Type
  geom_point(aes(color = CNV_Type), size = 5) +
  
  # Repelling Text Labels to prevent overlap
  geom_text_repel(aes(label = Label), 
                  size = 4.5, 
                  fontface = "bold", 
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  color = "black") +
  
  # Add the global correlation text to the top right corner
  annotate("text", x = Inf, y = Inf, label = global_label, 
           hjust = 1.1, vjust = 1.5, size = 6, fontface = "bold") +
  
  # Set specific point colors matching your screenshot
  scale_color_manual(values = c("DEL" = "red2", "DUP" = "royalblue")) +
  
  # Labels and Scaling
  labs(
    x = "N Sig. Genesets (UKBB)",
    y = "Effect-size Concordance (UKBB vs AoU)"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  coord_cartesian(ylim = c(-0.2, 1.05), clip = "off") +
  
  # Custom Prism-style Theme
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.8),
    axis.ticks = element_line(colour = "black", linewidth = 0.8),
    axis.text = element_text(size = 12, face = "bold", color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(20, 20, 20, 20)
  )

# Render the plot
print(p_summary)
```

![](Fig2_June2026_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

##### Panel E Stats

``` r
# ------------------------------------------------------------------------------
# 1. COMPUTE REGRESSION & CORRELATION STATS
# ------------------------------------------------------------------------------
# Using df_plot which was generated in the Panel E Figure chunk
formatted_stats_e <- df_plot %>%
  
  # 1. Clean the data: Keep only rows where both metrics exist
  filter(!is.na(N_Sig) & !is.na(Concordance)) %>%
  
  # 2. Add a global grouping label since we are calculating the overall trend
  mutate(Scope = "4 Mutual Traits") %>%
  group_by(Scope) %>%
  
  # 3. SAFETY NET: cor.test and lm require at least 3 points 
  filter(n() >= 3) %>%
  
  # 4. Calculate metrics on the fly
  summarize(
    
    # Correlation Metrics
    `Pearson r` = round(cor.test(N_Sig, Concordance)$estimate, 3),
    
    # Extract EXACT p-value (using scientific notation for readability)
    `Corr P-Value` = format(cor.test(N_Sig, Concordance)$p.value, scientific = TRUE, digits = 4),
    
    # Linear Regression Metrics (Y ~ X)
    `Slope (β)` = round(coef(lm(Concordance ~ N_Sig))[2], 4),
    `Intercept` = round(coef(lm(Concordance ~ N_Sig))[1], 4),
    `R-Squared` = round(summary(lm(Concordance ~ N_Sig))$r.squared, 3),
    
    # Add an N column so reviewers know how many data points formed the line
    `N Data Points` = n(),
    
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 2. RENDER THE TABLE
# ------------------------------------------------------------------------------
kable(
  formatted_stats_e, 
  align = "c",
  caption = "Linear Regression & Correlation Stats: Discovery Power vs. Reproducibility"
)
```

|      Scope      | Pearson r | Corr P-Value | Slope (β) | Intercept | R-Squared | N Data Points |
|:---------------:|:---------:|:------------:|:---------:|:---------:|:---------:|:-------------:|
| 4 Mutual Traits |   0.868   |  5.239e-03   |  0.0252   |  0.0394   |   0.753   |       8       |

Linear Regression & Correlation Stats: Discovery Power
vs. Reproducibility
