Fig2
================

## Fig. 2: Heatmap of effect sizes for whole body tissue, cell type gene sets across traits.

## Markdown for generating panel figures and statistics

## ——— Figure legend ———

Legend: Heatmap displays a representative set of the most significant
associations between 5 categories of traits and deletions and
duplications aggregated across (A) tissues and (B) whole body cell
types. Traits are categorized and shown along the x-axis, while the
y-axis lists gene sets, tissues, and cell types. The blue and red
intensity color scale reflects the negative and positive effect sizes.
Black asterisks indicate statistically significant associations between
traits and genes (FDR correction across 172 ✕ 43 ✕ 2 =14,792 tests).
C-F) Bar plots (with standard error) summarizing the differences in the
level of association/enrichments between the type of variants for C, E)
all gene-sets, and D, F) brain and non-brain functional gene sets. Black
asterisks (\*) indicate statistically significant proportion differences
(FDR-adjusted). Abbreviations: HDL: high-density lipoprotein; HbA1c:
glycated haemoglobin; BMD: bone mineral density; BMI: body mass index;
EA: educational attainment; FI: fluid intelligence; TDI: townsend
deprivation index; VSWM: visuospatial working memory; FVC: Forced vital
capacity; Del: deletion; Dup: duplication; SNP: single nucleotide
polymorphism.

``` r
#### Libraries for Heatmap & Point Range Plots

library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(ggplot2)
library(ggprism)
library(grid)
library(circlize)
```

    ## ========================================
    ## circlize version 0.4.16
    ## CRAN page: https://cran.r-project.org/package=circlize
    ## Github page: https://github.com/jokergoo/circlize
    ## Documentation: https://jokergoo.github.io/circlize_book/book/
    ## 
    ## If you use it in published research, please cite:
    ## Gu, Z. circlize implements and enhances circular visualization
    ##   in R. Bioinformatics 2014.
    ## 
    ## This message can be suppressed by:
    ##   suppressPackageStartupMessages(library(circlize))
    ## ========================================

``` r
library(ggpubr)
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

``` r
library(ComplexHeatmap)
```

    ## ========================================
    ## ComplexHeatmap version 2.22.0
    ## Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
    ## Github page: https://github.com/jokergoo/ComplexHeatmap
    ## Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
    ## 
    ## If you use it in published research, please cite either one:
    ## - Gu, Z. Complex Heatmap Visualization. iMeta 2022.
    ## - Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
    ##     genomic data. Bioinformatics 2016.
    ## 
    ## 
    ## The new InteractiveComplexHeatmap package can directly export static 
    ## complex heatmaps into an interactive Shiny app with zero effort. Have a try!
    ## 
    ## This message can be suppressed by:
    ##   suppressPackageStartupMessages(library(ComplexHeatmap))
    ## ========================================

#### Display 43\*(60+81+31) estimated effect sizes (B1) and FDR-adjusted p-values - For Deletion & Duplication Seperately

``` r
Path="/Users/sayekazem/Desktop/New_Project_2023_Summer/PaperMultitrait2024_Data_FinalVersion/EffectSizes_Euro/"
# Load necessary data and set row names
DUP_Fantom_43 <- read.csv(paste0(Path, 'DUP_Est_Fantom_43Traits_AncestryEuro_Jan2025.csv'))
rownames(DUP_Fantom_43) <- DUP_Fantom_43$X

DEL_Fantom_43 <- read.csv(paste0(Path, 'DEL_Est_Fantom_43Traits_AncestryEuro_Jan2025.csv'))
rownames(DEL_Fantom_43) <- DEL_Fantom_43$X

DUP_HPACell_43 <- read.csv(paste0(Path, 'DUP_Est_HPA_43Traits_AncestryEuro_Jan2025.csv'))
rownames(DUP_HPACell_43) <- DUP_HPACell_43$X

DEL_HPACell_43 <- read.csv(paste0(Path, 'DEL_Est_HPA_43Traits_AncestryEuro_Jan2025.csv'))
rownames(DEL_HPACell_43) <- DEL_HPACell_43$X

DUP_SilettiSuper_43 <- read.csv(paste0(Path, 'DUP_Est_SilettiSuper_43Traits_AncestryEuro_Jan2025.csv'))
rownames(DUP_SilettiSuper_43) <- DUP_SilettiSuper_43$X

DEL_SilettiSuper_43 <- read.csv(paste0(Path, 'DEL_Est_SilettiSuper_43Traits_AncestryEuro_Jan2025.csv'))
rownames(DEL_SilettiSuper_43) <- DEL_SilettiSuper_43$X

# Combine data for DUP and DEL
ALL_Ests_Final_DUP <- rbind(
  DUP_Fantom_43[2:length(DUP_Fantom_43)],
  DUP_HPACell_43[2:length(DUP_HPACell_43)],
  DUP_SilettiSuper_43[2:length(DUP_SilettiSuper_43)]
)

ALL_Ests_Final_DEL <- rbind(
  DEL_Fantom_43[2:length(DEL_Fantom_43)],
  DEL_HPACell_43[2:length(DEL_HPACell_43)],
  DEL_SilettiSuper_43[2:length(DEL_SilettiSuper_43)]
)

# Assign combined data to final variables
DUP_43 <- ALL_Ests_Final_DUP
DEL_43 <- ALL_Ests_Final_DEL

# Define column names
Column_name <- c(
  'Albumin', 'Birth Weight', 'BMD', 'BMI', 'Calcium', 'Cholesterol', 'FEV1', 'FVC', 'Gamma', 'Glucose', 'HbA1c', 'HDL', 'IGF1', 'LDL', 'Neuroticism', 'Platelet', 'Reaction', 'Redblood', 'Standing Height', 'TDI', 'Triglycerides', 'EA', 'FI', 'HeartRate', 'MET', 'Number of Sexual Partners', 'Symbol Digit Substitution', 'Trail Making 1', 'Trail Making 2', 'VSWM1', 'VSWM2', 'Menarche', 'Menopause', 'Depression (self-report)', 'GP Visit for Anxiety or Depression', 'Mood Swings', 'Risk Taking', 'Loneliness', 'Irritability', 'Miserableness', 'Hypertention', 'Guilty Feeling','MoodAnxiety'
)


# Set column names for final data
colnames(DEL_43) <- Column_name
colnames(DUP_43) <- Column_name

# Load p-value data and set row names
DUP_Fantom_43_p <- read.csv(paste0(Path, 'DUP_pval_Fantom_43Traits_AncestryEuro_Jan2025.csv'))
rownames(DUP_Fantom_43_p) <- DUP_Fantom_43_p$X

DEL_Fantom_43_p <- read.csv(paste0(Path, 'DEL_pval_Fantom_43Traits_AncestryEuro_Jan2025.csv'))
rownames(DEL_Fantom_43_p) <- DEL_Fantom_43_p$X

DUP_HPACell_43_p <- read.csv(paste0(Path, 'DUP_pval_HPA_43Traits_AncestryEuro_Jan2025.csv'))
rownames(DUP_HPACell_43_p) <- DUP_HPACell_43_p$X

DEL_HPACell_43_p <- read.csv(paste0(Path, 'DEL_pval_HPA_43Traits_AncestryEuro_Jan2025.csv'))
rownames(DEL_HPACell_43_p) <- DEL_HPACell_43_p$X

DUP_SilettiSuper_43_p <- read.csv(paste0(Path, 'DUP_pval_SilettiSuper_43Traits_AncestryEuro_Jan2025.csv'))
rownames(DUP_SilettiSuper_43_p) <- DUP_SilettiSuper_43_p$X

DEL_SilettiSuper_43_p <- read.csv(paste0(Path, 'DEL_pval_SilettiSuper_43Traits_AncestryEuro_Jan2025.csv'))
rownames(DEL_SilettiSuper_43_p) <- DEL_SilettiSuper_43_p$X

# Combine p-value data for DUP and DEL
ALL_pval_Final_DUP <- rbind(
  DUP_Fantom_43_p[2:length(DUP_Fantom_43_p)],
  DUP_HPACell_43_p[2:length(DUP_HPACell_43_p)],
  DUP_SilettiSuper_43_p[2:length(DUP_SilettiSuper_43_p)]
)

ALL_pval_Final_DEL <- rbind(
  DEL_Fantom_43_p[2:length(DEL_Fantom_43_p)],
  DEL_HPACell_43_p[2:length(DEL_HPACell_43_p)],
  DEL_SilettiSuper_43_p[2:length(DEL_SilettiSuper_43_p)]
)

# FDR Correction
ALL_pval_Final_DEL_vector <- as.vector(as.matrix(ALL_pval_Final_DEL))
ALL_pval_Final_DUP_vector <- as.vector(as.matrix(ALL_pval_Final_DUP))
combined_pvalues <- c(ALL_pval_Final_DEL_vector, ALL_pval_Final_DUP_vector)
fdr_corrected <- p.adjust(combined_pvalues, method = "fdr")

# Split corrected p-values
n_del <- length(ALL_pval_Final_DEL_vector)
del_fdr_corrected <- fdr_corrected[1:n_del]
dup_fdr_corrected <- fdr_corrected[(n_del + 1):length(fdr_corrected)]

# Reshape corrected p-values
dim(del_fdr_corrected) <- dim(as.matrix(ALL_pval_Final_DEL))
dim(dup_fdr_corrected) <- dim(as.matrix(ALL_pval_Final_DUP))
p_values_fdr_DEL <- as.data.frame(matrix(del_fdr_corrected, nrow = nrow(ALL_pval_Final_DEL), dimnames = dimnames(ALL_pval_Final_DEL)))
p_values_fdr_DUP <- as.data.frame(matrix(dup_fdr_corrected, nrow = nrow(ALL_pval_Final_DUP), dimnames = dimnames(ALL_pval_Final_DUP)))

# Set column names for corrected p-values
colnames(p_values_fdr_DEL) <- Column_name
colnames(p_values_fdr_DUP) <- Column_name


# Subset data for different categories
DEL_Fantom_43 <- DEL_43[1:60,]
DUP_Fantom_43 <- DUP_43[1:60,]

DEL_HPA_43 <- DEL_43[61:141,]
DUP_HPA_43 <- DUP_43[61:141,]

DEL_Siletti_43 <- DEL_43[142:172,]
DUP_Siletti_43 <- DUP_43[142:172,]

# Subset p-values for different categories
p_values_fdr_DEL_Fantom <- p_values_fdr_DEL[1:60,]
p_values_fdr_DUP_Fantom <- p_values_fdr_DUP[1:60,]

p_values_fdr_DEL_HPA <- p_values_fdr_DEL[61:141,]
p_values_fdr_DUP_HPA <- p_values_fdr_DUP[61:141,]

p_values_fdr_DEL_Siletti <- p_values_fdr_DEL[142:172,]
p_values_fdr_DUP_Siletti <- p_values_fdr_DUP[142:172,]
# Display a subset of the Effectsizes & FDR p-values

DEL_43[38:43,36:40]
```

    ##                   Mood Swings Risk Taking  Loneliness  Irritability
    ## salivary_gland    0.005624954  0.07718937 0.038568844  2.122400e-02
    ## occipital_cortex  0.109301140 -0.14114981 0.161640534 -3.237108e-03
    ## parietal_lobe     0.091364203 -0.11799118 0.107559454  3.544292e-05
    ## thalamus          0.084343451 -0.09349208 0.114267783  3.128459e-02
    ## urinary_bladder   0.026896205 -0.04163033 0.039638043  2.131925e-02
    ## prostate         -0.004818662 -0.04008795 0.008230059  2.327601e-02
    ##                  Miserableness
    ## salivary_gland     0.030000639
    ## occipital_cortex   0.100153042
    ## parietal_lobe      0.083059202
    ## thalamus           0.134323120
    ## urinary_bladder   -0.026481138
    ## prostate          -0.002541693

``` r
p_values_fdr_DEL[38:43,36:40]
```

    ##                  Mood Swings  Risk Taking   Loneliness Irritability
    ## salivary_gland   0.934473843 0.1135202748 5.320306e-01    0.7487988
    ## occipital_cortex 0.002474387 0.0008943507 6.125466e-05    0.9659713
    ## parietal_lobe    0.006602165 0.0029414220 4.568794e-03    0.9978481
    ## thalamus         0.002901150 0.0065470487 5.849230e-04    0.4606269
    ## urinary_bladder  0.617788774 0.4853068307 5.077773e-01    0.7453382
    ## prostate         0.925252273 0.3486757132 8.918103e-01    0.6135162
    ##                  Miserableness
    ## salivary_gland    5.776831e-01
    ## occipital_cortex  6.662418e-03
    ## parietal_lobe     1.502521e-02
    ## thalamus          1.853520e-07
    ## urinary_bladder   6.263662e-01
    ## prostate          9.619162e-01

#### Prepating Subset of Data for the Heatmap plot

##### Subsetting Data At Tissue Level - Brain vs. Non-Brain Traits and Brain vs. Non-Brain Tissue Gene Sets - Deletion

``` r
# Fig 1 A : Tissue level (Fantom)
# Convert data to matrices
DEL_Fantom_43_mat <- as.matrix(DEL_Fantom_43)
rownames(DEL_Fantom_43_mat) <- rownames(DEL_Fantom_43)
colnames(DEL_Fantom_43_mat) <- colnames(DEL_Fantom_43)

DUP_Fantom_43_mat <- as.matrix(DUP_Fantom_43)
rownames(DUP_Fantom_43_mat) <- rownames(DUP_Fantom_43)
colnames(DUP_Fantom_43_mat) <- colnames(DUP_Fantom_43)

# Order matrix by tissues
tissue_order <- c("olfactory_bulb", "frontal_lobe", "medial_frontal_gyrus", "insular_cortex",
                  "paracentral_gyrus", "parietal_lobe", "temporal_cortex", "medial_temporal_gyrus", "postcentral_gyrus",
                  "hippocampus", "occipital_lobe", "occipital_cortex", "occipital_pole",
                  "amygdala", "caudate", "putamen", "substantia_nigra", "globus_pallidus", "thalamus",
                  "nucleus_accumbens", "pituitary_gland", "locus_coeruleus", "corpus_callosum",
                  "cerebellum", "pons", "medulla_oblongata", "spinal_cord", "retina",
                  "salivary_gland", "tongue", "esophagus", "thyroid_gland", "thymus",
                  "tonsil", "heart_muscle", "skeletal_muscle", "smooth_muscle", "lung", "liver", "gallbladder",
                  "spleen", "lymph_node", "pancreas", "adipose_tissue", "colon", "small_intestine",
                  "appendix", "kidney", "urinary_bladder",
                  "breast", "endometrium", "cervix", "vagina", "ovary", "placenta", "prostate", "seminal_vesicle",
                  "ductus_deferens", "epididymis", "testis")

LabRow_Names <- c("Olfactory Bulb", "Frontal Lobe", "Medial Frontal Gyrus", "Insular Cortex",
                  "Paracentral Gyrus", "Parietal Lobe", "Temporal Cortex", "Medial Temporal Gyrus", "Postcentral Gyrus",
                  "Hippocampus", "Occipital Lobe", "Occipital Cortex", "Occipital Pole",
                  "Amygdala", "Caudate", "Putamen", "Substantia Nigra", "Globus Pallidus", "Thalamus",
                  "Nucleus Accumbens", "Pituitary Gland", "Locus Coeruleus", "Corpus Callosum",
                  "Cerebellum", "Pons", "Medulla Oblongata", "Spinal Cord", "Retina",
                  "Salivary Gland", "Tongue", "Esophagus", "Thyroid Gland", "Thymus",
                  "Tonsil", "Heart Muscle", "Skeletal Muscle", "Smooth Muscle", "Lung", "Liver", "Gallbladder",
                  "Spleen", "Lymph Node", "Pancreas", "Adipose Tissue", "Colon", "Small Intestine",
                  "Appendix", "Kidney", "Urinary Bladder",
                  "Breast", "Endometrium", "Cervix", "Vagina", "Ovary", "Placenta", "Prostate", "Seminal Vesicle",
                  "Ductus Deferens", "Epididymis", "Testis")

trait_name_all <- c('Albumin', 'IGF1', 'Cholesterol', 'HDL', 'LDL',
                'Calcium', 'Platelet', 'Redblood', 'Gamma', 'HbA1c', 'Triglycerides', 'Glucose',
                'Menarche', 'Menopause', 'Number of Sexual Partners', 'MET',
                'FVC', 'FEV1', 'Standing Height', 'Birth Weight', 'BMI', 'BMD', 'Hypertention', 'HeartRate',
                'EA', 'FI', 'Symbol Digit Substitution', 'Trail Making 1', 'Trail Making 2', 'Reaction', 'VSWM1', 'VSWM2', 'TDI',
                'Risk Taking', 'Depression (self-report)', 'GP Visit for Anxiety or Depression', 'Neuroticism', 'Irritability', 'Mood Swings', 'Miserableness', 'Guilty Feeling', 'Loneliness','MoodAnxiety')

trait_order <- trait_name_all

# Order rows and columns for DEL matrix
order_index <- match(rownames(DEL_Fantom_43_mat), tissue_order)
DEL_Fantom_43_mat_ordered_all <- DEL_Fantom_43_mat[order(order_index), ]
DEL_Fantom_43_mat_ordered_all <- DEL_Fantom_43_mat_ordered_all[, match(trait_order, colnames(DEL_Fantom_43_mat_ordered_all))]

# Order rows and columns for p-values
order_index <- match(rownames(p_values_fdr_DEL_Fantom), tissue_order)
p_values_fdr_DEL_Fantom_star_ordered_all <- p_values_fdr_DEL_Fantom[order(order_index), ]
p_values_fdr_DEL_Fantom_star_ordered_all <- p_values_fdr_DEL_Fantom_star_ordered_all[, match(trait_order, colnames(p_values_fdr_DEL_Fantom_star_ordered_all))]
p_values_fdr_DEL_Fantom_star_ordered_all <- ifelse(p_values_fdr_DEL_Fantom_star_ordered_all < 0.05, "*", "")

# Extract subset of datasets for heatmap plot
Tissues_Subset <- c("placenta", "adipose_tissue", "spleen", "liver", "skeletal_muscle", "heart_muscle", "thalamus", "occipital_cortex", "frontal_lobe", "temporal_cortex")
Traits_Subset <- c("Cholesterol", "HDL", "Platelet", "Redblood", "HbA1c",
                   "Menarche", "Menopause", "FVC", "BMI", "HeartRate", "BMD", "EA", "FI", "TDI", "VSWM1",
                   "Risk Taking", "Depression (self-report)", "Neuroticism", "Loneliness",'MoodAnxiety')

# Subset DEL_Fantom_43 matrix
DEL_Fantom_43_mat_ordered <- subset(DEL_Fantom_43_mat_ordered_all, select = Traits_Subset)
DEL_Fantom_43_mat_ordered <- DEL_Fantom_43_mat_ordered[rownames(DEL_Fantom_43_mat_ordered) %in% Tissues_Subset, ]

# Subset p-values
p_values_fdr_DEL_Fantom_star_ordered <- subset(p_values_fdr_DEL_Fantom_star_ordered_all, select = Traits_Subset)
p_values_fdr_DEL_Fantom_star_ordered <- p_values_fdr_DEL_Fantom_star_ordered[rownames(p_values_fdr_DEL_Fantom_star_ordered) %in% Tissues_Subset, ]

LabRow_Names <- c("Frontal Lobe", "Temporal Cortex", "Occipital Cortex", "Thalamus",
                  "Heart Muscle", "Skeletal Muscle", "Liver",
                  "Spleen", "Adipose Tissue", "Placenta")

trait_name <- c("Cholesterol", "HDL", "Platelet", "Redblood", "HbA1c",
                "Menarche", "Menopause", "FVC", "BMI", "HeartRate", "BMD", "EA", "FI", "TDI", "VSWM1",
                "Risk Taking", "Depression", "Neuroticism", "Loneliness",'MoodAnxiety')
```

#### Heatmap plot : Tissues & Traits

``` r
# Heatmap settings
colors <- colorRampPalette(c("#386cb0", "white", "#a50f15"))(100)
min_val <- -0.2
num_colors <- length(colors)
breaks <- seq(min_val, -min_val, length.out = num_colors)
color_function <- colorRamp2(breaks, colors)

# Save the heatmap as a PNG file
# png("/Heatmap_DEL_Fantom_FullVersion_FDR_Categorized_NOV2024.png", width = 14000, height = 15000, res = 800)

# Create row annotations
row_annotation <- data.frame(
  Tissue = c(rep("Brain", 4), rep("Non-Brain", 6))
)

row_annotation$Tissue <- factor(row_annotation$Tissue, levels = c("Brain", "Non-Brain"))

# Create column annotations
column_annotation <- data.frame(
  Trait = c(rep("Blood Assays",5),rep("Reproductive\n Factors",2),rep("Physical\n Measure",4),rep("Cognitive\n Metrics",4),rep("Mental\n Health",5))
)

column_annotation$Trait <- factor(column_annotation$Trait, levels = c("Blood Assays", "Reproductive\n Factors", "Physical\n Measure", "Cognitive\n Metrics", "Mental\n Health"))

# Define column heatmap annotation
column_ha <- HeatmapAnnotation(
  Trait = column_annotation$Trait,
  col = list(Trait = c("Blood Assays" = "#e78ac3", "Reproductive\n Factors" = "#8da0cb", "Physical\n Measure" = "#a6d854", "Cognitive\n Metrics" = "#66c2a5", "Mental\n Health" = "#fc8d62")),
  annotation_name = NULL,  # Remove both 'Supercluster' and 'Trait' labels
  show_legend = FALSE,
  show_annotation_name = FALSE
)

# Set row and column names for the ordered matrix
rownames(DEL_Fantom_43_mat_ordered) <- LabRow_Names
colnames(DEL_Fantom_43_mat_ordered) <- trait_name

# Create the heatmap without dendrograms
ht1 <- Heatmap(
  DEL_Fantom_43_mat_ordered, name = "Effect size",
  cluster_rows = FALSE,  # Disable row clustering
  cluster_columns = FALSE,  # Disable column clustering
  show_row_dend = FALSE,  # Do not show row dendrogram
  show_column_dend = FALSE,  # Do not show column dendrogram
  top_annotation = column_ha,
  row_split = row_annotation$Tissue,  # Split rows into 2 groups
  column_split = column_annotation$Trait,  # Split columns into 5 groups
  row_title = NULL,
  row_title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "sans"),
  row_names_gp = gpar(col = "black", fontsize = 10, fontfamily = "sans", fontface = "bold"),
  column_title_gp = gpar(col = c('#e78ac3', '#8da0cb', '#a6d854', '#66c2a5', '#fc8d62'), fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(col = "black", fontsize = 10, rot = 45, fontfamily = "sans", fontface = "bold"),
  column_names_rot = 45,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(p_values_fdr_DEL_Fantom_star_ordered[i, j], x = x, y = y - convertHeight(grobHeight(textGrob('*')), "mm"), gp = gpar(fontsize = 32, col = 'grey40', fontfamily = "sans", just = "center", fontface = "bold", hjust = 30))
  },
  col = color_function, row_gap = unit(0.5, "cm"),  # Adjust row margin (increase or decrease as needed)
  column_gap = unit(0.5, "cm"),
  show_heatmap_legend = FALSE,
  rect_gp = gpar(col = "lightgrey", lwd = 0.5)
)

# Draw the heatmap with padding
draw(ht1, padding = unit(c(0.5, 2, 0.5, 0.5), "cm"))
```

![](Fig2_files/figure-gfm/fig2_3-1.png)<!-- -->

#### Prepating Subset of Data for the Heatmap plot

##### Subsetting Data At Whole Body Cell Type Level - Brain vs. Non-Brain Traits and Brain vs. Non-Brain Cell Type Gene Sets - Deletion

``` r
# Fig 1 B : Whole Body Cell Type level (HPA)

# Convert DEL and DUP HPA data to matrices and set row and column names
DEL_HPA_43_mat <- as.matrix(DEL_HPA_43)
rownames(DEL_HPA_43_mat) <- rownames(DEL_HPA_43)
colnames(DEL_HPA_43_mat) <- colnames(DEL_HPA_43)

DUP_HPA_43_mat <- as.matrix(DUP_HPA_43)
rownames(DUP_HPA_43_mat) <- rownames(DUP_HPA_43)
colnames(DUP_HPA_43_mat) <- colnames(DUP_HPA_43)

# Define the order of cells for the heatmap
cell_order <- c("inhibitory_neurons", "excitatory_neurons", "astrocytes", "oligodendrocytes", "oligodendrocyte_precursor_cells", "microglial_cells",
                "muller_glia_cells", "melanocytes", "smooth_muscle_cells", "fibroblasts", "endometrial_stromal_cells", "mesothelial_cells",
                "endothelial_cells", "adipocytes", "lymphatic_endothelial_cells", "ovarian_stromal_cells", "granulosa_cells",
                "peritubular_cells", "leydig_cells", "sertoli_cells", "skeletal_myocytes", "cardiomyocytes", "t_cells", "nk_cells",
                "dendritic_cells", "plasma_cells", "b_cells", "erythroid_cells", "schwann_cells", "granulocytes", "monocytes",
                "kupffer_cells", "macrophages", "hofbauer_cells", "langerhans_cells", "proximal_tubular_cells", "collecting_duct_cells",
                "hepatocytes", "cholangiocytes", "exocrine_glandular_cells", "ductal_cells", "pancreatic_endocrine_cells", "distal_tubular_cells",
                "undifferentiated_cells", "intestinal_goblet_cells", "enteroendocrine_cells", "paneth_cells", "distal_enterocytes",
                "proximal_enterocytes", "gastric_mucus_secreting_cells", "prostatic_glandular_cells", "basal_prostatic_cells",
                "breast_myoepithelial_cells", "breast_glandular_cells", "salivary_duct_cells", "mucus_glandular_cells",
                "serous_glandular_cells", "alveolar_cells_type_2", "alveolar_cells_type_1", "secretory_cells", "glandular_and_luminal_cells",
                "ciliated_cells", "ionocytes", "basal_respiratory_cells", "club_cells", "squamous_epithelial_cells", "basal_squamous_epithelial_cells",
                "suprabasal_keratinocytes", "basal_keratinocytes", "syncytiotrophoblasts", "cytotrophoblasts", "extravillous_trophoblasts",
                "rod_photoreceptor_cells", "cone_photoreceptor_cells", "horizontal_cells", "bipolar_cells", "spermatogonia", "spermatocytes",
                "oocytes", "late_spermatids", "early_spermatids")

LabRow_Names <- cell_order

# Define the order of traits for the heatmap
trait_name <- trait_name_all

# Order rows and columns for DEL matrix
order_index <- match(rownames(DEL_HPA_43_mat), cell_order)
DEL_HPA_43_mat_ordered_all <- DEL_HPA_43_mat[order(order_index), ]
DEL_HPA_43_mat_ordered_all <- DEL_HPA_43_mat_ordered_all[, match(trait_order, colnames(DEL_HPA_43_mat_ordered_all))]

# Order rows and columns for p-values
order_index <- match(rownames(p_values_fdr_DEL_HPA), cell_order)
p_values_fdr_DEL_HPA_star_ordered_all <- p_values_fdr_DEL_HPA[order(order_index), ]
p_values_fdr_DEL_HPA_star_ordered_all <- p_values_fdr_DEL_HPA_star_ordered_all[, match(trait_order, colnames(p_values_fdr_DEL_HPA_star_ordered_all))]
p_values_fdr_DEL_HPA_star_ordered_all <- ifelse(p_values_fdr_DEL_HPA_star_ordered_all < 0.05, "*", "")

# Extract the subset of the datasets for heatmap plot (22*20 (Tissue*Trait))
Cells_Subset <- c("inhibitory_neurons", "excitatory_neurons", "astrocytes", "microglial_cells",
                  "mesothelial_cells", "adipocytes",
                  "skeletal_myocytes", "b_cells", "cholangiocytes",
                  "intestinal_goblet_cells")

Traits_Subset <- c("Cholesterol", "HDL", "Platelet", "Redblood", "HbA1c",
                   "Menarche", "Menopause", "FVC", "BMI", "HeartRate", "BMD", "EA", "FI", "TDI", "VSWM1",
                   "Risk Taking", "Depression (self-report)", "Neuroticism", "Loneliness",'MoodAnxiety')




# Subset DEL_HPA_43 matrix
DEL_HPA_43_mat_ordered <- subset(DEL_HPA_43_mat_ordered_all, select = Traits_Subset)
DEL_HPA_43_mat_ordered <- DEL_HPA_43_mat_ordered[rownames(DEL_HPA_43_mat_ordered) %in% Cells_Subset, ]

# Subset p-values
p_values_fdr_DEL_HPA_star_ordered <- subset(p_values_fdr_DEL_HPA_star_ordered_all, select = Traits_Subset)
p_values_fdr_DEL_HPA_star_ordered <- p_values_fdr_DEL_HPA_star_ordered[rownames(p_values_fdr_DEL_HPA_star_ordered) %in% Cells_Subset, ]

# Define the row names for the heatmap
LabRow_Names <- c("Inhibitory neurons", "Excitatory neurons", "Astrocytes", "Microglial cells",
                  "Mesothelial cells", "Adipocytes",
                  "Skeletal myocytes", "B-cells", "Cholangiocytes",
                  "Intestinal goblet cells")

# Define the trait names for the heatmap
trait_name <- c("Cholesterol", "HDL", "Platelet", "Redblood", "HbA1c",
                "Menarche", "Menopause", "FVC", "BMI", "HeartRate", "BMD", "EA", "FI", "TDI", "VSWM1",
                "Risk Taking", "Depression", "Neuroticism", "Loneliness",'MoodAnxiety')
```

#### Heatmap: Effect Sizes of Subset of Whole Body Cell Types on Traits For Deletion

``` r
row_annotation <- data.frame(
  Tissue = c(rep("Brain", 4), rep("Non-Brain", 6))
)

row_annotation$Tissue <- factor(row_annotation$Tissue, levels = c("Brain", "Non-Brain"))

# Create column annotations
column_annotation <- data.frame(
  Trait = c(rep("Blood Assays",5),rep("Reproductive\n Factors",2),rep("Physical\n Measure",4),rep("Cognitive\n Metrics",4),rep("Mental\n Health",5))
)

column_annotation$Trait <- factor(column_annotation$Trait, levels = c("Blood Assays", "Reproductive\n Factors", "Physical\n Measure", "Cognitive\n Metrics", "Mental\n Health"))

# Define column heatmap annotation
column_ha <- HeatmapAnnotation(
  Trait = column_annotation$Trait,
  col = list(Trait = c("Blood Assays" = "#e78ac3", "Reproductive\n Factors" = "#8da0cb", "Physical\n Measure" = "#a6d854", "Cognitive\n Metrics" = "#66c2a5", "Mental\n Health" = "#fc8d62")),
  annotation_name = NULL,  # Remove both 'Supercluster' and 'Trait' labels
  show_legend = FALSE,
  show_annotation_name = FALSE
)

# Set row and column names for the ordered matrix
rownames(DEL_HPA_43_mat_ordered) <- LabRow_Names
colnames(DEL_HPA_43_mat_ordered) <- trait_name

# Create the heatmap without dendrograms
ht2 <- Heatmap(
  DEL_HPA_43_mat_ordered, name = "Effect size",
  cluster_rows = FALSE,  # Disable row clustering
  cluster_columns = FALSE,  # Disable column clustering
  show_row_dend = FALSE,  # Do not show row dendrogram
  show_column_dend = FALSE,  # Do not show column dendrogram
  top_annotation = column_ha,
  row_split = row_annotation$Tissue,  # Split rows into 2 groups
  column_split = column_annotation$Trait,  # Split columns into 5 groups
  row_title = NULL,
  row_title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "sans"),
  row_names_gp = gpar(col = "black", fontsize = 10, fontfamily = "sans", fontface = "bold"),
  column_title_gp = gpar(col = c('#e78ac3', '#8da0cb', '#a6d854', '#66c2a5', '#fc8d62'), fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(col = "black", fontsize = 10, rot = 45, fontfamily = "sans", fontface = "bold"),
  column_names_rot = 45,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(p_values_fdr_DEL_HPA_star_ordered[i, j], x = x, y = y - convertHeight(grobHeight(textGrob('*')), "mm"), gp = gpar(fontsize = 32, col = 'grey40', fontfamily = "sans", just = "center", fontface = "bold", hjust = 30))
  },
  col = color_function, row_gap = unit(0.5, "cm"),  # Adjust row margin (increase or decrease as needed)
  column_gap = unit(0.5, "cm"),
  show_heatmap_legend = FALSE,
  rect_gp = gpar(col = "lightgrey", lwd = 0.5)
)

# Draw the heatmap with padding
draw(ht2, padding = unit(c(0.5, 2, 0.5, 0.5), "cm"))
```

![](Fig2_files/figure-gfm/fig2_5%20-1.png)<!-- -->

#### Prepating Subset of Data for the Heatmap plot

##### Subsetting Data At Brain Cell Type Level - Brain vs. Non-Brain Traits and Neuronal vs. Non-Neuronal Brain Cell Type Gene Sets - Deletion

``` r
# Convert DEL and DUP Siletti data to matrices and set row and column names
DEL_Siletti_43_mat <- as.matrix(DEL_Siletti_43)
rownames(DEL_Siletti_43_mat) <- rownames(DEL_Siletti_43)
colnames(DEL_Siletti_43_mat) <- colnames(DEL_Siletti_43)

DUP_Siletti_43_mat <- as.matrix(DUP_Siletti_43)
rownames(DUP_Siletti_43_mat) <- rownames(DUP_Siletti_43)
colnames(DUP_Siletti_43_mat) <- colnames(DUP_Siletti_43)

# Read and process the ordering of cells
#Ordering_Cell <- read.csv('/Users/sayekazem/Desktop/New_Project_2023_Summer/PaperMultitrait2024_Data/Pleiotropy_Data/Siletti_SuperClust/df_n_gene_TDEP_Z_siletti_superclusters_cell_type_group.csv')
#Ordering_Cell <- subset(Ordering_Cell, select = c('Supercluster_ID', 'Neuron'))
Ordering_Cell=read.csv(paste0(Path,"BrainCells_Categories.csv"))
Neuron <- Ordering_Cell$Supercluster_ID[Ordering_Cell$Neuron == 'TRUE']
NonNeuron <- Ordering_Cell$Supercluster_ID[Ordering_Cell$Neuron == 'FALSE']

# Define the order of cells for the heatmap
cell_order <- c(Neuron, NonNeuron)

# Set row names for the heatmap
LabRow_Names <- gsub("_", " ", cell_order)

# Define the order of traits for the heatmap
trait_order <- trait_name_all

# Order rows and columns for DEL matrix
order_index <- match(rownames(DEL_Siletti_43_mat), cell_order)
DEL_Siletti_43_mat_ordered_all <- DEL_Siletti_43_mat[order(order_index), ]
DEL_Siletti_43_mat_ordered_all <- DEL_Siletti_43_mat_ordered_all[, match(trait_order, colnames(DEL_Siletti_43_mat_ordered_all))]

# Order rows and columns for p-values
order_index <- match(rownames(p_values_fdr_DEL_Siletti), cell_order)
p_values_fdr_DEL_Siletti_star_ordered_all <- p_values_fdr_DEL_Siletti[order(order_index), ]
p_values_fdr_DEL_Siletti_star_ordered_all <- p_values_fdr_DEL_Siletti_star_ordered_all[, match(trait_order, colnames(p_values_fdr_DEL_Siletti_star_ordered_all))]
p_values_fdr_DEL_Siletti_star_ordered_all <- ifelse(p_values_fdr_DEL_Siletti_star_ordered_all < 0.05, "*", "")

# Extract the subset of the datasets for heatmap plot (22*20 (Tissue*Trait))

# Define the subset of cells for the heatmap
Cells_Subset <- c("amygdala_excitatory", "cge_interneuron", "upper_layer_intratelencephalic", "hippocampal_ca1_3",
                  "mge_interneuron", "upper_rhombic_lip", "microglia", "astrocyte", "choroid_plexus", "committed_oligodendrocyte_precursor")

# Define the subset of traits for the heatmap
Traits_Subset <- c("Cholesterol", "HDL", "Platelet", "Redblood", "HbA1c",
                   "Menarche", "Menopause", "FVC", "BMI", "HeartRate", "BMD", "EA", "FI", "TDI", "VSWM1",
                   "Risk Taking", "Depression (self-report)", "Neuroticism", "Loneliness",'MoodAnxiety')

# Subset DEL_Siletti_43 matrix
DEL_Siletti_43_mat_ordered <- subset(DEL_Siletti_43_mat_ordered_all, select = Traits_Subset)
DEL_Siletti_43_mat_ordered <- DEL_Siletti_43_mat_ordered[rownames(DEL_Siletti_43_mat_ordered) %in% Cells_Subset, ]

# Subset p-values
p_values_fdr_DEL_Siletti_star_ordered <- subset(p_values_fdr_DEL_Siletti_star_ordered_all, select = Traits_Subset)
p_values_fdr_DEL_Siletti_star_ordered <- p_values_fdr_DEL_Siletti_star_ordered[rownames(p_values_fdr_DEL_Siletti_star_ordered) %in% Cells_Subset, ]

# Define the row names for the heatmap
LabRow_Names <- c("Amygdala excitatory", "CGE interneuron", "Upper-layer IT", "Hippocampus CA1-3",
                  "MGE interneurons", "Upper rhombic lip", "Microglia", "Astrocyte", "Choroid plexus", "OPC")

# Define the trait names for the heatmap
trait_name <- c("Cholesterol", "HDL", "Platelet", "Redblood", "HbA1c",
                "Menarche", "Menopause", "FVC", "BMI", "HeartRate", "BMD", "EA", "FI", "TDI", "VSWM1",
                "Risk Taking", "Depression", "Neuroticism", "Loneliness",'MoodAnxiety')
```

#### Heatmap: Effect Sizes of Subset of Whole Brain Cell Types on Traits For Deletion

``` r
row_annotation <- data.frame(
  Tissue = c(rep("Neuronal", 6), rep("Non-Neuronal", 4))
)

row_annotation$Tissue <- factor(row_annotation$Tissue, levels = c("Neuronal", "Non-Neuronal"))

# Create column annotations
column_annotation <- data.frame(
  Trait = c(rep("Blood Assays",5),rep("Reproductive\n Factors",2),rep("Physical\n Measure",4),rep("Cognitive\n Metrics",4),rep("Mental\n Health",5))
)

column_annotation$Trait <- factor(column_annotation$Trait, levels = c("Blood Assays", "Reproductive\n Factors", "Physical\n Measure", "Cognitive\n Metrics", "Mental\n Health"))

# Define column heatmap annotation
column_ha <- HeatmapAnnotation(
  Trait = column_annotation$Trait,
  col = list(Trait = c("Blood Assays" = "#e78ac3", "Reproductive\n Factors" = "#8da0cb", "Physical\n Measure" = "#a6d854", "Cognitive\n Metrics" = "#66c2a5", "Mental\n Health" = "#fc8d62")),
  annotation_name = NULL,  # Remove both 'Supercluster' and 'Trait' labels
  show_legend = FALSE,
  show_annotation_name = FALSE
)

# Set row and column names for the ordered matrix
rownames(DEL_Siletti_43_mat_ordered) <- LabRow_Names
colnames(DEL_Siletti_43_mat_ordered) <- trait_name

# Create the heatmap without dendrograms
ht3 <- Heatmap(
  DEL_Siletti_43_mat_ordered, name = "Effect size",
  cluster_rows = FALSE,  # Disable row clustering
  cluster_columns = FALSE,  # Disable column clustering
  show_row_dend = FALSE,  # Do not show row dendrogram
  show_column_dend = FALSE,  # Do not show column dendrogram
  top_annotation = column_ha,
  row_split = row_annotation$Tissue,  # Split rows into 2 groups
  column_split = column_annotation$Trait,  # Split columns into 5 groups
  row_title = NULL,
  row_title_gp = gpar(fontsize = 8, fontface = "bold", fontfamily = "sans"),
  row_names_gp = gpar(col = "black", fontsize = 10, fontfamily = "sans", fontface = "bold"),
  column_title_gp = gpar(col = c('#e78ac3', '#8da0cb', '#a6d854', '#66c2a5', '#fc8d62'), fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(col = "black", fontsize = 10, rot = 45, fontfamily = "sans", fontface = "bold"),
  column_names_rot = 45,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(p_values_fdr_DEL_Siletti_star_ordered[i, j], x = x, y = y - convertHeight(grobHeight(textGrob('*')), "mm"), gp = gpar(fontsize = 32, col = 'grey40', fontfamily = "sans", just = "center", fontface = "bold", hjust = 30))
  },
  col = color_function, row_gap = unit(0.5, "cm"),  # Adjust row margin (increase or decrease as needed)
  column_gap = unit(0.5, "cm"),
  show_heatmap_legend = FALSE,
  rect_gp = gpar(col = "lightgrey", lwd = 0.5)
)

# Draw the heatmap with padding
draw(ht3, padding = unit(c(0.5, 2, 0.5, 0.5), "cm"))
```

![](Fig2_files/figure-gfm/fig2_7-1.png)<!-- -->

#### Merge All 3 Plots in a single plot

``` r
ht1_grob <- grid.grabExpr(draw(ht1, padding = unit(c(0, 1.5, 0, 0), "cm")))
ht2_grob <- grid.grabExpr(draw(ht2, padding = unit(c(0, 1.5, 0, 0), "cm")))
ht3_grob <- grid.grabExpr(draw(ht3, padding = unit(c(0, 1.5, 0, 0), "cm")))

p_stacked <- ggarrange(ggplotify::as.ggplot(ht1_grob), 
                       ggplotify::as.ggplot(ht2_grob), 
                       ggplotify::as.ggplot(ht3_grob),
                       ncol = 1, nrow = 3, align = "hv", widths = c(0.5, 0.5,0.5), heights = c(5, 5,5))

(p_stacked)
```

![](Fig2_files/figure-gfm/fig2_8-1.png)<!-- -->

#### The same Analysis is done for Duplication

<br><br>

#### Point Range Plot for Full Data

##### Percentage of Significant associations across Traits and Genesets Categories

``` r
###------------ ## Pleiotropy Across Brain & Non-Brain Tissues (Fantom) for Brain & Non-Brain Traits ------------------#############
#----------------------------------------------------------------------------------------------------------
# Define Brain Traits
BrainTraits <- c("Neuroticism", "Depression (self-report)", "GP Visit for Anxiety or Depression", "Mood Swings", "Risk Taking", "Loneliness", "Irritability", "Miserableness", "Guilty Feeling",
                 "Reaction", "TDI", "EA", "FI", "VSWM1", "VSWM2", "Symbol Digit Substitution",
                 "Trail Making 1", "Trail Making 2",'MoodAnxiety')

## Tissue Level (Fantom)####

# Subset p-values for Brain and Non-Brain traits
DEL_pval_Brain_Fantom <- subset(p_values_fdr_DEL_Fantom_star_ordered_all, select = BrainTraits)
DEL_pval_NBrain_Fantom <- p_values_fdr_DEL_Fantom_star_ordered_all[, setdiff(colnames(p_values_fdr_DEL_Fantom_star_ordered_all), BrainTraits)]

# Compute proportion in general (for bar plot)
BrainTrait_BrainTissue_DEL_Fantom <- mean(DEL_pval_Brain_Fantom[1:27, ] == '*') * 100
NonBrainTrait_BrainTissue_DEL_Fantom <- mean(DEL_pval_NBrain_Fantom[1:27, ] == '*') * 100

BrainTrait_NonBrainTissue_DEL_Fantom <- mean(DEL_pval_Brain_Fantom[28:60, ] == '*') * 100
NonBrainTrait_NonBrainTissue_DEL_Fantom <- mean(DEL_pval_NBrain_Fantom[28:60, ] == '*') * 100

#### Extract Mean and SD From Data (for point range plot) ########--------------------------------------------------------------------------------
# Brain Traits in Brain Tissue
BrainTrait_BrainTissue_DEL_Fantom_Mean_SD <- as.data.frame(DEL_pval_Brain_Fantom[1:27, ]) %>% mutate(across(everything(), ~ ifelse(. == "*", 1, 0)))
BrainTrait_BrainTissue_DEL_Fantom_Mean_SD <- BrainTrait_BrainTissue_DEL_Fantom_Mean_SD %>%
  rowwise() %>%
  summarise(
    Mean = mean(c_across(everything()))
    # SD = sd(c_across(everything()))
  )
BrainTrait_BrainTissue_DEL_Fantom_Mean_SD <- as.data.frame(BrainTrait_BrainTissue_DEL_Fantom_Mean_SD)
rownames(BrainTrait_BrainTissue_DEL_Fantom_Mean_SD) <- rownames(DEL_pval_Brain_Fantom)[1:27]
BrainTrait_BrainTissue_DEL_Fantom_Mean_SD$Tissue_Type <- "Brain"
BrainTrait_BrainTissue_DEL_Fantom_Mean_SD$Trait_Type <- "Brain"

BrainTrait_BrainTissue_DEL_Fantom_Mean_SD$SD <- sd(BrainTrait_BrainTissue_DEL_Fantom_Mean_SD$Mean)
BrainTrait_BrainTissue_DEL_Fantom_Mean_SD$Mean <- mean(BrainTrait_BrainTissue_DEL_Fantom_Mean_SD$Mean)
BrainTrait_BrainTissue_DEL_Fantom_Mean_SD <- BrainTrait_BrainTissue_DEL_Fantom_Mean_SD[1, ]
rownames(BrainTrait_BrainTissue_DEL_Fantom_Mean_SD) <- "stat1"

# Non-Brain Traits in Brain Tissue
NonBrainTrait_BrainTissue_DEL_Fantom_Mean_SD <- as.data.frame(DEL_pval_NBrain_Fantom[1:27, ]) %>% mutate(across(everything(), ~ ifelse(. == "*", 1, 0)))
NonBrainTrait_BrainTissue_DEL_Fantom_Mean_SD <- NonBrainTrait_BrainTissue_DEL_Fantom_Mean_SD %>%
  rowwise() %>%
  summarise(
    Mean = mean(c_across(everything()))
    # SD = sd(c_across(everything()))
  )
NonBrainTrait_BrainTissue_DEL_Fantom_Mean_SD <- as.data.frame(NonBrainTrait_BrainTissue_DEL_Fantom_Mean_SD)
rownames(NonBrainTrait_BrainTissue_DEL_Fantom_Mean_SD) <- rownames(DEL_pval_NBrain_Fantom)[1:27]
NonBrainTrait_BrainTissue_DEL_Fantom_Mean_SD$Tissue_Type <- "Brain"
NonBrainTrait_BrainTissue_DEL_Fantom_Mean_SD$Trait_Type <- "NonBrain"

NonBrainTrait_BrainTissue_DEL_Fantom_Mean_SD$SD <- sd(NonBrainTrait_BrainTissue_DEL_Fantom_Mean_SD$Mean)
NonBrainTrait_BrainTissue_DEL_Fantom_Mean_SD$Mean <- mean(NonBrainTrait_BrainTissue_DEL_Fantom_Mean_SD$Mean)
NonBrainTrait_BrainTissue_DEL_Fantom_Mean_SD <- NonBrainTrait_BrainTissue_DEL_Fantom_Mean_SD[1, ]
rownames(NonBrainTrait_BrainTissue_DEL_Fantom_Mean_SD) <- "stat2"

# Brain Traits in Non-Brain Tissue
BrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD <- as.data.frame(DEL_pval_Brain_Fantom[28:60, ]) %>% mutate(across(everything(), ~ ifelse(. == "*", 1, 0)))
BrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD <- BrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD %>%
  rowwise() %>%
  summarise(
    Mean = mean(c_across(everything()))
    # SD = sd(c_across(everything()))
  )
BrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD <- as.data.frame(BrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD)
rownames(BrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD) <- rownames(DEL_pval_Brain_Fantom)[28:60]
BrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD$Tissue_Type <- "NonBrain"
BrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD$Trait_Type <- "Brain"

BrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD$SD <- sd(BrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD$Mean)
BrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD$Mean <- mean(BrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD$Mean)
BrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD <- BrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD[1, ]
rownames(BrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD) <- "stat3"

# Non-Brain Traits in Non-Brain Tissue
NonBrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD <- as.data.frame(DEL_pval_NBrain_Fantom[28:60, ]) %>% mutate(across(everything(), ~ ifelse(. == "*", 1, 0)))
NonBrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD <- NonBrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD %>%
  rowwise() %>%
  summarise(
    Mean = mean(c_across(everything()))
    # SD = sd(c_across(everything()))
  )
NonBrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD <- as.data.frame(NonBrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD)
rownames(NonBrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD) <- rownames(DEL_pval_NBrain_Fantom)[28:60]
NonBrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD$Tissue_Type <- "NonBrain"
NonBrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD$Trait_Type <- "NonBrain"

NonBrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD$SD <- sd(NonBrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD$Mean)
NonBrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD$Mean <- mean(NonBrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD$Mean)
NonBrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD <- NonBrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD[1, ]
rownames(NonBrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD) <- "stat4"

# Combine data for point range plot
Point_Range_Data <- rbind(BrainTrait_BrainTissue_DEL_Fantom_Mean_SD, NonBrainTrait_BrainTissue_DEL_Fantom_Mean_SD, BrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD, NonBrainTrait_NonBrainTissue_DEL_Fantom_Mean_SD)

################# ---------- point range plot ###################
# Assuming 'data' contains columns Mean, SD, Tissue_Type, and Trait as described in your example.
# png("/Users/sayekazem/Desktop/New_Project_2023_Summer/PaperMultitrait2024_Figs/Fig2/PointRange_DEL_Fantom_Nov14.png", width = 2200, height = 1500, res = 450)
h1=ggplot(Point_Range_Data, aes(x = Trait_Type, y = Mean, color = Tissue_Type)) +
  geom_point(position = position_dodge(width = 0.3), size = 5, alpha = 0.8) + 
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), 
                  position = position_dodge(width = 0.3), 
                  linewidth = 0.8, alpha = 0.8, width = 0.3) +  # Adjust point and line sizes and transparency
  labs(x = "Trait Category", y = "Proportion of Associations",color = "Tissue Category") + 
  ylim(-0.06, 0.65) +
  scale_color_manual(values = c('Brain' = '#01665e', 'NonBrain' = 'seagreen3'),name = "Tissue Category") +  # Set custom colors
  theme_prism() +
  theme( 
        axis.text.y = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12) )


#####################################################

### Whole Body Cell types (HPA)


# Subset p-values for Brain and Non-Brain traits
DEL_pval_Brain_HPA <- subset(p_values_fdr_DEL_HPA_star_ordered_all, select = BrainTraits)
DEL_pval_NBrain_HPA <- p_values_fdr_DEL_HPA_star_ordered_all[, setdiff(colnames(p_values_fdr_DEL_HPA_star_ordered_all), BrainTraits)]

# Compute proportion in general (for bar plot)
BrainTrait_BrainCell_DEL_HPA <- mean(DEL_pval_Brain_HPA[1:7, ] == '*') * 100
NonBrainTrait_BrainCell_DEL_HPA <- mean(DEL_pval_NBrain_HPA[1:7, ] == '*') * 100

BrainTrait_NonBrainCell_DEL_HPA <- mean(DEL_pval_Brain_HPA[8:81, ] == '*') * 100
NonBrainTrait_NonBrainCell_DEL_HPA <- mean(DEL_pval_NBrain_HPA[8:81, ] == '*') * 100

#####################################################
#### Extract Mean and SD From Data (for point range plot) ########--------------------------------------------------------------------------------
# Brain Traits in Brain Tissue
BrainTrait_BrainTissue_DEL_HPA_Mean_SD <- as.data.frame(DEL_pval_Brain_HPA[1:7, ]) %>% mutate(across(everything(), ~ ifelse(. == "*", 1, 0)))
BrainTrait_BrainTissue_DEL_HPA_Mean_SD <- BrainTrait_BrainTissue_DEL_HPA_Mean_SD %>%
  rowwise() %>%
  summarise(
    Mean = mean(c_across(everything()))
    # SD = sd(c_across(everything()))
  )
BrainTrait_BrainTissue_DEL_HPA_Mean_SD <- as.data.frame(BrainTrait_BrainTissue_DEL_HPA_Mean_SD)
rownames(BrainTrait_BrainTissue_DEL_HPA_Mean_SD) <- rownames(DEL_pval_Brain_HPA)[1:7]
BrainTrait_BrainTissue_DEL_HPA_Mean_SD$Tissue_Type <- "Brain"
BrainTrait_BrainTissue_DEL_HPA_Mean_SD$Trait_Type <- "Brain"

BrainTrait_BrainTissue_DEL_HPA_Mean_SD$SD <- sd(BrainTrait_BrainTissue_DEL_HPA_Mean_SD$Mean)
BrainTrait_BrainTissue_DEL_HPA_Mean_SD$Mean <- mean(BrainTrait_BrainTissue_DEL_HPA_Mean_SD$Mean)
BrainTrait_BrainTissue_DEL_HPA_Mean_SD <- BrainTrait_BrainTissue_DEL_HPA_Mean_SD[1, ]
rownames(BrainTrait_BrainTissue_DEL_HPA_Mean_SD) <- "stat1"

# Non-Brain Traits in Brain Tissue
NonBrainTrait_BrainTissue_DEL_HPA_Mean_SD <- as.data.frame(DEL_pval_NBrain_HPA[1:7, ]) %>% mutate(across(everything(), ~ ifelse(. == "*", 1, 0)))
NonBrainTrait_BrainTissue_DEL_HPA_Mean_SD <- NonBrainTrait_BrainTissue_DEL_HPA_Mean_SD %>%
  rowwise() %>%
  summarise(
    Mean = mean(c_across(everything()))
    # SD = sd(c_across(everything()))
  )
NonBrainTrait_BrainTissue_DEL_HPA_Mean_SD <- as.data.frame(NonBrainTrait_BrainTissue_DEL_HPA_Mean_SD)
rownames(NonBrainTrait_BrainTissue_DEL_HPA_Mean_SD) <- rownames(DEL_pval_NBrain_HPA)[1:7]
NonBrainTrait_BrainTissue_DEL_HPA_Mean_SD$Tissue_Type <- "Brain"
NonBrainTrait_BrainTissue_DEL_HPA_Mean_SD$Trait_Type <- "NonBrain"

NonBrainTrait_BrainTissue_DEL_HPA_Mean_SD$SD <- sd(NonBrainTrait_BrainTissue_DEL_HPA_Mean_SD$Mean)
NonBrainTrait_BrainTissue_DEL_HPA_Mean_SD$Mean <- mean(NonBrainTrait_BrainTissue_DEL_HPA_Mean_SD$Mean)
NonBrainTrait_BrainTissue_DEL_HPA_Mean_SD <- NonBrainTrait_BrainTissue_DEL_HPA_Mean_SD[1, ]
rownames(NonBrainTrait_BrainTissue_DEL_HPA_Mean_SD) <- "stat2"

# Brain Traits in Non-Brain Tissue
BrainTrait_NonBrainTissue_DEL_HPA_Mean_SD <- as.data.frame(DEL_pval_Brain_HPA[8:81, ]) %>% mutate(across(everything(), ~ ifelse(. == "*", 1, 0)))
BrainTrait_NonBrainTissue_DEL_HPA_Mean_SD <- BrainTrait_NonBrainTissue_DEL_HPA_Mean_SD %>%
  rowwise() %>%
  summarise(
    Mean = mean(c_across(everything()))
    # SD = sd(c_across(everything()))
  )
BrainTrait_NonBrainTissue_DEL_HPA_Mean_SD <- as.data.frame(BrainTrait_NonBrainTissue_DEL_HPA_Mean_SD)
rownames(BrainTrait_NonBrainTissue_DEL_HPA_Mean_SD) <- rownames(DEL_pval_Brain_HPA)[8:81]
BrainTrait_NonBrainTissue_DEL_HPA_Mean_SD$Tissue_Type <- "NonBrain"
BrainTrait_NonBrainTissue_DEL_HPA_Mean_SD$Trait_Type <- "Brain"

BrainTrait_NonBrainTissue_DEL_HPA_Mean_SD$SD <- sd(BrainTrait_NonBrainTissue_DEL_HPA_Mean_SD$Mean)
BrainTrait_NonBrainTissue_DEL_HPA_Mean_SD$Mean <- mean(BrainTrait_NonBrainTissue_DEL_HPA_Mean_SD$Mean)
BrainTrait_NonBrainTissue_DEL_HPA_Mean_SD <- BrainTrait_NonBrainTissue_DEL_HPA_Mean_SD[1, ]
rownames(BrainTrait_NonBrainTissue_DEL_HPA_Mean_SD) <- "stat3"

# Non-Brain Traits in Non-Brain Tissue
NonBrainTrait_NonBrainTissue_DEL_HPA_Mean_SD <- as.data.frame(DEL_pval_NBrain_HPA[8:81, ]) %>% mutate(across(everything(), ~ ifelse(. == "*", 1, 0)))
NonBrainTrait_NonBrainTissue_DEL_HPA_Mean_SD <- NonBrainTrait_NonBrainTissue_DEL_HPA_Mean_SD %>%
  rowwise() %>%
  summarise(
    Mean = mean(c_across(everything()))
    # SD = sd(c_across(everything()))
  )
NonBrainTrait_NonBrainTissue_DEL_HPA_Mean_SD <- as.data.frame(NonBrainTrait_NonBrainTissue_DEL_HPA_Mean_SD)
rownames(NonBrainTrait_NonBrainTissue_DEL_HPA_Mean_SD) <- rownames(DEL_pval_NBrain_HPA)[8:81]
NonBrainTrait_NonBrainTissue_DEL_HPA_Mean_SD$Tissue_Type <- "NonBrain"
NonBrainTrait_NonBrainTissue_DEL_HPA_Mean_SD$Trait_Type <- "NonBrain"

NonBrainTrait_NonBrainTissue_DEL_HPA_Mean_SD$SD <- sd(NonBrainTrait_NonBrainTissue_DEL_HPA_Mean_SD$Mean)
NonBrainTrait_NonBrainTissue_DEL_HPA_Mean_SD$Mean <- mean(NonBrainTrait_NonBrainTissue_DEL_HPA_Mean_SD$Mean)
NonBrainTrait_NonBrainTissue_DEL_HPA_Mean_SD <- NonBrainTrait_NonBrainTissue_DEL_HPA_Mean_SD[1, ]
rownames(NonBrainTrait_NonBrainTissue_DEL_HPA_Mean_SD) <- "stat4"

# Combine data for point range plot
Point_Range_Data <- rbind(BrainTrait_BrainTissue_DEL_HPA_Mean_SD, NonBrainTrait_BrainTissue_DEL_HPA_Mean_SD, BrainTrait_NonBrainTissue_DEL_HPA_Mean_SD, NonBrainTrait_NonBrainTissue_DEL_HPA_Mean_SD)

################# ---------- point range plot ###################
# Assuming 'data' contains columns Mean, SD, Tissue_Type, and Trait as described in your example.
#png("/Users/sayekazem/Desktop/New_Project_2023_Summer/PaperMultitrait2024_Figs/Fig2/PointRange_DEL_HPA_Nov14.png", width = 2200, height = 1500, res = 450)

h2=ggplot(Point_Range_Data, aes(x = Trait_Type, y = Mean, color = Tissue_Type)) +
  geom_point(position = position_dodge(width = 0.3), size = 5, alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                position = position_dodge(width = 0.3),
                linewidth = 0.8, alpha = 0.8, width = 0.3) +  # Adjust point and line sizes and transparency
  labs(x = "Trait Category", y = "Proportion of Associations",color = "Cell Type Category") +
  ylim(-0.06, 0.65) +
  scale_color_manual(values = c('Brain' = '#7b3294', 'NonBrain' = '#c2a5cf',name = "Cell Type Category")) +  # Set custom colors
  theme_prism() +
  theme(axis.line = element_line(linewidth = 0.8),
        axis.text.y = element_text(size = 15),
        legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12) )



#####################
# Whole Brain Cell Types (Siletti) ###

DEL_pval_Brain_SilettiSuper=subset(p_values_fdr_DEL_Siletti_star_ordered_all,select=BrainTraits)

DEL_pval_NBrain_SilettiSuper <- p_values_fdr_DEL_Siletti_star_ordered_all[ , setdiff(colnames(p_values_fdr_DEL_Siletti_star_ordered_all), BrainTraits)]
#
BrainTrait_BrainCell_DEL_SilettiSuper=mean(DEL_pval_Brain_SilettiSuper[1:21,]=='*')*100
NonBrainTrait_BrainCell_DEL_SilettiSuper=mean(DEL_pval_NBrain_SilettiSuper[1:21,]=='*')*100

BrainTrait_NonBrainCell_DEL_SilettiSuper=mean(DEL_pval_Brain_SilettiSuper[22:31,]=='*')*100
NonBrainTrait_NonBrainCell_DEL_SilettiSuper=mean(DEL_pval_NBrain_SilettiSuper[22:31,]=='*')*100


#################################################################
#### Extract Mean and Sd From Data (for point range plot) ########--------------------------------------------------------------------------------
BrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD <- as.data.frame(DEL_pval_Brain_SilettiSuper[1:21,]) %>% mutate(across(everything(), ~ ifelse(. == "*", 1, 0)))
BrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD <- BrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD %>%
  rowwise() %>%
  summarise(
    Mean = mean(c_across(everything()))
    #SD = sd(c_across(everything()))
  )
BrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD=as.data.frame(BrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD)
rownames(BrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD) <- rownames(DEL_pval_Brain_SilettiSuper)[1:21]
BrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD$Tissue_Type="Neuronal"
BrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD$Trait_Type="Brain"

BrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD$SD=sd(BrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD$Mean)
BrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD$Mean=mean(BrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD$Mean)
BrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD=BrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD[1,]
rownames(BrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD)="stat1"

##
NonBrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD <- as.data.frame(DEL_pval_NBrain_SilettiSuper[1:21,]) %>% mutate(across(everything(), ~ ifelse(. == "*", 1, 0)))
NonBrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD <- NonBrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD %>%
  rowwise() %>%
  summarise(
    Mean = mean(c_across(everything()))
    #SD = sd(c_across(everything()))
  )
NonBrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD=as.data.frame(NonBrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD)
rownames(NonBrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD) <- rownames(DEL_pval_NBrain_SilettiSuper)[1:21]
NonBrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD$Tissue_Type="Neuronal"
NonBrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD$Trait_Type="NonBrain"

NonBrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD$SD=sd(NonBrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD$Mean)
NonBrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD$Mean=mean(NonBrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD$Mean)
NonBrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD=NonBrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD[1,]
rownames(NonBrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD)="stat2"


###

BrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD <- as.data.frame(DEL_pval_Brain_SilettiSuper[22:31,]) %>% mutate(across(everything(), ~ ifelse(. == "*", 1, 0)))
BrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD <- BrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD %>%
  rowwise() %>%
  summarise(
    Mean = mean(c_across(everything()))
    #SD = sd(c_across(everything()))
  )
BrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD=as.data.frame(BrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD)
rownames(BrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD) <- rownames(DEL_pval_Brain_SilettiSuper)[22:31]
BrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD$Tissue_Type="NonNeuronal"
BrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD$Trait_Type="Brain"

BrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD$SD=sd(BrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD$Mean)
BrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD$Mean=mean(BrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD$Mean)
BrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD=BrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD[1,]
rownames(BrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD)="stat3"

##
NonBrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD <- as.data.frame(DEL_pval_NBrain_SilettiSuper[22:31,]) %>% mutate(across(everything(), ~ ifelse(. == "*", 1, 0)))
NonBrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD <- NonBrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD %>%
  rowwise() %>%
  summarise(
    Mean = mean(c_across(everything()))
    #SD = sd(c_across(everything()))
  )
NonBrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD=as.data.frame(NonBrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD)
rownames(NonBrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD) <- rownames(DEL_pval_NBrain_SilettiSuper)[22:31]
NonBrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD$Tissue_Type="NonNeuronal"
NonBrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD$Trait_Type="NonBrain"

NonBrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD$SD=sd(NonBrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD$Mean)
NonBrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD$Mean=mean(NonBrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD$Mean)
NonBrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD=NonBrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD[1,]
rownames(NonBrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD)="stat4"
###
# Point_Range_Data For Fig2#
Point_Range_Data=(rbind(BrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD,NonBrainTrait_BrainTissue_DEL_SilettiSuper_Mean_SD,BrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD,NonBrainTrait_NonBrainTissue_DEL_SilettiSuper_Mean_SD))

################# ---------- point range plot ###################
# Assuming 'data' contains columns Mean, SD, Tissue_Type, and Trait as described in your example.
#png("/Users/sayekazem/Desktop/New_Project_2023_Summer/PaperMultitrait2024_Figs/Fig2/PointRange_DEL_SilettiSuper_Nov14.png", width = 2200 , height = 1500,res = 450)

h3=ggplot(Point_Range_Data, aes(x = Trait_Type, y = Mean, color = Tissue_Type)) +
  geom_point(position = position_dodge(width = 0.3), size = 5, alpha = 0.8) + 
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), 
                position = position_dodge(width = 0.3), 
                linewidth = 0.8, alpha = 0.8, width = 0.3) +  # Adjust point and line sizes and transparency
  labs(x = "Trait Category", y = "Proportion of Associations",color="Brain Cell Category") + 
  ylim(-0.06, 0.65) +
  scale_color_manual(values = c('Neuronal' = '#e66101', 'NonNeuronal' = '#fee08b'),name="Brain Cell Category") +  # Set custom colors
  theme_prism() +
  theme(axis.line = element_line(linewidth = 0.8), 
        axis.text.y = element_text(size = 15), 
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12) )

#dev.off()
p_stacked_2 <- ggarrange(h1, 
                       h2, 
                       h3,
                       ncol = 1, nrow = 3, align = "hv", widths = c(1, 1,1), heights = c(1, 1,1))

(p_stacked_2)
```

![](Fig2_files/figure-gfm/fig2_9-1.png)<!-- -->

<br><br>

#### The same Analysis is done for Duplication

<br><br>

#### Stats for different compariosns using Porportion Test

##### Test1 : “Del has more significant associations that Dup at Whole Body Tissue Level”

``` r
# pvalue Datasets

#All Together
DELs=p_values_fdr_DEL
DUPs=p_values_fdr_DUP
threshold <- 0.05
#p_values_fdr_DEL_Fantom
#p_values_fdr_DEL_HPA
#p_values_fdr_DEL_Siletti
#p_values_fdr_DUP_Fantom
#p_values_fdr_DUP_HPA
#p_values_fdr_DUP_Siletti
#################################

group_A_successes <- sum(p_values_fdr_DEL_Fantom<0.05)
group_A_total <- 60*43 
group_B_successes <- sum(p_values_fdr_DUP_Fantom<0.05)
group_B_total <- 60*43 
## One-tailed p-value: 0.003099545 

# Calculate proportions
prop_A <- group_A_successes / group_A_total
prop_B <- group_B_successes / group_B_total

# Calculate the pooled proportion
pooled_prop <- (group_A_successes + group_B_successes) / (group_A_total + group_B_total)

# Calculate the standard error
se <- sqrt(pooled_prop * (1 - pooled_prop) * (1 / group_A_total + 1 / group_B_total))

# Calculate the Z-score
z_observed <- (prop_A - prop_B) / se

# One-tailed p-value (assuming prop_A > prop_B)
p_value <- 1 - pnorm(z_observed)
```

``` r
# Print the results
paste0("Number of Significant associations: ",group_A_successes," & ",group_B_successes)
```

    ## [1] "Number of Significant associations: 568 & 352"

``` r
paste0("Percentage of Significant associations: ",prop_A," & ",prop_B)
```

    ## [1] "Percentage of Significant associations: 0.22015503875969 & 0.136434108527132"

``` r
paste0("Observed Difference: ", prop_A - prop_B)
```

    ## [1] "Observed Difference: 0.0837209302325581"

``` r
paste0("Z-score:", z_observed)
```

    ## [1] "Z-score:7.85600515931389"

``` r
paste0("One-tailed p-value: ", p_value)
```

    ## [1] "One-tailed p-value: 1.99840144432528e-15"

##### Test2 : “Del has more significant associations that Dup at Whole Body Cell Type Level”

``` r
# pvalue Datasets

#################################

group_A_successes <- sum(p_values_fdr_DEL_HPA<0.05)
group_A_total <- 81*43 
group_B_successes <- sum(p_values_fdr_DUP_HPA<0.05)
group_B_total <- 81*43 
## One-tailed p-value: 0.003099545 

# Calculate proportions
prop_A <- group_A_successes / group_A_total
prop_B <- group_B_successes / group_B_total

# Calculate the pooled proportion
pooled_prop <- (group_A_successes + group_B_successes) / (group_A_total + group_B_total)

# Calculate the standard error
se <- sqrt(pooled_prop * (1 - pooled_prop) * (1 / group_A_total + 1 / group_B_total))

# Calculate the Z-score
z_observed <- (prop_A - prop_B) / se

# One-tailed p-value (assuming prop_A > prop_B)
p_value <- 1 - pnorm(z_observed)
```

``` r
# Print the results
paste0("Number of Significant associations: ",group_A_successes," & ",group_B_successes)
```

    ## [1] "Number of Significant associations: 494 & 421"

``` r
paste0("Percentage of Significant associations: ",prop_A," & ",prop_B)
```

    ## [1] "Percentage of Significant associations: 0.141831754234855 & 0.120872810795291"

``` r
paste0("Observed Difference: ", prop_A - prop_B)
```

    ## [1] "Observed Difference: 0.0209589434395636"

``` r
paste0("Z-score:", z_observed)
```

    ## [1] "Z-score:2.58934829542889"

``` r
paste0("One-tailed p-value: ", p_value)
```

    ## [1] "One-tailed p-value: 0.00480788906913654"

##### Test3 : “Del has more significant associations that Dup at Brain Cell Type Level”

``` r
# pvalue Datasets

#################################

group_A_successes <- sum(p_values_fdr_DEL_Siletti<0.05)
group_A_total <- 31*43 
group_B_successes <- sum(p_values_fdr_DUP_Siletti<0.05)
group_B_total <- 31*43 
## One-tailed p-value: 0.003099545 

# Calculate proportions
prop_A <- group_A_successes / group_A_total
prop_B <- group_B_successes / group_B_total

# Calculate the pooled proportion
pooled_prop <- (group_A_successes + group_B_successes) / (group_A_total + group_B_total)

# Calculate the standard error
se <- sqrt(pooled_prop * (1 - pooled_prop) * (1 / group_A_total + 1 / group_B_total))

# Calculate the Z-score
z_observed <- (prop_A - prop_B) / se

# One-tailed p-value (assuming prop_A > prop_B)
p_value <- 1 - pnorm(z_observed)
```

``` r
# Print the results
paste0("Number of Significant associations: ",group_A_successes," & ",group_B_successes)
```

    ## [1] "Number of Significant associations: 196 & 162"

``` r
paste0("Percentage of Significant associations: ",prop_A," & ",prop_B)
```

    ## [1] "Percentage of Significant associations: 0.147036759189797 & 0.121530382595649"

``` r
paste0("Observed Difference: ", prop_A - prop_B)
```

    ## [1] "Observed Difference: 0.0255063765941485"

``` r
paste0("Z-score:", z_observed)
```

    ## [1] "Z-score:1.93129925968404"

``` r
paste0("One-tailed p-value: ", p_value)
```

    ## [1] "One-tailed p-value: 0.0267230274477795"

##### Test4 : “There is a similar proportion of significant associations with whole body cell type compared to tissue gene sets”

``` r
Fanotm=(rbind(p_values_fdr_DEL_Fantom,p_values_fdr_DUP_Fantom))
HPA=(rbind(p_values_fdr_DEL_HPA,p_values_fdr_DUP_HPA))

# pvalue Datasets

#################################

group_A_successes <- sum(Fanotm<0.05)
group_A_total <- 60*43*2 
group_B_successes <- sum(HPA<0.05)
group_B_total <- 81*43*2 
## One-tailed p-value: 0.003099545 

# Calculate proportions
prop_A <- group_A_successes / group_A_total
prop_B <- group_B_successes / group_B_total

# Calculate the pooled proportion
pooled_prop <- (group_A_successes + group_B_successes) / (group_A_total + group_B_total)

# Calculate the standard error
se <- sqrt(pooled_prop * (1 - pooled_prop) * (1 / group_A_total + 1 / group_B_total))

# Calculate the Z-score
z_observed <- (prop_A - prop_B) / se

# One-tailed p-value (assuming prop_A > prop_B)
p_value <- 1 - pnorm(z_observed)
```

``` r
# Print the results
paste0("Number of Significant associations: ",group_A_successes," & ",group_B_successes)
```

    ## [1] "Number of Significant associations: 920 & 915"

``` r
paste0("Percentage of Significant associations: ",prop_A," & ",prop_B)
```

    ## [1] "Percentage of Significant associations: 0.178294573643411 & 0.131352282515073"

``` r
paste0("Observed Difference: ", prop_A - prop_B)
```

    ## [1] "Observed Difference: 0.0469422911283376"

``` r
paste0("Z-score:", z_observed)
```

    ## [1] "Z-score:7.13169607504534"

``` r
paste0("One-tailed p-value: ", p_value)
```

    ## [1] "One-tailed p-value: 4.95714580495132e-13"

##### Test5 : “Brain tissue gene sets showed a higher proportion of significant associations with complex traits than NonBrain”

``` r
################ Brain Vs NonBrain ########################
#Fantom
FanAll=rbind(p_values_fdr_DEL_Fantom,p_values_fdr_DUP_Fantom)

Brain_tissues=c("olfactory_bulb", "frontal_lobe", "medial_frontal_gyrus", "insular_cortex", "paracentral_gyrus",
                "parietal_lobe", "temporal_cortex", "medial_temporal_gyrus", "postcentral_gyrus", "hippocampus",
                "occipital_lobe", "occipital_cortex", "occipital_pole", "amygdala", "caudate",
                "putamen", "substantia_nigra", "globus_pallidus", "thalamus", "nucleus_accumbens",
                "pituitary_gland", "locus_coeruleus", "corpus_callosum", "cerebellum", "pons",
                "medulla_oblongata", "spinal_cord",
                
                "olfactory_bulb1", "frontal_lobe1", "medial_frontal_gyrus1", "insular_cortex1", "paracentral_gyrus1",
                "parietal_lobe1", "temporal_cortex1", "medial_temporal_gyrus1", "postcentral_gyrus1", "hippocampus1",
                "occipital_lobe1", "occipital_cortex1", "occipital_pole1", "amygdala1", "caudate1",
                "putamen1", "substantia_nigra1", "globus_pallidus1", "thalamus1", "nucleus_accumbens1",
                "pituitary_gland1", "locus_coeruleus1", "corpus_callosum1", "cerebellum1", "pons1",
                "medulla_oblongata1", "spinal_cord1")
NonBrain_tissues=rownames(FanAll[!rownames(FanAll)%in%Brain_tissues,])


FanAll_Brain_tissue=subset(FanAll,rownames(FanAll)%in%Brain_tissues)
FanAll_NonBrain_tissue=subset(FanAll,rownames(FanAll)%in%NonBrain_tissues)


group_A_successes <- sum(FanAll_Brain_tissue<0.05)
group_A_total <- 27*43*2 
group_B_successes <- sum(FanAll_NonBrain_tissue<0.05)
group_B_total <- (33)*43*2 
## One-tailed p-value: 0.003099545 

# Calculate proportions
prop_A <- group_A_successes / group_A_total
prop_B <- group_B_successes / group_B_total

# Calculate the pooled proportion
pooled_prop <- (group_A_successes + group_B_successes) / (group_A_total + group_B_total)

# Calculate the standard error
se <- sqrt(pooled_prop * (1 - pooled_prop) * (1 / group_A_total + 1 / group_B_total))

# Calculate the Z-score
z_observed <- (prop_A - prop_B) / se

# One-tailed p-value (assuming prop_A > prop_B)
p_value <- 1 - pnorm(z_observed)
```

``` r
# Print the results
paste0("Number of Significant associations: ",group_A_successes," & ",group_B_successes)
```

    ## [1] "Number of Significant associations: 648 & 272"

``` r
paste0("Percentage of Significant associations: ",prop_A," & ",prop_B)
```

    ## [1] "Percentage of Significant associations: 0.27906976744186 & 0.0958421423537703"

``` r
paste0("Observed Difference: ", prop_A - prop_B)
```

    ## [1] "Observed Difference: 0.18322762508809"

``` r
paste0("Z-score:", z_observed)
```

    ## [1] "Z-score:17.1070949009609"

``` r
paste0("One-tailed p-value: ", p_value)
```

    ## [1] "One-tailed p-value: 0"

##### Test6 : “Brain celltypes gene sets showed a higher proportion of significant associations with complex traits than NonBrain”

``` r
################ Brain Vs NonBrain ########################
#HPA
HPAAll=rbind(p_values_fdr_DEL_HPA,p_values_fdr_DUP_HPA)

Brain_cells=c("inhibitory_neurons", "excitatory_neurons","astrocytes","oligodendrocytes","oligodendrocyte_precursor_cells","microglial_cells",
              'muller_glia_cells',
              "inhibitory_neurons1", "excitatory_neurons1","astrocytes1","oligodendrocytes1","oligodendrocyte_precursor_cells1","microglial_cells1",
              'muller_glia_cells1')
NonBrain_cells=rownames(HPAAll[!rownames(HPAAll)%in%Brain_cells,])
#dim(HPAAll_NonBrain_cell)

HPAAll_Brain_cell=subset(HPAAll,rownames(HPAAll)%in%Brain_cells)
HPAAll_NonBrain_cell=subset(HPAAll,rownames(HPAAll)%in%NonBrain_cells)


group_A_successes <- sum(HPAAll_Brain_cell<0.05)
group_A_total <- 7*43*2 
group_B_successes <- sum(HPAAll_NonBrain_cell<0.05)
group_B_total <- (81-7)*43*2 
## One-tailed p-value: 0.003099545 

# Calculate proportions
prop_A <- group_A_successes / group_A_total
prop_B <- group_B_successes / group_B_total

# Calculate the pooled proportion
pooled_prop <- (group_A_successes + group_B_successes) / (group_A_total + group_B_total)

# Calculate the standard error
se <- sqrt(pooled_prop * (1 - pooled_prop) * (1 / group_A_total + 1 / group_B_total))

# Calculate the Z-score
z_observed <- (prop_A - prop_B) / se

# One-tailed p-value (assuming prop_A > prop_B)
p_value <- 1 - pnorm(z_observed)
```

``` r
# Print the results
paste0("Number of Significant associations: ",group_A_successes," & ",group_B_successes)
```

    ## [1] "Number of Significant associations: 103 & 812"

``` r
paste0("Percentage of Significant associations: ",prop_A," & ",prop_B)
```

    ## [1] "Percentage of Significant associations: 0.17109634551495 & 0.127592708988058"

``` r
paste0("Observed Difference: ", prop_A - prop_B)
```

    ## [1] "Observed Difference: 0.0435036365268923"

``` r
paste0("Z-score:", z_observed)
```

    ## [1] "Z-score:3.02034466039955"

``` r
paste0("One-tailed p-value: ", p_value)
```

    ## [1] "One-tailed p-value: 0.00126243593465225"

##### Test7 : “Brain cell types gene sets showed higher associations with brain traits than non-brain traits”

``` r
################ Brain Vs NonBrain ########################
#HPA
#By Trait
Brain_Traits = c(
  "EA", "FI", "Symbol Digit Substitution", "Trail Making 1", "Trail Making 2", "Reaction",
  "VSWM1", "VSWM2", "TDI", "Risk Taking", "Depression (self-report)", "GP Visit for Anxiety or Depression",
  "Neuroticism", "Irritability", "Mood Swings", "Miserableness", "Guilty Feeling", "Loneliness",'MoodAnxiety'
)
HPAAll_Brain_cell_BrainTrait=subset(HPAAll,rownames(HPAAll)%in%Brain_cells , colnames(HPAAll)%in%Brain_Traits)
#dim(HPAAll_Brain_cell_BrainTrait)
HPAAll_Brain_cell_NBrainTrait=subset(HPAAll,rownames(HPAAll)%in%Brain_cells , !colnames(HPAAll)%in%Brain_Traits)
#

group_A_successes <- sum(HPAAll_Brain_cell_BrainTrait<0.05)
group_A_total <- 7*19*2 
group_B_successes <- sum(HPAAll_Brain_cell_NBrainTrait<0.05)
group_B_total <- 7*(43-19)*2 
## One-tailed p-value: 0.003099545 

# Calculate proportions
prop_A <- group_A_successes / group_A_total
prop_B <- group_B_successes / group_B_total

# Calculate the pooled proportion
pooled_prop <- (group_A_successes + group_B_successes) / (group_A_total + group_B_total)

# Calculate the standard error
se <- sqrt(pooled_prop * (1 - pooled_prop) * (1 / group_A_total + 1 / group_B_total))

# Calculate the Z-score
z_observed <- (prop_A - prop_B) / se

# One-tailed p-value (assuming prop_A > prop_B)
p_value <- 1 - pnorm(z_observed)
```

``` r
# Print the results
paste0("Number of Significant associations: ",group_A_successes," & ",group_B_successes)
```

    ## [1] "Number of Significant associations: 59 & 44"

``` r
paste0("Percentage of Significant associations: ",prop_A," & ",prop_B)
```

    ## [1] "Percentage of Significant associations: 0.221804511278195 & 0.130952380952381"

``` r
paste0("Observed Difference: ", prop_A - prop_B)
```

    ## [1] "Observed Difference: 0.0908521303258145"

``` r
paste0("Z-score:", z_observed)
```

    ## [1] "Z-score:2.93950993725208"

``` r
paste0("One-tailed p-value: ", p_value)
```

    ## [1] "One-tailed p-value: 0.00164365864558735"

##### Test8 : “NonBrain cell types gene sets showed higher associations with Nonbrain traits than brain traits”

``` r
################ Brain Vs NonBrain ########################
#HPA
#By Trait
HPAAll_NonBrain_cell_BrainTrait=subset(HPAAll,rownames(HPAAll)%in%NonBrain_cells , colnames(HPAAll)%in%Brain_Traits)
#dim(HPAAll_Brain_cell_BrainTrait)
HPAAll_NonBrain_cell_NBrainTrait=subset(HPAAll,rownames(HPAAll)%in%NonBrain_cells , !colnames(HPAAll)%in%Brain_Traits)

group_A_successes <- sum(HPAAll_NonBrain_cell_NBrainTrait<0.05)
group_A_total <- (81-7)*(43-19)*2 
group_B_successes <- sum(HPAAll_NonBrain_cell_BrainTrait<0.05)
group_B_total <- (81-7)*(19)*2 
## One-tailed p-value: 0.003099545 

# Calculate proportions
prop_A <- group_A_successes / group_A_total
prop_B <- group_B_successes / group_B_total

# Calculate the pooled proportion
pooled_prop <- (group_A_successes + group_B_successes) / (group_A_total + group_B_total)

# Calculate the standard error
se <- sqrt(pooled_prop * (1 - pooled_prop) * (1 / group_A_total + 1 / group_B_total))

# Calculate the Z-score
z_observed <- (prop_A - prop_B) / se

# One-tailed p-value (assuming prop_A > prop_B)
p_value <- 1 - pnorm(z_observed)
```

``` r
# Print the results
paste0("Number of Significant associations: ",group_A_successes," & ",group_B_successes)
```

    ## [1] "Number of Significant associations: 493 & 319"

``` r
paste0("Percentage of Significant associations: ",prop_A," & ",prop_B)
```

    ## [1] "Percentage of Significant associations: 0.138795045045045 & 0.113442389758179"

``` r
paste0("Observed Difference: ", prop_A - prop_B)
```

    ## [1] "Observed Difference: 0.0253526552868658"

``` r
paste0("Z-score:", z_observed)
```

    ## [1] "Z-score:3.01043968453622"

``` r
paste0("One-tailed p-value: ", p_value)
```

    ## [1] "One-tailed p-value: 0.00130434876917307"

##### Test9 : “For Neuronal gene sets, duplications had a higher proportion of associations with brain-related traits”

``` r
################ Brain Vs NonBrain ########################
#Siletti
#Neuornal: 
Neuornal=c("amygdala_excitatory", "cerebellar_inhibitory", "cge_interneuron",
           "deep_layer_corticothalamic_and_6b", "deep_layer_intratelencephalic", "deep_layer_near_projecting",
           "eccentric_medium_spiny_neuron", "hippocampal_ca1_3", "hippocampal_ca4",
           "hippocampal_dentate_gyrus", "lamp5_lhx6_and_chandelier", "lower_rhombic_lip",
           "mammillary_body", "medium_spiny_neuron", "mge_interneuron",
           "midbrain_derived_inhibitory", "miscellaneous", "splatter",
           "thalamic_excitatory", "upper_layer_intratelencephalic", "upper_rhombic_lip"
)
NonNeuron=rownames(p_values_fdr_DEL_Siletti[!rownames(p_values_fdr_DEL_Siletti)%in%Neuornal,])

#Brain_Traits

#############
DUP_Nuron_BrainTrait= subset(p_values_fdr_DUP_Siletti,rownames(p_values_fdr_DUP_Siletti)%in%Neuornal , colnames(p_values_fdr_DUP_Siletti)%in%Brain_Traits)
  
DUP_Nuron_NonBrainTrait=subset(p_values_fdr_DUP_Siletti,rownames(p_values_fdr_DUP_Siletti)%in%Neuornal , !colnames(p_values_fdr_DUP_Siletti)%in%Brain_Traits)
##

group_A_successes <- sum(DUP_Nuron_BrainTrait<0.05)
group_A_total <- (21)*19 
group_B_successes <- sum(DUP_Nuron_NonBrainTrait<0.05)
group_B_total <- (21)*(43-19)
## One-tailed p-value: 0.003099545 

# Calculate proportions
prop_A <- group_A_successes / group_A_total
prop_B <- group_B_successes / group_B_total

# Calculate the pooled proportion
pooled_prop <- (group_A_successes + group_B_successes) / (group_A_total + group_B_total)

# Calculate the standard error
se <- sqrt(pooled_prop * (1 - pooled_prop) * (1 / group_A_total + 1 / group_B_total))

# Calculate the Z-score
z_observed <- (prop_A - prop_B) / se

# One-tailed p-value (assuming prop_A > prop_B)
p_value <- 1 - pnorm(z_observed)
```

``` r
# Print the results
paste0("Number of Significant associations: ",group_A_successes," & ",group_B_successes)
```

    ## [1] "Number of Significant associations: 93 & 45"

``` r
paste0("Percentage of Significant associations: ",prop_A," & ",prop_B)
```

    ## [1] "Percentage of Significant associations: 0.233082706766917 & 0.0892857142857143"

``` r
paste0("Observed Difference: ", prop_A - prop_B)
```

    ## [1] "Observed Difference: 0.143796992481203"

``` r
paste0("Z-score:", z_observed)
```

    ## [1] "Z-score:5.96382832447019"

``` r
paste0("One-tailed p-value: ", p_value)
```

    ## [1] "One-tailed p-value: 1.23197807511133e-09"

##### Test10 : “For Neuronal gene sets, deletion had a higher proportion of associations with Nonbrain-related traits”

``` r
################ Brain Vs NonBrain ########################
#Siletti
#Neuornal: 
DEL_Nuron_BrainTrait= subset(p_values_fdr_DEL_Siletti,rownames(p_values_fdr_DEL_Siletti)%in%Neuornal , colnames(p_values_fdr_DEL_Siletti)%in%Brain_Traits)

DEL_Nuron_NonBrainTrait=subset(p_values_fdr_DEL_Siletti,rownames(p_values_fdr_DEL_Siletti)%in%Neuornal , !colnames(p_values_fdr_DEL_Siletti)%in%Brain_Traits)
##

group_A_successes <- sum(DEL_Nuron_NonBrainTrait<0.05)
group_A_total <- (21)*(43-19)
group_B_successes <- sum(DEL_Nuron_BrainTrait<0.05)
group_B_total <- (21)*19
## One-tailed p-value: 0.003099545 

# Calculate proportions
prop_A <- group_A_successes / group_A_total
prop_B <- group_B_successes / group_B_total

# Calculate the pooled proportion
pooled_prop <- (group_A_successes + group_B_successes) / (group_A_total + group_B_total)

# Calculate the standard error
se <- sqrt(pooled_prop * (1 - pooled_prop) * (1 / group_A_total + 1 / group_B_total))

# Calculate the Z-score
z_observed <- (prop_A - prop_B) / se

# One-tailed p-value (assuming prop_A > prop_B)
p_value <- 1 - pnorm(z_observed)
```

``` r
# Print the results
paste0("Number of Significant associations: ",group_A_successes," & ",group_B_successes)
```

    ## [1] "Number of Significant associations: 100 & 30"

``` r
paste0("Percentage of Significant associations: ",prop_A," & ",prop_B)
```

    ## [1] "Percentage of Significant associations: 0.198412698412698 & 0.075187969924812"

``` r
paste0("Observed Difference: ", prop_A - prop_B)
```

    ## [1] "Observed Difference: 0.123224728487886"

``` r
paste0("Z-score:", z_observed)
```

    ## [1] "Z-score:5.23819953914587"

``` r
paste0("One-tailed p-value: ", p_value)
```

    ## [1] "One-tailed p-value: 8.10753474356218e-08"

``` r
####### Microgilia, Astrocytes and oligo finding are concordant using two different resources ###############
DEL_Sill_4=DEL_Siletti_43_mat_ordered_all[c("astrocyte","microglia","oligodendrocyte","oligodendrocyte_precursor"),]
DEL_HPA_4=DEL_HPA_43_mat_ordered_all[c("astrocytes","microglial_cells","oligodendrocytes","oligodendrocyte_precursor_cells"),]
rownames(DEL_HPA_4)=rownames(DEL_Sill_4)

#Corr
## DEL
common_rows=rownames(DEL_HPA_4)
correlations_DEL <- sapply(common_rows, function(row) {
  cor(as.numeric(DEL_Sill_4[row, ]), as.numeric(DEL_HPA_4[row, ]))
})
############################
print(correlations_DEL)
```

    ##                 astrocyte                 microglia           oligodendrocyte 
    ##                 0.5952922                 0.7790018                 0.5756332 
    ## oligodendrocyte_precursor 
    ##                 0.3759430
