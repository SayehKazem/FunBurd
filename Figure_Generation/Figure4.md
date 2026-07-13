Fig4
================
Kuldeep Kumar

## Fig. 4: Dissecting pleiotropy, gene function, and genetic constraint.

#### – Figure legend –

**Legend:** (A) Correlation between the fraction of constraint genes
(LOEUF top-decile) within a gene set and the number of traits showing
significant associations with that gene set for deletions and
duplications. Each data point is a gene set. Y-axis: the percentage of
significantly associated traits for a variant (functional pleiotropy);
X-axis: the percentage of top decile LOEUF genes within the gene set.
(B) Box plots show the distribution of constraint gene percentages for
all Brain and Non-Brain gene sets (172). Each point represents a gene
set, and asterisks (\*) indicate statistically significant differences
between groups. (C) Shows the same information on correlations in panel
(A) for different constraint metrics. (D) Functional pleiotropy,
normalized by the fraction of genetic constraint at the tissue level, is
shown for 27 brain and 33 non-brain gene sets – deletions are presented
on the left and duplications on the right. X-axis: represents the
proportion of intolerant genes for the different gene sets. Y-axis:
functional pleiotropy normalized by genetic constraint (centile). I.e.,
the 50th centile shows median functional pleiotropy computed across 100
randomly sampled gene sets. Circles and triangles represent brain and
non-brain gene sets, respectively. Gray shaded ribbons indicating
25th–75th (ribbon 1), 10th–90th (ribbon 2), and 5th–95th (ribbon 3)
centiles. This is followed by violin plots showing the distribution of
normalized functional pleiotropy across brain and non-brain traits.
Orange and green asterisks demonstrate significantly (FDR-corrected, q
\< 0.05) increased functional pleiotropy compared to what is expected
for a gene set with a comparable fraction of genetic constraint. The
grey star shows a significant difference between the brain and non-brain
gene sets, normalized functional pleiotropy. (E) The same analyses at
the whole-body cell type level, including 7 brain and 74 non-brain gene
sets.

#### Libraries

``` r
library(ggplot2)
library(ggprism)
library(ggpubr)
library(dplyr)
library(tidyr)
library(patchwork)
library(knitr)
library(openxlsx)
library(svglite)

# Define Core Colors and Sizes
brain_colorcode  <- "#d95f02"      
nonbrain_colorcode <- "#66a61e"  
in_base_size_ggprism <- 18
```

#### Panel A: Correlation between constraint and Del/Dup/GWAS functional pleiotropy

Description: Scatter plots evaluating the relationship between the
fraction of constraint genes (LOEUF) and the percentage of significantly
associated traits. \####

##### Panel A Data

``` r
##### Panel A Data #####
# 1. Load the Pre-computed Correlation & P-Jaccard Statistics
load("df_constraint_pleiotropy_corr.RData") 

stats_panel_a <- df_constraint_pleiotropy_corr %>%
  # Filter for the specific constraint measure used in Panel A
  filter(constraint == "measure_LEOUF_topdecile") %>%
  # Extract BOTH the correlation and the adjusted Jaccard p-value
  select(stat, cor, adj_pvalue_Jaccard)

# 2. Load the Geneset-level plotting coordinates
load("df_geneset_level_stats.RData") 

df_fp <- df_geneset_level_stats %>%
  mutate(
    # Convert constraint to percentage for the X-axis
    constraint_pct = constraint_measure_LEOUF_topdecile * 100
  ) %>%
  select(Geneset, Geneset_Type, Geneset_Cat, constraint_pct, 
         functional_pleiotropy_deletion, functional_pleiotropy_duplication, functional_pleiotropy_GWASenrichment)

kable(head(df_fp), format = "markdown")
```

| Geneset               | Geneset_Type    | Geneset_Cat | constraint_pct | functional_pleiotropy_deletion | functional_pleiotropy_duplication | functional_pleiotropy_GWASenrichment |
|:----------------------|:----------------|:------------|---------------:|-------------------------------:|----------------------------------:|-------------------------------------:|
| adipocytes            | SC_HPA81        | NonBrain    |      10.185934 |                       9.302326 |                         16.279070 |                             0.000000 |
| adipose_tissue        | Tissue_Fantom60 | NonBrain    |       7.137491 |                       4.651163 |                          4.651163 |                             0.000000 |
| alveolar_cells_type_1 | SC_HPA81        | NonBrain    |      11.168831 |                      34.883721 |                         11.627907 |                             5.263158 |
| alveolar_cells_type_2 | SC_HPA81        | NonBrain    |      10.783488 |                      13.953488 |                          6.976744 |                             0.000000 |
| amygdala              | Tissue_Fantom60 | Brain       |      17.880317 |                      30.232558 |                         13.953488 |                            18.421053 |
| amygdala_excitatory   | SC_Siletti31    | Brain       |      13.460076 |                      11.627907 |                         18.604651 |                            47.368421 |

``` r
kable(stats_panel_a, align = "c", caption = "Exact Stats Extracted for Panel A")
```

|                   stat                   |    cor    | adj_pvalue_Jaccard |
|:----------------------------------------:|:---------:|:------------------:|
|  functional_burden_pleiotropy_deletion   | 0.3898868 |     0.0229091      |
| functional_burden_pleiotropy_duplication | 0.4384765 |     0.0026250      |
|   functional_pleiotropy_GWASenrichment   | 0.5182434 |     0.0420000      |

Exact Stats Extracted for Panel A

##### Panel A Figure

``` r
# Function to plot scatter and map the exact 'cor' and 'p-jaccard' from the stats dataframe
plot_pleio_scatter_jaccard <- function(data, y_var, color_hex, stat_label, stats_df) {
  
  # Map BOTH correlation and p-value from the stats dataframe
  row_match <- stats_df %>% filter(stat == stat_label)
  cor_val <- row_match$cor
  pval_adj <- row_match$adj_pvalue_Jaccard
  
  # Format the label (add asterisk if FDR < 0.05)
  cor_label <- paste0("r=", round(cor_val, 2), ifelse(pval_adj < 0.05, "*", ""))
  
  ggplot(data, aes(x = constraint_pct, y = .data[[y_var]])) +
    geom_point(color = color_hex, size = 2, alpha = 1) + 
    # Add the extracted correlation label to the plot
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, label = cor_label, size = 6, fontface="bold") +
    geom_smooth(method = "lm", formula=y ~ x, se = FALSE, linewidth = 1.5, fullrange = TRUE, color = "black") + 
    scale_y_continuous(limits = c(0, 75)) +
    theme_prism(axis_text_angle = 90, base_size = 14) +
    labs(x = "% constraint genes", y = "% functional pleiotropy") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

# Generate the three plots using their exact stat labels
p_scatter_del <- plot_pleio_scatter_jaccard(
  df_fp, "functional_pleiotropy_deletion", "red3", 
  "functional_burden_pleiotropy_deletion", stats_panel_a
)

p_scatter_dup <- plot_pleio_scatter_jaccard(
  df_fp, "functional_pleiotropy_duplication", "royalblue", 
  "functional_burden_pleiotropy_duplication", stats_panel_a
) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

p_scatter_gwas <- plot_pleio_scatter_jaccard(
  df_fp, "functional_pleiotropy_GWASenrichment", "darkgreen", 
  "functional_pleiotropy_GWASenrichment", stats_panel_a
) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Layout the final figure
p_scatter_Fig4A <- (p_scatter_del + p_scatter_dup + p_scatter_gwas) + plot_layout(widths = c(1, 1, 1))

print(p_scatter_Fig4A)
```

    ## Warning: Removed 7 rows containing missing values or values outside the scale range
    ## (`geom_smooth()`).

![](Fig4_June2026_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

##### Panel A Stats

``` r
kable(head(stats_panel_a), format = "markdown")
```

| stat                                     |       cor | adj_pvalue_Jaccard |
|:-----------------------------------------|----------:|-------------------:|
| functional_burden_pleiotropy_deletion    | 0.3898868 |          0.0229091 |
| functional_burden_pleiotropy_duplication | 0.4384765 |          0.0026250 |
| functional_pleiotropy_GWASenrichment     | 0.5182434 |          0.0420000 |

#### Panel B: Distribution of constraint across Brain and Non-Brain gene sets

Description: Box plots comparing the percentage of constraint genes
between functional categories.

##### Panel B Data

``` r
df_box <- df_fp %>%
  select(Geneset, Geneset_Cat, constraint_pct) %>%
  mutate(Geneset_Cat = factor(Geneset_Cat, levels = c("Brain", "NonBrain")))

kable(head(df_box), format = "markdown")
```

| Geneset               | Geneset_Cat | constraint_pct |
|:----------------------|:------------|---------------:|
| adipocytes            | NonBrain    |      10.185934 |
| adipose_tissue        | NonBrain    |       7.137491 |
| alveolar_cells_type_1 | NonBrain    |      11.168831 |
| alveolar_cells_type_2 | NonBrain    |      10.783488 |
| amygdala              | Brain       |      17.880317 |
| amygdala_excitatory   | Brain       |      13.460076 |

##### Panel B Figure

``` r
p_box_Fig4B <- ggplot(df_box, aes(x = Geneset_Cat, y = constraint_pct, fill = Geneset_Cat)) +
  # Boxplot with black outline
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.5, linewidth = 0.8, color = "black") + 
  # Jitter points: shape 21 enables 'fill' aesthetic
  geom_jitter(aes(color=Geneset_Cat,shape = Geneset_Cat), stroke = 0.4, 
              width = 0.15, size = 2, alpha = 0.9) + 
  scale_fill_manual(values = c("Brain" = brain_colorcode, "NonBrain" = nonbrain_colorcode)) +
  scale_color_manual(values = c("Brain" = brain_colorcode, "NonBrain" = nonbrain_colorcode)) +
  scale_shape_manual(values = c("Brain" = 16, "NonBrain" = 17)) +
  theme_prism(base_size = in_base_size_ggprism) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold")
  ) +
  labs(x = NULL, y = "% constraint genes") +
  coord_flip()

print(p_box_Fig4B)
```

![](Fig4_June2026_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

##### Panel B Stats

``` r
res_wilcox_B <- wilcox.test(
  df_box %>% filter(Geneset_Cat == "Brain") %>% pull(constraint_pct),
  df_box %>% filter(Geneset_Cat == "NonBrain") %>% pull(constraint_pct)
)

stats_panel_b <- data.frame(
  Test = "Wilcoxon Rank Sum Test",
  Comparison = "Brain vs Non-Brain Constraint (%)",
  W_Statistic = res_wilcox_B$statistic,
  P_Value = format(res_wilcox_B$p.value, digits = 4)
)

kable(stats_panel_b, align = "c", caption = "Statistical Comparison of Constraint Percentages")
```

|     |          Test          |            Comparison             | W_Statistic |  P_Value  |
|:----|:----------------------:|:---------------------------------:|:-----------:|:---------:|
| W   | Wilcoxon Rank Sum Test | Brain vs Non-Brain Constraint (%) |    5673     | 4.152e-12 |

Statistical Comparison of Constraint Percentages

#### Panel C: Constraint vs Functional pleiotropy across metrics

Description: Evaluates functional pleiotropy correlations across seven
different established metrics of genetic constraint.

##### Panel C Data & Stats

``` r
load("df_constraint_pleiotropy_corr.RData")

# Map the pre-computed dataframe directly to plotting coordinates
df_metrics_cor <- df_constraint_pleiotropy_corr %>%
  mutate(
    # Map the raw constraint strings to clean metric labels
    Metric = case_when(
      grepl("LEOUF", constraint) ~ "LOEUF",
      grepl("missense_Z", constraint) ~ "missense Z",
      grepl("s_het", constraint) ~ "s-het",
      grepl("CDS", constraint) ~ "constraint CDS",
      grepl("pHaplo", constraint) ~ "pHaplo",
      grepl("pTriplo", constraint) ~ "pTriplo",
      grepl("GeneLenBP", constraint) ~ "gene length",
      TRUE ~ constraint
    ),
    # Map the stat column to Del, Dup, GWAS
    Type = case_when(
      grepl("deletion", stat) ~ "Del",
      grepl("duplication", stat) ~ "Dup",
      grepl("GWASenrichment", stat) ~ "GWAS",
      TRUE ~ stat
    ),
    # Assign significance using the Jaccard adjusted p-value directly
    sig = ifelse(adj_pvalue_Jaccard < 0.05, "q<0.05", "n.s.")
  ) %>%
  # Set factors for proper plot ordering
  mutate(
    Type = factor(Type, levels = c("Del", "Dup", "GWAS")),
    sig = factor(sig, levels = c("q<0.05", "n.s."))
  ) %>%
  select(Metric, Type, cor, pvalue_Jaccard, adj_pvalue_Jaccard, sig)

# Order Y-axis metrics based on Deletion correlation strength
metric_order <- df_metrics_cor %>% filter(Type == "Del") %>% arrange(cor) %>% pull(Metric)
df_metrics_cor$Metric <- factor(df_metrics_cor$Metric, levels = metric_order)

kable(head(df_metrics_cor), format = "markdown")
```

| Metric         | Type |       cor | pvalue_Jaccard | adj_pvalue_Jaccard | sig     |
|:---------------|:-----|----------:|---------------:|-------------------:|:--------|
| LOEUF          | Del  | 0.3898868 |          0.012 |          0.0229091 | q\<0.05 |
| missense Z     | Del  | 0.5020228 |          0.001 |          0.0026250 | q\<0.05 |
| s-het          | Del  | 0.2636889 |          0.112 |          0.1470000 | n.s.    |
| constraint CDS | Del  | 0.4633309 |          0.001 |          0.0026250 | q\<0.05 |
| pHaplo         | Del  | 0.1100412 |          0.321 |          0.3210000 | n.s.    |
| pTriplo        | Del  | 0.4063301 |          0.008 |          0.0168000 | q\<0.05 |

##### Panel C Figure

``` r
p_Fig4C <- ggplot(df_metrics_cor, aes(y = Metric, x = cor, group = Type)) +
  geom_point(aes(color = Type, fill = Type, shape = sig), size = 5, position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = c("Del" = "red3", "Dup" = "royalblue", "GWAS" = "darkgreen")) +
  scale_color_manual(values = c("Del" = "red3", "Dup" = "royalblue", "GWAS" = "darkgreen")) +
  scale_shape_manual(values = c("q<0.05" = 21, "n.s." = 4)) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "black") +
  geom_vline(xintercept = 0.4, linetype = "dotted", color = "black") +
  theme_prism(base_size = in_base_size_ggprism) +
  theme(legend.position = c(0.85, 0.2)) +
  labs(x = "Correlation", y = NULL) +
  coord_cartesian(xlim = c(0, 0.7))

print(p_Fig4C)
```

![](Fig4_June2026_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

#### Panel D & E: Normative Null Models (Tissue & Cell Type)

Description: Functional pleiotropy normalized by the fraction of genetic
constraint, plotted against null distributions for Tissue (Fantom60) and
Cell Type (HPA81).

##### Panel D & E Data

``` r
# 1. Extract and format the exact plotting coordinates from geneset level stats
# No external null distribution data is required.
df_normative <- df_geneset_level_stats %>%
  mutate(
    # Convert constraint fraction to percentage for the X-axis
    constraint_pct = constraint_measure_LEOUF_topdecile * 100
  ) %>%
  select(Geneset, Geneset_Type, Geneset_Cat, constraint_pct, 
         normative_model_centile_del_pleiotropy, normative_model_centile_dup_pleiotropy)

# Ensure Geneset_Cat is a factor for consistent color mapping
df_normative$Geneset_Cat <- factor(df_normative$Geneset_Cat, levels = c("Brain", "NonBrain"))

kable(head(df_normative), format = "markdown")
```

| Geneset               | Geneset_Type    | Geneset_Cat | constraint_pct | normative_model_centile_del_pleiotropy | normative_model_centile_dup_pleiotropy |
|:----------------------|:----------------|:------------|---------------:|---------------------------------------:|---------------------------------------:|
| adipocytes            | SC_HPA81        | NonBrain    |      10.185934 |                                     54 |                                     86 |
| adipose_tissue        | Tissue_Fantom60 | NonBrain    |       7.137491 |                                     33 |                                     31 |
| alveolar_cells_type_1 | SC_HPA81        | NonBrain    |      11.168831 |                                     94 |                                     71 |
| alveolar_cells_type_2 | SC_HPA81        | NonBrain    |      10.783488 |                                     68 |                                     46 |
| amygdala              | Tissue_Fantom60 | Brain       |      17.880317 |                                     90 |                                     79 |
| amygdala_excitatory   | SC_Siletti31    | Brain       |      13.460076 |                                     62 |                                     91 |

##### Panel D & E Figure

``` r
# --- 1. Compute Stats First (Needed for Dynamic Annotations) ---
compute_normative_stats <- function(data, gset_type, metric_col, variant_label) {
  sub_df <- data %>% filter(Geneset_Type == gset_type)
  
  brain_vals <- sub_df %>% filter(Geneset_Cat == "Brain") %>% pull(.data[[metric_col]])
  nbrain_vals <- sub_df %>% filter(Geneset_Cat == "NonBrain") %>% pull(.data[[metric_col]])
  
  # Unified Wilcoxon Rank-sum Tests
  p_brain_50 <- wilcox.test(brain_vals, mu = 50,exact = FALSE)$p.value
  p_nbrain_50 <- wilcox.test(nbrain_vals, mu = 50,exact = FALSE)$p.value
  p_b_vs_nb <- wilcox.test(brain_vals, nbrain_vals,exact = FALSE)$p.value
  
  return(data.frame(Geneset_Type = gset_type, Variant = variant_label,
                    P_Brain_vs_50 = p_brain_50, P_NonBrain_vs_50 = p_nbrain_50, P_Brain_vs_NonBrain = p_b_vs_nb))
}

# Run across all categories to get raw p-values
df_stats_raw <- bind_rows(
  compute_normative_stats(df_normative, "Tissue_Fantom60", "normative_model_centile_del_pleiotropy", "Del"),
  compute_normative_stats(df_normative, "Tissue_Fantom60", "normative_model_centile_dup_pleiotropy", "Dup"),
  compute_normative_stats(df_normative, "SC_HPA81", "normative_model_centile_del_pleiotropy", "Del"),
  compute_normative_stats(df_normative, "SC_HPA81", "normative_model_centile_dup_pleiotropy", "Dup")
)

# 1. Pool all 8 tests for 'shift from 50%' and adjust FDR
pooled_vs_50_pvals <- c(df_stats_raw$P_Brain_vs_50, df_stats_raw$P_NonBrain_vs_50)
pooled_vs_50_fdr <- p.adjust(pooled_vs_50_pvals, method = "fdr")

# 2. Pool all 4 tests for 'brain vs non-brain' and adjust FDR
pooled_b_vs_nb_fdr <- p.adjust(df_stats_raw$P_Brain_vs_NonBrain, method = "fdr")

# Map the correctly pooled FDR values back
df_stats_DE_final <- df_stats_raw %>%
  mutate(
    # Indices 1:4 belong to Brain, indices 5:8 belong to NonBrain
    FDR_Brain_vs_50 = pooled_vs_50_fdr[1:4],
    FDR_NonBrain_vs_50 = pooled_vs_50_fdr[5:8],
    FDR_Brain_vs_NonBrain = pooled_b_vs_nb_fdr
  )

# --- 2. Helper Function to Build Individual Components ---
build_normative_plots <- function(data_obs, gset_type, var_y, variant_type, stats_df, y_label = NULL, hide_y = FALSE) {
  
  df_sub <- data_obs %>% filter(Geneset_Type == gset_type)
  plot_stats <- stats_df %>% filter(Geneset_Type == gset_type, Variant == variant_type)
  
  # Scatter Plot (Using manual annotations for the Centile bands instead of a null dataset)
  p_scatter <- ggplot(df_sub, aes(x = constraint_pct, y = .data[[var_y]])) +
    # Background Centile Ribbons
    annotate("rect", ymin = 5, ymax = 95, xmin = -Inf, xmax = Inf, fill = "grey80", alpha = 0.4) +
    annotate("rect", ymin = 10, ymax = 90, xmin = -Inf, xmax = Inf, fill = "grey70", alpha = 0.4) +
    annotate("rect", ymin = 25, ymax = 75, xmin = -Inf, xmax = Inf, fill = "grey60", alpha = 0.4) +
    geom_hline(yintercept = 50, color = "black", linewidth = 1) +
    
    # Scatter Points
    geom_point(aes(color = Geneset_Cat, shape = Geneset_Cat), size = 3) +
    scale_color_manual(values = c("Brain" = brain_colorcode, "NonBrain" = nonbrain_colorcode)) +
    scale_shape_manual(values = c("Brain" = 16, "NonBrain" = 17)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 22), breaks = seq(0, 22, by = 5)) +
    scale_y_continuous(limits = c(0, 125), breaks = c(10, 25, 50, 75, 90)) + 
    theme_prism(base_size = in_base_size_ggprism, base_line_size = in_base_size_ggprism/24) +
    theme(legend.position = "none") +
    labs(x = "% constraint genes", y = y_label)

  if (hide_y) {
    p_scatter <- p_scatter + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  
  # Violin Plot
  p_violin <- ggplot(df_sub, aes(x = Geneset_Cat, y = .data[[var_y]], color = Geneset_Cat,shape=Geneset_Cat)) +
    geom_hline(yintercept = 50, linetype = "dotted", color = "black", linewidth = 2) +
    geom_violin(fill = "white", alpha = 1, trim = TRUE, linewidth = 0.8) +
    geom_point(position = position_jitter(seed = 1, width = 0.2), size = 1.5, alpha = 0.8) +
    stat_summary(color = "black", fun.data = "mean_cl_boot", geom = "pointrange", size = 0.7) +
    scale_color_manual(values = c("Brain" = brain_colorcode, "NonBrain" = nonbrain_colorcode)) +
     scale_shape_manual(values = c("Brain" = 16, "NonBrain" = 17)) +
    scale_x_discrete(labels = c("Brain", "NonBrain")) +
    scale_y_continuous(limits = c(0, 125), breaks = c(10, 25, 50, 75, 90)) +
    theme_prism(base_size = in_base_size_ggprism, base_line_size = in_base_size_ggprism/24) +
    theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

  # Add Significance Annotations dynamically based on the stats
  sig_size <- 10
  
  if (plot_stats$FDR_Brain_vs_50 < 0.05) {
    p_violin <- p_violin + annotate("text", x = 1, y = 105, label = "*", size = sig_size, color = brain_colorcode, fontface = "bold")
  }
  if (plot_stats$FDR_NonBrain_vs_50 < 0.05) {
    p_violin <- p_violin + annotate("text", x = 2, y = 105, label = "*", size = sig_size, color = nonbrain_colorcode, fontface = "bold")
  }
  if (plot_stats$FDR_Brain_vs_NonBrain < 0.05) {
    p_violin <- p_violin + 
      annotate("segment", x = 1, xend = 2, y = 115, yend = 115, color = "grey50", linewidth = 1) +
      annotate("segment", x = 1, xend = 1, y = 110, yend = 115, color = "grey50", linewidth = 1) +
      annotate("segment", x = 2, xend = 2, y = 110, yend = 115, color = "grey50", linewidth = 1) +
      annotate("text", x = 1.5, y = 118, label = "*", size = sig_size, color = "grey50", fontface = "bold")
  }

  return(p_scatter + p_violin + plot_layout(widths = c(0.72, 0.28)))
}

# --- 3. Assemble Final Composite Panels ---
y_label_text <- "functional pleiotropy\nnormalized by constraint (centiles)"

pD_Del <- build_normative_plots(df_normative, "Tissue_Fantom60", "normative_model_centile_del_pleiotropy", "Del", df_stats_DE_final, y_label = y_label_text)
pD_Dup <- build_normative_plots(df_normative, "Tissue_Fantom60", "normative_model_centile_dup_pleiotropy", "Dup", df_stats_DE_final, hide_y = TRUE)
row_D <- pD_Del | pD_Dup

pE_Del <- build_normative_plots(df_normative, "SC_HPA81", "normative_model_centile_del_pleiotropy", "Del", df_stats_DE_final, y_label = y_label_text)
pE_Dup <- build_normative_plots(df_normative, "SC_HPA81", "normative_model_centile_dup_pleiotropy", "Dup", df_stats_DE_final, hide_y = TRUE)
row_E <- pE_Del | pE_Dup

p_single_stacked <- row_D / row_E
print(p_single_stacked)
```

![](Fig4_June2026_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

##### Panel D & E Stats

``` r
# --- Simplified and Robust Stats Calculation ---
# Display the stats dataframe that was computed and utilized in the Figure chunk
kable(df_stats_DE_final %>% select(Geneset_Type, Variant, starts_with("FDR"), starts_with("P")), align = "c", caption = "FDR-Adjusted Wilcoxon Tests for Normative Centiles (Panels D & E)")
```

|  Geneset_Type   | Variant | FDR_Brain_vs_50 | FDR_NonBrain_vs_50 | FDR_Brain_vs_NonBrain | P_Brain_vs_50 | P_NonBrain_vs_50 | P_Brain_vs_NonBrain |
|:---------------:|:-------:|:---------------:|:------------------:|:---------------------:|:-------------:|:----------------:|:-------------------:|
| Tissue_Fantom60 |   Del   |    0.0000647    |     0.3421268      |       0.0000001       |   0.0000081   |    0.2138293     |      0.0000000      |
| Tissue_Fantom60 |   Dup   |    0.0414258    |     0.6104321      |       0.1062800       |   0.0111889   |    0.6104321     |      0.0531400      |
|    SC_HPA81     |   Del   |    0.3625316    |     0.4015900      |       0.2446723       |   0.2718987   |    0.3513912     |      0.1835043      |
|    SC_HPA81     |   Dup   |    0.0685759    |     0.0414258      |       0.2631414       |   0.0342880   |    0.0155347     |      0.2631414      |

FDR-Adjusted Wilcoxon Tests for Normative Centiles (Panels D & E)
