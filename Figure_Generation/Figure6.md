Fig6
================
Sayeh Kazem

## Fig. 6: Gene dosage responses across traits and functional gene sets.

#### – Figure legend –

**Legend:** (A) Proportion of gene dosage responses significant for
deletions-only, duplication-only, and both deletion-duplication, across
brain and non-brain traits. Asterisks (*) denote significant overlap.
(B) Brain and non-brain traits showed similar deletion-duplication
association overlap. (C) Proportions of monotonic (dark orange) and
non-monotonic (grey) responses across brain and non-brain traits. Total
responses include only true gene dosage responses (FDR-significant
effects for both deletions and duplications); proportions are calculated
within this subset. (D) Box plots depict the distribution of monotonic
responses for brain and non-brain traits; asterisks (*) indicate
statistically significant differences. (E) Gene dosage responses across
Educational Attainment (EA), Fluid Intelligence (FI), Body Mass Index
(BMI), and Platelet. Each line connects deletion (left) and duplication
(right) association effect sizes for a functional gene set. Dark orange
lines represent monotonic responses (significant effect sizes for both
dosages with opposing directionality). Grey lines represent
non-monotonic responses. (F-G) Summary of deletion-duplication effect
size correlations across traits.

#### Libraries

``` r
library(ggplot2)
library(ggprism)
library(ggpubr)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(knitr)
```

#### Panel A & B: Gene Dosage Response Functional Overlap

Description: Proportions of significant associations isolated to
Del-Only, Dup-Only, and Both, plotted alongside the overall Del-Dup
Overlap percentage distribution.

##### Panel A & B Data

``` r
# Load the pre-computed trait-level statistics
load('df_traitlevel_stats.RData')

# Define Custom Ordering matching the figure
custom_order <- c(
  # Brain Traits
  "Neuroticism", "MoodAnxiety", "FI", "Loneliness", "Reaction", "EA", "TDI",
  # Non-Brain Traits
  "Platelet", "Menarche", "Gamma", "IGF1", "HDL", "StandingHeight", "FEV1", "BMI", "FVC"
)

# Prepare Overlap Data
df_overlap <- df_traitlevel_stats %>%
  filter(Trait %in% custom_order) %>%
  mutate(
    Trait = factor(Trait, levels = custom_order),
    Trait_Cat = ifelse(Trait_Cat == "Brain", "Brain Traits", "Non-Brain Traits")
  ) %>%
  select(Trait, Trait_Cat, 
         Del_Dup = proportion_del_dup_bothsig, 
         Del_only = proportion_del_sig_only, 
         Dup_only = proportion_dup_sig_only) %>%
  pivot_longer(cols = c(Dup_only, Del_Dup, Del_only), names_to = "Category", values_to = "Percentage") %>%
  mutate(
    Percentage = Percentage * 100,
    Category = factor(Category, levels = c("Dup_only", "Del_Dup", "Del_only"))
  )

kable(head(df_overlap), format = "markdown")
```

| Trait | Trait_Cat        | Category | Percentage |
|:------|:-----------------|:---------|-----------:|
| BMI   | Non-Brain Traits | Dup_only |   18.66667 |
| BMI   | Non-Brain Traits | Del_Dup  |   21.33333 |
| BMI   | Non-Brain Traits | Del_only |   60.00000 |
| EA    | Brain Traits     | Dup_only |   35.10638 |
| EA    | Brain Traits     | Del_Dup  |   26.59574 |
| EA    | Brain Traits     | Del_only |   38.29787 |

##### Panel A Figure

``` r
p_overlap_bars <- ggplot(df_overlap, aes(fill = Category, y = Percentage, x = Trait)) + 
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap(~Trait_Cat, scales = "free_x") +
  scale_fill_manual(values = c("Del_only" = "red3", "Del_Dup" = "grey60", "Dup_only" = "royalblue")) +
  scale_y_continuous(limits = c(0, 101), expand = c(0, 0)) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.8), 
    axis.ticks = element_line(colour = "black", linewidth = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 14, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold")
  ) +
  labs(y = "Proportion of\nAssociations", x = NULL)

print(p_overlap_bars)
```

![](Fig6_June2026_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

##### Panel A Stats

``` r
#Panel A Stats: Overlap Significance per Trait
df_overlap_sig <- df_traitlevel_stats %>%
  filter(Trait %in% custom_order) %>%
  mutate(Trait_Cat = ifelse(Trait_Cat == "Brain", "Brain Traits", "Non-Brain Traits")) %>%
  select(Trait, Trait_Cat, adj_pvalue_count_del_dup_overlap) %>%
  arrange(Trait_Cat, adj_pvalue_count_del_dup_overlap)

kable(head(df_overlap_sig), align = "c", caption = "FDR-Adjusted P-Values for Del-Dup Overlap Significance")
```

|    Trait    |  Trait_Cat   | adj_pvalue_count_del_dup_overlap |
|:-----------:|:------------:|:--------------------------------:|
|     TDI     | Brain Traits |            0.0215000             |
|  Reaction   | Brain Traits |            0.1863333             |
|     EA      | Brain Traits |            0.4204444             |
| Neuroticism | Brain Traits |            0.8321765             |
| Loneliness  | Brain Traits |            0.9340556             |
|     FI      | Brain Traits |            1.0000000             |

FDR-Adjusted P-Values for Del-Dup Overlap Significance

##### Panel B Figure

``` r
df_overlap_box <- df_overlap %>%
  filter(Category == "Del_Dup")

p_overlap_box <- ggplot(df_overlap_box, aes(x = Trait_Cat, y = Percentage, fill = Trait_Cat)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 1, linewidth = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7, color = 'grey30') +
  scale_fill_manual(values = c("Brain Traits" = "white", "Non-Brain Traits" = "white")) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.8),
    axis.ticks = element_line(colour = "black", linewidth = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(y = "Del-Dup Overlap (%)") +
  coord_cartesian(ylim = c(0, 35), clip = "off")

print(p_overlap_box)
```

![](Fig6_June2026_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

##### Panel B Stats

``` r
# Panel B Stats: Wilcoxon Test (Brain vs Non-Brain Overlap Percentage)
overlap_brain <- df_overlap_box %>% filter(Trait_Cat == "Brain Traits") %>% pull(Percentage)
overlap_nonbrain <- df_overlap_box %>% filter(Trait_Cat == "Non-Brain Traits") %>% pull(Percentage)

wilcox_overlap <- wilcox.test(overlap_brain, overlap_nonbrain)

overlap_test_stats <- data.frame(
  Test = "Wilcoxon Rank Sum Test",
  Comparison = "Brain vs Non-Brain Del-Dup Overlap (%)",
  W_Statistic = wilcox_overlap$statistic,
  P_Value = format(wilcox_overlap$p.value, digits = 4)
)

kable(overlap_test_stats, align = "c", caption = "Statistical Comparison of Del-Dup Overlap Percentages")
```

|     |          Test          |               Comparison               | W_Statistic | P_Value |
|:----|:----------------------:|:--------------------------------------:|:-----------:|:-------:|
| W   | Wilcoxon Rank Sum Test | Brain vs Non-Brain Del-Dup Overlap (%) |     31      |    1    |

Statistical Comparison of Del-Dup Overlap Percentages

#### Panel C & D : True Monotonicity Proportions

##### Description: Evaluates the proportion of true gene dosage responses that are monotonic vs non-monotonic.

##### Panel C & D Data

``` r
df_mono <- df_traitlevel_stats %>%
  filter(Trait %in% custom_order) %>%
  # Filter to traits with sufficient true responses for robust statistics
  filter(true_gene_dosage_responses_del_and_dup_sig > 5) %>%
  mutate(Trait = factor(Trait, levels = custom_order)) %>%
  select(Trait, Trait_Cat, 
         rel_mono = true_gene_dosage_responses_proportion_monotonic,
         rel_ushape = true_gene_dosage_responses_proportion_nonmonotonic) %>%
  pivot_longer(cols = c(rel_mono, rel_ushape), names_to = "Type", values_to = "Prop") %>%
  mutate(
    StackGroup = ifelse(Type == "rel_mono", 1, 2),
    Trait_Cat = ifelse(Trait_Cat == "Brain", "Brain Traits", "Non-Brain Traits")
  )

kable(head(df_mono), format = "markdown")
```

| Trait | Trait_Cat        | Type       |      Prop | StackGroup |
|:------|:-----------------|:-----------|----------:|-----------:|
| BMI   | Non-Brain Traits | rel_mono   | 0.5000000 |          1 |
| BMI   | Non-Brain Traits | rel_ushape | 0.5000000 |          2 |
| EA    | Brain Traits     | rel_mono   | 0.0800000 |          1 |
| EA    | Brain Traits     | rel_ushape | 0.9200000 |          2 |
| FEV1  | Non-Brain Traits | rel_mono   | 0.1153846 |          1 |
| FEV1  | Non-Brain Traits | rel_ushape | 0.8846154 |          2 |

##### Panel C Figure

``` r
p_mono_bars <- ggplot(df_mono, aes(x = Trait, y = Prop, fill = interaction(Trait_Cat, Type), group = StackGroup)) +
  geom_col(color = "white", linewidth = 0.5, position = position_stack(reverse = TRUE)) + 
  facet_wrap(~Trait_Cat, scales = "free_x") +
  scale_fill_manual(values = c(
    "Brain Traits.rel_mono" = "darkorange2", 
    "Non-Brain Traits.rel_mono" = "darkorange2",
    "Brain Traits.rel_ushape" = "grey75", 
    "Non-Brain Traits.rel_ushape" = "grey75"
  )) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.8), 
    axis.ticks = element_line(colour = "black", linewidth = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size = 14,face = "bold", color = "black"),
    axis.text.y = element_text(size = 14, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none", 
    strip.background = element_blank(),
    strip.text = element_blank() # Hide to flush cleanly under Panel A
  ) +
  labs(y = "Proportion of\nMonotonic/non-monotonic", x = NULL)

print(p_mono_bars)
```

![](Fig6_June2026_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

##### Panel D Figure

``` r
df_mono_box <- df_traitlevel_stats %>%
  filter(true_gene_dosage_responses_del_and_dup_sig > 5) %>%
  mutate(Trait_Cat = ifelse(Trait_Cat == "Brain", "Brain Traits", "Non-Brain Traits"))

p_mono_box <- ggplot(df_mono_box, aes(x = Trait_Cat, y = true_gene_dosage_responses_proportion_monotonic, fill = Trait_Cat)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.8, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.8, size = 2, color = "darkorange2") +
  scale_fill_manual(values = c("Brain Traits" = "white", "Non-Brain Traits" = "white")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(y = "Monotonic Proportion")

print(p_mono_box)
```

    ## Warning: Removed 2 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](Fig6_June2026_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

##### Panel D Stats

``` r
# Panel D Stats: Wilcoxon Test (Brain vs Non-Brain Monotonicity Proportion)
# Using the exact subset() methodology:
wilcox_mono <- wilcox.test(
  subset(df_mono_box, Trait_Cat == 'Brain Traits')$true_gene_dosage_responses_proportion_monotonic,
  subset(df_mono_box, Trait_Cat == 'Non-Brain Traits')$true_gene_dosage_responses_proportion_monotonic
)
```

    ## Warning in wilcox.test.default(subset(df_mono_box, Trait_Cat == "Brain
    ## Traits")$true_gene_dosage_responses_proportion_monotonic, : cannot compute
    ## exact p-value with ties

``` r
mono_test_stats <- data.frame(
  Test = "Wilcoxon Rank Sum Test",
  Comparison = "Brain vs Non-Brain Monotonic Proportion",
  W_Statistic = wilcox_mono$statistic,
  P_Value = format(wilcox_mono$p.value, digits = 4)
)

kable(mono_test_stats, align = "c", caption = "Statistical Comparison of Monotonicity (Panel D)")
```

|     |          Test          |               Comparison                | W_Statistic | P_Value |
|:----|:----------------------:|:---------------------------------------:|:-----------:|:-------:|
| W   | Wilcoxon Rank Sum Test | Brain vs Non-Brain Monotonic Proportion |     8.5     | 0.01634 |

Statistical Comparison of Monotonicity (Panel D)

#### Panel E : Examples of True Gene Dosage Responses

##### Description: Illustration of true gene dosage responses across 2 brain (EA, FI) and 2 non-brain traits (BMI, Platelet). Evaluates only genesets significant for BOTH DEL and DUP.

#### Panel E Data

``` r
load('df_funburd_ukbb.RData')
target_traits_E <- c('EA', 'FI', 'BMI', 'Platelet')

# Pivot wider to evaluate DEL and DUP significance simultaneously 
df_true_responses <- df_funburd_ukbb %>%
  filter(Trait %in% target_traits_E) %>%
  select(Trait, Trait_Cat, Geneset, CNV_Type, Effectsize, FDR_p_value) %>%
  pivot_wider(
    names_from = CNV_Type, 
    values_from = c(Effectsize, FDR_p_value),
    names_glue = "{.value}_{CNV_Type}"
  ) %>%
  filter(FDR_p_value_DEL < 0.05 & FDR_p_value_DUP < 0.05) %>%
  mutate(
    is_monotonic = ifelse(sign(Effectsize_DEL) != sign(Effectsize_DUP), "Monotonic", "Non-Monotonic"),
    concordance = ifelse(is_monotonic == "Monotonic", "darkorange2", "grey")
  )

# Pivot longer for plotting lines between DEL and DUP
df_lines <- df_true_responses %>%
  pivot_longer(
    cols = c(Effectsize_DEL, Effectsize_DUP),
    names_to = "Dosage_Type",
    values_to = "Effect_Size"
  ) %>%
  mutate(
    Dosage_Type = ifelse(Dosage_Type == "Effectsize_DEL", "Del", "Dup"),
    concordance = factor(concordance, levels = c("darkorange2", "grey")),
    Trait = factor(Trait, levels = target_traits_E)
  ) %>%
  arrange(desc(concordance))

kable(head(df_lines), format = "markdown")
```

| Trait | Trait_Cat | Geneset               | FDR_p_value_DUP | FDR_p_value_DEL | is_monotonic  | concordance | Dosage_Type | Effect_Size |
|:------|:----------|:----------------------|----------------:|----------------:|:--------------|:------------|:------------|------------:|
| EA    | Brain     | alveolar_cells_type_1 |       0.0001770 |       0.0005808 | Non-Monotonic | grey        | Del         |  -0.0493210 |
| EA    | Brain     | alveolar_cells_type_1 |       0.0001770 |       0.0005808 | Non-Monotonic | grey        | Dup         |  -0.0391073 |
| FI    | Brain     | alveolar_cells_type_1 |       0.0450629 |       0.0001088 | Non-Monotonic | grey        | Del         |  -0.0702922 |
| FI    | Brain     | alveolar_cells_type_1 |       0.0450629 |       0.0001088 | Non-Monotonic | grey        | Dup         |  -0.0301889 |
| EA    | Brain     | appendix              |       0.0079526 |       0.0237238 | Non-Monotonic | grey        | Del         |  -0.0337634 |
| EA    | Brain     | appendix              |       0.0079526 |       0.0237238 | Non-Monotonic | grey        | Dup         |  -0.0232749 |

#### Panel E Figure

``` r
p_lines <- ggplot(df_lines, aes(x = Dosage_Type, y = Effect_Size, group = Geneset, color = concordance)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_line(alpha = 0.4, linewidth = 0.8) +
  facet_wrap(~Trait, scales = "fixed", nrow = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(-0.15, 0.15)) +
  scale_color_manual(values = c("darkorange2" = "darkorange2", "grey" = "grey80")) +
  theme_minimal() +
  labs(y = "Effect-Sizes", x = NULL) +
  theme(
    legend.position = "none",
    text = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "bold", size = 12),
    strip.text = element_text(face = "bold", size = 14)
  )

print(p_lines)
```

![](Fig6_June2026_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

#### Panel F & G : Del-Dup Effect Size Correlations

##### Description: Distribution of effect size correlations between deletions and duplications across brain and non-brain traits

#### Panel F & G Data

``` r
corr_brain_traits <- c("RiskTaking", "Depression", "SymbolDigit", "TrailMaking2", 
                       "Neuroticism", "VSWM2", "MoodAnxiety", "FI", "Miserableness", 
                       "Loneliness", "MoodSwings", "VSWM1", "Reaction", "EA", "TDI")

corr_nbrain_traits <- c("BirthWeight", "Platelet", "Menarche", "HbA1c", "Cholesterol", 
                        "Triglycerides", "Albumin", "SexPartners", "Gamma", "IGF1", 
                        "HDL", "StandingHeight", "FEV1", "BMI", "FVC", "HeartRate", "Redblood", "Menopause")

# Filter and format correlation data
df_corr <- df_traitlevel_stats %>%
  filter(Trait %in% c(corr_brain_traits, corr_nbrain_traits), !is.na(cor_del_dup)) %>%
  mutate(
    Trait_Cat = ifelse(Trait %in% corr_brain_traits, "Brain Traits", "Non-Brain Traits")
  ) %>%
  # Arrange traits by correlation magnitude to match the plot order
  arrange(Trait_Cat, cor_del_dup) %>%
  mutate(Trait = factor(Trait, levels = unique(Trait)))

# Calculate specific means for Brain and Non-Brain to plot the dashed lines per facet
df_corr_means <- df_corr %>%
  group_by(Trait_Cat) %>%
  summarize(mean_cor = mean(cor_del_dup, na.rm = TRUE))

kable(head(df_corr %>% select(Trait, Trait_Cat, cor_del_dup, adj_JaccardPvalue_cor_del_dup)), format = "markdown")
```

| Trait        | Trait_Cat    | cor_del_dup | adj_JaccardPvalue_cor_del_dup |
|:-------------|:-------------|------------:|------------------------------:|
| RiskTaking   | Brain Traits |  -0.7892288 |                     0.0000000 |
| Depression   | Brain Traits |  -0.7705149 |                     0.0000000 |
| SymbolDigit  | Brain Traits |  -0.7232680 |                     0.0000000 |
| TrailMaking2 | Brain Traits |  -0.4618233 |                     0.0226667 |
| Neuroticism  | Brain Traits |  -0.4475877 |                     0.0370909 |
| VSWM2        | Brain Traits |  -0.4059057 |                     0.0575385 |

#### Panel F Figure

``` r
p_corr_bars <- ggplot(df_corr, aes(x = Trait, y = cor_del_dup, fill = Trait_Cat)) +
  geom_col(color = "black", linewidth = 0.5) +
  facet_wrap(~Trait_Cat, scales = "free_x") +
  scale_fill_manual(values = c("Brain Traits" = "purple4", "Non-Brain Traits" = "purple4")) +  
  # Category-specific dashed mean lines
  geom_hline(data = df_corr_means, aes(yintercept = mean_cor), linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1) + 
  # Significance asterisks placed just below the 0 line inside the negative bars
  geom_text(aes(y = -0.05, label = ifelse(adj_JaccardPvalue_cor_del_dup < 0.05, "*", "")), 
            color = "white", size = 15, fontface = "bold", vjust = 1) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-0.8, 0.45)) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.8), 
    axis.ticks = element_line(colour = "black", linewidth = 0.8),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16, face = "bold", color = "black"),
    axis.text.y = element_text(size = 16, face = "bold", color = "black"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 18, face = "bold"), 
    panel.spacing = unit(2.5, "lines") # <-- This creates the space between Brain and Non-Brain facets
  ) +
  labs(x = "", y = "Del-Dup Correlation")

print(p_corr_bars)
```

![](Fig6_June2026_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

#### Panel G Figure

``` r
#### Panel G Figure #####
p_corr_box <- ggplot(df_corr, aes(x = Trait_Cat, y = cor_del_dup, fill = Trait_Cat)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 1, linewidth = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7, color = 'purple4') +
  scale_fill_manual(values = c("Brain Traits" = "white", "Non-Brain Traits" = "white")) +
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.8),
    axis.ticks = element_line(colour = "black", linewidth = 0.8),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(y = "Del-Dup Correlation") +
  coord_cartesian(ylim = c(-0.8, 0.45), clip = "off")

print(p_corr_box)
```

![](Fig6_June2026_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

#### Panel F & G Stats

``` r
# Display Wilcoxon Rank Sum Test for Correlation
# 1. Panel F Stats: Correlation Significance per Trait
df_corr_sig <- df_corr %>%
  select(Trait, Trait_Cat, adj_JaccardPvalue_cor_del_dup) %>%
  arrange(Trait_Cat, adj_JaccardPvalue_cor_del_dup)

kable(head(df_corr_sig), align = "c", caption = "FDR-Adjusted P-Values for Del-Dup Correlation (Panel F)")
```

|    Trait     |  Trait_Cat   | adj_JaccardPvalue_cor_del_dup |
|:------------:|:------------:|:-----------------------------:|
|  RiskTaking  | Brain Traits |           0.0000000           |
|  Depression  | Brain Traits |           0.0000000           |
| SymbolDigit  | Brain Traits |           0.0000000           |
| TrailMaking2 | Brain Traits |           0.0226667           |
|      FI      | Brain Traits |           0.0226667           |
| Neuroticism  | Brain Traits |           0.0370909           |

FDR-Adjusted P-Values for Del-Dup Correlation (Panel F)

``` r
# 2. Panel G Stats: Wilcoxon Test (Brain vs Non-Brain Correlation)
wilcox_corr <- wilcox.test(
  subset(df_corr, Trait_Cat == 'Brain Traits')$cor_del_dup,
  subset(df_corr, Trait_Cat == 'Non-Brain Traits')$cor_del_dup
)

corr_test_stats <- data.frame(
  Test = "Wilcoxon Rank Sum Test",
  Comparison = "Brain vs Non-Brain Del-Dup Correlation",
  W_Statistic = wilcox_corr$statistic,
  P_Value = format(wilcox_corr$p.value, digits = 4)
)

kable(corr_test_stats, align = "c", caption = "Statistical Comparison of Del-Dup Correlations (Panel G)")
```

|     |          Test          |               Comparison               | W_Statistic | P_Value |
|:----|:----------------------:|:--------------------------------------:|:-----------:|:-------:|
| W   | Wilcoxon Rank Sum Test | Brain vs Non-Brain Del-Dup Correlation |     126     | 0.7618  |

Statistical Comparison of Del-Dup Correlations (Panel G)
