#Download & load the data
options(timeout = 600)  # 10 minutes

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "tidyverse", "glmnet", "survival", "survminer"))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)


# Clinical & survival data
query_clinical <- GDCquery(
  project = "TCGA-OV",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)
GDCdownload(query_clinical)
clinical_data <- GDCprepare(query_clinical)
names(clinical_data)

#### Get best survival data by Extracting patient and follow‑up tables
clinical_patient <- clinical_data$clinical_patient_ov
clinical_followup <- clinical_data$clinical_follow_up_v1.0_ov

library(TCGAbiolinks)
library(dplyr)
# cleaned clinical_patient
clinical_patient <- clinical_patient %>%
  mutate(
    death_time = as.numeric(na_if(death_days_to, "--")),
    follow_up_time = as.numeric(na_if(last_contact_days_to, "--")),
    vital_status = tolower(vital_status)
  ) %>%
  select(
    bcr_patient_barcode,
    death_time,
    follow_up_time,
    vital_status,
    everything()
  )

# cleaned clinical_followup
clinical_followup <- clinical_followup %>%
  mutate(
    death_time = as.numeric(na_if(death_days_to, "--")),
    follow_up_time = as.numeric(na_if(last_contact_days_to, "--")),
    vital_status = tolower(vital_status)
  ) %>%
  group_by(bcr_patient_barcode) %>%
  arrange(desc(follow_up_time)) %>%  # keep the latest follow-up per patient
  slice(1) %>%
  ungroup() %>%
  select(
    bcr_patient_barcode,
    death_time,
    follow_up_time,
    vital_status
  )

# Combine patient & follow‑up and take the most updated info: use follow-up if available, else patient
clinical_combined <- full_join(
  clinical_patient,
  clinical_followup,
  by = "bcr_patient_barcode",
  suffix = c("_pat", "_fu")
) %>%
  mutate(
    # take follow-up if present, else patient
    death_time = coalesce(death_time_fu, death_time_pat),
    follow_up_time = coalesce(follow_up_time_fu, follow_up_time_pat),
    vital_status = coalesce(vital_status_fu, vital_status_pat),
    # survival time and event
    os_time = ifelse(!is.na(death_time), death_time, follow_up_time),
    os_event = ifelse(vital_status == "dead", 1, 0)
  ) %>%
  select(bcr_patient_barcode, os_time, os_event)

# Filter to complete cases
survival_data <- clinical_combined %>%
  filter(!is.na(os_time), !is.na(os_event))

# inspect
head(survival_data)
summary(survival_data$os_time)
table(survival_data$os_event)


# Merge survival_data with clinical_patient
library(dplyr)
clinical_patient_up <- clinical_patient %>%
  inner_join(survival_data)

write.csv(clinical_patient_up, "clinical_patient_up.csv", row.names = FALSE)




#±§§§§
# 📦 Load packages
library(dplyr)
library(survival)
library(glmnet)
library(caret)

#Check survival times
table(clinical_patient_up$os_time <= 0, useNA = "ifany")
clinical_patient_up <- clinical_patient_up %>%
  filter(os_time > 0)


# Define columns to exclude
exclude_cols <- c(
  "patient_id", "bcr_patient_barcode", "bcr_patient_uuid",
  "form_completion_date", "prospective_collection", "retrospective_collection",
  "last_contact_days_to", "death_days_to", "death_time", "follow_up_time",
  "os_time", "os_event", "project_code", "system_version", "informed_consent_verified"
)

# Select predictors and outcome
predictors_df <- clinical_patient_up %>% select(-all_of(exclude_cols))
outcome_time <- clinical_patient_up$os_time
outcome_event <- clinical_patient_up$os_event


# Remove constant columns
nzv <- nearZeroVar(predictors_df, saveMetrics = TRUE)

# Show which columns are constant
cat("Dropping constant columns:\n")
print(rownames(nzv[nzv$zeroVar == TRUE, ]))

predictors_df_clean <- predictors_df[, !nzv$zeroVar, drop = FALSE]

# Make sure character columns become factors
predictors_df_clean[] <- lapply(predictors_df_clean, function(x) {
  if (is.character(x)) x <- as.factor(x)
  if (is.factor(x)) x <- droplevels(x) # drop unused levels
  return(x)
})

# 📦 Step 4: Prepare dummy-encoded predictors with model.matrix
# Add a dummy outcome just to use model.matrix
mm <- model.matrix(~ ., data = predictors_df_clean)[, -1, drop = FALSE]
X <- as.matrix(mm)

# 📦 Step 5: Survival outcome
y <- Surv(time = outcome_time, event = outcome_event)

# 📈 Step 6: Run LASSO Cox
set.seed(123)
fit_lasso <- cv.glmnet(
  x = X,
  y = y,
  family = "cox",
  alpha = 1  # LASSO
)

# Plot CV curve
plot(fit_lasso)

# Get selected features at best lambda
lasso_coef <- coef(fit_lasso, s = "lambda.min")
lasso_features <- rownames(lasso_coef)[which(lasso_coef != 0)]
cat("📊 LASSO selected features:\n")
print(lasso_features)



# Gene expression data
library(tidyverse)
library(biomaRt)
library(survival)
library(survminer)

# 1. Load TPM data (example: TSV with Ensembl transcript IDs in first column)
tpm <- read_tsv("TCGA-OV.star_tpm.tsv.gz.tsv")

# Check colnames
# First col is Ensembl transcript ID with version suffix, e.g. ENSG00000000003.15
# Samples columns after that

# 2. Clean Ensembl IDs to gene-level by removing version suffix
tpm <- tpm %>%
  rename(Ensembl_transcript = 1) %>%
  mutate(Ensembl_ID_clean = str_remove(Ensembl_transcript, "\\..*$"))

# 3. Use biomaRt to map Ensembl gene IDs to HGNC symbol & get gene biotype
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_annotation <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
  filters = "ensembl_gene_id",
  values = unique(tpm$Ensembl_ID_clean),
  mart = mart
)

# 4. Filter protein-coding genes and join with TPM data
tpm_annotated <- tpm %>%
  inner_join(gene_annotation %>% filter(gene_biotype == "protein_coding"),
             by = c("Ensembl_ID_clean" = "ensembl_gene_id"))

# 5. Aggregate TPM expression by gene symbol (sum across transcripts)
tpm_gene <- tpm_annotated %>%
  dplyr::select(-Ensembl_transcript, -Ensembl_ID_clean, -gene_biotype) %>%
  group_by(hgnc_symbol) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup()

# Ensure sample IDs match: you may need to clean sample IDs in both datasets
# For example, trim to first 12 characters of sample barcodes (TCGA-XX-XXXX)
colnames(tpm_gene)[-1] <- substr(colnames(tpm_gene)[-1], 1, 12)

# Subset clinical data to samples present in TPM data
surv_data_filtered <- survival_data %>%
  filter(bcr_patient_barcode %in% colnames(tpm_gene)[-1])

# Reorder TPM columns to match survival data samples
tpm_gene <- tpm_gene %>%
  dplyr::select(hgnc_symbol, all_of(surv_data_filtered$bcr_patient_barcode))

# 7. Run univariate Cox models for each gene
cox_results <- map_dfr(tpm_gene$hgnc_symbol, function(gene) {
  expr <- as.numeric(tpm_gene %>% filter(hgnc_symbol == gene) %>% dplyr::select(-hgnc_symbol))
  # Build dataframe for survival
  df <- surv_data_filtered %>%
    mutate(expr = expr)
  
  cox_model <- coxph(Surv(os_time, os_event) ~ expr, data = df)
  cox_summary <- summary(cox_model)
  
  tibble(
    gene = gene,
    coef = cox_summary$coefficients[1],
    exp_coef = cox_summary$conf.int[1],
    p_value = cox_summary$coefficients[5]
  )
})

# 8. Adjust p-values using BH method
cox_results <- cox_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

# 9. Output significant genes (FDR < 0.05)
significant_genes <- cox_results %>% filter(p_adj < 0.05)

#@##@#@#@#@#@#
library(dplyr)


library(dplyr)

# Check actual gene symbol column name
colnames(tpm_gene)[1]

# Compute variance per gene
gene_variance <- tpm_gene %>%
  rowwise() %>%
  mutate(variance = var(c_across(-1), na.rm = TRUE)) %>%
  ungroup()

# Keep top 5000 by variance
top_genes <- gene_variance %>%
  arrange(desc(variance)) %>%
  slice_head(n = 5000) %>%
  dplyr::select(-variance)

# Result
tpm_gene_filtered <- top_genes


# 7. Run univariate Cox models for each gene
cox_results <- map_dfr(tpm_gene_filtered$hgnc_symbol, function(gene) {
  expr <- as.numeric(tpm_gene_filtered %>% filter(hgnc_symbol == gene) %>% dplyr::select(-hgnc_symbol))
  # Build dataframe for survival
  df <- surv_data_filtered %>%
    mutate(expr = expr)
  
  cox_model <- coxph(Surv(os_time, os_event) ~ expr, data = df)
  cox_summary <- summary(cox_model)
  
  tibble(
    gene = gene,
    coef = cox_summary$coefficients[1],
    exp_coef = cox_summary$conf.int[1],
    p_value = cox_summary$coefficients[5]
  )
})


# 8. Adjust p-values using BH method
cox_results <- cox_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

# 9. Output significant genes (FDR < 0.05)
significant_genes <- cox_results %>% filter(p_adj < 0.1)
#####result -> failed no gene was statistically significat for predicting survival

####continue with relevant gene list based on current literature.
library(survival)
library(dplyr)

# Define your gene list
genes_of_interest <- c(
  "BRCA1", "BRCA2", "RAD51", "PALB2", "ATM", "ATR", "CHEK1", "CHEK2", "FANCA", "FANCD2",
  "TP53", "BCL2", "BAX", "MDM2",
  "VEGFA", "FLT1", "KDR",
  "CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2",
  "MLH1", "MSH2", "MSH6", "PMS2",
  "PARP1", "HIF1A", "CCNE1", "MYC", "AKT1", "PIK3CA", "KRAS", "NRAS"
)

# Filter your gene expression matrix
tpm_selected <- tpm_gene %>% filter(hgnc_symbol %in% genes_of_interest)


tpm_selected$hgnc_symbol
library(tibble)
library(dplyr)

# Make rows = samples, columns = genes # Transpose & prepare for survival analysis
expr_matrix <- tpm_selected %>%
  tibble::column_to_rownames("hgnc_symbol") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample")

head(expr_matrix)
dim(expr_matrix)  # should be ~400 × 36


# Merge with survival
surv_data <- survival_data %>%
  rename(sample = bcr_patient_barcode)

merged_data <- inner_join(surv_data, expr_matrix, by = "sample")

# Run univariate Cox
results <- lapply(genes_of_interest, function(gene) {
  cox <- coxph(Surv(os_time, os_event) ~ merged_data[[gene]], data = merged_data)
  summary_cox <- summary(cox)
  tibble(
    gene = gene,
    HR = summary_cox$coef[2],
    lower95 = summary_cox$conf.int[,"lower .95"],
    upper95 = summary_cox$conf.int[,"upper .95"],
    pval = summary_cox$coef[5]
  )
})

results_df <- bind_rows(results) %>%
  mutate(padj_bonferroni = p.adjust(pval, method = "bonferroni"))

# View significant ones
results_df %>% filter(padj_bonferroni < 0.05)
# since nothing survives Bonferroni correction, you can also look at raw p-values < 0.05 to get candidates
results_df %>% filter(pval < 0.05)
results_df %>% filter(pval < 0.1)

#or check if they survive BH correction
results_df <- results_df %>%
  mutate(padj_bh = p.adjust(pval, method = "BH"))
results_df %>% filter(padj_bh < 0.05)
####$$$$ no they didnot. done now we only have CTLA4, MYC as predictor and normal p-value 
## testins our hypothesis for biologically relevant genes no need for corecction it was not large number of hypothesis testing

#KM plots: MYC & CTLA4
library(survival)
library(survminer)
library(dplyr)

# assume: 
# - survival_data: your clinical data (bcr_patient_barcode, os_time, os_event)
# - tpm_selected: your TPM matrix of 34 genes (rownames = hgnc_symbol, columns = samples)

# transpose and merge survival with gene expression
expr_df <- as.data.frame(t(tpm_selected[ , -1]))  # drop hgnc_symbol column if present
colnames(expr_df) <- tpm_selected$hgnc_symbol
expr_df$sample <- colnames(tpm_selected)[-1]

surv_df <- survival_data %>% 
  rename(sample = bcr_patient_barcode)

merged_df <- inner_join(surv_df, expr_df, by = "sample")

# Function to plot KM by median
library(survminer)
plot_km_gene <- function(gene) {
  merged_df <- merged_df %>%
    mutate(gene_group = ifelse(.data[[gene]] > median(.data[[gene]], na.rm=TRUE), 
                               "High", "Low"))
  
  fit <- survfit(Surv(os_time, os_event) ~ gene_group, data = merged_df)
  
  ggsurvplot(
    fit,
    data = merged_df,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    legend.title = gene,
    legend.labs = c("High", "Low"),
    xlab = "Time (days)",
    ylab = "Survival Probability",
    title = paste("Kaplan-Meier Curve for", gene),
    
    # Style options
    palette = c("red2", "pink2"),  # golden and blue
    font.main = c(16, "bold", "gray17"), 
    font.x = c(14, "bold"),
    font.y = c(14, "bold"),
    font.tickslab = c(12, "plain"),
    
    risk.table.title = "Number at risk",
    risk.table.height = 0.25,
    risk.table.fontsize = 4,
    
    legend = "top",
    legend.box = "horizontal",
    
    ggtheme = theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top"
      )
  )
}

# plot for MYC
plot_km_gene("MYC")

# plot for CTLA4
plot_km_gene("CTLA4")


# Forest plot using ggplot2
plot_forest <- function(results_df) {
  results_df <- results_df %>%
    mutate(significant = ifelse(pval < 0.05, "Yes", "No"))
  
  ggplot(results_df, aes(x = HR, y = reorder(gene, HR), xmin = lower95, xmax = upper95, color = significant)) +
    geom_point(size = 3) +
    geom_errorbarh(height = 0.25, size = 1) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    scale_x_log10(
      breaks = c(0.5, 1, 2, 5),
      labels = c("0.5", "1", "2", "5")
    ) +
    scale_color_manual(values = c("Yes" = "#D55E00", "No" = "gray60")) +
    labs(
      title = "Hazard Ratios for Selected Genes",
      x = "Hazard Ratio (log scale)",
      y = NULL,
      color = "Significant (p < 0.05)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(face = "italic"),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    )
}

# Example usage:
plot_forest(results_df)



#####Build the single dataset#####
library(dplyr)

clinical_clean <- clinical_patient_up %>%
    dplyr::select(
    bcr_patient_barcode,
    initial_pathologic_dx_year,
    age_at_initial_pathologic_diagnosis,
    tumor_status,
    tumor_grade,
    residual_disease_largest_nodule,
    lymphovascular_invasion_indicator,
    karnofsky_score,
    clinical_stage,
    tissue_source_site,
    os_time,
    os_event
  )

# Take a quick look
glimpse(clinical_clean)
head(clinical_data)

#reshape tpm_selected
library(dplyr)
library(tidyr)

# remove gene column temporarily and transpose
tpm_long <- tpm_selected %>%
  pivot_longer(
    cols = -hgnc_symbol,
    names_to = "bcr_patient_barcode",
    values_to = "expression"
  )

# reshape to wide: one row per patient
tpm_wide <- tpm_long %>%
  pivot_wider(
    names_from = hgnc_symbol,
    values_from = expression
  )

# look at result
glimpse(tpm_wide)

#merge clinical and gene expression
merged_data <- clinical_clean %>%
  inner_join(tpm_wide, by = "bcr_patient_barcode")



#Analysing
#Descriptive Statistics 
merged_data <- read.csv("merged_data.csv")



#Continuous variables:
merged_data <- merged_data %>% mutate(
  age_at_initial_pathologic_diagnosis = as.numeric(age_at_initial_pathologic_diagnosis),
  karnofsky_score = as.numeric(karnofsky_score),
  residual_disease_largest_nodule = as.numeric(residual_disease_largest_nodule),
  os_time = as.numeric(os_time),
  os_event = as.numeric(os_event)
)

merged_data %>% 
  summarise(
    age_mean = mean(age_at_initial_pathologic_diagnosis, na.rm = TRUE),
    age_sd = sd(age_at_initial_pathologic_diagnosis, na.rm = TRUE),
    karnofsky_mean = mean(karnofsky_score, na.rm = TRUE),
    karnofsky_sd = sd(karnofsky_score, na.rm = TRUE),
    residual_disease_mean = mean(residual_disease_largest_nodule, na.rm = TRUE),
    residual_disease_sd = sd(residual_disease_largest_nodule, na.rm = TRUE),
    os_time_mean = mean(os_time, na.rm = TRUE),
    os_time_sd = sd(os_time, na.rm = TRUE)
  )


library(psych)
merged_data %>% 
  dplyr::select(age_at_initial_pathologic_diagnosis, karnofsky_score, 
         residual_disease_largest_nodule, os_time) %>%
  describe()


#Categorical variables:
categorical_vars <- c("tumor_grade", "clinical_stage", "tumor_status",
                      "lymphovascular_invasion_indicator", "tissue_source_site", "os_event")

for (var in categorical_vars) {
  cat("===== ", var, " =====\n")
  print(merged_data %>% count(.data[[var]]))
  cat("\n")
}
# repeat for other categorical variables
library(janitor)
merged_data %>%
  tabyl(tumor_status) %>%
  adorn_totals()


#Continuous vars — Boxplots:Age
library(ggplot2)

ggplot(merged_data, aes(x = "", 
                        y = age_at_initial_pathologic_diagnosis)) +
  geom_boxplot(fill = "#69b3a2", color = "#1b4d3e", width = 0.3, outlier.color = "red", outlier.shape = 16, outlier.size = 2) +
  geom_jitter(width = 0.1, alpha = 0.3, color = "gray40", size = 1) +
  labs(title = "Age at Initial Pathologic Diagnosis",
       subtitle = "Distribution of age in years (n = 421)",
       y = "Age (years)",
       x = "") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

#Tumor Status (barplot)
merged_data %>%
  count(tumor_status) %>%
  mutate(perc = n / sum(n) * 100,
         label = paste0(tumor_status, "\n", round(perc, 1), "%")) %>%
  ggplot(aes(x = "", y = perc, fill = tumor_status)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
  scale_fill_brewer(palette = "Pastel1") +
  labs(title = "Tumor Status Distribution") +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "none"
  )



#Tumor Grade (barplot)
ggplot(merged_data, aes(x = tumor_grade)) +
  geom_bar(fill = "#d35400", color = "black") +
  labs(
    title = "Tumor Grade",
    subtitle = "Distribution of histologic tumor grades",
    x = "Tumor Grade", y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )


#Lymphovascular Invasion (pie chart)
library(ggplot2)
# Clean up (optional: drop NAs or relabel)
lv_data <- merged_data %>%
  filter(!is.na(lymphovascular_invasion_indicator)) %>%
  mutate(lymphovascular_invasion_indicator = factor(
    lymphovascular_invasion_indicator,
    levels = c("NO", "YES"),
    labels = c("No Invasion", "Invasion")
  ))
# Summary counts
lv_counts <- lv_data %>%
  count(lymphovascular_invasion_indicator) %>%
  mutate(perc = n / sum(n) * 100,
         label = paste0(lymphovascular_invasion_indicator, "\n", round(perc, 1), "%"))

ggplot(lv_counts, aes(x = "", y = perc, fill = lymphovascular_invasion_indicator)) +
  geom_col(width = 1, color = "gray 14") +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), color = "maroon 1", size = 5) +
  scale_fill_viridis_d(option = "inferno") +
  labs(
    title = "Lymphovascular Invasion in Patients",
    fill = "Invasion"
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )
#or Bar Chart by Tumor Grade
ggplot(merged_data, aes(x = tumor_grade, fill = lymphovascular_invasion_indicator)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(option = "plasma") +
  labs(
    title = "Proportion of Lymphovascular Invasion by Tumor Grade",
    x = "Tumor Grade", y = "Proportion",
    fill = "Invasion"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )


#Clinical Stage (barplot)
ggplot(merged_data, aes(x = clinical_stage)) +
  geom_bar(fill = "#16a085", color = "black") +
  labs(
    title = "Clinical Stage at Diagnosis",
    subtitle = "Stage distribution",
    x = "Clinical Stage", y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5)
  )


#Tissue Source Site (barplot)
ggplot(merged_data, aes(x = tissue_source_site)) +
  geom_bar(fill = "#c0392b", color = "black") +
  labs(
    title = "Tissue Source Site",
    subtitle = "Distribution of source sites",
    x = "Tissue Site", y = "Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


#age vs clinical_stage
#This shows if younger or older patients tend to present at certain stages
ggplot(merged_data, aes(x = clinical_stage, y = age_at_initial_pathologic_diagnosis, fill = clinical_stage)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Age at Diagnosis by Clinical Stage",
    x = "Clinical Stage", y = "Age (years)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    legend.position = "none"
  )

#Violin plot: Age vs Tumor Grade
library(ggplot2)
library(viridis)
# Convert tumor_grade to factor
merged_data$tumor_grade <- factor(merged_data$tumor_grade)

ggplot(merged_data, aes(x = tumor_grade, 
                        y = age_at_initial_pathologic_diagnosis, 
                        fill = tumor_grade)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = "gray30") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_viridis(discrete = TRUE, option = "plasma") +
  labs(
    title = "Age at Diagnosis by Tumor Grade",
    x = "Tumor Grade",
    y = "Age (years)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "none"
  )


# Survival Summary:
table(merged_data$os_event)
survfit(Surv(os_time, os_event) ~ 1, data = merged_data) %>% summary()

#simple Kaplan-Meier plot:
library(survival)
library(survminer)
km_fit <- survfit(Surv(os_time, os_event) ~ 1, data = merged_data)
ggsurvplot(km_fit, data = merged_data,
           conf.int = TRUE, risk.table = TRUE,
           title = "Overall Survival (All Patients)")
######end of discriptive#######


#vis and presentation
merged_data <- read_csv('merged_data.csv')
library(dplyr)

cleaned_data <- merged_data %>%
  dplyr::select(
    bcr_patient_barcode,
    age_at_initial_pathologic_diagnosis,
    tumor_status,
    tumor_grade,
    residual_disease_largest_nodule,
    clinical_stage,
    tissue_source_site,
    os_time,
    os_event,
    CTLA4,
    MYC
  )


library(survival)
library(survminer)

# Make gene expression standardized or scaled if desired
cleaned_data$CTLA4 <- scale(cleaned_data$CTLA4)
cleaned_data$MYC   <- scale(cleaned_data$MYC)

colnames(merged_data)
# Convert character columns to factors
library(dplyr)

cleaned_data <- cleaned_data %>%
  mutate(
    tumor_status = as.factor(tumor_status),
    tumor_grade = as.factor(tumor_grade),
    residual_disease_largest_nodule = as.factor(residual_disease_largest_nodule),
    clinical_stage = as.factor(clinical_stage)
  )


# Fit Cox model
cox_model <- coxph(Surv(os_time, os_event) ~ 
                     CTLA4 + MYC +
                     age_at_initial_pathologic_diagnosis +
                     tumor_status +
                     tumor_grade +
                     residual_disease_largest_nodule +
                     clinical_stage +
                     tissue_source_site,
                   data = cleaned_data)

summary(cox_model)

colSums(is.na(cleaned_data))
cleaned_data_complete <- cleaned_data %>%
  drop_na(os_time, os_event, age_at_initial_pathologic_diagnosis,
          tumor_status, tumor_grade, residual_disease_largest_nodule,
          clinical_stage, CTLA4, MYC)

cx<- coxph(Surv(os_time, os_event) ~ age_at_initial_pathologic_diagnosis +
        tumor_status + tumor_grade + residual_disease_largest_nodule + 
        clinical_stage + CTLA4 + MYC, data = cleaned_data_complete)

# Forest plot
ggforest(cox_model, data = cleaned_data)

#AGE normal distribution check
summary(cleaned_data$agećat_initial_pathologic_diagnosis)
hist(as.numeric(cleaned_data$age_at_initial_pathologic_diagnosis), 
     breaks = 20, main = "Age at diagnosis", xlab = "Age")
#AGE quantile
quantile(as.numeric(cleaned_data$age_at_initial_pathologic_diagnosis), 
         probs = c(0.25, 0.5, 0.75))

#AGE Linear vs categorical
library(survival)
cox_age_cont <- coxph(Surv(os_time, os_event) ~ as.numeric(age_at_initial_pathologic_diagnosis), data = cleaned_data)
cox_age_cat <- coxph(Surv(os_time, os_event) ~ cut(as.numeric(age_at_initial_pathologic_diagnosis), 
                                                   breaks = c(-Inf, 51, 68, Inf), 
                                                   labels = c("<51", "51–68", ">68")), 
                     data = cleaned_data)
AIC(cox_age_cont, cox_age_cat)


#modify 
library(dplyr)

# Check the actual quantiles of age
quantile(cleaned_data$age_at_initial_pathologic_diagnosis, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

# Use these to cut into 4 groups
cleaned_data <- cleaned_data %>%
  mutate(
    age_group = cut(
      age_at_initial_pathologic_diagnosis,
      breaks = c(-Inf, 51, 59, 68, Inf),
      labels = c("<=51", "52-59", "60-68", ">68")
    )
  )

# See the distribution
table(cleaned_data$age_group, useNA = "ifany")


# Create the Surv object
y_surv <- with(cleaned_data, Surv(os_time, os_event))

# Full Cox with all predictors
cox_full <- coxph(y_surv ~ . - bcr_patient_barcode, data = cleaned_data)
summary(cox_full)

# C-index for full model
lp_full <- predict(cox_full, type = "lp")    # linear predictor
cindex_full <- rcorr.cens(lp_full, y_surv)["C Index"]
cat("C-index (full Cox model):", round(cindex_full, 3), "\n")



#@#@## final model 
library(survival)

cox_fit <- coxph(
  Surv(os_time, os_event) ~ 
    age_group + tumor_status + tumor_grade +
    residual_disease_largest_nodule + clinical_stage +
    tissue_source_site + CTLA4 + MYC,
  data = cleaned_data
)

summary(cox_fit)
#test the ph assumption
cox_zph <- cox.zph(cox_fit)
plot(cox_zph)

cox_strat <- coxph(
  Surv(os_time, os_event) ~ tumor_status + residual_disease_largest_nodule + tissue_source_site + CTLA4 + MYC + strata(age_group) + strata(tumor_grade) + strata(clinical_stage),
  data = cleaned_data
)
summary(cox_strat)
#to confirm the PH assumption now holds:
zph_strat <- cox.zph(cox_strat)
print(zph_strat)


#eror!
plot(cox_zph, var = "age_at_initial_pathologic_diagnosis")
plot(ph_test, var = "CTLA4")
plot(ph_test, var = "MYC")



ggforest(
  cox_fit,
  data = cleaned_data,
  main = "Hazard ratios of prognostic factors",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize = 1,
  refLabel = "Reference",
  noDigits = 2
)

fit <- survfit(cox_fit)
plot(fit, xlab = "Time", ylab = "Survival Probability")


#martingle
martingale_resid <- residuals(cox_fit, type = "martingale")
str(martingale_resid)

# Plot against a continuous predictor (e.g., age)
plot(cleaned_data$age_group, martingale_resid,
     xlab = "Age", ylab = "Martingale residuals",
     main = "Martingale residuals vs Age")
abline(h = 0, lty = 2)

# Plot Schoenfeld residuals for tissue_source (optional) # plots residuals vs time for each covariate
plot(cox_fit)
str(cox_fit)

####@#@#@#@#@#±±±±§§§§§
library(survival)
library(survminer)   # for pretty plots if you like

# Fit your multivariable Cox model
cox_model <- coxph(
  Surv(os_time, os_event) ~ 
    age_at_initial_pathologic_diagnosis + 
    tumor_status + tumor_grade +
    residual_disease_largest_nodule + 
    clinical_stage + tissue_source_site + 
    CTLA4 + MYC,
  data = cleaned_data
)

summary(cox_model)

# ----------------------------------------------------
# 1️⃣ Martingale residuals vs continuous predictors
# ----------------------------------------------------

martingale_resid <- residuals(cox_model, type = "martingale")

# Plot vs Age
plot(
  cleaned_data$age_at_initial_pathologic_diagnosis,
  martingale_resid,
  xlab = "Age at diagnosis",
  ylab = "Martingale residuals",
  main = "Martingale residuals vs Age"
)
abline(h = 0, lty = 2)

# You can repeat for gene expressions:
plot(
  cleaned_data$CTLA4,
  martingale_resid,
  xlab = "CTLA4 (scaled)",
  ylab = "Martingale residuals",
  main = "Martingale residuals vs CTLA4"
)
abline(h = 0, lty = 2)

plot(
  cleaned_data$MYC,
  martingale_resid,
  xlab = "MYC (scaled)",
  ylab = "Martingale residuals",
  main = "Martingale residuals vs MYC"
)
abline(h = 0, lty = 2)

# ----------------------------------------------------
# 2️⃣ Schoenfeld residuals + PH assumption test
# ----------------------------------------------------

ph_test <- cox.zph(cox_model)
print(ph_test)

# Global test p < 0.05 suggests PH violated.
# Plot Schoenfeld residuals for each covariate
plot(ph_test, var = "age_at_initial_pathologic_diagnosis")
plot(ph_test, var = "CTLA4")
plot(ph_test, var = "MYC")
# You can loop through other categorical ones if desired

# Or all at once
plot(ph_test)

# ----------------------------------------------------
# 3️⃣ Deviance residuals QQ plot
# ----------------------------------------------------
deviance_resid <- residuals(cox_fit, type = "deviance")

# Remove bad values
good_resid <- deviance_resid[is.finite(deviance_resid)]

# Plot QQ
qqnorm(good_resid, main = "QQ plot of deviance residuals")
qqline(good_resid)


deviance_resid <- residuals(cox_model, type = "deviance")

qqnorm(deviance_resid, main = "QQ plot of deviance residuals")
qqline(deviance_resid)

# ----------------------------------------------------
# 4️⃣ Concordance & global fit
# ----------------------------------------------------

cat("Concordance (c-index):", summary(cox_model)$concordance[1], "\n")

anova(cox_model, test = "Chisq")

# You can also use:
# survminer::ggcoxdiagnostics(cox_model, type = "dfbeta") etc.


###$%#%


library(survival)

# Assuming your data is in a dataframe called df
# with columns: time, status, tissue_source

# Table of tissue_source counts
table(cleaned_data$tissue_source)

# Univariate Cox model for tissue_source
cox_uni <- coxph(Surv(os_time, os_event) ~ tissue_source_site, data = cleaned_data)
summary(cox_uni)

# Test PH assumption
ph_test <- cox.zph(cox_uni)
print(ph_test)

# Plot Schoenfeld residuals for tissue_source (optional)
plot(ph_test)


#
library(survival)

# Fit Cox model with stratification on tissue_source
cox_model <- coxph(Surv(os_time, os_event) ~  age_group + tumor_status + tumor_grade +
                     residual_disease_largest_nodule + clinical_stage +
                     strata(tissue_source_site) + CTLA4 + MYC,
                   data = cleaned_data)

# View model summary
summary(cox_model)

# Check proportional hazards assumption again (optional)
ph_test <- cox.zph(cox_model)
print(ph_test)
plot(ph_test)  # visual check for residuals


table(cleaned_data$tumor_status)
cleaned_data$tumor_status_bin <- ifelse(cleaned_data$tumor_status == "TUMOR FREE", 1, 0)
table(cleaned_data$tumor_status_bin)
logit_fit <- glm(
  tumor_status_bin ~ age_group + tumor_grade + clinical_stage +
    residual_disease_largest_nodule + tissue_source_site + CTLA4 + MYC,
  data = cleaned_data,
  family = binomial
)

summary(logit_fit)


logit_fit <- glm(
  os_event ~ age_group + tumor_status + tumor_grade + 
    residual_disease_largest_nodule + clinical_stage + 
    tissue_source_site + CTLA4 + MYC,
  data = cleaned_data,
  family = binomial()
)

summary(logit_fit)

exp(cbind(OR = coef(logit_fit), confint(logit_fit)))



install.packages("forestplot")

library(forestplot)

# Extract estimates
coefs <- summary(cox_strat)$coefficients
ci <- confint(cox_strat)

# Create a table
hr_table <- data.frame(
  Variable = rownames(coefs),
  HR = exp(coefs[,1]),
  Lower = exp(ci[,1]),
  Upper = exp(ci[,2]),
  p = coefs[,5]
)

# Filter only significant variables or those of interest
hr_table <- hr_table[hr_table$Variable %in% c("CTLA4", "MYC"),]

# Make forest plot
forestplot::forestplot(
  labeltext = hr_table$Variable,
  mean = hr_table$HR,
  lower = hr_table$Lower,
  upper = hr_table$Upper,
  xlog = TRUE,
  col = forestplot::fpColors(box="royalblue", line="darkblue", summary="royalblue"),
  title = "Hazard Ratios with 95% CI"
)

#@#@
summary(cleaned_data)
sapply(cleaned_data, function(x) length(unique(x)))


library(survival)

cox_fit <- coxph(
  Surv(os_time, os_event) ~ .,
  data = cleaned_data
)

summary(cox_fit)


cox.zph(cox_fit)


logit_fit <- glm(
  os_event ~ tumor_status + residual_disease_largest_nodule +
    tissue_source_site + CTLA4 + MYC +
    age_group + tumor_grade + clinical_stage,
  data = cleaned_data,
  family = binomial()
)

summary(logit_fit)

# Odds ratios
exp(coef(logit_fit))

# OR with 95% CI
exp(cbind(OR = coef(logit_fit), confint(logit_fit)))


OR_CI <- exp(cbind(
  OR = coef(logit_fit),
  confint.default(logit_fit)
))
round(OR_CI, 3)

summary(logit_fit)
table(cleaned_data$tissue_source_site, cleaned_data$os_event)

# Look at total counts per level
table_sites <- table(cleaned_data$tissue_source_site)
rare_levels <- names(table_sites[table_sites <= 5])

# Recode into "Other"
cleaned_data$tissue_source_site_collapsed <- as.character(cleaned_data$tissue_source_site)
cleaned_data$tissue_source_site_collapsed[cleaned_data$tissue_source_site_collapsed %in% rare_levels] <- "Other"

# Make it a factor again
cleaned_data$tissue_source_site_collapsed <- factor(cleaned_data$tissue_source_site_collapsed)


logit_fit <- glm(os_event ~ age_group + tumor_status + tumor_grade + 
                   residual_disease_largest_nodule + clinical_stage + 
                   tissue_source_site_collapsed + CTLA4 + MYC,
                 data = cleaned_data,
                 family = binomial)

summary(logit_fit)

# OR and CI
exp(cbind(OR = coef(logit_fit), confint.default(logit_fit)))

library(broom)
tidy(logit_fit, conf.int = TRUE, exponentiate = TRUE)

logit_fit <- glm(os_event ~ age_group + tumor_status + tumor_grade + 
                   residual_disease_largest_nodule + clinical_stage + 
                   tissue_source_site_collapsed + CTLA4 + MYC,
                 data = cleaned_data,
                 family = binomial)

library(broom)

logit_results <- tidy(logit_fit, exponentiate = TRUE, conf.int = TRUE)

sig_vars <- logit_results[logit_results$p.value < 0.05, ]
print(sig_vars)


library(survival)
library(survminer)
surv_obj <- Surv(time = cleaned_data$os_time, event = cleaned_data$os_event)
km_fit <- survfit(surv_obj ~ tumor_status, data = cleaned_data)

ggsurvplot(
  km_fit,
  data = cleaned_data,
  risk.table = TRUE,            # show number at risk table
  pval = TRUE,                  # show log-rank test p-value
  conf.int = TRUE,              # show confidence intervals
  xlab = "Time (days)",         # or months, depending on your units
  ylab = "Survival probability",
  legend.title = "Tumor Status",
  legend.labs = levels(cleaned_data$tumor_status),
  palette = c("#E69F00", "#56B4E9", "#009E73")
)


cleaned_km_data <- subset(
  cleaned_data,
  tumor_status != "Not Available" & !is.na(tumor_status)
)

library(dplyr)

cleaned_km_data <- cleaned_data %>%
  filter(tumor_status != "Not Available", !is.na(tumor_status))

surv_obj <- Surv(time = cleaned_km_data$os_time, event = cleaned_km_data$os_event)
km_fit <- survfit(surv_obj ~ tumor_status, data = cleaned_km_data)

ggsurvplot(
  km_fit,
  data = cleaned_km_data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  xlab = "Time (days)",
  ylab = "Survival probability",
  legend.title = "Tumor Status",
  legend.labs = levels(cleaned_km_data$tumor_status),
  palette = c("#E69F00", "#56B4E9", "#009E73")
)


km_data_rd <- subset(
  cleaned_data,
  residual_disease_largest_nodule != "Not Available" & !is.na(residual_disease_largest_nodule)
)

# Optional: drop unused factor levels
km_data_rd$residual_disease_largest_nodule <- droplevels(merged_data$residual_disease_largest_nodule)
surv_obj_rd <- Surv(time = merged_data$os_time, event = merged_data$os_event)
km_fit_rd <- survfit(surv_obj_rd ~ residual_disease_largest_nodule, data = merged_data)
ggsurvplot(
  km_fit_rd,
  data = km_data_rd,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  xlab = "Time (days)",
  ylab = "Survival probability",
  legend.title = "Residual Disease (mm)",
  legend.labs = levels(km_data_rd$residual_disease_largest_nodule),
  palette = "Dark2"
)



