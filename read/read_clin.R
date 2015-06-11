library(data.table)
library(yaml)    # Read config
library(utils)   # Progress bar
library(R.utils) # Verbose
library(stringr) # String parsing

# Define verbose and verbose level
verbose <- Verbose(threshold = -1)

# Read base directories from config.yml, define cancer and data type -----------
ruler(verbose)
enter(verbose, "Preparation")
cat(verbose, "Reading config.yml")
config <- yaml.load_file("/share/scratch/arj32/tcga-analysis/config.yml")
raw_data_dir    <- config$data_dirs$raw_data
parsed_data_dir <- config$data_dirs$parsed_data
cancer_type     <- config$cancer_type
data_type       <- "clin"

# Define absolute directories --------------------------------------------------
cat(verbose, "Defining directories to read from")
src_file_dir <- paste( raw_data_dir, cancer_type, "/", data_type, "/", sep = "")
output_dir <- paste(parsed_data_dir, cancer_type, "/", data_type, "/",  sep="")
filename <- paste(
  "nationwidechildrens.org_clinical_patient_", cancer_type, ".txt", sep = ""
)
exit(verbose)

# Read from only data file  ----------------------------------------------------
enter(verbose, "Reading source files")
cat(verbose, "Reading from ", src_file_dir)

clinical <- fread(
  paste(src_file_dir, filename, sep = ""), 
  select = c(2, 7:8, 11:16, 21, 41, 44, 50, 56, 99)
)[-(1:2),]
cat(verbose, "Cleaning up data")
clinical[,menopause_status := menopause_status %>%
  str_replace(" \\(.*\\)", "") %>%
  str_replace("\\[.*\\]", "Unknown")]
clinical[,gender := str_to_title(gender)]
clinical[,tumor_status := tumor_status %>%
  str_replace("TUMOR FREE", "Free") %>%
  str_replace("WITH TUMOR", "With") %>%
  str_replace("\\[.*\\]", "Unknown")]
clinical[,last_contact_days_to := as.numeric(last_contact_days_to)]
clinical[,death_days_to := as.numeric(death_days_to)]
clinical[,ajcc_pathologic_tumor_stage := ajcc_pathologic_tumor_stage %>%
  str_replace("Stage ", "") %>%
  str_replace("\\[.*\\]", "Unknown")]
clinical[,er_status_by_ihc := er_status_by_ihc %>%
  str_replace("\\[.*\\]", "Unknown")]
clinical[,pr_status_by_ihc := pr_status_by_ihc %>%
  str_replace("\\[.*\\]", "Unknown")]
clinical[,her2_status_by_ihc := her2_status_by_ihc %>%
  str_replace("\\[.*\\]", "Unknown")]
clinical[,histological_type := histological_type %>%
  str_replace("Infiltrating Carcinoma NOS", "Infiltrating Ductal Carcinoma") %>%
  str_replace("Mixed Histology (please specify)", "Mixed") %>%
  str_replace("Other  specify", "Other") %>%
  str_replace("\\[.*\\]", "Unknown")]
exit(verbose)

# Save to .RDS file ------------------------------------------------------------
enter(verbose, "Saving files")
cat(
  verbose, "Saving to cancer file to ",
  paste(output_dir, cancer_type, "_", data_type, "_cancer.Rds", sep = "")
)
saveRDS(
  clinical,
  file = paste(
    output_dir, cancer_type, "_", data_type, "_cancer.Rds", sep = ""
  )
)
exit(verbose)
