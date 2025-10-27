# --- Setup --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(TwoSampleMR)
  library(ieugwasr)
  library(stringr)
  library(qvalue)
  # library(genetics.binaRies) # used implicitly for PLINK binary
})

# Check OpenGWAS access (update your token in ~/.Renviron if needed)
# Visit https://api.opengwas.io/profile/ -> set IEU token, then restart R.
user()  # should print your user info

# --- Paths (edit these) --------------------------------------

REF_BFILE <- "/Users/nagrodzkij/cam/IMT/MS/lipidomics/LD/working/eur_cleaned_10feb"

OUTCOME_TSV  <- "/Users/nagrodzkij/cam/IMT/MS/lipidomics/extended/severity/imsgc_mssev_discovery_hg38.tsv"
EXPO_DIR     <- "/Users/nagrodzkij/cam/IMT/MS/lipidomics/extended/exposures"
WORK_DIR     <- "/Users/nagrodzkij/cam/IMT/MS/lipidomics/extended/severity/working/new"
RESULTS_DIR  <- "/Users/nagrodzkij/cam/IMT/MS/lipidomics/extended/severity/results/new"

# Combined results output
COMBINED_CSV <- file.path(RESULTS_DIR, "_combined.csv")

# Ensure required paths exist
stopifnot(
  file.exists(REF_BFILE),
  file.exists(OUTCOME_TSV),
  dir.exists(EXPO_DIR),
  dir.exists(WORK_DIR),
  dir.exists(RESULTS_DIR)
)

# --- Outcome data -------------------------------------------------------------

# Read MS severity GWAS; rename to TwoSampleMR conventions; create chr:pos "SNP"
ms_dat <- readr::read_tsv(
  OUTCOME_TSV,
  col_types = "dcccddddddcc"
) %>%
  dplyr::rename(
    effect_allele = A1,
    other_allele  = A2,
    beta          = BETA,
    se            = SE,
    rsID          = SNP,
    eaf           = AF1,
    pval          = P
  ) %>%
  tibble::add_column(Phenotype = "MS severity") %>%
  dplyr::mutate(SNP = stringr::str_sub(chrpos, 4, -1))

# Format outcome for TwoSampleMR
outcome_dat <- TwoSampleMR::format_data(
  ms_dat %>%
    dplyr::select(-SNP) %>%
    dplyr::rename(SNP = rsID),
  type = "outcome"
)

# --- Exposure files -----------------------------------------------------------

gwas_raw_files <- list.files(EXPO_DIR, pattern = "filtered", full.names = TRUE)
n_files <- length(gwas_raw_files)

if (n_files == 0) {
  stop("No exposure files found in: ", EXPO_DIR)
}

combinedresults_list <- vector("list", length = n_files)

# --- Main loop ----------------------------------------------------------------

i <- 0
for (file in gwas_raw_files) {
  i <- i + 1
  message("── Processing [", i, "/", n_files, "]: ", file)
  
  # Exposure name
  # NOTE: adjust if filename pattern changes.
  exposure_name <- stringr::str_sub(file, 60, -5)
  
  # Output file paths for this exposure
  fexp_name        <- file.path(EXPO_DIR, paste0(exposure_name, ".csv"))
  fharm_name       <- file.path(WORK_DIR, paste0(exposure_name, "_harmonise.csv"))
  fressingle_name  <- file.path(RESULTS_DIR, paste0(exposure_name, "_resultssingle.csv"))
  fresmrpresso_name<- file.path(RESULTS_DIR, paste0(exposure_name, "_mrpresso.csv"))
  
  # Skip if single-SNP results already exist
  if (file.exists(fressingle_name)) {
    message("   → Results single already done for ", exposure_name)
    next
  }
  
  # --- Exposure QC & formatting  -----------------------------
  
  dat <- readr::read_tsv(file)
  message("   • filtering ", basename(file))
  
  dat <- dat %>%
    dplyr::filter(nchar(effect_allele) == 1 & nchar(other_allele) == 1) %>%
    dplyr::filter(chromosome %in% c(1:22)) %>%
    dplyr::filter(p_val <= 5e-8) %>%
    dplyr::filter(effect_allele_frequency >= 0.05 & effect_allele_frequency <= 0.95) %>%
    dplyr::filter(INFO_UKBB_EUR > 0.9 & INFO_EstBB > 0.9) %>%
    dplyr::filter(
      !(effect_allele == "A" & other_allele == "T"),
      !(effect_allele == "T" & other_allele == "A"),
      !(effect_allele == "C" & other_allele == "G"),
      !(effect_allele == "G" & other_allele == "C")
    ) %>%
    dplyr::mutate(SNP = paste0(chromosome, ":", base_pair_location)) %>%
    dplyr::select(
      chromosome, base_pair_location, effect_allele, other_allele,
      beta, standard_error, effect_allele_frequency, SNP, p_val
    ) %>%
    tibble::add_column(Phenotype = exposure_name)
  
  message("   • filtering done")
  
  # Keep only SNPs present in outcome GWAS & with matching alleles; bring rsID
  dat <- dat %>%
    dplyr::filter(SNP %in% ms_dat$SNP) %>%
    dplyr::inner_join(
      ms_dat %>% dplyr::select(SNP, rsID, effect_allele, other_allele),
      by = c("SNP", "effect_allele", "other_allele")
    )
  
  message("   • matching on chr:pos done")
  
  # Format exposure for TwoSampleMR
  exposure_dat <- TwoSampleMR::format_data(
    dat %>%
      dplyr::select(-SNP) %>%
      dplyr::rename(
        se   = standard_error,
        pval = p_val,
        eaf  = effect_allele_frequency,
        SNP  = rsID
      ) %>%
      dplyr::select(SNP, dplyr::everything()),
    type = "exposure"
  )
  
  # --- LD clumping   -------------------------------------
  clumped_dat1 <- ieugwasr::ld_clump(
    dplyr::tibble(rsid = dat$rsID, pval = dat$p_val),
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = REF_BFILE,
    clump_p = 5e-8
  )
  
  exposure_dat <- exposure_dat %>%
    dplyr::filter(SNP %in% clumped_dat1$rsid)
  
  # --- Harmonise & MR ---------------------------------------------------------
  
  combo_dat <- TwoSampleMR::harmonise_data(exposure_dat, outcome_dat)
  write.csv(combo_dat, fharm_name, row.names = FALSE)
  
  res <- TwoSampleMR::mr(combo_dat, method_list = c("mr_ivw_mre"))
  
  # Keep IVW-MRE row 
  singleresults <- res[1, , drop = FALSE]
  
  # Heterogeneity (IVW-MRE)
  heterogeneity <- TwoSampleMR::mr_heterogeneity(combo_dat, method_list = c("mr_ivw_mre"))
  singleresults <- merge(
    singleresults, heterogeneity,
    by = c("id.exposure", "id.outcome", "outcome", "exposure")
  )
  
  # Pleiotropy (Egger intercept)
  pleiotropy <- TwoSampleMR::mr_pleiotropy_test(combo_dat)
  singleresults <- merge(
    singleresults, pleiotropy,
    by = c("id.exposure", "id.outcome", "outcome", "exposure")
  )
  
  # MR-PRESSO 
  mr_presso <- run_mr_presso(combo_dat)
  mr_presso_global <- mr_presso[[1]][[2]][[1]][["Pvalue"]]
  mr_presso_raw <- mr_presso[[1]][[1]][1, , drop = FALSE]
  mr_presso_outlier <- mr_presso[[1]][[1]][2, , drop = FALSE]
  rownames(mr_presso_outlier) <- 1
  mr_presso_combined <- merge(mr_presso_global, mr_presso_raw)
  mr_presso_combined <- merge(mr_presso_combined, mr_presso_outlier, by="row.names")
  mr_presso_combined <- dplyr::select(mr_presso_combined, -c('Row.names'))
  
  write.csv(mr_presso_combined, fresmrpresso_name, row.names = FALSE)
  
  
  # Merge PRESSO outlier row
  singleresults <- merge(singleresults, mr_presso_combined)
  
  
  # Standardize column names
  newcolumns <- c(
    "id.exposure","id.outcome","outcome","exposure",
    "method1","nsnp.method1","b.method1","se.method1","pval.method1",
    "method_heterogen","Q.heterogen","Q_df.heterogen","Q_pval.heterogen",
    "egger_intercept","se.egger_intercept","pval.egger_intercept",
    "pval.MR_PRESSO_global",
    "Exposure","MR_PRESSO_analysis","b.MR_PRESSO","sd.MR_PRESSO","T-stat.MR_presso","pval.MR_PRESSO",
    "Exposure","MR_PRESSO_analysis","b.MR_PRESSO","sd.MR_PRESSO","T-stat.MR_presso","pval.MR_PRESSO"
  )
  colnames(singleresults) <- newcolumns
  
  # Collect results
  combinedresults_list[[i]] <- singleresults
  
  # Single-SNP MR results (saved per-exposure)
  single_snp_res <- TwoSampleMR::mr_singlesnp(combo_dat)
  write.csv(single_snp_res, fressingle_name, row.names = FALSE)
}

# --- Combine & save all results ----------------------------------------------

# Keep only non-null entries
combinedresults_list <- Filter(Negate(is.null), combinedresults_list)
combinedresults_list <- lapply(combinedresults_list, function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)
  if ("pval.MR_PRESSO_global" %in% names(df)) {
    df$pval.MR_PRESSO_global <- suppressWarnings(as.character(df$pval.MR_PRESSO_global))
  }
  df
})
combinedresults_list <- Filter(Negate(is.null), combinedresults_list)

combinedresults_df <- dplyr::bind_rows(combinedresults_list)

if ("pval.method1" %in% names(combinedresults_df)) {
  combinedresults_df <- combinedresults_df[order(combinedresults_df$pval.method1), ]
}

write.csv(combinedresults_df, COMBINED_CSV, row.names = FALSE)
message("Combined results written to: ", COMBINED_CSV)
