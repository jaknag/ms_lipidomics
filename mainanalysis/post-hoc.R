
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
exp1 <- 'filtered_GCST90451306_Triglycerides_in_VLDL.tsv'
exp2 <- 'filtered_GCST90451354_Triglycerides_in_chylomicrons_and_extremely_large_VLDL.tsv'

exp1_path <- file.path(EXPO_DIR, exp1)
exp2_path <- file.path(EXPO_DIR, exp2)


# --- Main loop -----------------------------------------------------------

dat <- readr::read_tsv(exp2_path)

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
  ) 

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

res <- TwoSampleMR::mr(combo_dat, method_list = c("mr_ivw_mre", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression"))

mr_presso <- TwoSampleMR::run_mr_presso(combo_dat, SignifThreshold = 0.1)
