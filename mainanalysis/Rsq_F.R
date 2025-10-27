# Load necessary library
library(dplyr)

ms_dat = read_tsv("/Users/nagrodzkij/cam/IMT/MS/lipidomics/extended/severity/imsgc_mssev_discovery_hg38.tsv",col_types = "dcccddddddcc") %>%
  dplyr::rename("effect_allele" = A1,"other_allele"=A2,"beta"=BETA,"se"=SE,"rsID"=SNP,"eaf"=AF1,"pval"=P) %>%
  add_column(Phenotype="MS severity") %>%
  dplyr::mutate(SNP=str_sub(chrpos,4,-1))

# ms_dat = read_tsv("/Users/nagrodzkij/cam/IMT/MS/lipidomics/extended/risk/outcomes/ms_risk.tsv",col_types = "dcccdddddddc") %>% 
#   dplyr::rename("effect_allele" = A1,"other_allele"=A2,"beta"=BETA,"se"=SE) %>%
#   add_column(Phenotype="MS")

bfile_ref = "/Users/nagrodzkij/cam/IMT/MS/lipidomics/LD/working/backups/eur_cleaned_10feb"

outcome_dat = format_data(ms_dat %>% 
                            dplyr::select(-SNP) %>% 
                            dplyr::rename("SNP"=rsID, "pval"=P),
                          type="outcome")

save_location = "/Users/nagrodzkij/cam/IMT/MS/lipidomics/extended/severity/results/Rsq/"

# navigate to folder with all GWAS stats
setwd("/Users/nagrodzkij/cam/IMT/MS/lipidomics/extended/exposures")

# read in exposure data
gwas_raw_files = list.files(pattern="filtered",full.names = T)
count = length(gwas_raw_files)
results_df <- data.frame(exposure = character(), R2_total = numeric(), F_stat = numeric(), stringsAsFactors = FALSE)

i=0
for(file in gwas_raw_files){
  i=i+1
  message("iteration ", i, " out of ", length(gwas_raw_files))
  #Obtain exposure name from filename
  exposure_name = str_sub(file,12,-5)
  
  df = read_tsv(file)
  df = df %>% 
    filter(nchar(effect_allele)==1 & nchar(other_allele)==1) %>% 
    filter(chromosome %in% c(1:22)) %>% 
    filter(p_val <= 5e-8) %>% 
    filter(effect_allele_frequency >= 0.05 & effect_allele_frequency <= 0.95) %>%
    filter(INFO_UKBB_EUR > 0.9 & INFO_EstBB > 0.9) %>% 
    filter(
      !(effect_allele == "A" & other_allele == "T") &
        !(effect_allele == "T" & other_allele == "A") &
        !(effect_allele == "C" & other_allele == "G") &
        !(effect_allele == "G" & other_allele == "C")
    ) %>% 
    mutate(SNP = paste0(chromosome,":",base_pair_location)) %>% 
    add_column("Phenotype"=exposure_name)
  
  message("filtering done")
  
  # filter exposure SNPs to those in outcome GWAS based on chrom:pos & matching effect/non-effect alleles
  df = df %>% 
    filter(SNP %in% ms_dat$SNP) %>% 
    inner_join(ms_dat %>% 
                 dplyr::select(SNP,rsID,effect_allele,other_allele),
               by=c("SNP","effect_allele","other_allele"))
  
  message("matching on chr:pos done")
  
  # # format
  # exposure_dat = format_data(df %>% 
  #                              dplyr::select(-SNP) %>%
  #                              dplyr::rename("se"= standard_error,"pval"=p_val,"eaf"=effect_allele_frequency, "SNP"=rsID) %>% 
  #                              dplyr::select(SNP,everything()),
  #                            type = "exposure")
  
  clumped_dat1 = ld_clump(
    dplyr::tibble(rsid=df$rsID, pval=df$p_val),
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = bfile_ref,
    clump_p = 5e-8, 
    #clump_kb=50000, clump_r2 = 0.00001
  )
  
  df_clumped = df %>%
    filter(rsID %in% clumped_dat1$rsid)
  
  df_clumped_name = paste0(save_location,exposure_name,"_clumped",".csv")
  write.csv(df_clumped,df_clumped_name,row.names = FALSE)
  
  df_clumped = df_clumped %>%
    dplyr::rename("EAF"=effect_allele_frequency, "Beta"=beta, "Sample_Size"=n)
  
  # Compute per-SNP variance explained (R² for each SNP)
  # Assume variance = 1 as metabolites had inverse normal transformation applied 
  df_clumped$R2_SNP <- 2 * df_clumped$EAF * (1 - df_clumped$EAF) * (df_clumped$Beta^2)
  
  # Compute total R² (sum across all SNPs)
  R2_total <- sum(df_clumped$R2_SNP, na.rm = TRUE)
  
  # Get sample size (assuming the largest sample size is the correct one)
  N <- max(df_clumped$Sample_Size, na.rm = TRUE)
  
  # Count number of SNPs used as instruments
  k <- nrow(df_clumped)
  
  # Compute F-statistic
  F_stat <- (R2_total * (N - k - 1)) / (k * (1 - R2_total))
  
  results_df <- rbind(results_df, data.frame(exposure = exposure_name, R2_total = R2_total, F_stat = F_stat))
  
}


R2_Fstat_file = paste0(save_location,"_R2_Fstat",".csv")
write.csv(results_df, R2_Fstat_file, row.names=FALSE)

