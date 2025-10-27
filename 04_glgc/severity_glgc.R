library(tidyverse)
library(TwoSampleMR)
library(ieugwasr)
library(stringr)

user() #Ensure this returns your details; if not, need to update access key 

#Ensure access key to OpenGWAS API is correct
#Go to https://api.opengwas.io/profile/
#Get the token
#Open terminal and type in 'open $HOME/.Renviron' 
#Update the key 

# Point to LD reference panel locally
bfile_ref = "/Users/nagrodzkij/cam/IMT/MS/lipidomics/LD/working/eur_cleaned_10feb"

# read in outcome GWAS, rename columns, add column 'SNP' with chr:pos 
ms_dat = read_tsv("/Users/nagrodzkij/cam/IMT/MS/lipidomics/extended/severity/imsgc_mssev_discovery_hg38.tsv",col_types = "dcccddddddcc") %>% 
  dplyr::rename("effect_allele" = A1,"other_allele"=A2,"beta"=BETA,"se"=SE,"rsID"=SNP,"eaf"=AF1,"pval"=P) %>%
  add_column(Phenotype="MS severity") %>%
  dplyr::mutate(SNP=str_sub(chrpos,4,-1))

# format outcome data
outcome_dat = format_data(ms_dat %>% 
                            dplyr::select(-SNP) %>% 
                            dplyr::rename("SNP"=rsID),type="outcome")

glgc_path = "/Users/nagrodzkij/cam/IMT/MS/lipidomics/extended/glgc/"

# navigate to folder with all GWAS stats
setwd(paste0(glgc_path, "exposures/"))

# read in exposure data
gwas_raw_files = list.files(pattern="filtered",full.names = T)

count = length(gwas_raw_files)
combinedresults_list = vector("list",length = count)

i=0
for(file in gwas_raw_files){
  i=i+1
  message("Iteration ", i)
  
  #Obtain exposure name from filename
  exposure_name = str_sub(file,12,-5)
  
  #Set filenames for saving csv files
  fexp_name = paste0(glgc_path,"exposures/",exposure_name,".csv")
  fharm_name = paste0(glgc_path,"severity/working/",exposure_name,"_harmonise",".csv")
  fres_name = paste0(glgc_path, "severity/results/",exposure_name,"_results",".csv")
  fressingle_name = paste0(glgc_path, "severity/results/",exposure_name,"_resultssingle",".csv")
  fresmrpresso_name = paste0(glgc_path, "severity/results/",exposure_name,"_mrpresso",".csv")
  combinedresults = paste0(glgc_path, "severity/results/","_combined",".csv")
  
  #Check if already done
  if (file.exists(fressingle_name)){
    print(paste('Results single already done for ', exposure_name, sep = ' '))} 
  else {
    
    # filter exposure summary stas
    dat = read_tsv(file)
    message("filtering ",file)
    # basic QC 
    dat = dat %>% 
      filter(nchar(effect_allele)==1 & nchar(other_allele)==1) %>% 
      filter(chromosome %in% c(1:22)) %>% 
      filter(p_value <= 5e-8) %>% 
      filter(effect_allele_frequency >= 0.05 & effect_allele_frequency <= 0.95) %>%
      # filter(INFO_UKBB_EUR > 0.9 & INFO_EstBB > 0.9) %>% 
      filter(
        !(effect_allele == "A" & other_allele == "T") &
          !(effect_allele == "T" & other_allele == "A") &
          !(effect_allele == "C" & other_allele == "G") &
          !(effect_allele == "G" & other_allele == "C")
      ) %>% 
      mutate(SNP = paste0(chromosome,":",base_pair_location)) %>% 
      dplyr::select(chromosome,base_pair_location,effect_allele,other_allele,beta,standard_error,effect_allele_frequency,SNP,
                    p_value) %>%
      add_column("Phenotype"=exposure_name)
    
    message("Filtering done")
    
    # filter exposure SNPs to those in outcome GWAS based on chrom:pos & matching effect/non-effect alleles
    dat = dat %>% 
      filter(SNP %in% ms_dat$SNP) %>% 
      inner_join(ms_dat %>% 
                   dplyr::select(SNP,rsID,effect_allele,other_allele),
                 by=c("SNP","effect_allele","other_allele"))
    
    message("Matching on chr:pos done")
    
    # format
    exposure_dat = format_data(dat %>% 
                                 dplyr::select(-SNP) %>%
                                 dplyr::rename("se"= standard_error,"pval"=p_value,"eaf"=effect_allele_frequency,"SNP"=rsID) %>% 
                                 dplyr::select(SNP,everything()),
                               type = "exposure")
    
    # clump 
    clumped_dat1 = ld_clump(
      dplyr::tibble(rsid=dat$rsID, pval=dat$p_value),
      plink_bin = genetics.binaRies::get_plink_binary(),
      bfile = bfile_ref,
      clump_p = 5e-8, 
      #clump_kb=50000, clump_r2 = 0.00001
    )
    
    exposure_dat = exposure_dat %>%
      filter(SNP %in% clumped_dat1$rsid)
    
    # harmonise 
    combo_dat = harmonise_data(exposure_dat,outcome_dat)
    write.csv(combo_dat,fharm_name, row.names = FALSE)
    message("Harmonisation done")
    
    # run basic MR 
    res = mr(combo_dat, method_list=c("mr_ivw_mre"))
    write.csv(res,fres_name, row.names = FALSE)
    
    # pull results from this exposure to later add to combined dataframe
    #row = which(res == "Inverse variance weighted (multiplicative random effects)", arr.ind=TRUE)[1,'row'] #Which row is the inverse variance weighted in? 
    #singleresults = res[row, ]
    singleresults = res[1, ]
    
    # plot basic results 
    # ggplot(res,
    #      aes(b,method))+
    #   geom_point(size=3)+
    #   geom_errorbarh(mapping = aes(xmin = b - 1.96*se,xmax = b+1.96*se,y=method),height=0.3)+
    #   theme_bw()+
    #   geom_vline(xintercept=0)+
    #   labs(x="MR estimate of exposure on MS risk")
    # 
    # heterogeneity 
    heterogeneity = mr_heterogeneity(combo_dat, method_list=c("mr_ivw_mre"))
    # pull heterogeneity results from this exposure
    singleresults = merge(singleresults,heterogeneity, by = c("id.exposure","id.outcome","outcome","exposure"))
    
    # pleitropy 
    pleiotropy = mr_pleiotropy_test(combo_dat)
    # pull pleiotropy from this exposure
    singleresults = merge(singleresults,pleiotropy, by = c("id.exposure","id.outcome","outcome","exposure"))
    
    # MR-PRESSO (pleiotropy robust method)
    mr_presso = run_mr_presso(combo_dat)
    mr_presso = mr_presso[[1]]
    mr_presso = mr_presso[[1]]
    write.csv(mr_presso,fresmrpresso_name, row.names = FALSE)
    mr_presso_raw = mr_presso[1,]
    mr_presso_outlier = mr_presso[2,] 
    # pull mr_presso results from this exposure
    singleresults = merge(singleresults, mr_presso_raw)
    
    #Change column names in combined dataframe to avoid repeating column names
    newcolumns = c("id.exposure","id.outcome","outcome","exposure","method1","nsnp.method1","b.method1","se.method1",
                   "pval.method1","method_heterogen","Q.heterogen","Q_df.heterogen","Q_pval.heterogen",
                   "egger_intercept","se.egger_intercept","pval.egger_intercept",
                   "Exposure","MR_PRESSO_analysis","b.MR_PRESSO","sd.MR_PRESSO","T-stat.MR_presso","pval.MR_PRESSO")
    
    newcolumns -> colnames(singleresults)
    
    
    singleresults = merge(singleresults, mr_presso_outlier)
    
    
    #Change column names in combined dataframe to avoid repeating column names
    newcolumns2 = c("id.exposure","id.outcome","outcome","outcome2","exposure",
                    "method1","nsnp.method1","b.method1","se.method1",
                    "pval.method1","method_heterogen","Q.heterogen","Q_df.heterogen","Q_pval.heterogen",
                    "egger_intercept","se.egger_intercept","pval.egger_intercept",
                    "MR_PRESSO_analysis","b.MR_PRESSO","sd.MR_PRESSO","T-stat.MR_presso","pval.MR_PRESSO",
                    "MR_PRESSO_analysis2", "b.MR_PRESSO2","sd.MR_PRESSO2", "T-stat.MR_presso2","pval.MR_PRESSO2") 
    
    newcolumns2 -> colnames(singleresults)
    
    # add results from this exposure to combined dataframe 
    combinedresults_list [[i]] = singleresults
    
    
    # single snp 
    single_snp_res = mr_singlesnp(combo_dat)
    write.csv(single_snp_res, fressingle_name, row.names=FALSE)
    
    # plot single snp results 
    # ggplot(combo_dat,
    #        aes(beta.exposure,beta.outcome))+
    #   geom_point()+
    #   theme_bw()+
    #   geom_vline(xintercept = 0)+
    #   geom_hline(yintercept=0)+
    #   geom_abline(slope = res[res$method=="Inverse variance weighted (multiplicative random effects)",]$b)
  }
}

columns = colnames(singleresults)
combinedresults_df = data.frame(matrix(ncol=length(columns), nrow=0)) #Create new df for combined results
colnames(combinedresults_df) <- columns

combinedresults_df <- dplyr::bind_rows(combinedresults_list)
combinedresults_df <- combinedresults_df[order(combinedresults_df$pval.method1),]

write.csv(combinedresults_df, combinedresults, row.names=FALSE)
