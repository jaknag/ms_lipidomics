library(tidyverse)
library(TwoSampleMR)
library(ieugwasr)
library(stringr)
library(MVMR)

rawdat = read_tsv("/Users/nagrodzkij/cam/IMT/MS/lipidomics/extended/mvmr/mvmr_bmi_vldl_raw.tsv")

# X1 = BMI; X2 = VLDL-TG; Y = severity MS

F.data <- format_mvmr(
  BXGs = rawdat[, c(4, 6)],
  BYG = rawdat[, 2],
  seBXGs = rawdat[, c(5, 7)],
  seBYG = rawdat[, 3],
  RSID = rawdat[, 1]
)
head(F.data)

sres <- strength_mvmr(r_input = F.data, gencov = 0)

pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)

res <- ivw_mvmr(r_input = F.data)

res1 <- qhet_mvmr(F.data, CI = TRUE, iterations = 100)
