library(tidyverse)
library(TwoSampleMR)
library(ieugwasr)
library(stringr)
library(ggplot2)

user() #Ensure this returns your details; if not, need to update access key 


#Ensure access key to OpenGWAS API is correct
#Go to https://api.opengwas.io/profile/
#Get the token
#Open terminal and type in 'open $HOME/.Renviron' 
#Update the key 


combinedresults = paste0("/Users/nagrodzkij/cam/IMT/MS/lipidomics/extended/severity/results/","_severity_joined",".csv")
results = read.csv(combinedresults)

#Obtain exposure name to plot for
exposure_name_extended = results$exposure[3] #Choose the row with the desired exposure

fharm_name = paste0("/Users/nagrodzkij/cam/IMT/MS/lipidomics/extended/severity/working/",exposure_name_extended,"_harmonise",".csv")
combo_dat = read.csv(fharm_name)

res_loo <- mr_leaveoneout(combo_dat, method= mr_ivw_mre)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

