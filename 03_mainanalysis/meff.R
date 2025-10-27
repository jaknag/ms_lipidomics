library(tidyverse)
library(poolr)

cor_matrix <- read.csv("/Users/nagrodzkij/cam/IMT/MS/lipidomics/literature/correlation_matrix_lipids2.csv")
cor_matrix <- subset(cor_matrix,select=-X)

max(cor_matrix)
min(cor_matrix)

#Adjust values (if >1, make it 1; if <1, make it -1)
cor_matrix[cor_matrix > 1] <- 1
cor_matrix[cor_matrix < -1] <- -1

eff_tests <- meff(cor_matrix, method="liji")

print(eff_tests)
