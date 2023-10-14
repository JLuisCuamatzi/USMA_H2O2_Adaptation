library(data.table)
library(dplyr)
library(ggplot2)
rm(list = ls())
setwd("C:/Users/SilverStone/Documents/PaperToMBE/GitLab_Repository/02_MappingAndCoveragePlotting/CNV_by_CNVnator/cnvFiles/")

list.files()

df.01 <- fread("2021EE01.cnvnator")

df.18 <- fread("2021EE18.cnvnator")

names.df <- c("CNV_Type", "Coordinate", "Size_bp", "RD", "e-value")

df.01 <- df.01[,c(1:5)]
df.18 <- df.18[,c(1:5)]

names(df.01) <- names.df
names(df.18) <- names.df

library(writexl)

writexl::write_xlsx(x = df.01, path = "CNVnator.Results.SG200.xlsx", col_names = T)
writexl::write_xlsx(x = df.18, path = "CNVnator.Results.UmH2O2R.xlsx", col_names = T)

getwd()

# create ID
df.01$ID <- paste0(df.01$V1, df.01$V2)
df.18$ID <- paste0(df.18$V1, df.18$V2)

# remove CNV from SG200
df.18 <- df.18[ !(df.18$ID %in% df.01$ID), ]

# el numer de obs redce de 96 a 54
df.18
# graficar distribución de tamaño, col V
ggplot(data = df.18)+
  geom_density(aes(x = V3)) +
  scale_x_log10()

ggplot(data = df.18)+
  geom_density(aes(x = V4))
rm(list = ls())



