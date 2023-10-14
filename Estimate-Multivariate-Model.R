#  Initialise  ----

rm(list=ls())
library(MTVGARCH)
library(knitr)

setwd("C:/Source/Repos/tv_betas")

# Get Data  ----
# Get all the final estimated univariate models & compile them into a multivariate ntvgarch Object

tvg_list = list()
tvg_list[[1]] = readRDS("Results/ANZ_Final_TVG_model.RDS")
tvg_list[[2]] = readRDS("Results/CBA_Final_TVG_model.RDS")
tvg_list[[3]] = readRDS("Results/NAB_Final_TVG_model.RDS")
tvg_list[[4]] = readRDS("Results/WBC_Final_TVG_model.RDS")
tvg_list[[5]] = readRDS("Results/STW_Final_TVG_model.RDS")
tvg_list[[6]] = readRDS("Results/SPY_Final_TVG_model.RDS")
tvg_list[[7]] = readRDS("Results/PR_Final_TVG_model.RDS")
tvg_list[[8]] = readRDS("Results/CRAU_Final_TVG_model.RDS")
tvg_list[[9]] = readRDS("Results/CRUS_Final_TVG_model.RDS")

# Make a multivaraite model object ----

tvg_names = c("ANZ","CBA","NAB","WBC","STW","SPY","PR","CRAU","CRUS")
tvBetas <- ntvgarch(tvg_list,tvg_names)


# Estimate the multivariate model ----

??
    