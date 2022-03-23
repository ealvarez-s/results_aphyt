

if (Sys.info()['sysname']=="Windows") {OS<-"C:/"} else {OS<-paste("/Users/",Sys.info()['user'],"/", sep="")}

## Required packages
library(RNetCDF)
library(chron)
library(abind)
library(Hmisc)
library(Metrics)
library(akima)

library(plotrix)
library(plot3D)

library(unikn)      
library(viridisLite)
library(wesanderson)


# Custom functions
source("functions_misc.R")

# Working directory
global_path<-"/Users/echo_user/GIT_locales/results_aphyt/"



## 01_Generate_OASIM_climatology.R


## 02_Insitu_Satellite_2008_2012.R
## 03_Insitu_HPLC_Chla.R
## 04_Insitu_BO_Valente_Casey_Chase_CSIRO_Astrid_modbands.R
## 05_Insitu_BO_Valente_Casey_Chase_CSIRO_Astrid_make_grids.R


## 10_Optics1_phyto_REcoM_carbon.R
  # Figure 2
## 10_Optics2_Marshall_sensitivity_carbon.R

## 10_Optics3_NOMAD_PPC.R
  # Figure 3



## 11_Interpolate1_2Dsurface_Chla_Rrs.R

## 11_Interpolate2_2Dderived_Zeu_NPP.R

## 11_Interpolate3_3Dglobal_Chla_IOP.R


## 12_Figures1_Sensitivity_analysis.R
  # Figure S1
  # Figure S2
  # Figure S3

## 12_Figures2_Surface_Chla_Rrs_IOPs_PPC.R
  # Figure S4
  # Figure S5
  # Figure S6
  # Figure 4
  # Figure 5
  # Figure 7
  # Figure 8

## 12_Figures3_Depth_Chla_IOPs_PPC.R
  # Figure 6

## 13_Check_BioOp_Rships.R
  # Figure 9
  # Figure 10

## 14_Check_Spectral_Shapes.R
  # Figure 11
  # Figure 12












