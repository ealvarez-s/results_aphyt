
# Working directory
global_path<-"/Users/echo_user/GIT_locales/results_aphyt/"
setwd(global_path)

# Custom functions
source("Scripts/functions_misc.R")


## Required packages
library(abind)
library(akima)
library(chron)
library(Hmisc)
library(lmodel2)
library(Metrics)
library(RNetCDF)
library(plotrix)
library(plot3D)
library(unikn)      
library(viridisLite)
library(wesanderson)


# Script to format OASIM data as surface forcing,
# its output is available in Dat_WB_OASIM/binary/
#   01_Generate_OASIM_climatology.R

# Scripts to format observations, in situ collected and from satellite,
# their output is available in Dat_observations/ and Dat_satellite/
#   02_Insitu_Satellite_2008_2012.R
#   03_Insitu_HPLC_Chla.R
#   04_Insitu_BO_Valente_Casey_Chase_CSIRO_Astrid_modbands.R
#   05_Insitu_BO_Valente_Casey_Chase_CSIRO_Astrid_make_grids.R

# Scripts to format MITgcm-REcoM2 output,
# their output is available in Res_model/interpolated/
#   11_Interpolate1_2Dsurface_Chla_Rrs.R
#   11_Interpolate2_2Dderived_Zeu_NPP.R
#   11_Interpolate3_3Dglobal_Chla_IOP.R



# The next scripts reproduce all the results in the manuscript,
# the required input is fully available in the repository,
# and they should reproduce all figures available in Figures/

#   10_Optics1_phyto_REcoM_carbon.R
#     Figure 2

#   10_Optics2_Marshall_sensitivity_carbon.R
#     phytoplankton absorption .dat files for SA

#   10_Optics3_NOMAD_PPC.R
#     Figure 3

#   12_Figures1_Sensitivity_analysis.R
#     Figure S1
#     Figure S2
#     Figure S3

#   12_Figures2_Surface_Chla_Rrs_IOPs_PPC.R
#     Figure S4
#     Figure S5
#     Figure S6
#     Figure 4
#     Figure 5
#     Figure 7
#     Figure 8

#   12_Figures3_Depth_Chla_IOPs_PPC.R
#     Figure 6

#   13_Check_BioOp_Rships.R
#     Figure 9
#     Figure 10

#   14_Check_Spectral_Shapes.R
#     Figure 11
#     Figure 12
