source("00_Summary_EDITME.R")
source("functions_misc.R")

## Sensitivity analysis and metrics
## Chla, Reflectance, IOP's, PPC

############################################
## METRICS CORRELATION with SATELITE/In situ
############################################

o_dir<-paste(global_path,"Res_model/",sep="")        
p_dir<-paste(global_path,"Res_model/interpolated/",sep="")
path_figures<-paste(global_path,"Figures/",sep="")
# List of simulations
experimentos<-read.csv(paste(o_dir,"run_log_marshall_PPC_PS_SA_G100.csv",sep=""), sep=",")

        modelnames<-experimentos$name[c(10:29,9)]
        ruta<-o_dir   #paste(OS,substring(experimentos$path[c(10:29,9)],17,100), sep="")
        puntos<-rep(21,length(modelnames))
        nombre <- substring(modelnames,19,100)
        simbolos<-letters[seq( from = 1, to = 20 )] 
        #nombre <- c(expression(paste("EXP-1: ",{{"a*"}}[PH],"(", lambda, ") constant",sep="")),
        #expression(paste("EXP-2: ",{{"a*"}}[PH],"(", lambda, ") variable",sep="")))
    
        
        
## Reference values for global(aph) and global(aps)
###################        
                pdir<-paste(global_path,"Phyto_optics/",sep="") 
                wd <- c(12.5,rep(25,11),12.5)
                apps<-rep(NA, length=length(modelnames))
                apph<-rep(NA, length=length(modelnames))
                
                MEDIA_apph_dia <-rep(NA, length=length(modelnames))
                MEDIA_apps_dia <-rep(NA, length=length(modelnames))
                MEDIA_qy_dia   <-rep(NA, length=length(modelnames))
                PEAK_apph_dia <-rep(NA, length=length(modelnames))
                PEAK_apps_dia <-rep(NA, length=length(modelnames))                

                MEDIA_apph_phy <-rep(NA, length=length(modelnames))
                MEDIA_apps_phy <-rep(NA, length=length(modelnames))
                MEDIA_qy_phy   <-rep(NA, length=length(modelnames))
                PEAK_apph_phy <-rep(NA, length=length(modelnames))
                PEAK_apps_phy <-rep(NA, length=length(modelnames))                 
                                
                for (w in 1:length(modelnames)){ 
                      
             ### Diatoms
                      alpha <- 0.19/86400
                      if (w==21) {
                      filename<-paste(pdir,"optics_phyto_recom_carbon.dat",sep="")  
                      } else{
                      filename<-paste(pdir,"sensitivity/optics_phyto_recom_carbon_",w,".dat",sep="")}
                      datos <- read.table(filename, skip=21)
                      lambda <- datos[,1]
                      # AP_PS
                      dat <- datos[,3]
                      media_apps_dia <- sum(wd*dat)/(700-400)
                      media_qy_dia <- alpha / media_apps_dia
                          MEDIA_apps_dia[w] <- media_apps_dia 
                          PEAK_apps_dia[w] <- dat[lambda=="450"]
                          MEDIA_qy_dia[w] <- media_qy_dia
                        
                      # AP
                      dat <- datos[,2]
                      media_apph_dia <- sum(wd*dat)/(700-400)
                          MEDIA_apph_dia[w] <- media_apph_dia 
                          PEAK_apph_dia[w] <- dat[lambda=="450"]
                      
             ### Small Phyto
                      alpha <- 0.14/86400
                      if (w==21) {
                      filename<-paste(pdir,"optics_phyto_recom_carbon.dat",sep="")  
                      } else{
                      filename<-paste(pdir,"sensitivity/optics_phyto_recom_carbon_",w,".dat",sep="")
                      }
                      datos <- read.table(filename, skip=7, nrows=13)
                      lambda <- datos[,1]
                      # AP_PS
                      dat <- datos[,3]         
                      media_apps_phy <- sum(wd*dat)/(700-400)
                      media_qy_phy <- alpha /  media_apps_phy
                      
                          MEDIA_apps_phy[w] <- media_apps_phy 
                          PEAK_apps_phy[w] <- dat[lambda=="450"]
                          MEDIA_qy_phy[w] <- media_qy_phy                      
                      
                      # AP
                      dat <- datos[,2]
                      media_apph_phy <- sum(wd*dat)/(700-400)
                          MEDIA_apph_phy[w] <- media_apph_phy 
                          PEAK_apph_phy[w] <- dat[lambda=="450"]
                          
                      # MODEL
                      modelname<-modelnames[w]
                      #o_dir<-ruta[w]              
                      load(paste(p_dir,"2cloro_dia/",modelname,".RData",sep=""))
                      diatoms <- res
                      #dim(diatoms)
                      load(paste(p_dir,"2cloro_phy/",modelname,".RData",sep=""))
                      smallphy <- res
                      #dim(smallphy)                  
                      m<-(media_apps_dia*diatoms+media_apps_phy*smallphy)/(diatoms+smallphy)
                      apps[w]<-round(mean(m,na.rm=T),3)
                      m<-(media_apph_dia*diatoms+media_apph_phy*smallphy)/(diatoms+smallphy)
                      apph[w]<-round(mean(m,na.rm=T),3)
                      }
#######################
        
      tablaS1<-data.frame(c(simbolos,"21"),PEAK_apph_phy,MEDIA_apph_phy,PEAK_apps_phy,MEDIA_apps_phy,MEDIA_qy_phy,
                                   PEAK_apph_dia,MEDIA_apph_dia,PEAK_apps_dia,MEDIA_apps_dia,MEDIA_qy_dia, apps,  apph)    
      write.csv(tablaS1, file=paste(pdir,"table_S1.csv",sep=""))       
                
      
      
#############
#### CHLA log
#############        
        R_chlal_ist   <-rep(NA,length=length(modelnames))
        AE_chlal_ist  <-rep(NA,length=length(modelnames))
        n_chlal_ist  <-rep(NA,length=length(modelnames))
        R_chlal_sat   <-rep(NA,length=length(modelnames))
        AE_chlal_sat  <-rep(NA,length=length(modelnames))
        n_chlal_sat  <-rep(NA,length=length(modelnames))

########################         
        ##### In situ 
        #### MAREDAT and others
        is_dir <- paste(global_path,"Dat_observations/HPLC/grids_pigments",sep="")  
        load(paste(is_dir,"/TChla_insitu_new2.RData",sep=""))                           
        longitude<-as.numeric(dimnames(TChla_ALL)$x)
        latitude<-as.numeric(dimnames(TChla_ALL)$y)        
        matriz1  <- (apply(TChla_ALL[,,1:2], MARGIN=c("x","y"), mean, na.rm=TRUE))/1e+3 # ng L-1 to ug L-1
        ##### In situ Valente and others
        is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="") 
        load(paste(is_dir,"/media_annual_Chla.Rdata",sep=""))
        matriz2  <- apply(res[,,1:2], MARGIN=c("x","y"), mean, na.rm=TRUE) # ug L-1
        nueva<-abind(matriz1, matriz2, along=0.5)
        matriz<-apply(nueva,MARGIN=c(2,3), mean, na.rm=TRUE)        
        
        # SAT
        o_dir <- paste(global_path,"Dat_observations/satellite/",sep="")
        load(paste(o_dir,"climatologies_2012_2018/media_anual_log_chl_OCCCI_2012_2018.RData", sep=""))
        media <- apply(res, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
   
    for (w in 1:length(modelnames)){ 
        modelname<-modelnames[w]
        #o_dir<-ruta[w]              
        # MODEL
        load(paste(p_dir,"1cloro_log/",modelname,".RData",sep=""))
        mod <- res
        media[is.na(mod)] <- NA
        
        # Correlation In situ
        x<-cbind(c(mod), c(log10(matriz)))
        x[x=="-Inf"] <- NA
        x<-x[complete.cases(x)==TRUE,]
        R_chlal_ist[w]<- rcorr(x, type="pearson")$r[1,2]
        AE_chlal_ist[w]<- (sum(x[,1]-x[,2]))/nrow(x)        
        n_chlal_ist[w]<-nrow(x)
          
        # Correlation Sat
        x<-cbind(c(mod), c(media))
        x[x=="-Inf"] <- NA
        x<-x[complete.cases(x)==TRUE,]
        R_chlal_sat[w]<- rcorr(x, type="pearson")$r[1,2]
        AE_chlal_sat[w]<- (sum(x[,1]-x[,2]))/nrow(x)
        n_chlal_sat[w]<- nrow(x)
    }  # end loop model
########################    

  
     
###############    
## REFLECTANCE
###############
        R_rrs400   <-rep(NA,length=length(modelnames))
        R_rrs450   <-rep(NA,length=length(modelnames))
        R_rrs475   <-rep(NA,length=length(modelnames))
        R_rrs500   <-rep(NA,length=length(modelnames))
        R_rrs550   <-rep(NA,length=length(modelnames))
        R_rrs675   <-rep(NA,length=length(modelnames))
            R_rrs412   <-rep(NA,length=length(modelnames))
            R_rrs443   <-rep(NA,length=length(modelnames))
            R_rrs490   <-rep(NA,length=length(modelnames))
            R_rrs510   <-rep(NA,length=length(modelnames))
            R_rrs555   <-rep(NA,length=length(modelnames))
            R_rrs670   <-rep(NA,length=length(modelnames))
        n_rrs400   <-rep(NA,length=length(modelnames))
        n_rrs450   <-rep(NA,length=length(modelnames))
        n_rrs475   <-rep(NA,length=length(modelnames))
        n_rrs500   <-rep(NA,length=length(modelnames))
        n_rrs550   <-rep(NA,length=length(modelnames))
        n_rrs675   <-rep(NA,length=length(modelnames))
            n_rrs412   <-rep(NA,length=length(modelnames))
            n_rrs443   <-rep(NA,length=length(modelnames))
            n_rrs490   <-rep(NA,length=length(modelnames))
            n_rrs510   <-rep(NA,length=length(modelnames))
            n_rrs555   <-rep(NA,length=length(modelnames))
            n_rrs670   <-rep(NA,length=length(modelnames))           
#####################            
             
        # MASCARA
        modelname<-modelnames[3]
        load(paste(p_dir,"1reflectance/",modelname,".RData",sep=""))
        mod <- res
        z <- mod[,,12]

          # SAT
          o_dir <- paste(global_path,"Dat_observations/satellite/",sep="")
          landas_sat <- c("Rrs_412","Rrs_443","Rrs_490","Rrs_510","Rrs_555","Rrs_670")
          load(paste(o_dir,"climatologies_2012_2018/media_seasonal_Rrs_2deg_OCCCI_2012_2018.RData",sep=""))
          res <- med
          
          media412 <- apply(res[,,,1], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
          media412[is.na(z)] <- NA
          
          media443 <- apply(res[,,,2], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
          media443[is.na(z)] <- NA
          
          media490 <- apply(res[,,,3], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
          media490[is.na(z)] <- NA
          
          media510 <- apply(res[,,,4], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
          media510[is.na(z)] <- NA
          
          media555 <- apply(res[,,,5], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
          media555[is.na(z)] <- NA
          
          media670 <- apply(res[,,,6], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
          media670[is.na(z)] <- NA
        
   
       ##### In situ Rrs
       is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
       load(paste(is_dir,"/media_annual_Rrs.Rdata",sep=""))           
       #dim(matriz_annual)
       #lambda

        for (w in 1:length(modelnames)){  
                    modelname<-modelnames[w]
                            # MODEL
                            load(paste(p_dir,"1reflectance/",modelname,".RData",sep=""))
                            mod <- res
                            z <- mod[,,12]

                            # Correlation 443 vs 450
                            x<-cbind(c(mod[,,1]), c(media412))
                            x[x=="-Inf"] <- NA
                            x<-x[complete.cases(x)==TRUE,]
                            R_rrs412[w]<- rcorr(x, type="pearson")$r[1,2  ]
                            n_rrs412[w]<-nrow(x) 
                              
                            # Correlation 443 vs 450
                            x<-cbind(c(mod[,,3]), c(media443))
                            x[x=="-Inf"] <- NA
                            x<-x[complete.cases(x)==TRUE,]
                            R_rrs443[w]<- rcorr(x, type="pearson")$r[1,2  ]
                            n_rrs443[w]<-nrow(x) 
                            
                            # Correlation 490 vs 475
                            x<-cbind(c(mod[,,4]), c(media490))
                            x[x=="-Inf"] <- NA
                            x<-x[complete.cases(x)==TRUE,]
                            R_rrs490[w]<- rcorr(x, type="pearson")$r[1,2  ]
                            n_rrs490[w]<-nrow(x) 
                            
                            # Correlation 510 vs 500
                            x<-cbind(c(mod[,,5]), c(media510))
                            x[x=="-Inf"] <- NA
                            x<-x[complete.cases(x)==TRUE,]
                            R_rrs510[w]<- rcorr(x, type="pearson")$r[1,2  ]
                            n_rrs510[w]<-nrow(x) 
                            
                            # Correlation 555 vs 550
                            x<-cbind(c(mod[,,7]), c(media555))
                            x[x=="-Inf"] <- NA
                            x<-x[complete.cases(x)==TRUE,]
                            R_rrs555[w]<- rcorr(x, type="pearson")$r[1,2]
                            n_rrs555[w]<-nrow(x) 
                            
                            # Correlation 670 vs 675
                            x<-cbind(c(mod[,,12]), c(media670))
                            x[x=="-Inf"] <- NA
                            x<-x[complete.cases(x)==TRUE,]
                            R_rrs670[w]<- rcorr(x, type="pearson")$r[1,2]
                            n_rrs670[w]<-nrow(x) 

                            # Correlation 450 vs 450
                            x<-cbind(c(mod[,,1]), c(matriz_annual[,,1]))
                            x[x=="-Inf"] <- NA
                            x<-x[complete.cases(x)==TRUE,]
                            R_rrs400[w]<- rcorr(x, type="pearson")$r[1,2  ]
                            n_rrs400[w]<-nrow(x) 
                            
                            # Correlation 443 vs 450
                            x<-cbind(c(mod[,,3]), c(matriz_annual[,,3]))
                            x[x=="-Inf"] <- NA
                            x<-x[complete.cases(x)==TRUE,]
                            R_rrs450[w]<- rcorr(x, type="pearson")$r[1,2  ]
                            n_rrs450[w]<-nrow(x) 
                            
                            # Correlation 490 vs 475
                            x<-cbind(c(mod[,,4]), c(matriz_annual[,,4]))
                            x[x=="-Inf"] <- NA
                            x<-x[complete.cases(x)==TRUE,]
                            R_rrs475[w]<- rcorr(x, type="pearson")$r[1,2  ]
                            n_rrs475[w]<-nrow(x) 
                            
                            # Correlation 510 vs 500
                            x<-cbind(c(mod[,,5]), c(matriz_annual[,,5]))
                            x[x=="-Inf"] <- NA
                            x<-x[complete.cases(x)==TRUE,]
                            R_rrs500[w]<- rcorr(x, type="pearson")$r[1,2  ]
                            n_rrs500[w]<-nrow(x) 
                            
                            # Correlation 555 vs 550
                            x<-cbind(c(mod[,,7]), c(matriz_annual[,,7]))
                            x[x=="-Inf"] <- NA
                            x<-x[complete.cases(x)==TRUE,]
                            R_rrs550[w]<- rcorr(x, type="pearson")$r[1,2]
                            n_rrs550[w]<-nrow(x) 
                            
                            # Correlation 670 vs 675
                            x<-cbind(c(mod[,,12]), c(matriz_annual[,,12]))
                            x[x=="-Inf"] <- NA
                            x<-x[complete.cases(x)==TRUE,]
                            R_rrs675[w]<- rcorr(x, type="pearson")$r[1,2] 
                            n_rrs675[w]<-nrow(x) 
                            }
#####################        
        

    
#################
## IOP's surface
################ 
       
       R_ap     <-rep(NA,length=length(modelnames))
       R_aphs   <-rep(NA,length=length(modelnames))
       R_anap   <-rep(NA,length=length(modelnames))
       R_acdom  <-rep(NA,length=length(modelnames))
         n_ap     <-rep(NA,length=length(modelnames))
         n_aphs   <-rep(NA,length=length(modelnames))
         n_anap   <-rep(NA,length=length(modelnames))
         n_acdom  <-rep(NA,length=length(modelnames))         
## Sat       
        R_tot   <-rep(NA,length=length(modelnames))
        R_aph   <-rep(NA,length=length(modelnames))
        R_adg   <-rep(NA,length=length(modelnames))
        R_bbp   <-rep(NA,length=length(modelnames))        
          n_tot   <-rep(NA,length=length(modelnames))
          n_aph   <-rep(NA,length=length(modelnames))
          n_adg   <-rep(NA,length=length(modelnames))
          n_bbp   <-rep(NA,length=length(modelnames)) 
        
### Atot
########        
     # SAT
     o_dir <- paste(global_path,"Dat_observations/satellite/",sep="")
     load(paste(o_dir,"climatologies_2012_2018/media_anual_Atot_2deg_OCCCI_2012_2018.RData",sep=""))
     landa_sat <- dimnames(res)$l
     sat <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)        
     media <- sat[,,2]
        
     for (w in 1:length(modelnames)){ 
        modelname<-modelnames[w]
        load(paste(p_dir,"2iop_atot/",modelname,".RData",sep=""))
        mod <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
        #dim(mod)
        z <- mod
        z[z<0.01] <- NA
        media[is.na(z)] <- NA
        # Correlation Sat
        x<-cbind(c(z), c(media))
        x[x=="-Inf"] <- NA
        x<-x[complete.cases(x)==TRUE,]
        R_tot[w]<- rcorr(x, type="pearson")$r[1,2]
        n_tot[w]<- nrow(x)
        }
#######        
      
### Aph
########
       # In situ
       is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="") 
       load(paste(is_dir,"/media_annual_Aph.Rdata",sep=""))           
     
       # SAT
       o_dir <- paste(global_path,"Dat_observations/satellite/",sep="")
       load(paste(o_dir,"climatologies_2012_2018/media_anual_Aph_2deg_OCCCI_2012_2018.RData",sep=""))
       landa_sat <- dimnames(res)$l
       sat <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)        
       media <- sat[,,2]

        for (w in 1:length(modelnames)){ 
        modelname<-modelnames[w]
        load(paste(p_dir,"2iop_aph/",modelname,".RData",sep=""))
        mod <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
        z <- mod
        z[z<1e-5] <- NA
        media[is.na(z)] <- NA
        
        # Correlation Sat
        x<-cbind(c(z), c(media))
        x[x=="-Inf"] <- NA
        x<-x[complete.cases(x)==TRUE,]
        R_aph[w]<- rcorr(x, type="pearson")$r[1,2]
        n_aph[w]<- nrow(x)    
        
        # Correlation in situ
        x<-cbind(c(z), c(matriz_annual[,,1,3]))
        x[x=="-Inf"] <- NA
        x<-x[complete.cases(x)==TRUE,]
        R_aphs[w]<- rcorr(x, type="pearson")$r[1,2]
        n_aphs[w]<- nrow(x)
        }        
#########     
       
### Anap
########       
       # In situ
       is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="") 
       load(paste(is_dir,"/media_annual_Anap.Rdata",sep=""))           

       for (w in 1:length(modelnames)){ 
         modelname<-modelnames[w]
         load(paste(p_dir,"2iop_apt/",modelname,".RData",sep=""))
         mod <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
         z <- mod
         z[z<1e-5] <- NA
         # Correlation in situ
         x<-cbind(c(z), c(matriz_annual[,,1,3]))
         x[x=="-Inf"] <- NA
         x<-x[complete.cases(x)==TRUE,]
         R_anap[w]<- rcorr(x, type="pearson")$r[1,2]    
         n_anap[w]<- nrow(x)
       }               
#########       
       
### Acdom
#########       
       # In situ
       is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="") 
       load(paste(is_dir,"/media_annual_Acdom.Rdata",sep=""))           
       
       for (w in 1:length(modelnames)){ 
         modelname<-modelnames[w]
         load(paste(p_dir,"2iop_acd/",modelname,".RData",sep=""))
         mod <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
         z <- mod
         z[z<1e-5] <- NA
         # Correlation in situ
         x<-cbind(c(z), c(matriz_annual[,,1,3]))
         x[x=="-Inf"] <- NA
         x<-x[complete.cases(x)==TRUE,]
         R_acdom[w]<- rcorr(x, type="pearson")$r[1,2]  
         n_acdom[w]<- nrow(x)
       }          
##########
       
### Ap
##########       
       # In situ
       is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="") 
       load(paste(is_dir,"/media_annual_Ap.Rdata",sep=""))           
       
       for (w in 1:length(modelnames)){ 
         modelname<-modelnames[w]
         load(paste(p_dir,"2iop_app/",modelname,".RData",sep=""))
         mod <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
         z <- mod
         z[z<1e-5] <- NA
         # Correlation in situ
         x<-cbind(c(z), c(matriz_annual[,,1,3]))
         x[x=="-Inf"] <- NA
         x<-x[complete.cases(x)==TRUE,]
         R_ap[w]<- rcorr(x, type="pearson")$r[1,2]  
         n_ap[w]<- nrow(x)
       }        
###########         
        
### Adg
##########       
         # SAT
         o_dir <- paste(global_path,"Dat_observations/satellite/",sep="")
         load(paste(o_dir,"climatologies_2012_2018/media_anual_Adg_2deg_OCCCI_2012_2018.RData",sep=""))
         landa_sat <- dimnames(res)$l
         sat <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)        
         media <- sat[,,2]
            for (w in 1:length(modelnames)){ 
            modelname<-modelnames[w]
            load(paste(p_dir,"2iop_adg/",modelname,".RData",sep=""))
            mod <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
            #dim(mod)
            z <- mod
            z[z<1e-5] <- NA
            media[is.na(z)] <- NA
            # Correlation
            x<-cbind(c(z), c(media))
            x[x=="-Inf"] <- NA
            x<-x[complete.cases(x)==TRUE,]
            R_adg[w]<- rcorr(x, type="pearson")$r[1,2]
            n_adg[w]<- nrow(x)
            }        
##########
         
### Bbp
###########         
           # SAT
           o_dir <- paste(global_path,"Dat_observations/satellite/",sep="")
           load(paste(o_dir,"climatologies_2012_2018/media_anual_Bbp_2deg_OCCCI_2012_2018.RData",sep=""))
           landa_sat <- dimnames(res)$l
           sat <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)
           media <- sat[,,2]
           media[is.na(z)] <- NA
                  for (w in 1:length(modelnames)){ 
                  modelname<-modelnames[w]
                  load(paste(p_dir,"2iop_bbp/",modelname,".RData",sep=""))
                  mod <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
                  z <- mod
                  media[is.na(z)] <- NA
                  # Correlation
                  x<-cbind(c(z), c(media))
                  x[x=="-Inf"] <- NA
                  x<-x[complete.cases(x)==TRUE,]
                  R_bbp[w]<- rcorr(x, type="pearson")$r[1,2]
                  n_bbp[w]<- nrow(x)
                  }        
##########
           
##################        
        
  

########        
#### PPC        
########        

        # In situ 
        is_dir <- paste(global_path,"Dat_observations/HPLC/grids_pigments",sep="")  
        load(paste(is_dir,"/PPCinsitu_new3.RData",sep=""))
        matriz1  <- apply(PPCinsitu_ALL[,,1:2], MARGIN=c("x","y"), mean, na.rm=TRUE)

        ## CORRELATIONS
        R_PPC    <- rep(NA,length(modelnames)) 
        AE_PPC   <- rep(NA,length(modelnames)) 
        R_PPC_s  <- rep(NA,length(modelnames)) 
        AE_PPC_s <- rep(NA,length(modelnames)) 
        n_PPC   <- rep(NA,length(modelnames)) 
########################        
  for (w in 1:length(modelnames)){ 
          modelname<-modelnames[w]
          o_dir<-ruta[w]              
          # MODEL
          load(paste(p_dir,"1ppc_total/",modelname,".RData",sep=""))
          mod <- res
          # Correlation
          matriz1[matriz1=="Inf" | matriz1=="-Inf"] <- NA
          x<-cbind(c(mod), c(matriz1))
          x[x=="-Inf"] <- NA
          x<-x[complete.cases(x)==TRUE,]
            R_PPC[w]<- rcorr(x, type="pearson")$r[1,2]
            AE_PPC[w]<- (sum(x[,1]-x[,2]))/nrow(x)
            n_PPC[w]<- nrow(x)
      }  # end loop model

########################    
        
        
     
        
        
#################################
###  FIGURE S1  Spectra for EXP-S
#################################          

    png(file=paste(path_figures,"FigureS1_sensitivity_analysis.png", sep=""), width = 1200, height = 420, pointsize=18)
        #layout(matrix(c(1:5,5), ncol=3, byrow=FALSE))
        par(mfrow=c(1,3))    
        par(mar=c(5,5,1,1)) 
        par(oma=c(0,0,2,2.5))
        layout.show(n=3)
        lettersize=2
        sizelab=1.4
        sizeaxis=1.2
          oscuros<-c(usecol(pal_petrol,n=20))
          claros <-c(usecol(pal_bordeaux,n=20))      
          mycol<-usecol(pal = c(rev(pal_karpfenblau), pal_peach),n=20, alpha = 0.75)    
          mycol<-usecol(pal = pal_grau, n=20, alpha = 0.75)  
          simbolos<-letters[seq( from = 1, to = 20 )]         
        pdir<-paste(OS,"Datos/Res_C20_radtrans/",sep="")
        wd <- c(12.5,rep(25,11),12.5)
        
### A Small Phyto
#################        
        plot(1,1, las=1, type="b", xlim=c(400,700), ylim=c(0,0.1), pch=19, lty=1, lwd=2,
             col="cyan2", ylab=expression(paste("a*"," (",m^2, " mg ", Chl^-1,")")),
             xlab=expression(lambda), mgp=c(3,1,0), yaxs="i",
             cex.lab=sizelab, cex.axis=sizeaxis) 
        
        for (i in c(1:20)){
          filename<-paste(global_path,"Phyto_optics/sensitivity/optics_phyto_recom_carbon_",i,".dat",sep="")
          datos <- read.table(filename, skip=7, nrows=13)
          lambda <- datos[,1]
          # AP
          dat <- datos[,2]
          points(lambda,dat, las=1, type="l", pch=19, lty=1, lwd=2, col=oscuros[i])
          if (i==12) {points(lambda,dat, las=1, type="l", pch=19, lty=1, lwd=6, col="black")}
          # AP_PS
          dat <- datos[,3]
          points(lambda,dat, las=1, type="l", pch=19, lty=2, lwd=2, col=oscuros[i])
          if (i==12) {points(lambda,dat, las=1, type="l", pch=19, lty=3, lwd=6, col="black")}
          }
        mtext(3,text="a)",at=400,line=0.5,cex=lettersize, font=1)
        
        legend(x=550,y=0.06,legend=c(expression("a*"[PH](lambda)),
                                     expression("a*"[PS](lambda))),lty=c(1,2), col=oscuros[10], bty="n", cex=sizelab)                  
        
        ##ESCALA
        par(new=T)
        par(mar=c(25.5,7,1.5,2))
        mat <- matrix(c(1:20), nrow=20)
        image(mat, col=oscuros, las=1, xaxt="n", yaxt="n", main="", cex.main=0.8) 
        axis(1,at=seq(0,1, length=20), labels=simbolos, cex.axis=0.7, las=1)
#################        
        
### B Diatoms
#############        
        par(mar=c(5,5,1,1)) 
        plot(1,1, las=1, type="b", xlim=c(400,700), ylim=c(0,0.1), pch=19, lty=1, lwd=2,
             col="cyan2", ylab=expression(paste("a*"," (",m^2, " mg ", Chl^-1,")")),
             xlab=expression(lambda), mgp=c(3,1,0), yaxs="i",
             cex.lab=sizelab, cex.axis=sizeaxis) 
        
        for (i in c(1:20)){
          filename<-paste(global_path,"Phyto_optics/sensitivity/optics_phyto_recom_carbon_",i,".dat",sep="")
          datos <- read.table(filename, skip=21)
          lambda <- datos[,1]
          # AP
          dat <- datos[,2]
          points(lambda,dat, las=1, type="l", pch=19, lty=1, lwd=2, col=claros[i])
          if (i==12){points(lambda,dat, las=1, type="l", pch=19, lty=1, lwd=6, col="black")}
          # AP_PS
          dat <- datos[,3]
          points(lambda,dat, las=1, type="l", pch=19, lty=2, lwd=2, col=claros[i])
          if (i==12){points(lambda,dat, las=1, type="l", pch=19, lty=3, lwd=6, col="black")}
          }
        mtext(3,text="b)",at=400,line=0.5, cex=lettersize, font=1)
        legend(x=550,y=0.06,legend=c(expression("a*"[PH](lambda)),
                                     expression("a*"[PS](lambda))),lty=c(1,2), col=claros[10], bty="n", cex=sizelab)                  
        
        #legend(x="topleft",legend=c(paste("<aPS> = ", media22,sep=""),paste("<aPS> = ",round(media2,3), sep="")),text.col=c("red2","cyan2"), bty = "n")
        #legend(x="topright",legend=c(paste("QY = ", formatC(QY22[1],format="e",digits=2),sep=""),paste("QY = ", formatC(QY22[2],format="e",digits=2),sep="")),text.col=c("red2","cyan2"), bty = "n")    
        
        ##ESCALA
        par(new=T)
        par(mar=c(25.5,7,1.5,2))
        mat <- matrix(c(1:20), nrow=20)
        image(mat, col=claros, las=1, xaxt="n", yaxt="n", main="", cex.main=0.8) 
        axis(1,at=seq(0,1, length=20), labels=simbolos, cex.axis=0.7, las=1)
#############        
        
        
#### C ALL combinations alpha
#############################        
        par(mar=c(5,5,1,1)) 
        ab <- seq(0.0030,0.041, length=100)    
        qy <- seq(2.1e-5, 4.0e-4, length=100)  
        res <- matrix(c(rep(ab, times=length(qy))*rep(qy,each=length(ab))*86400),ncol=length(ab), byrow=TRUE)
        colnames(res) <- ab
        rownames(res) <- qy
        #range(res,na.rm=T)
        res[res>=0.9]<-0.9
        res[res<=0.009]<-0.009
        
        image2D(x=qy, y=ab, z=res, las=1,
                xlab=expression(paste(phi," (mmolC ", J^-1,")")),
                ylab=expression(paste(bar("a*")," (",m^2, " mg ", Chl^-1,")")),
                zlim=c(0.009,0.9), mgp=c(3.1,0.8,0), col=mycol, colkey=T,
                cex.lab=sizelab, cex.axis=sizeaxis)
        contour(x=qy, y=ab, z=res, levels=c(0.14), col=oscuros[10], add=TRUE, lty=1, lwd=3, labcex=2.0, font=2)
        contour(x=qy, y=ab, z=res, levels=c(0.19), col=claros[10], add=TRUE, lty=1, lwd=3, labcex=2.0, font=2)
        mtext(3,text="c)",at=3e-5,line=0.5, cex=lettersize, font=1)
        mtext(4,text=expression(paste(alpha," (",m^2," mmolC ", " mg", Chl^-1, J^-1," x 86400s/d)")),
              line=2,at=0.022,cex=0.9)
        
        load(paste(global_path,"Phyto_optics/tabla_total_sensitivity.Rdata",sep=""))     
        #tabla_total

        for (i in c(1:20)){
          ### Diatoms
          alpha <- 0.19/86400
          filename<-paste(global_path,"Phyto_optics/sensitivity/optics_phyto_recom_carbon_",i,".dat",sep="")
          datos <- read.table(filename, skip=21)
          lambda <- datos[,1]
          # AP_PS
          dat <- datos[,3]
          media1 <- sum(wd*dat)/(700-400)
          QYi <- alpha / media1
          points(x=QYi, y=media1, pch=21, cex=0.75,  bg=claros[i], col=claros[i])          
          # AP
          dat <- datos[,2]
          media2 <- sum(wd*dat)/(700-400)
          points(x=QYi, y=media2, pch=22, cex=0.75,  bg=claros[i], col=claros[i])          
          ### Small Phyto
          alpha <- 0.14/86400
          filename<-paste(global_path,"Phyto_optics/sensitivity/optics_phyto_recom_carbon_",i,".dat",sep="")
          datos <- read.table(filename, skip=7, nrows=13)
          lambda <- datos[,1]
          # AP_PS
          dat <- datos[,3]         
          media1 <- sum(wd*dat)/(700-400)
          QYi <- alpha / media1
          points(x=QYi, y=media1, pch=21, cex=0.75, bg=oscuros[i], col=oscuros[i])            
          # AP
          dat <- datos[,2]
          media2 <- sum(wd*dat)/(700-400)
          points(x=QYi, y=media2, pch=22, cex=0.75,  bg=oscuros[i], col=oscuros[i])      
        }
             # Initial spectra: 12, l 
                  i=12
                  # Diatoms
                  alpha <- 0.19/86400
                  filename<-paste(global_path,"Phyto_optics/sensitivity/optics_phyto_recom_carbon_",i,".dat",sep="")
                  datos <- read.table(filename, skip=21)
                  lambda <- datos[,1]
                  # AP_PS
                  dat <- datos[,3]
                  media1 <- sum(wd*dat)/(700-400)
                  QYi <- alpha / media1
                  points(x=QYi, y=media1, pch=21, cex=2.0,         bg=claros[10])          
                  # AP
                  dat <- datos[,2]
                  media2 <- sum(wd*dat)/(700-400)
                  points(x=QYi, y=media2, pch=22, cex=2.0,         bg=claros[10])          
         
                           
                  ### Small Phyto
                  alpha <- 0.14/86400
                  filename<-paste(global_path,"Phyto_optics/sensitivity/optics_phyto_recom_carbon_",i,".dat",sep="")
                  datos <- read.table(filename, skip=7, nrows=13)
                  lambda <- datos[,1]
                  # AP_PS
                  dat <- datos[,3]         
                  media1 <- sum(wd*dat)/(700-400)
                  QYi <- alpha / media1
                  points(x=QYi, y=media1, pch=21, cex=2.0,         bg=oscuros[10])            
                  # AP
                  dat <- datos[,2]
                  media2 <- sum(wd*dat)/(700-400)
                  points(x=QYi, y=media2, pch=22, cex=2.0,         bg=oscuros[10])          
                  
            # Spectra scaled to mean values in (Alvarez et al 2019 GBC): 3, c 
                  i=3 
                  ### Diatoms
                  alpha <- 0.19/86400
                  filename<-paste(global_path,"Phyto_optics/sensitivity/optics_phyto_recom_carbon_",i,".dat",sep="")
                  datos <- read.table(filename, skip=21)
                  lambda <- datos[,1]
                  # AP_PS
                  dat <- datos[,3]
                  media1 <- sum(wd*dat)/(700-400)
                  QYi <- alpha / media1
                  points(x=QYi, y=media1, pch=21, cex=2.0,         bg=claros[10])          
                  # AP
                  dat <- datos[,2]
                  media2 <- sum(wd*dat)/(700-400)
                  points(x=QYi, y=media2, pch=22, cex=2.0,         bg=claros[10])          
                  
                  ### Small Phyto
                  alpha <- 0.14/86400
                  filename<-paste(global_path,"Phyto_optics/sensitivity/optics_phyto_recom_carbon_",i,".dat",sep="")
                  datos <- read.table(filename, skip=7, nrows=13)
                  lambda <- datos[,1]
                  # AP_PS
                  dat <- datos[,3]         
                  media1 <- sum(wd*dat)/(700-400)
                  QYi <- alpha / media1
                  points(x=QYi, y=media1, pch=21, cex=2.0,         bg=oscuros[10])            
                  # AP
                  dat <- datos[,2]
                  media2 <- sum(wd*dat)/(700-400)
                  points(x=QYi, y=media2, pch=22, cex=2.0,         bg=oscuros[10])          
                  
        ## Segments 
        segments(x0=tabla_total[3,8],x1=tabla_total[3,8],y0=tabla_total[3,4],y1=tabla_total[3,6])   
        segments(x0=tabla_total[3,8],x1=tabla_total[3,9],y0=tabla_total[3,4],y1=tabla_total[3,4])
        segments(x0=tabla_total[12,8],x1=tabla_total[12,9],y0=tabla_total[12,4],y1=tabla_total[12,5])  
        segments(x0=tabla_total[12,8],x1=tabla_total[12,8],y0=tabla_total[12,4],y1=tabla_total[12,6])
        
        legend(x="topright",legend=c(expression(bar("a*")[PH]),
                                     expression(bar("a*")[PS])),pch=c(22,21), col="black", bty="n", cex=sizelab)
        
#############################
        
        dev.off()     
        
        
        
        
        
        
     
##########################################     
### FIGURE S2 Metrics Sensitivity Analysis (EXP-S)    
##########################################   
      png(file=paste(path_figures,"FigureS2_SA_metrics_vertical.png",sep=""),width=750, height=1000, units = "px", pointsize = 18, bg = "white") 

##########################################                 
      par(mfrow=c(4,2))
      par(mar=c(3,2,2,2))
      par(oma=c(3,3,3,0))
      layout.show(8)
      lettersize=2
      titlesize=1.6
      simbolos<-letters[seq( from = 1, to = 20 )] 

#### Chla
#### A
      plot(apph[1:20],R_chlal_sat[1:20], type="p", las=1, xlim=c(0.005,0.04), ylim=c(0.66,0.70),
                  xlab="m2 mg Chla-1", col="grey30",pch=simbolos, font=2, cex=1.5, yaxs="i", xaxs="i") 
                  regre<-lm(R_chlal_sat[1:20]~apph[1:20])
                  equis<-seq(0.007,0.039,length=10)
                  ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                  points(equis, ies, type="l", col="grey30")
      points(apps[1:20],R_chlal_sat[1:20], type="p", las=1,col="grey70",pch=simbolos, font=2, cex=1.5)
                  regre<-lm(R_chlal_sat[1:20]~apps[1:20])
                  ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                  points(equis, ies, type="l", col="grey70")
                  segments(x0=apps[1:20], x1=apph[1:20], y0=R_chlal_sat[1:20], lwd=c(rep(1,2),3,rep(1,8),3,rep(1,8)))  
       mtext(3,text="a)",line=0.5, at=0.007,cex=lettersize, font=1)
       legend(x="topright",legend=c(expression(bar("a*")[PH]),
                                    expression(bar("a*")[PS])),pch=19, col=c("grey30","grey70"), bty="n")   
       legend(x="topleft", legend=paste("n=",unique(n_chlal_sat)), bty="n")
       
#### B
       plot(apph[1:20],R_chlal_ist[1:20], type="p", las=1, xlim=c(0.005,0.04), ylim=c(0.51,0.54),
                       xlab="m2 mg Chla-1", col="grey30",pch=simbolos, font=2, cex=1.5, yaxs="i", xaxs="i") 
                  regre<-lm(R_chlal_ist[1:20]~apph[1:20])
                  ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                  points(equis, ies, type="l", col="grey30")
       points(apps[1:20],R_chlal_ist[1:20], type="p", las=1, col="grey70",pch=simbolos, font=2, cex=1.5)
                  regre<-lm(R_chlal_ist[1:20]~apps[1:20])
                  ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                  points(equis, ies, type="l", col="grey70")
                  segments(x0=apps[1:20], x1=apph[1:20], y0=R_chlal_ist[1:20], lwd=c(rep(1,2),3,rep(1,8),3,rep(1,8)))             
       mtext(3,text="b)",line=0.5,at=0.007,cex=lettersize, font=1)
       legend(x="topright",legend=c(expression(bar("a*")[PH]),
                                    expression(bar("a*")[PS])),pch=19, col=c("grey30","grey70"), bty="n")                  
       legend(x="topleft", legend=paste("n=",unique(n_chlal_ist)), bty="n")
                  
#### Rrs
#### C
             plot(apph[1:20],R_rrs443[1:20], type="p", las=1, xlim=c(0.005,0.04),  ylim=c(0.0,0.60), xlab="m2 mg Chla-1",pch=simbolos,
                      font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="blue2", yaxs="i", xaxs="i") 
                  regre<-lm(R_rrs443[1:20]~apph[1:20])
                  ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                  points(equis, ies, type="l", col="blue2")
             points(apph[1:20],R_rrs412[1:20], type="p", las=1,pch=simbolos,font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="purple")
                 regre<-lm(R_rrs412[1:20]~apph[1:20])
                 ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                 points(equis, ies, type="l", col="purple")
             points(apph[1:20],R_rrs490[1:20], type="p", las=1,pch=simbolos,font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="lightblue4")
                 regre<-lm(R_rrs490[1:20]~apph[1:20])
                 ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                 points(equis, ies, type="l", col="lightblue4")
             points(apph[1:20],R_rrs555[1:20], type="p", las=1,pch=simbolos,font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="grey70")
                 regre<-lm(R_rrs555[1:20]~apph[1:20])
                 ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                 points(equis, ies, type="l", col="grey70")
             points(apph[1:20],R_rrs670[1:20], type="p", las=1,pch=simbolos,font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="grey30")        
                 regre<-lm(R_rrs670[1:20]~apph[1:20])
                 ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                 points(equis, ies, type="l", col="grey30")
           mtext(3,text="c)",line=0.5,at=0.007,cex=lettersize, font=1)
           
           legend(x=0.005,y=0.42,
                  legend=c(expression(R[RS]("412nm")),expression(R[RS]("443nm")),expression(R[RS]("490nm"))),
                  pch=19,col=c("purple","blue2","lightblue4"), bty="n", cex=0.9)                  
           #legend(x=0.018,y=0.37,
           legend(x=0.005,y=0.15,
                  legend=c(expression(R[RS]("555nm")),expression(R[RS]("670nm"))),pch=19,
                  col=c("grey70","grey30"), bty="n", cex=0.9)                  
        
           legend(x=0.015,y=0.1, legend=paste("n=",unique(n_rrs670)), bty="n")
                                
#### D      
            plot(apph[1:20],R_rrs450[1:20], type="p", las=1, xlim=c(0.005,0.04),  ylim=c(0.0,0.60), xlab="m2 mg Chla-1",pch=simbolos,
                      font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="blue2", yaxs="i", xaxs="i") 
                 regre<-lm(R_rrs450[1:20]~apph[1:20])
                 ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                 points(equis, ies, type="l", col="blue2")
            points(apph[1:20],R_rrs400[1:20], type="p", las=1,pch=simbolos,font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="purple")
                 regre<-lm(R_rrs400[1:20]~apph[1:20])
                 ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                 points(equis, ies, type="l", col="purple")
            points(apph[1:20],R_rrs475[1:20], type="p", las=1,pch=simbolos,font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="lightblue4")
                 regre<-lm(R_rrs475[1:20]~apph[1:20])
                 ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                 points(equis, ies, type="l", col="lightblue4")
            points(apph[1:20],R_rrs500[1:20], type="p", las=1,pch=simbolos,font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="lightblue")
                 regre<-lm(R_rrs500[1:20]~apph[1:20])
                 ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                 points(equis, ies, type="l", col="lightblue")                 
            points(apph[1:20],R_rrs550[1:20], type="p", las=1,pch=simbolos, font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="grey70")
                 regre<-lm(R_rrs550[1:20]~apph[1:20])
                 ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                 points(equis, ies, type="l", col="grey70")
            points(apph[1:20],R_rrs675[1:20], type="p", las=1,pch=simbolos,font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="grey30")        
                 regre<-lm(R_rrs675[1:20]~apph[1:20])
                 ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                 points(equis, ies, type="l", col="grey30")               
           mtext(3,text="d)",line=0.5,at=0.007,cex=lettersize, font=1)
                 
           legend(x=0.005,y=0.42,legend=c(expression(R[RS]("400nm")),expression(R[RS]("450nm")),expression(R[RS]("475nm"))),pch=19,
                 col=c("purple","blue2","lightblue4"), bty="n", cex=0.9)                  

           legend(x=0.018,y=0.42,legend=c(expression(R[RS]("500nm")),expression(R[RS]("550nm")),expression(R[RS]("675nm"))),pch=19,
                  col=c("lightblue","grey70","grey30"), bty="n", cex=0.9)                  
           
           #text(x=0.015,y=0.3, pos=4, labels=paste("n=",unique(n_rrs675)))                  
           legend(x=0.030,y=0.41, legend=paste("n=",unique(n_rrs675)), bty="n")       
                   
#### IOPs sat
#### E
        plot(apph[1:20],R_tot[1:20], type="p", las=1, xlim=c(0.005,0.04),  ylim=c(0.2,0.5), xlab="m2 mg Chla-1",pch=simbolos,
                          font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="grey30", yaxs="i", xaxs="i") 
                          regre<-lm(R_tot[1:20]~apph[1:20])
                          ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                          points(equis, ies, type="l", col="grey30")    
        points(apph[1:20],R_aph[1:20], type="p", las=1,pch=simbolos, font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="palegreen4")
                          regre<-lm(R_aph[1:20]~apph[1:20])
                          ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                          points(equis, ies, type="l", col="palegreen4")
        points(apph[1:20],R_adg[1:20], type="p", las=1,pch=simbolos,font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="palegreen2")
                          regre<-lm(R_adg[1:20]~apph[1:20])
                          ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                          points(equis, ies, type="l", col="palegreen2")

         mtext(3,text="e)",line=0.5,at=0.007,cex=lettersize, font=1)
                          
         legend(x="topright",legend=c(expression(a[TOT]("443nm")),expression(a[PH]("443nm")),expression(a[DG]("443nm"))),pch=19,
                                 col=c("grey30","palegreen4","palegreen2"), bty="n")                  
                  
         legend(x="topleft", legend=paste("n=",unique(n_adg)), bty="n")                   
                                           
#### IOPs
#### F
               plot(apph[1:20],R_ap[1:20], type="p", las=1, xlim=c(0.005,0.04),  ylim=c(0.0,0.50), xlab="m2 mg Chla-1",pch=simbolos,
                               font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="grey70", yaxs="i", xaxs="i") 
                          regre<-lm(R_ap[1:20]~apph[1:20])
                          ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                          points(equis, ies, type="l", col="grey70")     
               points(apph[1:20],R_aphs[1:20], type="p", las=1,pch=simbolos,font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="palegreen4")
                          regre<-lm(R_aphs[1:20]~apph[1:20])
                          ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                          points(equis, ies, type="l", col="palegreen4")
               points(apph[1:20],R_acdom[1:20], type="p", las=1,pch=simbolos,font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="palegreen2")
                          regre<-lm(R_acdom[1:20]~apph[1:20])
                          ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                          points(equis, ies, type="l", col="palegreen2")
               points(apph[1:20],R_anap[1:20], type="p", las=1,pch=simbolos,font=2, cex=c(rep(1.5,2),2,rep(1.5,8),2,rep(1.5,8)), col="khaki4")        
                          regre<-lm(R_anap[1:20]~apph[1:20])
                          ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                          points(equis, ies, type="l", col="khaki4")                        
               mtext(3,text="f)",line=0.5,at=0.007,cex=lettersize, font=1)
                          
               legend(x="topright",legend=c(expression(a[P]("450nm")),expression(a[PH]("450nm")),
                                   expression(a[CDOM]("450nm")),expression(a[NAP]("450nm"))),pch=19,
                                   col=c("grey70","palegreen4","palegreen2","khaki4"), bty="n")  
                             
               legend(x="topleft", legend=paste("n=",unique(n_aphs)), bty="n")                                                       
                          
#### PPC
#### G
              plot(1,1,type="n",yaxt="n",xaxt="n",bty="n")
                          
              plot(apph[1:20],R_PPC[1:20], type="p", las=1, xlim=c(0.005,0.04),   ylim=c(0.46,0.54), xlab="m2 mg Chla-1",pch=simbolos,
                    font=2, cex=1.5, col="grey30", yaxs="i", xaxs="i")
                    regre<-lm(R_PPC[1:20]~apph[1:20])
                    ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                    points(equis, ies, type="l", col="grey30")     
              points(apps[1:20],R_PPC[1:20], type="p", las=1,pch=simbolos,font=2, cex=1.5, col="grey70")
                    regre<-lm(R_PPC[1:20]~apps[1:20])
                    ies<-regre$coefficients[2]*equis+regre$coefficients[1]
                    points(equis, ies, type="l", col="grey70")    
                    #plot(apps,AE_PPC, type="p", las=1, xlim=c(0.005,0.025),  ylim=c(-0.2,0.2), xlab="m2 mg Chla-1", col="blue") 
                    segments(x0=apps[1:20], x1=apph[1:20], y0=R_PPC[1:20], lwd=c(rep(1,2),3,rep(1,8),3,rep(1,8)))  

              mtext(3,text="g)",line=0.5,at=0.007,cex=lettersize, font=1)
              legend(x="topright",legend=c(expression(bar("a*")[PH]),
                                           expression(bar("a*")[PS])),pch=19, col=c("grey30","grey70"), bty="n")                  
                    
              legend(x="topleft", legend=paste("n=",unique(n_PPC)), bty="n")  
              
                    mtext(1,at=0.5,line=0.8,outer=T,text=expression(paste(bar("a*")," (",m^2, " mg ", Chl^-1,")")))
                    mtext(2,at=0.5,line=1.2,outer=T,text="Pearson's R")
                    mtext(4,at=0.5,line=1,outer=T,text="")                    
                          
##########################################  
         mtext(3,text=c("Satellite",expression(italic("In situ"))),line=-1.0,
                 at=c(0.25,0.75), cex=titlesize, font=1, outer=T) 
         dev.off()
         
         
                    
                    
                    
#################################         
### FIGURE S3 Contours of metrics (EXP-S)       
#################################         
  png(file=paste(path_figures,"FigureS3_SA_contours.png",sep=""),width = 800, height = 1000, units = "px", pointsize = 18, bg = "white") 
         layout(matrix(c(1:5,5), ncol=2, byrow=TRUE), heights = c(1,1,1.5))
         simbolos<-letters[seq( from = 1, to = 20 )]        
         par(mar=c(4,4.5,2,2))  
         par(oma=c(2,2,1,1))  
         layout.show(5)
         lettersize=2
         sizelab=1.4
         labinner=1.6
         legendsize=1.0
##########################################           
         # Clases
         clases_apph  <- seq(0.005,0.040,by=0.0001)
         clases_apph_mid <- rep(NA, length=length(clases_apph)-1)
         for (i in 1:length(clases_apph_mid)) { clases_apph_mid[i] <- ((clases_apph[i+1]-clases_apph[i])/2)+clases_apph[i] }
         indice1 <- apph[1:20]
         factor1 <- cut(indice1, breaks=clases_apph,  include.lowest = TRUE,right = TRUE)
         niveles1<- levels(factor1)        
         
         clases_apps  <-seq(0.005,0.040,by=0.0001)
         clases_apps_mid <- rep(NA, length=length(clases_apps)-1)
         for (i in 1:length(clases_apps_mid)) { clases_apps_mid[i] <- ((clases_apps[i+1]-clases_apps[i])/2)+clases_apps[i] }
         indice2 <- apps[1:20]
         factor2 <- cut(indice2, breaks=clases_apps,  include.lowest = TRUE,right = TRUE)
         niveles2<- levels(factor2)      
      
#### Chla
         value   <- R_chlal_sat[1:20]
         experiment <- data.frame(value, factor1, factor2)
         res <- aggregate(value~factor1+factor2,experiment, mean, na.rm=TRUE, drop=FALSE)         
         equis<-rep(clases_apph_mid,times=length(clases_apps_mid))
         ies  <-rep(clases_apps_mid,each=length(clases_apph_mid))
         zetas<-res$value
         tabla<-cbind(equis, ies, zetas)
         tabla<-tabla[complete.cases(tabla),]
         resul<-interp(x=tabla[,1],y=tabla[,2],z=tabla[,3],xo=clases_apph_mid, yo=clases_apps_mid)
         #RESUL3 <- matrix(res$value, nrow=length(niveles1), ncol=length(niveles2), byrow=FALSE)
         rownames(resul$z) <- clases_apph_mid   #niveles1    # range(RESUL, na.rm=TRUE)
         colnames(resul$z) <- clases_apps_mid   #niveles2 
         resul_chla<-resul$z
         #seecol("grad_all")
         mycol<-usecol(pal_grau,n=25)[c(1:20)]
         image2D(x=clases_apph_mid, y=clases_apps_mid, z=resul_chla, las=1, xlim=c(0.010,0.036), ylim=c(0.005,0.026),xlab="",ylab="",col=mycol)
         contour(x=clases_apph_mid, y=clases_apps_mid, z=resul_chla, levels=c(0.68), add=T,
                 lwd=2, col="black",labcex=sizelab)
         text(x=0.011, y=0.024, labels=expression(paste(log[10],"(TChla)",sep="")),pos=4, font=2, cex=labinner)
         mtext(3,text="a)",line=0.5, at=0.011,cex=lettersize, font=1)
         mtext(3,text="R",line=0.2, at=0.0373,cex=legendsize, font=1)
         
#### PPC
         value   <- R_PPC[1:20]
         experiment <- data.frame(value, factor1, factor2)
         res <- aggregate(value~factor1+factor2,experiment, mean, na.rm=TRUE, drop=FALSE)         
         equis<-rep(clases_apph_mid,times=length(clases_apps_mid))
         ies  <-rep(clases_apps_mid,each=length(clases_apph_mid))
         zetas<-res$value
         tabla<-cbind(equis, ies, zetas)
         tabla<-tabla[complete.cases(tabla),]
         resul<-interp(x=tabla[,1],y=tabla[,2],z=tabla[,3],xo=clases_apph_mid, yo=clases_apps_mid)
         #RESUL3 <- matrix(res$value, nrow=length(niveles1), ncol=length(niveles2), byrow=FALSE)
         rownames(resul$z) <- clases_apph_mid   #niveles1    # range(RESUL, na.rm=TRUE)
         colnames(resul$z) <- clases_apps_mid   #niveles2 
         resul_ppc<-resul$z
         mycol<-usecol(pal_peach,n=25)[c(1:20)]        
         image2D(x=clases_apph_mid, y=clases_apps_mid, z=resul_ppc, las=1, xlim=c(0.010,0.036), ylim=c(0.005,0.026),xlab="",ylab="",col=mycol)
         contour(x=clases_apph_mid, y=clases_apps_mid, z=resul_ppc, levels=c(0.50), add=T,
                 lwd=2, col="black",labcex=sizelab)
         text(x=0.011, y=0.024, labels="PPC:TChla",pos=4, font=1, cex=labinner)
         mtext(3,text="b)",line=0.5, at=0.011,cex=lettersize, font=1)
         mtext(3,text="R",line=0.2, at=0.0373,cex=legendsize, font=1)
         
#### Rrs 443
         value   <- R_rrs443[1:20]
         experiment <- data.frame(value, factor1, factor2)
         res <- aggregate(value~factor1+factor2,experiment, mean, na.rm=TRUE, drop=FALSE)         
         equis<-rep(clases_apph_mid,times=length(clases_apps_mid))
         ies  <-rep(clases_apps_mid,each=length(clases_apph_mid))
         zetas<-res$value
         tabla<-cbind(equis, ies, zetas)
         tabla<-tabla[complete.cases(tabla),]
         resul<-interp(x=tabla[,1],y=tabla[,2],z=tabla[,3],xo=clases_apph_mid, yo=clases_apps_mid)
         #RESUL3 <- matrix(res$value, nrow=length(niveles1), ncol=length(niveles2), byrow=FALSE)
         rownames(resul$z) <- clases_apph_mid   #niveles1    # range(RESUL, na.rm=TRUE)
         colnames(resul$z) <- clases_apps_mid   #niveles2 
         resul_rrs<-resul$z
         mycol<-usecol(pal_karpfenblau,n=25)[c(1:20)]
         image2D(x=clases_apph_mid, y=clases_apps_mid, z=resul_rrs, las=1, xlim=c(0.010,0.036), ylim=c(0.005,0.026),xlab="",ylab="", col=mycol)
         contour(x=clases_apph_mid, y=clases_apps_mid, z=resul_rrs, levels=c(0.55), add=T, lwd=2,
                 labcex=sizelab)
         text(x=0.011, y=0.024, labels=expression(R[RS]("443nm")),pos=4, font=2, cex=labinner)         
         mtext(3,text="c)",line=0.5, at=0.011,cex=lettersize, font=1)
         mtext(3,text="R",line=0.2, at=0.0373,cex=legendsize, font=1)
         
#### Aph 443
         value   <- R_aph[1:20]
         experiment <- data.frame(value, factor1, factor2)
         res <- aggregate(value~factor1+factor2,experiment, mean, na.rm=TRUE, drop=FALSE)         
         equis<-rep(clases_apph_mid,times=length(clases_apps_mid))
         ies  <-rep(clases_apps_mid,each=length(clases_apph_mid))
         zetas<-res$value
         tabla<-cbind(equis, ies, zetas)
         tabla<-tabla[complete.cases(tabla),]
         resul<-interp(x=tabla[,1],y=tabla[,2],z=tabla[,3],xo=clases_apph_mid, yo=clases_apps_mid)
         #RESUL3 <- matrix(res$value, nrow=length(niveles1), ncol=length(niveles2), byrow=FALSE)
         rownames(resul$z) <- clases_apph_mid   #niveles1    # range(RESUL, na.rm=TRUE)
         colnames(resul$z) <- clases_apps_mid   #niveles2 
         resul_aph<-resul$z
         mycol<-usecol(pal_seegruen,n=25)[c(1:20)]
         image2D(x=clases_apph_mid, y=clases_apps_mid, z=resul_aph, las=1, xlim=c(0.010,0.036), ylim=c(0.005,0.026),xlab="",ylab="", col=mycol)
         contour(x=clases_apph_mid, y=clases_apps_mid, z=resul_aph, levels=c(0.34), add=T, lwd=2,
                 labcex=sizelab)
         text(x=0.011, y=0.024, labels=expression(a[PH]("443nm")),pos=4, font=2, cex=labinner)           
         mtext(3,text="d)",line=0.5, at=0.011,cex=lettersize, font=1)
         mtext(3,text="R",line=0.2, at=0.0373,cex=legendsize, font=1)
         
#### Combined         
         plot(x=apph[1:20], y=apps[1:20], type="n",las=1, xlim=c(0.010,0.036), ylim=c(0.005,0.026),
              pch=simbolos, xlab="", ylab="", cex.axis=1.4)
         contour(x=clases_apph_mid, y=clases_apps_mid, z=resul_chla, levels=c(0.68),      add=T, lwd=2, col="grey30")
         contour(x=clases_apph_mid, y=clases_apps_mid, z=resul_ppc,  levels=c(0.50),      add=T, lwd=2, col="orange2")
         contour(x=clases_apph_mid, y=clases_apps_mid, z=resul_rrs,  levels=c(0.55),      add=T, lwd=2, col="blue")
         contour(x=clases_apph_mid, y=clases_apps_mid, z=resul_aph,  levels=c(0.34),      add=T, lwd=2, col="green2")
         points(x=apph[1:20], y=apps[1:20], las=1, xlim=c(0.007,0.036), ylim=c(0.005,0.026),pch=simbolos,font=2, cex=1.5)
         mtext(3,text="e)",at=0.0095, line=0.5,cex=lettersize, font=1)
         
         mtext(1,at=0.5,line=0.1,outer=T,text=expression(paste(bar("a*")[PH]," (", m^2, " mg ", Chl^-1,")")))
         mtext(2,at=0.5,line=0.1,outer=T,text=expression(paste(bar("a*")[PS]," (", m^2, " mg ", Chl^-1,")")))
         #mtext(2,at=0.5,line=1.2,outer=T,text="Pearson's R")
         #mtext(4,at=0.5,line=1,outer=T,text="")           
         
##########################################           
      dev.off()   
         