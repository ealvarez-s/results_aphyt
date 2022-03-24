source("00_Summary_EDITME.R")
source("functions_misc.R")

########################################
##  In situ bio-optical data  GLOBAL  ##
########################################
###          make 3D grids           ###
########################################
###  All datasets combined:  Valente + Casey + Chase + CSIRO + Astrid  ###
##########################################################################

p_dir<-paste(global_path,"Dat_observations/optics/total_tables/",sep="")
s_dir<-paste(global_path,"Dat_observations/optics/grids_optics/",sep="")
# Input: .Rdata matrices with rrs and iop data averaged in the wavebands (25nm) of the model and collocated (in /total_tables/)
# Output:  grids with 2x2 degrees and model depth layers (in /grids_optics/)


### Define grid
###############
### lon=2 x lat=2 x depth=model, time monthly

    o_dir <-paste(global_path,"Res_model/", sep="")
    archivo<-"grid.nc"
    filename <- paste(o_dir,archivo, sep="")
    filen <- open.nc(filename)
    filedepth <- read.nc(filen)       
    profun<-abs(filedepth$Z)
    prof_model <- profun 
    profunU<-abs(filedepth$Zu)        
    profunL<-abs(filedepth$Zl)    
    #cbind(profunL, prof_model, profunU)
    latitud  <- cbind(seq(-90,88,by=2),seq(-89,89,by=2),seq(-88,90,by=2))        #filedepth$Y
    longitud <- cbind(seq(-180,178,by=2),seq(-179,179,by=2),seq(-178,180,by=2))  #filedepth$X
    clases_lat  <-seq(-90,90,by=2)
    clases_log  <-seq(-180,180,by=2)
    clases_pro  <-c(0,profunU)
    clases_time <-c(0:12)    
    
    
    

########
### Chla: Valente, Chase, Astrid
########                  
  
            # Tablas 12nm
            load(paste(p_dir,"Valente2_ChlaRssIOP_modbands12.RData", sep=""))
            #colnames(tabla_total)
                  #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                  #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                  #dim(tabla_total)
            tabla1<-tabla_total[,c(1:5,8)] 
            colnames(tabla1)<-c("lista","fecha3","latitude","longitude","depth","MTChla")
            tabla1[,5]<-rep(5,length=nrow(tabla1))
            tabla1[,c(6)][tabla1[,c(6)]=="NaN"]<-NA
            flag1<-tabla_total[,which(colnames(tabla_total)=="QFChla")] 
            tabla1[,c(6)][flag1==1]<-NA
                  #range(tabla_total$fecha3[!is.na(tabla1[,c(6)])])
            load(paste(p_dir,"Chase_RssHplc_modbands12.RData", sep=""))
            #colnames(tabla_total)
                  #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                  #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                  #dim(tabla_total)
                  #tabla_total$depth
            tabla2<-tabla_total[,c(1:5,6)] 
            colnames(tabla2)<-c("lista","fecha3","latitude","longitude","depth","MTChla")
            tabla2[,c(6)][tabla2[,c(6)]=="NaN"]<-NA
            
            load(paste(p_dir,"Astrid5_IopHplc2nm_modbands12.RData", sep=""))
            #colnames(tabla_total)
                  #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                  #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                  #dim(tabla_total)
                  #tabla_total$depth
            tabla3<-tabla_total[,c(1:5,75)] 
            colnames(tabla3)<-c("lista","fecha3","latitude","longitude","depth","MTChla")
            tabla3[,c(6)][tabla3[,c(6)]=="NaN"]<-NA
             flag3<-tabla_total[,which(colnames(tabla_total)=="MYFLAG3")]
             flag5<-tabla_total[,which(colnames(tabla_total)=="MYFLAG5")]
             tabla1[,c(6)][flag3==3]<-NA            
 
            tabla_total<-rbind(tabla1,tabla2,tabla3)
            #colnames(tabla_total)
            #sum(!is.na(tabla_total[,6]))
            #range(tabla_total$fecha3[!is.na(tabla_total$MTChla)],na.rm=T)
            
       # Fill grid with ALL data
            lati<-tabla_total$latitude    # plot(long, lati)
            long<-tabla_total$longitude
            depth<-abs(tabla_total$depth)      # CUIDADO: lo lee como factor
            #range(depth, na.rm=TRUE)
            #sort(unique(depth, na.rm=TRUE))
            datetime <- as.Date(tabla_total$fecha3, format="%Y-%m-%d")  
            #range(datetime, na.rm=TRUE)
            YEAR  <- years(datetime)   #    range(YEAR, na.rm=TRUE)
            times <- as.numeric(format(datetime,"%m"))    
            month<-sort(unique(times))
            DAY <- days(datetime)
            JULIAN <- julian(datetime, origin = as.Date("1997-01-01"))
            #range(JULIAN, na.rm=TRUE)
        
            # Create grid levels
            indice1 <- c(long)   # range(indice1, na.rm=TRUE)
            indice2 <- c(lati)   # range(indice2, na.rm=TRUE)
            indice3 <- c(depth)  # range(indice3, na.rm=TRUE)
            indice4 <- c(times)  # range(indice4, na.rm=TRUE)
            factor1 <-cut(indice1, breaks=clases_log,   include.lowest = TRUE,right = TRUE)
            factor2 <-cut(indice2, breaks=clases_lat,   include.lowest = TRUE,right = TRUE)
            factor3 <-cut(indice3, breaks=clases_pro,   include.lowest = TRUE,right = TRUE)
            factor4 <-cut(indice4, breaks=clases_time,  include.lowest = FALSE,right = TRUE)
            niveles1<-levels(factor1)
            niveles2<-levels(factor2)
            niveles3<-levels(factor3)
            niveles4<-levels(factor4)           
        
            # Fill grid cells 
            valor1 <- tabla_total$MTChla
            experiment <- data.frame(valor1, factor1, factor2, factor3, factor4)
            #names(experiment)   
            res <- aggregate(valor1~factor1+factor2+factor3+factor4,experiment,
                             mean, na.rm=TRUE, drop=F, simplify=FALSE)
            #dim(res)
            #head(res)
            valor <- res$valor1
            #length(valor)
            valor[valor=="NULL"]<-NA
            valor2<-unlist(valor)
            #length(valor2)
            
            chla_diat <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(niveles4)),
                                   dimnames=list(x=niveles1,y=niveles2,z=niveles3,t=niveles4))
            for (i in 1:nrow(res)){
                chla_diat[which(niveles1==res[i,1]),which(niveles2==res[i,2]),which(niveles3==res[i,3]),which(niveles4==res[i,4])] <- valor2[i]} # end loop i
            
            chla_diat[chla_diat==0] <-NA
            res <- chla_diat        
            dimnames(res)<- list(x=longitud[,2],y=latitud[,2], z=profun,t=month)
            
#############################    
            
            # SEASONAL 
            save(longitud,latitud,profun,month,res,file=paste(s_dir,"media_seasonal_Chla.Rdata",sep=""))
            # ANNUAL
            res[res<=0] <- NA
            res <- apply(res, MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)
            save(longitud,latitud,profun,res,file=paste(s_dir,"media_annual_Chla.Rdata",sep=""))
            
                ## FIGURE Chla
                png(file=paste(s_dir,"MTChla.png", sep=""),width=600, height=350, pointsize=8) 
                ##############################                 
                par(mfrow=c(1,1))
                par(mar=c(3,3,1,2))
                image2D(x=longitud[,2], y=latitud[,2], z=log10(res[,,1]), xlim=c(-180,180), ylim=c(-90,90),ylab="", xlab="",las=1, main="",xaxs="i", yaxs="i")
                load(paste(global_path,"Misc/batimetria_world_degree6.Rdata",sep=""))
                par(mar=c(3,3,1,4))
                contour(x=lon, y=lat, z=prof, levels=c(0), lwd=0.4, ylab="", xlab="",drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",bty="n", add=TRUE)
                dev.off()         
                ##############################



                
                
########
#### Rrs: Valente, Chase, Casey, Astrid
########              
  
                # Tablas 12nm
                load(paste(p_dir,"Valente_ChlaRssIOP_modbands12.RData", sep="")) 
                #colnames(tabla_total)
                      #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                      #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                      #dim(tabla_total)                
                lambda <- as.numeric(substring(colnames(tabla_total),5,8))[which(substring(colnames(tabla_total),1,3)=="rss")]                   
                tabla1<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="rss"))] 
                colnames(tabla1)<-c("lista","fecha3","latitude","longitude","depth",lambda)
                tabla1[,5]<-rep(5,length=nrow(tabla1))
                tabla1[,c(6:18)][tabla1[,c(6:18)]=="NaN"]<-NA
                
                load(paste(p_dir,"Chase_RssHplc_modbands12.RData", sep=""))   
                      #dim(tabla_total)
                      #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                      #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                      #dim(tabla_total)                
                tabla2<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="rss"))] 
                colnames(tabla2)<-c("lista","fecha3","latitude","longitude","depth",lambda)
                tabla2[,c(6:18)][tabla2[,c(6:18)]=="NaN"]<-NA
                
                load(paste(p_dir,"Casey_RssIop_modbands12.RData", sep=""))  
                      #dim(tabla_total)
                      #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                      #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                      #dim(tabla_total)                
                tabla3<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="rss"))] 
                colnames(tabla3)<-c("lista","fecha3","latitude","longitude","depth",lambda)
                tabla3[,c(6:18)][tabla3[,c(6:18)]=="NaN"]<-NA
                
                load(paste(p_dir,"Astrid5_IopHplc2nm_modbands12.RData", sep=""))  
                #colnames(tabla_total)
                    #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                    #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                    #dim(tabla_total)                
                tabla4<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="rrs"))] 
                colnames(tabla4)<-c("lista","fecha3","latitude","longitude","depth",lambda)
                tabla4[,c(6:18)][tabla4[,c(6:18)]=="NaN"]<-NA                
                
                tabla_total<-rbind(tabla1,tabla2,tabla3,tabla4)
                #colnames(tabla_total)                
                sum(!is.na(tabla_total[,10]))
                range(tabla_total$fecha3[!is.na(tabla_total$"450")],na.rm=T)
                
                
          # Fill grid with ALL data
                lati<-tabla_total$latitude    # plot(long, lati)
                long<-tabla_total$longitude
                depth<-abs(tabla_total$depth)      # CUIDADO: lo lee como factor
                depth[depth==999]<-NA
                #range(depth, na.rm=TRUE)
                #sort(unique(depth, na.rm=TRUE))
                datetime <- as.Date(tabla_total$fecha3, format="%Y-%m-%d")  
                #range(datetime, na.rm=TRUE)
                YEAR  <- years(datetime)   #    range(YEAR, na.rm=TRUE)
                times <- as.numeric(format(datetime,"%m"))    
                month<-sort(unique(times))
                DAY <- days(datetime)
                JULIAN <- julian(datetime, origin = as.Date("1997-01-01"))
                #range(JULIAN, na.rm=TRUE)
  
                # Create grid levels
                indice1 <- c(long)   # range(indice1, na.rm=TRUE)
                indice2 <- c(lati)   # range(indice2, na.rm=TRUE)
                indice3 <- c(depth)  # range(indice3, na.rm=TRUE)
                indice4 <- c(times)  # range(indice4, na.rm=TRUE)
                factor1 <-cut(indice1, breaks=clases_log,   include.lowest = TRUE,right = TRUE)
                factor2 <-cut(indice2, breaks=clases_lat,   include.lowest = TRUE,right = TRUE)
                factor3 <-cut(indice3, breaks=clases_pro,   include.lowest = TRUE,right = TRUE)
                factor4 <-cut(indice4, breaks=clases_time,  include.lowest = FALSE,right = TRUE)
                niveles1<-levels(factor1)
                niveles2<-levels(factor2)
                niveles3<-levels(factor3)
                niveles4<-levels(factor4) 
        
                ## Fill the grid
                matriz_seasonal <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles4),length(lambda)),
                                                   dimnames=list(x=niveles1,y=niveles2,t=niveles4,l=lambda))
                matriz_annual <- array(NA, dim=c(length(niveles1),length(niveles2),length(lambda)),
                                                 dimnames=list(x=niveles1, y=niveles2, l=lambda))
        
                for (w in c(1:length(lambda))){        
                  valor1 <- tabla_total[,5+w]
                  chla_diat <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles4)),
                                         dimnames=list(x=niveles1,y=niveles2,t=niveles4))
                  experiment <- data.frame(valor1, factor1, factor2, factor4)
                  if(sum(complete.cases(experiment))>0){
                  res <- aggregate(valor1~factor1+factor2+factor4,experiment,
                                   mean, na.rm=TRUE, drop=TRUE, simplify=FALSE, na.action = na.omit)
                  valor <- unlist(res$valor1)
                  for (i in 1:nrow(res)){
                    chla_diat[which(niveles1==res[i,1]),which(niveles2==res[i,2]),which(niveles4==res[i,3])] <- valor[i]} # end loop i
                  chla_diat[chla_diat==0] <-NA
                  }
                  dimnames(chla_diat) <- list(x=longitud[,2], y=latitud[,2], t=month)
                  res <- chla_diat
                  # SEASONAL
                  if (sum(!is.na(res))>=10){ matriz_seasonal[,,,w] <- res }
                  # ANNUAL
                  res[res<=0] <- NA
                  res <- apply(res, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)
                  if (sum(!is.na(res))>=10){ matriz_annual[,,w] <- res}
                }  # end loop w lambda
#########################
        
              save(longitud,latitud,month,lambda,  matriz_seasonal,  file=paste(s_dir,"media_seasonal_Rrs.Rdata",sep=""))
              save(longitud,latitud,lambda,        matriz_annual,    file=paste(s_dir,"media_annual_Rrs.Rdata",sep=""))


                          
########
##### ACDOM: Valente, CSIRO, Casey, Astrid
########                         
               # Tablas 12nm
               load(paste(p_dir,"Valente_ChlaRssIOP_modbands12.RData", sep="")) 
               #colnames(tabla_total)
                          #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                          #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                          #dim(tabla_total)       
                          lambda <- as.numeric(substring(colnames(tabla_total),5,8))[which(substring(colnames(tabla_total),1,3)=="adg")]                   
                          tabla1<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="adg"))] 
                          colnames(tabla1)<-c("lista","fecha3","latitude","longitude","depth",lambda)
                          tabla1[,5]<-rep(5,length=nrow(tabla1))
                          tabla1[,c(6:18)][tabla1[,c(6:18)]=="NaN"]<-NA
                          
                          load(paste(p_dir,"Casey_RssIop_modbands12.RData", sep=""))  
                          #dim(tabla_total)
                          #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                          #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                          #dim(tabla_total)        
                          tabla2<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="acd"))] 
                          colnames(tabla2)<-c("lista","fecha3","latitude","longitude","depth",lambda)
                          tabla2[,c(6:18)][tabla2[,c(6:18)]=="NaN"]<-NA
                          
                          load(paste(p_dir,"CSIRO_Iop_modbands12.RData", sep=""))   
                          #dim(tabla_total)
                          #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                          #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                          #dim(tabla_total)        
                          tabla3<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="acd"))] 
                          colnames(tabla3)<-c("lista","fecha3","latitude","longitude","depth",lambda)
                          tabla3[,c(6:18)][tabla3[,c(6:18)]=="NaN"]<-NA
                          
                 load(paste(p_dir,"Astrid5_IopHplc2nm_modbands12.RData", sep=""))
                 #colnames(tabla_total)
                          #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                          #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                          #dim(tabla_total)        
                          tabla4<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,5)=="acdom"))] 
                          colnames(tabla4)<-c("lista","fecha3","latitude","longitude","depth",lambda)
                          tabla4[,c(6:18)][tabla4[,c(6:18)]=="NaN"]<-NA                          
                          
                          
                          tabla_total<-rbind(tabla1,tabla2,tabla3,tabla4)
                          #dim(tabla_total)
                          #colnames(tabla_total)                
                          #sum(!is.na(tabla_total[,9]))
                          #dim(tabla_total[!is.na(tabla_total[,9]),])
                            #plot(lambda,tabla_total[!is.na(tabla_total[,9]),][1,6:18], type="b")
                          #points(lambda,tabla_total[!is.na(tabla_total[,9]),][2,6:18], type="b")
                          #range(tabla_total$fecha3[!is.na(tabla_total$"450")],na.rm=T)
                          
                          # Fill grid with ALL data
                          lati<-tabla_total$latitude    # plot(long, lati)
                          long<-tabla_total$longitude
                          depth<-abs(tabla_total$depth)      # CUIDADO: lo lee como factor
                          depth[depth==999]<-NA
                          #range(depth, na.rm=TRUE)
                          #sort(unique(depth, na.rm=TRUE))
                          datetime <- as.Date(tabla_total$fecha3, format="%Y-%m-%d")  
                          #range(datetime, na.rm=TRUE)
                          YEAR  <- years(datetime)   #    range(YEAR, na.rm=TRUE)
                          times <- as.numeric(format(datetime,"%m"))    
                          month<-sort(unique(times))
                          DAY <- days(datetime)
                          JULIAN <- julian(datetime, origin = as.Date("1997-01-01"))
                          #range(JULIAN, na.rm=TRUE)
                          
                          # Create grid levels
                          indice1 <- c(long)   # range(indice1, na.rm=TRUE)
                          indice2 <- c(lati)   # range(indice2, na.rm=TRUE)
                          indice3 <- c(depth)  # range(indice3, na.rm=TRUE)
                          indice4 <- c(times)  # range(indice4, na.rm=TRUE)
                          factor1 <-cut(indice1, breaks=clases_log,   include.lowest = TRUE,right = TRUE)
                          factor2 <-cut(indice2, breaks=clases_lat,   include.lowest = TRUE,right = TRUE)
                          factor3 <-cut(indice3, breaks=clases_pro,   include.lowest = TRUE,right = TRUE)
                          factor4 <-cut(indice4, breaks=clases_time,  include.lowest = FALSE,right = TRUE)
                          niveles1<-levels(factor1)
                          niveles2<-levels(factor2)
                          niveles3<-levels(factor3)
                          niveles4<-levels(factor4) 
                          
                          ## Fill the grid
                          matriz_seasonal <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(niveles4),length(lambda)),
                                                   dimnames=list(x=niveles1,y=niveles2,z=niveles3,t=niveles4,l=lambda))
                          matriz_annual <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(lambda)),
                                                 dimnames=list(x=niveles1, y=niveles2, z=niveles3, l=lambda))
                          #dim(matriz_seasonal)
                          
                          for (w in c(1:length(lambda))){
                            valor1 <- tabla_total[,5+w]
                            experiment <- data.frame(valor1, factor1, factor2, factor3, factor4)
                            #names(experiment)   
                            chla_diat <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(niveles4)),
                                               dimnames=list(x=niveles1,y=niveles2,z=niveles3,t=niveles4))                            
                            if(sum(complete.cases(experiment))>0){
                              res <- aggregate(valor1~factor1+factor2+factor3+factor4,experiment,mean, na.rm=TRUE, drop=F, simplify=FALSE)
                              #dim(res)
                              #head(res)
                              valor <- res$valor1
                              #length(valor)
                              valor[valor=="NULL"]<-NA
                              valor2<-unlist(valor)
                              #length(valor2)
                              
                              for (i in 1:nrow(res)){
                                chla_diat[which(niveles1==res[i,1]),which(niveles2==res[i,2]),which(niveles3==res[i,3]),which(niveles4==res[i,4])] <- valor2[i]} # end loop i
                              }
                            chla_diat[chla_diat==0] <-NA
                            res <- chla_diat        
                            dimnames(res)<- list(x=longitud[,2], y=latitud[,2], z=profun, t=month) 
                            
                            # SEASONAL
                            if (sum(!is.na(res))>=10){ matriz_seasonal[,,,,w] <- res }
                            # ANNUAL
                            res[res<=0] <- NA
                            res <- apply(res, MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)
                            #image2D(log10(res[,,1]), zlim=c(-3,0))
                            if (sum(!is.na(res))>=10){ matriz_annual[,,,w] <- res }
                            
                            print(paste(lambda[w]," ready"))
                          }  # end loop w lambda
######################### 
                          
             save(longitud,latitud,profun,month,lambda,  matriz_seasonal,  file=paste(s_dir,"media_seasonal_Acdom.Rdata",sep=""))
             save(longitud,latitud,profun,lambda,        matriz_annual,    file=paste(s_dir,"media_annual_Acdom.Rdata",sep=""))
                          


########
##### APH: Valente, Astrid, CSIRO, Casey
########              
                        
            # Tablas 12nm
            load(paste(p_dir,"Valente2_ChlaRssIOP_modbands12.RData", sep="")) 
            #colnames(tabla_total)
                              #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                              #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                              #dim(tabla_total)       
            lambda <- as.numeric(substring(colnames(tabla_total),5,8))[which(substring(colnames(tabla_total),1,3)=="aph")]                   
            tabla1<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="aph"))] 
            colnames(tabla1)<-c("lista","fecha3","latitude","longitude","depth",lambda)
            tabla1[,5]<-rep(5,length=nrow(tabla1))
            tabla1[,c(6:18)][tabla1[,c(6:18)]=="NaN"]<-NA
            
            load(paste(p_dir,"Casey_RssIop_modbands12.RData", sep=""))  
                              #dim(tabla_total)
                              #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                              #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                              #dim(tabla_total)        
            tabla2<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="aph"))] 
            colnames(tabla2)<-c("lista","fecha3","latitude","longitude","depth",lambda)
            tabla2[,c(6:18)][tabla2[,c(6:18)]=="NaN"]<-NA
                # Completar   
                tabla_ap<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ap_"))] 
                tabla_anap<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ana"))] 
                tablita<-(tabla_ap[,6:18]-tabla_anap[,6:18])
                tabla2[,6:18][is.na(tabla2[,6:18])]<-tablita[is.na(tabla2[,6:18])]
                
            load(paste(p_dir,"CSIRO_Iop_modbands12.RData", sep=""))   
                              #dim(tabla_total)
                              #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                              #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                              #dim(tabla_total)        
            tabla3<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="aph"))] 
            colnames(tabla3)<-c("lista","fecha3","latitude","longitude","depth",lambda)
            tabla3[,c(6:18)][tabla3[,c(6:18)]=="NaN"]<-NA
                # Completar   
                tabla_ap<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ap_"))] 
                tabla_anap<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ad_"))] 
                tablita<-(tabla_ap[,6:18]-tabla_anap[,6:18])
                tabla3[,6:18][is.na(tabla3[,6:18])]<-tablita[is.na(tabla3[,6:18])]
                
                
            load(paste(p_dir,"Astrid5_IopHplc2nm_modbands12.RData", sep=""))   
            #colnames(tabla_total)
                              #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                              #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                              #dim(tabla_total)        
            tabla4<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="aph"))] 
            colnames(tabla4)<-c("lista","fecha3","latitude","longitude","depth",lambda)
            tabla4[,c(6:18)][tabla4[,c(6:18)]=="NaN"]<-NA
                # Completar   
                tabla_ap<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ap_"))] 
                tabla_anap<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ana"))] 
                tablita<-(tabla_ap[,6:18]-tabla_anap[,6:18])
                tabla4[,6:18][is.na(tabla4[,6:18])]<-tablita[is.na(tabla4[,6:18])]

                
            tabla_total<-rbind(tabla1,tabla2,tabla3,tabla4)
            #names(tabla_total)
            #range(tabla_total$fecha3[!is.na(tabla_total$"450")],na.rm=T)
            #range(tabla_total$fecha3)
            #plot(tabla_total$latitude,tabla_total[,8],cex=0.4)
            #colnames(tabla_total)                
            #sum(!is.na(tabla_total[,9]))
            #dim(tabla_total[!is.na(tabla_total[,9]),])
              #plot(lambda,tabla_total[!is.na(tabla_total[,9]),][1,6:18], type="b", las=1, ylim=c(0,0.05))
            #points(lambda,tabla_total[!is.na(tabla_total[,9]),][10,6:18], type="b")
            #save(tabla_total,file=paste(s_dir,"tabla_total_Aph_test.Rdata",sep=""))
            
            # Fill grid with ALL data
                          lati<-tabla_total$latitude    # plot(long, lati)
                          long<-tabla_total$longitude
                          depth<-abs(tabla_total$depth)      # CUIDADO: lo lee como factor
                          depth[depth==999]<-NA
                          #range(depth, na.rm=TRUE)
                          #sort(unique(depth, na.rm=TRUE))
                          datetime <- as.Date(tabla_total$fecha3, format="%Y-%m-%d")  
                          #range(datetime, na.rm=TRUE)
                          YEAR  <- years(datetime)   #    range(YEAR, na.rm=TRUE)
                          times <- as.numeric(format(datetime,"%m"))    
                          month<-sort(unique(times))
                          DAY <- days(datetime)
                          JULIAN <- julian(datetime, origin = as.Date("1997-01-01"))
                          #range(JULIAN, na.rm=TRUE)
                          
                          #colnames(tabla_total)
                          
                          # Create grid levels
                          indice1 <- c(long)   # range(indice1, na.rm=TRUE)
                          indice2 <- c(lati)   # range(indice2, na.rm=TRUE)
                          indice3 <- c(depth)  # range(indice3, na.rm=TRUE)
                          indice4 <- c(times)  # range(indice4, na.rm=TRUE)
                          factor1 <-cut(indice1, breaks=clases_log,   include.lowest = TRUE,right = TRUE)
                          factor2 <-cut(indice2, breaks=clases_lat,   include.lowest = TRUE,right = TRUE)
                          factor3 <-cut(indice3, breaks=clases_pro,   include.lowest = TRUE,right = TRUE)
                          factor4 <-cut(indice4, breaks=clases_time,  include.lowest = FALSE,right = TRUE)
                          niveles1<-levels(factor1)
                          niveles2<-levels(factor2)
                          niveles3<-levels(factor3)
                          niveles4<-levels(factor4) 
                          
                          ## Fill the grid
                          matriz_seasonal <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(niveles4),length(lambda)),
                                                   dimnames=list(x=niveles1,y=niveles2,z=niveles3,t=niveles4,l=lambda))
                          matriz_annual <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(lambda)),
                                                 dimnames=list(x=niveles1, y=niveles2, z=niveles3, l=lambda))
                          #dim(matriz_seasonal)
                          
                          for (w in c(1:length(lambda))){
                            valor1 <- tabla_total[,5+w]
                            experiment <- data.frame(valor1, factor1, factor2, factor3, factor4)
                            #names(experiment)   
                            chla_diat <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(niveles4)),
                                               dimnames=list(x=niveles1,y=niveles2,z=niveles3,t=niveles4))                            
                            if(sum(complete.cases(experiment))>0){
                            res <- aggregate(valor1~factor1+factor2+factor3+factor4,experiment,mean, na.rm=TRUE, drop=F, simplify=FALSE)
                            #dim(res)
                            #head(res)
                            valor <- res$valor1
                            #length(valor)
                            valor[valor=="NULL"]<-NA
                            valor2<-unlist(valor)
                            #length(valor2)

                            for (i in 1:nrow(res)){
                              chla_diat[which(niveles1==res[i,1]),which(niveles2==res[i,2]),which(niveles3==res[i,3]),which(niveles4==res[i,4])] <- valor2[i]} # end loop i
                            }
                            chla_diat[chla_diat==0] <-NA
                            res <- chla_diat        
                            dimnames(res)<- list(x=longitud[,2], y=latitud[,2], z=profun, t=month) 
                            
                            # SEASONAL
                            if (sum(!is.na(res))>=10){ matriz_seasonal[,,,,w] <- res }
                            # ANNUAL
                            res[res<=0] <- NA
                            res <- apply(res, MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)
                            #image2D(log10(res[,,1]), zlim=c(-3,0))
                            if (sum(!is.na(res))>=10){ matriz_annual[,,,w] <- res }
                            
                            print(paste(lambda[w]," ready"))
                          }  # end loop w lambda
######################### 
                          
            save(longitud,latitud,profun,month,lambda,  matriz_seasonal,  file=paste(s_dir,"media_seasonal_Aph.Rdata",sep=""))
            save(longitud,latitud,profun,lambda,        matriz_annual,    file=paste(s_dir,"media_annual_Aph.Rdata",sep=""))
                       
                        
                        
                                          
########
##### ANAP: Astrid, CSIRO, Casey
########                          
                        # Tablas 12nm
                        load(paste(p_dir,"Astrid5_IopHplc2nm_modbands12.RData", sep=""))   
                        #dim(tabla_total)
                        #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                        #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                        #dim(tabla_total)        
                        lambda <- as.numeric(substring(colnames(tabla_total),6,9))[which(substring(colnames(tabla_total),1,3)=="ana")]                   
                        tabla1<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ana"))] 
                        colnames(tabla1)<-c("lista","fecha3","latitude","longitude","depth",lambda)
                        tabla1[,5]<-rep(5,length=nrow(tabla1))
                        tabla1[,c(6:18)][tabla1[,c(6:18)]=="NaN"]<-NA
                            # Completar   
                            tabla_ap<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ap_"))] 
                            tabla_aph<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="aph"))] 
                            tablita<-(tabla_ap[,6:18]-tabla_aph[,6:18])
                            tabla1[,6:18][is.na(tabla1[,6:18])]<-tablita[is.na(tabla1[,6:18])]

                            
                        load(paste(p_dir,"Casey_RssIop_modbands12.RData", sep=""))  
                        #dim(tabla_total)
                        #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                        #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                        #dim(tabla_total)        
                        tabla2<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ana"))] 
                        colnames(tabla2)<-c("lista","fecha3","latitude","longitude","depth",lambda)
                        tabla2[,c(6:18)][tabla2[,c(6:18)]=="NaN"]<-NA
                            # Completar   
                            tabla_ap<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ap_"))] 
                            tabla_aph<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="aph"))] 
                            tablita<-(tabla_ap[,6:18]-tabla_aph[,6:18])
                            tabla2[,6:18][is.na(tabla2[,6:18])]<-tablita[is.na(tabla2[,6:18])]
 
                            
                        load(paste(p_dir,"CSIRO_Iop_modbands12.RData", sep=""))   
                        #dim(tabla_total)
                        #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                        #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                        #dim(tabla_total)        
                        tabla3<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ad_"))] 
                        colnames(tabla3)<-c("lista","fecha3","latitude","longitude","depth",lambda)
                        tabla3[,c(6:18)][tabla3[,c(6:18)]=="NaN"]<-NA
                            # Completar   
                            tabla_ap<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ap_"))] 
                            tabla_aph<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="aph"))] 
                            tablita<-(tabla_ap[,6:18]-tabla_aph[,6:18])
                            tabla3[,6:18][is.na(tabla3[,6:18])]<-tablita[is.na(tabla3[,6:18])]
  
                            
                        tabla_total<-rbind(tabla1,tabla2,tabla3)
                        #dim(tabla_total)
                        #colnames(tabla_total)                
                        #sum(!is.na(tabla_total[,9]))
                        #dim(tabla_total[!is.na(tabla_total[,9]),])
                        #  plot(lambda,tabla_total[!is.na(tabla_total[,9]),][1,6:18], type="b", las=1, ylim=c(0,0.1))
                        #points(lambda,tabla_total[!is.na(tabla_total[,9]),][10,6:18], type="b")
                        
                        
                        # Fill grid with ALL data
                        lati<-tabla_total$latitude    # plot(long, lati)
                        long<-tabla_total$longitude
                        depth<-abs(tabla_total$depth)      # CUIDADO: lo lee como factor
                        depth[depth==999]<-NA
                        #range(depth, na.rm=TRUE)
                        #sort(unique(depth, na.rm=TRUE))
                        datetime <- as.Date(tabla_total$fecha3, format="%Y-%m-%d")  
                        #range(datetime, na.rm=TRUE)
                        YEAR  <- years(datetime)   #    range(YEAR, na.rm=TRUE)
                        times <- as.numeric(format(datetime,"%m"))    
                        month<-sort(unique(times))
                        DAY <- days(datetime)
                        JULIAN <- julian(datetime, origin = as.Date("1997-01-01"))
                        #range(JULIAN, na.rm=TRUE)
                        
                        #colnames(tabla_total)
                        
                        # Create grid levels
                        indice1 <- c(long)   # range(indice1, na.rm=TRUE)
                        indice2 <- c(lati)   # range(indice2, na.rm=TRUE)
                        indice3 <- c(depth)  # range(indice3, na.rm=TRUE)
                        indice4 <- c(times)  # range(indice4, na.rm=TRUE)
                        factor1 <-cut(indice1, breaks=clases_log,   include.lowest = TRUE,right = TRUE)
                        factor2 <-cut(indice2, breaks=clases_lat,   include.lowest = TRUE,right = TRUE)
                        factor3 <-cut(indice3, breaks=clases_pro,   include.lowest = TRUE,right = TRUE)
                        factor4 <-cut(indice4, breaks=clases_time,  include.lowest = FALSE,right = TRUE)
                        niveles1<-levels(factor1)
                        niveles2<-levels(factor2)
                        niveles3<-levels(factor3)
                        niveles4<-levels(factor4) 
                        
                        ## Fill the grid
                        matriz_seasonal <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(niveles4),length(lambda)),
                                                 dimnames=list(x=niveles1,y=niveles2,z=niveles3,t=niveles4,l=lambda))
                        matriz_annual <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(lambda)),
                                               dimnames=list(x=niveles1, y=niveles2, z=niveles3, l=lambda))
                        #dim(matriz_seasonal)
                        
                        for (w in c(1:length(lambda))){
                          valor1 <- tabla_total[,5+w]
                          experiment <- data.frame(valor1, factor1, factor2, factor3, factor4)
                          #names(experiment)   
                          chla_diat <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(niveles4)),
                                             dimnames=list(x=niveles1,y=niveles2,z=niveles3,t=niveles4))                            
                          if(sum(complete.cases(experiment))>0){
                            res <- aggregate(valor1~factor1+factor2+factor3+factor4,experiment,mean, na.rm=TRUE, drop=F, simplify=FALSE)
                            #dim(res)
                            #head(res)
                            valor <- res$valor1
                            #length(valor)
                            valor[valor=="NULL"]<-NA
                            valor2<-unlist(valor)
                            #length(valor2)
                            
                            for (i in 1:nrow(res)){
                              chla_diat[which(niveles1==res[i,1]),which(niveles2==res[i,2]),which(niveles3==res[i,3]),which(niveles4==res[i,4])] <- valor2[i]} # end loop i
                          }
                          chla_diat[chla_diat==0] <-NA
                          res <- chla_diat        
                          dimnames(res)<- list(x=longitud[,2], y=latitud[,2], z=profun, t=month) 
                          
                          # SEASONAL
                          if (sum(!is.na(res))>=10){ matriz_seasonal[,,,,w] <- res }
                          # ANNUAL
                          res[res<=0] <- NA
                          res <- apply(res, MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)
                          #image2D(log10(res[,,1]), zlim=c(-3,0))
                          if (sum(!is.na(res))>=10){ matriz_annual[,,,w] <- res }
                          
                          print(paste(lambda[w]," ready"))
                        }  # end loop w lambda
######################### 
                        
            save(longitud,latitud,profun,month,lambda,  matriz_seasonal,  file=paste(s_dir,"media_seasonal_Anap.Rdata",sep=""))
            save(longitud,latitud,profun,lambda,        matriz_annual,    file=paste(s_dir,"media_annual_Anap.Rdata",sep=""))
                        

                        
                      
                                  
########
##### AP: Astrid, CSIRO, Casey
########                         
                        # Tablas 12nm
                        load(paste(p_dir,"Astrid5_IopHplc2nm_modbands12.RData", sep=""))   
                        #dim(tabla_total)
                        #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                        #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                        #dim(tabla_total)        
                        lambda <- as.numeric(substring(colnames(tabla_total),4,6))[which(substring(colnames(tabla_total),1,3)=="ap_")]                   
                        tabla1<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ap_"))] 
                        colnames(tabla1)<-c("lista","fecha3","latitude","longitude","depth",lambda)
                        tabla1[,5]<-rep(5,length=nrow(tabla1))
                        tabla1[,c(6:18)][tabla1[,c(6:18)]=="NaN"]<-NA
                            # Completar   
                            tabla_anap<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ana"))] 
                            tabla_aph<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="aph"))] 
                            tablita<-(tabla_anap[,6:18]+tabla_aph[,6:18])
                            tabla1[,6:18][is.na(tabla1[,6:18])]<-tablita[is.na(tabla1[,6:18])]
 
                            
                        load(paste(p_dir,"Casey_RssIop_modbands12.RData", sep=""))  
                        #dim(tabla_total)
                        #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                        #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                        #dim(tabla_total)        
                        tabla2<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ap_"))] 
                        colnames(tabla2)<-c("lista","fecha3","latitude","longitude","depth",lambda)
                        tabla2[,c(6:18)][tabla2[,c(6:18)]=="NaN"]<-NA
                            # Completar   
                            tabla_anap<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ana"))] 
                            tabla_aph<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="aph"))] 
                            tablita<-(tabla_anap[,6:18]+tabla_aph[,6:18])
                            tabla2[,6:18][is.na(tabla2[,6:18])]<-tablita[is.na(tabla2[,6:18])]
 
                            
                        load(paste(p_dir,"CSIRO_Iop_modbands12.RData", sep=""))   
                        #dim(tabla_total)
                        #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                        #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                        #dim(tabla_total)        
                        tabla3<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ap_"))] 
                        colnames(tabla3)<-c("lista","fecha3","latitude","longitude","depth",lambda)
                        tabla3[,c(6:18)][tabla3[,c(6:18)]=="NaN"]<-NA
                            # Completar   
                            tabla_anap<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ad_"))] 
                            tabla_aph<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="aph"))] 
                            tablita<-(tabla_anap[,6:18]+tabla_aph[,6:18])
                            tabla3[,6:18][is.na(tabla3[,6:18])]<-tablita[is.na(tabla3[,6:18])]

                            
                        tabla_total<-rbind(tabla1,tabla2,tabla3)
                        #dim(tabla_total)
                        #colnames(tabla_total)                
                        #sum(!is.na(tabla_total[,9]))
                        #dim(tabla_total[!is.na(tabla_total[,9]),])
                        #  plot(lambda,tabla_total[!is.na(tabla_total[,9]),][1,6:18], type="b", las=1, ylim=c(0,0.1))
                        #points(lambda,tabla_total[!is.na(tabla_total[,9]),][10,6:18], type="b")              
                        
                        # Fill grid with ALL data
                        lati<-tabla_total$latitude    # plot(long, lati)
                        long<-tabla_total$longitude
                        depth<-abs(tabla_total$depth)      # CUIDADO: lo lee como factor
                        depth[depth==999]<-NA
                        #range(depth, na.rm=TRUE)
                        #sort(unique(depth, na.rm=TRUE))
                        datetime <- as.Date(tabla_total$fecha3, format="%Y-%m-%d")  
                        #range(datetime, na.rm=TRUE)
                        YEAR  <- years(datetime)   #    range(YEAR, na.rm=TRUE)
                        times <- as.numeric(format(datetime,"%m"))    
                        month<-sort(unique(times))
                        DAY <- days(datetime)
                        JULIAN <- julian(datetime, origin = as.Date("1997-01-01"))
                        #range(JULIAN, na.rm=TRUE)
                        
                        #colnames(tabla_total)
                        
                        # Create grid levels
                        indice1 <- c(long)   # range(indice1, na.rm=TRUE)
                        indice2 <- c(lati)   # range(indice2, na.rm=TRUE)
                        indice3 <- c(depth)  # range(indice3, na.rm=TRUE)
                        indice4 <- c(times)  # range(indice4, na.rm=TRUE)
                        factor1 <-cut(indice1, breaks=clases_log,   include.lowest = TRUE,right = TRUE)
                        factor2 <-cut(indice2, breaks=clases_lat,   include.lowest = TRUE,right = TRUE)
                        factor3 <-cut(indice3, breaks=clases_pro,   include.lowest = TRUE,right = TRUE)
                        factor4 <-cut(indice4, breaks=clases_time,  include.lowest = FALSE,right = TRUE)
                        niveles1<-levels(factor1)
                        niveles2<-levels(factor2)
                        niveles3<-levels(factor3)
                        niveles4<-levels(factor4) 
                        
                        ## Fill the grid
                        matriz_seasonal <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(niveles4),length(lambda)),
                                                 dimnames=list(x=niveles1,y=niveles2,z=niveles3,t=niveles4,l=lambda))
                        matriz_annual <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(lambda)),
                                               dimnames=list(x=niveles1, y=niveles2, z=niveles3, l=lambda))
                        #dim(matriz_seasonal)
                        
                        for (w in c(1:length(lambda))){
                          valor1 <- tabla_total[,5+w]
                          experiment <- data.frame(valor1, factor1, factor2, factor3, factor4)
                          #names(experiment)   
                          chla_diat <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(niveles4)),
                                             dimnames=list(x=niveles1,y=niveles2,z=niveles3,t=niveles4))                            
                          if(sum(complete.cases(experiment))>0){
                            res <- aggregate(valor1~factor1+factor2+factor3+factor4,experiment,mean, na.rm=TRUE, drop=F, simplify=FALSE)
                            #dim(res)
                            #head(res)
                            valor <- res$valor1
                            #length(valor)
                            valor[valor=="NULL"]<-NA
                            valor2<-unlist(valor)
                            #length(valor2)
                            
                            for (i in 1:nrow(res)){
                              chla_diat[which(niveles1==res[i,1]),which(niveles2==res[i,2]),which(niveles3==res[i,3]),which(niveles4==res[i,4])] <- valor2[i]} # end loop i
                          }
                          chla_diat[chla_diat==0] <-NA
                          res <- chla_diat        
                          dimnames(res)<- list(x=longitud[,2], y=latitud[,2], z=profun, t=month) 
                          
                          # SEASONAL
                          if (sum(!is.na(res))>=10){ matriz_seasonal[,,,,w] <- res }
                          # ANNUAL
                          res[res<=0] <- NA
                          res <- apply(res, MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)
                          #image2D(log10(res[,,1]), zlim=c(-3,0))
                          if (sum(!is.na(res))>=10){ matriz_annual[,,,w] <- res }
                          
                          print(paste(lambda[w]," ready"))
                        }  # end loop w lambda
######################### 
                        
            save(longitud,latitud,profun,month,lambda,  matriz_seasonal,  file=paste(s_dir,"media_seasonal_Ap.Rdata",sep=""))
            save(longitud,latitud,profun,lambda,        matriz_annual,    file=paste(s_dir,"media_annual_Ap.Rdata",sep=""))
                        
                        
                        
                        
## Chla-specific coefficients: Valente (surface), Astrid (depth)                        
##############################                               
                 # Tables 12nm
                 load(paste(p_dir,"Valente2_ChlaRssIOP_modbands12.RData", sep="")) 
                        #dim(tabla_total)
                        #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                        #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                        #dim(tabla_total)       
                        lambda <- as.numeric(substring(colnames(tabla_total),5,8))[which(substring(colnames(tabla_total),1,3)=="aph")]                   
                        tabla1<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="aph"))] 
                        colnames(tabla1)<-c("lista","fecha3","latitude","longitude","depth",lambda)
                        tabla1[,5]<-rep(5,length=nrow(tabla1))
                        tabla1[,c(6:18)][tabla1[,c(6:18)]=="NaN"]<-NA
                        cloro1<-tabla_total[,which(colnames(tabla_total)=="Chla_HPLC.mg.m..3.")] 
                        flag1<-tabla_total[,which(colnames(tabla_total)=="QFChla")] 
                        cloro1[flag1==1]<-NA
                        for (i in c(1:13)){
                          tabla1[,5+i]<-(tabla1[,5+i]/cloro1)}
                        
                 load(paste(p_dir,"Astrid5_IopHplc2nm_modbands12.RData", sep=""))   
                        #dim(tabla_total)
                        #ano<-as.numeric(as.character(years(tabla_total$fecha3)))
                        #tabla_total<-tabla_total[ano>=2007 & ano<=2018,]
                        #dim(tabla_total)        
                        tabla4<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="aph"))] 
                        colnames(tabla4)<-c("lista","fecha3","latitude","longitude","depth",lambda)
                        tabla4[,c(6:18)][tabla4[,c(6:18)]=="NaN"]<-NA
                            # Completar   
                            tabla_ap<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ap_"))] 
                            tabla_anap<-tabla_total[,c(1:5,which(substring(colnames(tabla_total),1,3)=="ana"))] 
                            tablita<-(tabla_ap[,6:18]-tabla_anap[,6:18])
                            tabla4[,6:18][is.na(tabla4[,6:18])]<-tablita[is.na(tabla4[,6:18])]
                        cloro4<-tabla_total[,which(colnames(tabla_total)=="MTChla")]  
                        flag3<-tabla_total[,which(colnames(tabla_total)=="MYFLAG3")]
                        cloro4[flag3==1]<-NA
                        for (i in c(1:13)){
                          tabla4[,5+i]<-(tabla4[,5+i]/cloro4)}
                        
                tabla_total<-rbind(tabla1,tabla4)
                        #dim(tabla_total)
                        #colnames(tabla_total)                
                        #sum(!is.na(tabla_total[,9]))
                        #dim(tabla_total[!is.na(tabla_total[,9]),])
                        #  plot(lambda,tabla_total[!is.na(tabla_total[,9]),][1,6:18], type="b", las=1, ylim=c(0,0.15))
                        #points(lambda,tabla_total[!is.na(tabla_total[,9]),][10,6:18], type="b")
                        #points(lambda,tabla_total[!is.na(tabla_total[,9]),][100,6:18], type="b")
                        #points(lambda,tabla_total[!is.na(tabla_total[,9]),][500,6:18], type="b")
                        
                        # Fill grid with ALL data
                        lati<-tabla_total$latitude    # plot(long, lati)
                        long<-tabla_total$longitude
                        depth<-abs(tabla_total$depth)      # CUIDADO: lo lee como factor
                        depth[depth==999]<-NA
                        #range(depth, na.rm=TRUE)
                        #sort(unique(depth, na.rm=TRUE))
                        datetime <- as.Date(tabla_total$fecha3, format="%Y-%m-%d")  
                        #range(datetime, na.rm=TRUE)
                        YEAR  <- years(datetime)   #    range(YEAR, na.rm=TRUE)
                        times <- as.numeric(format(datetime,"%m"))    
                        month<-sort(unique(times))
                        DAY <- days(datetime)
                        JULIAN <- julian(datetime, origin = as.Date("1997-01-01"))
                        #range(JULIAN, na.rm=TRUE)
                        
                        #colnames(tabla_total)
                        
                        # Create grid levels
                        indice1 <- c(long)   # range(indice1, na.rm=TRUE)
                        indice2 <- c(lati)   # range(indice2, na.rm=TRUE)
                        indice3 <- c(depth)  # range(indice3, na.rm=TRUE)
                        indice4 <- c(times)  # range(indice4, na.rm=TRUE)
                        factor1 <-cut(indice1, breaks=clases_log,   include.lowest = TRUE,right = TRUE)
                        factor2 <-cut(indice2, breaks=clases_lat,   include.lowest = TRUE,right = TRUE)
                        factor3 <-cut(indice3, breaks=clases_pro,   include.lowest = TRUE,right = TRUE)
                        factor4 <-cut(indice4, breaks=clases_time,  include.lowest = FALSE,right = TRUE)
                        niveles1<-levels(factor1)
                        niveles2<-levels(factor2)
                        niveles3<-levels(factor3)
                        niveles4<-levels(factor4) 
                        
                        ## Fill the grid
                        matriz_seasonal <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(niveles4),length(lambda)),
                                                 dimnames=list(x=niveles1,y=niveles2,z=niveles3,t=niveles4,l=lambda))
                        matriz_annual <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(lambda)),
                                               dimnames=list(x=niveles1, y=niveles2, z=niveles3, l=lambda))
                        #dim(matriz_seasonal)
                        
                        for (w in c(1:length(lambda))){
                          valor1 <- tabla_total[,5+w]
                          experiment <- data.frame(valor1, factor1, factor2, factor3, factor4)
                          #names(experiment)   
                          chla_diat <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(niveles4)),
                                             dimnames=list(x=niveles1,y=niveles2,z=niveles3,t=niveles4))                            
                          if(sum(complete.cases(experiment))>0){
                            res <- aggregate(valor1~factor1+factor2+factor3+factor4,experiment,mean, na.rm=TRUE, drop=F, simplify=FALSE)
                            #dim(res)
                            #head(res)
                            valor <- res$valor1
                            #length(valor)
                            valor[valor=="NULL"]<-NA
                            valor2<-unlist(valor)
                            #length(valor2)
                            
                            for (i in 1:nrow(res)){
                              chla_diat[which(niveles1==res[i,1]),which(niveles2==res[i,2]),which(niveles3==res[i,3]),which(niveles4==res[i,4])] <- valor2[i]} # end loop i
                          }
                          chla_diat[chla_diat==0] <-NA
                          res <- chla_diat        
                          dimnames(res)<- list(x=longitud[,2], y=latitud[,2], z=profun, t=month) 
                          
                          # SEASONAL
                          if (sum(!is.na(res))>=10){ matriz_seasonal[,,,,w] <- res }
                          # ANNUAL
                          res[res<=0] <- NA
                          res <- apply(res, MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)
                          #image2D(log10(res[,,1]), zlim=c(-3,0))
                          if (sum(!is.na(res))>=10){ matriz_annual[,,,w] <- res }
                          
                          print(paste(lambda[w]," ready"))
                        }  # end loop w lambda
############################## 
                        
             save(longitud,latitud,profun,month,lambda,  matriz_seasonal,  file=paste(s_dir,"media_seasonal_Aps.Rdata",sep=""))
             save(longitud,latitud,profun,lambda,        matriz_annual,    file=paste(s_dir,"media_annual_Aps.Rdata",sep=""))
                        
                        