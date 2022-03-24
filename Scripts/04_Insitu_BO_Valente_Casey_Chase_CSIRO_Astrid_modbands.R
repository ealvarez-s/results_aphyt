source("00_Summary_EDITME.R")
source("functions_misc.R")

########################################
##  In situ bio-optical data  GLOBAL  ##
########################################
###        average wavebands         ###
########################################
# Input: original .csv tables of rrs and iop samples (in /csvs/)
# Output:  .Rdata matrices with rrs and iop data averaged in the wavebands (25nm) of the model and collocated (in /total_tables/)
# one for Valente: "Valente2_ChlaRssIOP_modbands12.RData"
# one for CSIRO: "CSIRO_Iop_modbands12"    # Note: CSIRO IOPs and HPLC data are not collocated, there's no ID in the samples to match them
# other for Bracher : "Astrid6_IopHplc2nm_modbands12.RData"
# other for Chase: "Chase_RssHplc_modbands12.RData"
# other for Casey: "Casey_RssIop_modbands12.RDat"
# all are combined in the next script, it can be called with source at the end of this script

modelo_lambda <- c(400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0,600.0, 625.0, 650.0, 675.0, 700.0)
lambda_extremos<-sort(unique(c(modelo_lambda[-length(modelo_lambda)]+(diff(modelo_lambda)/2),400,700)))
#cbind(lambda_extremos[-length(lambda_extremos)],modelo_lambda,lambda_extremos[-1])
s_dir<-paste(global_path,"Dat_observations/optics/total_tables/",sep="")




########################
## Valente et al 2019 ##
########################     
p_dir<-paste(global_path,"Dat_observations/optics/csvs/datasets_Valente/",sep="")

  ## RRS
  reflectancia<-read.csv(paste(p_dir,"insitudb_rrs.csv", sep=""))
  lon<-reflectancia$Longitude
  lat<-reflectancia$Latitude
  #unique(reflectancia$Comment..dataset.)
  #sort(unique(reflectancia$Comment..subdataset.))
  #names(reflectancia)
  nombres_lambda <- names(reflectancia)[c(10:620)]
  numeros_lambda <- as.numeric(substring(nombres_lambda,16,nchar(nombres_lambda)-4))          
      reflectancia[c(10:620)][reflectancia[c(10:620)]<=0]<-NA
      reflectancia[c(10:620)][reflectancia[c(10:620)]>0.15]<-NA
  nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
  
           for (i in c(1:length(modelo_lambda))){
                  cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
                  if (length(cuales[!is.na(cuales)])>1) {
                    res<-rowMeans(reflectancia[,9+cuales],na.rm=T)
                    res[res==0]<-NA
                    nueva_reflec[,i]<-res
                  } else if (length(cuales[!is.na(cuales)])==1) {
                    res<-reflectancia[,9+cuales]
                    res[res==0]<-NA
                    nueva_reflec[,i]<-res
                  } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
           }

                # ## FIGURA
                # par(mfrow=c(1,2))
                # par(mar=c(4,4,1,1))
                # 
                # reflectancia[reflectancia==0]<-NA
                # matplot(numeros_lambda, t(reflectancia[,c(5:615)]), type="b", pch=19,
                #         col=rainbow(10), xlim=c(300,800), ylim=c(0,0.1), las=1) 
                # matplot(modelo_lambda, t(nueva_reflec), type="b", pch=19,
                #         col=rainbow(10), xlim=c(300,800), ylim=c(0,0.1), las=1)
                # 
                # 
            colnames(nueva_reflec) <- paste("rss_",modelo_lambda, sep="")
            # colnames(reflectancia)[c(1:9,621:624)]
            tabla_rss <- cbind(reflectancia[,c(1:9,621:624)],nueva_reflec) 
            # colnames(tabla_rss)
            # unique(tabla_rss$Comment..dataset.)
            # save(tabla_rss,modelo_lambda, paste(p_dir,"insitudb_rrs_modbands12.RData", sep=""))

  ## APH
  reflectancia<-read.csv(paste(p_dir,"insitudb_aph_withoutAWI&NOMAD.csv", sep=""))
  #colnames(reflectancia)
  lon<-reflectancia$Longitude
  lat<-reflectancia$Latitude
  #unique(reflectancia$Comment..dataset.)
  #sort(unique(reflectancia$Comment..subdataset.))
  #names(reflectancia)
  nombres_lambda <- names(reflectancia)[c(11:560)]
  numeros_lambda <- as.numeric(substring(nombres_lambda,nchar(nombres_lambda)-6,nchar(nombres_lambda)-4))          
    reflectancia[c(11:560)][reflectancia[c(11:560)]<=0.0001]<-NA
    reflectancia[c(11:560)][reflectancia[c(11:560)]>10]<-NA
  nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))

          for (i in c(1:length(modelo_lambda))){
            cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
            if (length(cuales[!is.na(cuales)])>1) {
              res<-rowMeans(reflectancia[,10+cuales],na.rm=T)
              res[res==0]<-NA
              nueva_reflec[,i]<-res
            } else if (length(cuales[!is.na(cuales)])==1) {
              res<-reflectancia[,10+cuales]
              res[res==0]<-NA
              nueva_reflec[,i]<-res
            } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
          }
          # ## FIGURA
          # par(mfrow=c(1,2))
          # par(mar=c(4,4,1,1))
          # reflectancia[reflectancia==0]<-NA
          # matplot(numeros_lambda, t(reflectancia[,c(5:554)]), type="b", pch=19,col=rainbow(10), xlim=c(300,800), ylim=c(0,0.1), las=1) 
          # matplot(modelo_lambda, t(nueva_reflec), type="b", pch=19,col=rainbow(10), xlim=c(300,800), ylim=c(0,0.3), las=1)
        colnames(nueva_reflec) <- paste("aph_",modelo_lambda, sep="")
        #colnames(reflectancia)[c(1:10,561:564)]
        tabla_aph <- cbind(reflectancia[,c(1:10,561:564)],nueva_reflec) 
        # colnames(tabla_aph)

          
  ## ACDOM
  reflectancia <- read.csv(paste(p_dir,"insitudb_acdom.csv", sep=""))
  lon<-reflectancia$Longitude
  lat<-reflectancia$Latitude
  #unique(reflectancia$Comment..dataset.)
  #sort(unique(reflectancia$Comment..subdataset.)) 
  #names(reflectancia)
  nombres_lambda <- names(reflectancia)[c(10:42)]
  numeros_lambda <- as.numeric(substring(nombres_lambda,nchar(nombres_lambda)-29,nchar(nombres_lambda)-27))          
     reflectancia[c(10:42)][reflectancia[c(10:42)]<=0.0001]<-NA
     reflectancia[c(10:42)][reflectancia[c(10:42)]>10]<-NA
  nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
  #dim(nueva_reflec)

          for (i in c(1:length(modelo_lambda))){
            cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
            if (length(cuales[!is.na(cuales)])>1) {
              res<-rowMeans(reflectancia[,9+cuales],na.rm=T)
              res[res==0]<-NA
              nueva_reflec[,i]<-res
            } else if (length(cuales[!is.na(cuales)])==1) {
              res<-reflectancia[,9+cuales]
              res[res==0]<-NA
              nueva_reflec[,i]<-res
            } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
          }
          
          # ## FIGURA
          # par(mfrow=c(1,2))
          # par(mar=c(4,4,1,1))
          # reflectancia[reflectancia==0]<-NA
          # matplot(numeros_lambda, t(reflectancia[,c(5:37)]), type="b", pch=19,col=rainbow(10), xlim=c(300,800), ylim=c(0,0.5), las=1) 
          # matplot(modelo_lambda, t(nueva_reflec), type="b", pch=19,col=rainbow(10), xlim=c(300,800), ylim=c(0,0.3), las=1)
        colnames(nueva_reflec) <- paste("adg_",modelo_lambda, sep="")
        # colnames(reflectancia)[c(1:9,43:46)]
        tabla_adg <- cbind(reflectancia[,c(1:9,43:46)],nueva_reflec) 
        # colnames(tabla_adg)

          
          
  ## BBP
  reflectancia <- read.csv(paste(p_dir,"insitudb_bbp.csv", sep=""))
  lon<-reflectancia$Longitude
  lat<-reflectancia$Latitude
  #unique(reflectancia$Comment..dataset.)
  #sort(unique(reflectancia$Comment..subdataset.))
  #names(reflectancia)
  nombres_lambda <- names(reflectancia)[c(10:41)]
  numeros_lambda <- as.numeric(substring(nombres_lambda,nchar(nombres_lambda)-6,nchar(nombres_lambda)-4))          
     reflectancia[c(10:41)][reflectancia[c(10:41)]<=0.0001]<-NA
     reflectancia[c(10:41)][reflectancia[c(10:41)]>10]<-NA
  nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
  #dim(nueva_reflec)

          for (i in c(1:length(modelo_lambda))){
            cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
            if (length(cuales[!is.na(cuales)])>1) {
              res<-rowMeans(reflectancia[,9+cuales],na.rm=T)
              res[res==0]<-NA
              nueva_reflec[,i]<-res
            } else if (length(cuales[!is.na(cuales)])==1) {
              res<-reflectancia[,9+cuales]
              res[res==0]<-NA
              nueva_reflec[,i]<-res
            } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
          }
          
          # ## FIGURA
          # par(mfrow=c(1,2))
          # par(mar=c(4,4,1,1))
          # reflectancia[reflectancia==0]<-NA
          # matplot(numeros_lambda, t(reflectancia[,c(5:36)]), type="b", pch=19, col=rainbow(10), xlim=c(300,800), ylim=c(0,0.1), las=1) 
          # matplot(modelo_lambda, t(nueva_reflec), type="b", pch=19, col=rainbow(10), xlim=c(300,800), ylim=c(0,0.1), las=1)
          # 
       colnames(nueva_reflec) <- paste("bbp_",modelo_lambda, sep="")
        # colnames(reflectancia)[c(1:9,42:45)]
       tabla_bbp <- cbind(reflectancia[,c(1:9,42:45)],nueva_reflec) 
        # colnames(tabla_bbp)

  ## Chla 
  clorofila<-read.csv(paste(p_dir,"insitudb_chla.csv", sep=""))
  #names(clorofila)
  #dim(clorofila)
  clorofila[c(10)][clorofila[c(10)]<=0.001]<-NA
  clorofila[c(10)][clorofila[c(10)]>100]<-NA  
  #colnames(clorofila)
  #lon<-clorofila$Longitude
  #lat<-clorofila$Latitude

########################  
  
##########################################
### Match all tables: Rrs, IOPs & Chla ###  
##########################################      

  lista <- sort(unique(c(clorofila$ID,tabla_rss$ID,tabla_aph$ID,tabla_adg$ID,tabla_bbp$ID)))
  #length(lista)

  cloro_ordered <- clorofila[match(lista,clorofila$ID),] 
  tabla_rss_ordered <- tabla_rss[match(lista,tabla_rss$ID),]
  tabla_aph_ordered <- tabla_aph[match(lista,tabla_aph$ID),] 
  tabla_adg_ordered <- tabla_adg[match(lista,tabla_adg$ID),] 
  tabla_bbp_ordered <- tabla_bbp[match(lista,tabla_bbp$ID),] 
  #tabla_kd_ordered  <- tabla_kd[match(lista,tabla_kd$ID),]
  
  fechas <- data.frame(cloro_ordered$Date.Time,      tabla_rss_ordered$Date.Time,
                       tabla_aph_ordered$Date.Time,  tabla_adg_ordered$Date.Time,
                       tabla_bbp_ordered$Date.Time)
  fecha <- apply(fechas, MARGIN=1, FUN=unique, na.rm=TRUE)
  fecha2 <- unlist(fecha)
  fecha3 <- fecha2[!is.na(fecha2)]
  #length(fecha3)
  #range(fecha3,na.rm=T)
  
  latitude <- apply(cbind(cloro_ordered$Latitude,      tabla_rss_ordered$Latitude,
                          tabla_aph_ordered$Latitude,  tabla_adg_ordered$Latitude,
                          tabla_bbp_ordered$Latitude),
                          MARGIN=1, FUN=mean, na.rm=TRUE)
  #sum(!is.na(latitude))
  #length(latitude)
  
  longitude <- apply(cbind(cloro_ordered$Longitude,      tabla_rss_ordered$Longitude,
                           tabla_aph_ordered$Longitude,  tabla_adg_ordered$Longitude,
                           tabla_bbp_ordered$Longitude),
                           MARGIN=1,FUN=mean, na.rm=TRUE)
  #sum(!is.na(longitude))
  #length(longitude)
  
  depth <- apply(cbind(cloro_ordered$Depth_water.m.,      tabla_rss_ordered$Depth.water..m.,
                       tabla_aph_ordered$Depth.water..m., tabla_adg_ordered$Depth.water..m.,
                       tabla_bbp_ordered$Depth.water..m.),
                       MARGIN=1,FUN=mean, na.rm=TRUE)
  #sum(!is.na(depth))
  #sort(unique(depth))
  depth<-rep(5,length(depth))
  
  
  QFTime <- apply(cbind(cloro_ordered$QF_time,      tabla_rss_ordered$QF.time,
                        tabla_aph_ordered$QF.time,  tabla_adg_ordered$QF.time,
                        tabla_bbp_ordered$QF.time),
                        MARGIN=1,FUN=mean, na.rm=TRUE)
  #length(QFTime)
  QFTime[is.nan(QFTime)] <- NA
  #unique(QFTime)
  
  QFChla <- apply(cbind(cloro_ordered$QF_Chl,      tabla_rss_ordered$QF_Chl,
                        tabla_aph_ordered$QF_Chl,  tabla_adg_ordered$QF_Chl,
                        tabla_bbp_ordered$QF_Chl),
                        MARGIN=1,FUN=mean, na.rm=TRUE)
  #length(QFChla)    
  QFChla[is.nan(QFChla)] <- NA
  #unique(QFChla)
  
  
  metadata <- data.frame(lista, fecha3, latitude, longitude, depth, QFChla, QFTime)
  
  tabla_total <- data.frame(metadata,cloro_ordered[,10:11],
                            tabla_rss_ordered[,14:ncol(tabla_rss_ordered)],
                            tabla_aph_ordered[,15:ncol(tabla_aph_ordered)],
                            tabla_adg_ordered[,14:ncol(tabla_adg_ordered)],
                            tabla_bbp_ordered[,14:ncol(tabla_bbp_ordered)])
  
  #dim(tabla_total)
  #length(QFTime)
  #length(QFChla)
  tabla_total<-tabla_total[is.na(QFTime) & is.na(QFChla),]
  #head(tabla_total)
  tabla_total<-tabla_total[!is.nan(tabla_total$longitude),]
  tabla_total<-tabla_total[!is.na(tabla_total$longitude),]
  #dim(tabla_total)  
  #colnames(tabla_total)
##########################################
save(tabla_total, file=paste(s_dir,"Valente2_ChlaRssIOP_modbands12.RData", sep=""))    


  
  

########################
##  Chase et al 2017  ##
########################         
p_dir<-paste(global_path,"Dat_observations/optics/csvs/",sep="")
  
#### RRS
          reflectancia<-read.csv(paste(p_dir,"TOTAL_Chase_rrs.csv", sep=""))
          #reflectancia$Depth
          #dim(reflectancia)
          #names(reflectancia)
          lon<-reflectancia$Longitude..degrees_east.
          lat<-reflectancia$Latitude..degrees_north.
          #unique(reflectancia$CRUISE)
          #unique(reflectancia$Provider)            
          #dim(reflectancia)            
          #names(reflectancia)
          nombres_lambda <- names(reflectancia)[c(19:155)]
          numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))
              reflectancia[c(19:155)][reflectancia[c(19:155)]==(-999)]<-NA
              reflectancia[c(19:155)][reflectancia[c(19:155)]<=0]<-NA
              reflectancia[c(19:155)][reflectancia[c(19:155)]>0.15]<-NA
          nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
          
          for (i in c(1:length(modelo_lambda))){
            cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
            if (length(cuales[!is.na(cuales)])>1) {
              res<-rowMeans(reflectancia[,18+cuales],na.rm=T)
              res[is.nan(res)]<-NA
              res[res==0]<-NA
              nueva_reflec[,i]<-res
            } else if (length(cuales[!is.na(cuales)])==1) {
              res<-reflectancia[,18+cuales]
              res[is.nan(res)]<-NA
              res[res==0]<-NA
              nueva_reflec[,i]<-res
            } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
          }

                      ## FIGURA
                      #par(mfrow=c(1,2))
                      #par(mar=c(4,4,1,1))
                      # reflectancia[reflectancia==0]<-NA
                      #matplot(numeros_lambda, t(reflectancia[,c(19:619)]), type="b", pch=19,
                      #        col=rainbow(10), xlim=c(300,800), ylim=c(0,0.02), las=1) 
                      #matplot(modelo_lambda, t(nueva_reflec), type="b", pch=19,
                      #        col=rainbow(10), xlim=c(300,800), ylim=c(0,0.02), las=1)
                      
                colnames(nueva_reflec) <- paste("rss_",modelo_lambda, sep="")
                # colnames(reflectancia)  reflectancia$Depth
                tabla_rss <- cbind(reflectancia[,c(1:18)],nueva_reflec) 
                # colnames(tabla_rss)
                # unique(tabla_rss$Provider)
                # save(tabla_rss,modelo_lambda, paste(p_dir,"insitudb_rrs_modbands12.RData", sep=""))

#### HPLC 
       clorofila<-read.csv(paste(p_dir,"TOTAL_Chase_hplc.csv", sep=""))
       #clorofila$Depth
       #names(clorofila)
       #dim(clorofila)
       clorofila[c(16)][clorofila[c(16)]<=0.001]<-NA
       clorofila[c(16)][clorofila[c(16)]>100]<-NA  
       #dim(clorofila)
       #colnames(clorofila)
       #lon<-clorofila$Longitude
       #lat<-clorofila$Latitude
########################

#### Match all tables  
#####################       
  #dim(tabla_rss)
  #dim(clorofila) 

  lista <- sort(unique(c(clorofila$ID,tabla_rss$ID)))
  #length(lista)
  
  cloro_ordered      <- clorofila[match(lista,clorofila$ID),] 
  tabla_rss_ordered  <- tabla_rss[match(lista,tabla_rss$ID),]

  #names(tabla_rss_ordered)  tabla_rss_ordered$Depth   tabla_rss$Depth
  #names(cloro_ordered)

  #colnames(cloro_ordered)
  #head(cloro_ordered)
  cloro_ordered_date<-paste(cloro_ordered[,5],"-",sprintf("%02d", cloro_ordered[,6]),"-",sprintf("%02d", cloro_ordered[,7]),sep="")
                      #,"T",sprintf("%02d", cloro_ordered[,9]),":",sprintf("%02d", cloro_ordered[,10]),sep="")
  cloro_ordered_date[is.na(cloro_ordered[,5])]<-NA
  tabla_rss_ordered_date<-paste(tabla_rss_ordered[,5],"-",sprintf("%02d", tabla_rss_ordered[,6]),"-",sprintf("%02d", tabla_rss_ordered[,7]),sep="")
                      #,"T",sprintf("%02d", tabla_rss_ordered[,9]),":",sprintf("%02d", tabla_rss_ordered[,10]),sep="")  
  tabla_rss_ordered_date[is.na(tabla_rss_ordered[,5])]<-NA
  
  fechas <- data.frame(cloro_ordered_date, tabla_rss_ordered_date)
  fecha  <- apply(fechas, MARGIN=1, function(x)unique(x[!is.na(x)]))
  fecha3<-as.character(fecha)
  fecha3[fecha3=="character(0)"]<-NA
  
  #sum(!is.na(fecha))
  #length(fecha)
  
  latitude <- apply(cbind(cloro_ordered$Latitude,    tabla_rss_ordered$Latitude),    MARGIN=1, FUN=mean, na.rm=TRUE)
  #sum(!is.na(latitude))
  #length(latitude)
  
  longitude <- apply(cbind(cloro_ordered$Longitude,  tabla_rss_ordered$Longitude),MARGIN=1,FUN=mean, na.rm=TRUE)
  #sum(!is.na(longitude))
  #length(longitude)
  
  depth <- apply(cbind(cloro_ordered$Depth,    tabla_rss_ordered$Depth),MARGIN=1, FUN=mean, na.rm=TRUE)
  #sum(!is.na(depth))
  #sort(unique(depth))
  #length(depth)

  metadata <- data.frame(lista, fecha3, latitude, longitude, depth)  
  #colnames(cloro_ordered)
  #colnames(tabla_rss_ordered)
  tabla_total <- data.frame(metadata, cloro_ordered[,16:ncol(cloro_ordered)],
                            tabla_rss_ordered[,19:ncol(tabla_rss_ordered)])
  #dim(tabla_total)
  #head(tabla_total)
  tabla_total<-tabla_total[!is.nan(tabla_total$longitude),]
  #dim(tabla_total)
#####################  
save(tabla_total, file=paste(s_dir,"Chase_RssHplc_modbands12.RData", sep=""))    
  

  

  

########################
##  Casey et al 2019  ##
########################     
p_dir<-paste(global_path,"Dat_observations/optics/csvs/",sep="")
  
#### RRS
    reflectancia<-read.csv(paste(p_dir,"TOTAL_Casey_rrs.csv", sep=""))
    #names(reflectancia)
    lon<-reflectancia$Longitude..degrees_east.
    lat<-reflectancia$Latitude..degrees_north.
    #unique(reflectancia$CRUISE)
    #unique(reflectancia$Provider)            

    nombres_lambda <- names(reflectancia)[c(20:620)]
    numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
    reflectancia[c(20:620)][reflectancia[c(20:620)]==(-999)]<-NA
      reflectancia[c(20:620)][reflectancia[c(20:620)]<=0]<-NA
      reflectancia[c(20:620)][reflectancia[c(20:620)]>0.15]<-NA
    nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
    
    for (i in c(1:length(modelo_lambda))){
      cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
      if (length(cuales[!is.na(cuales)])>1) {
        res<-rowMeans(reflectancia[,19+cuales],na.rm=T)
        res[is.nan(res)]<-NA
        res[res==0]<-NA
        nueva_reflec[,i]<-res
      } else if (length(cuales[!is.na(cuales)])==1) {
        res<-reflectancia[,19+cuales]
        res[is.nan(res)]<-NA
        res[res==0]<-NA
        nueva_reflec[,i]<-res
      } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
    }

            ## FIGURA
            #par(mfrow=c(1,2))
            #par(mar=c(4,4,1,1))
            # reflectancia[reflectancia==0]<-NA
            #matplot(numeros_lambda, t(reflectancia[,c(20:620)]), type="b", pch=19,
            #        col=rainbow(10), xlim=c(300,800), ylim=c(0,0.02), las=1) 
            #matplot(modelo_lambda, t(nueva_reflec), type="b", pch=19,
            #        col=rainbow(10), xlim=c(300,800), ylim=c(0,0.02), las=1)
            
         colnames(nueva_reflec) <- paste("rss_",modelo_lambda, sep="")
         #colnames(reflectancia)
         tabla_rss <- cbind(reflectancia[,c(1:19)],nueva_reflec) 
         #colnames(tabla_rss)
         #unique(tabla_rss$Provider)
         #save(tabla_rss,modelo_lambda, paste(p_dir,"insitudb_rrs_modbands12.RData", sep=""))
            
#### AP
    absorcion<-read.csv(paste(p_dir,"TOTAL_Casey_ap.csv", sep=""))
    #names(absorcion)
    #dim(absorcion)
    lon<-absorcion$Longitude..degrees_east.
    lat<-absorcion$Latitude..degrees_north.
    #unique(absorcion$CRUISE)
    #unique(absorcion$Provider)            
    
    nombres_lambda <- names(absorcion)[c(19:619)]
    numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
    absorcion[c(19:619)][absorcion[c(19:619)]==(-999)]<-NA
      absorcion[c(19:619)][absorcion[c(19:619)]<=0.0001]<-NA
      absorcion[c(19:619)][absorcion[c(19:619)]>10]<-NA
    nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(absorcion))
    
    for (i in c(1:length(modelo_lambda))){
      cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
      if (length(cuales[!is.na(cuales)])>1) {
        res<-rowMeans(absorcion[,18+cuales],na.rm=T)
        res[is.nan(res)]<-NA
        res[res==0]<-NA
        nueva_reflec[,i]<-res
      } else if (length(cuales[!is.na(cuales)])==1) {
        res<-absorcion[,18+cuales]
        res[is.nan(res)]<-NA
        res[res==0]<-NA
        nueva_reflec[,i]<-res
      } else {nueva_reflec[,i]<-rep(NA,length=nrow(absorcion))}
    }

          ## FIGURA
          #par(mfrow=c(1,2))
          #par(mar=c(4,4,1,1))
          
          # absorcion[absorcion==0]<-NA
          #matplot(numeros_lambda, t(absorcion[,c(19:619)]), type="b", pch=19,col=rainbow(10), xlim=c(300,800), ylim=c(0,0.04), las=1) 
          #matplot(modelo_lambda, t(nueva_reflec), type="b", pch=19,col=rainbow(10), xlim=c(300,800), ylim=c(0,0.04), las=1)
          
      colnames(nueva_reflec) <- paste("ap_",modelo_lambda, sep="")
      #colnames(absorcion)
      tabla_ap <- cbind(absorcion[,c(1:18)],nueva_reflec) 
      #colnames(tabla_ap)
      #unique(tabla_ap$Provider)
      #save(tabla_aph,modelo_lambda, paste(p_dir,"insitudb_rrs_modbands12.RData", sep=""))


#### APH
    absorcion<-read.csv(paste(p_dir,"TOTAL_Casey_aph.csv", sep=""))
    #names(absorcion)
    lon<-absorcion$Longitude..degrees_east.
    lat<-absorcion$Latitude..degrees_north.
    #unique(absorcion$CRUISE)
    #unique(absorcion$Provider)            
    
    nombres_lambda <- names(absorcion)[c(19:619)]
    numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
    absorcion[c(19:619)][absorcion[c(19:619)]==(-999)]<-NA
        absorcion[c(19:619)][absorcion[c(19:619)]<=0.0001]<-NA
        absorcion[c(19:619)][absorcion[c(19:619)]>10]<-NA
    nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(absorcion))
    
    for (i in c(1:length(modelo_lambda))){
      cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
      if (length(cuales[!is.na(cuales)])>1) {
        res<-rowMeans(absorcion[,18+cuales],na.rm=T)
        res[is.nan(res)]<-NA
        res[res==0]<-NA
        nueva_reflec[,i]<-res
      } else if (length(cuales[!is.na(cuales)])==1) {
        res<-absorcion[,18+cuales]
        res[is.nan(res)]<-NA
        res[res==0]<-NA
        nueva_reflec[,i]<-res
      } else {nueva_reflec[,i]<-rep(NA,length=nrow(absorcion))}
    }
  
            ## FIGURA
            #par(mfrow=c(1,2))
            #par(mar=c(4,4,1,1))
            # absorcion[absorcion==0]<-NA
            #matplot(numeros_lambda, t(absorcion[,c(19:619)]), type="b", pch=19,col=rainbow(10), xlim=c(300,800), ylim=c(0,0.5), las=1) 
            #matplot(modelo_lambda, t(nueva_reflec), type="b", pch=19,col=rainbow(10), xlim=c(300,800), ylim=c(0,0.5), las=1)
            
       colnames(nueva_reflec) <- paste("aph_",modelo_lambda, sep="")
       # colnames(absorcion)
       tabla_aph <- cbind(absorcion[,c(1:18)],nueva_reflec) 
       #colnames(tabla_aph)
       #unique(tabla_cp$Provider)
       #save(tabla_aph,modelo_lambda, paste(p_dir,"insitudb_rrs_modbands12.RData", sep=""))
            
       
#### ANAP
       absorcion<-read.csv(paste(p_dir,"TOTAL_Casey_anap.csv", sep=""))
       #names(absorcion)
       lon<-absorcion$Longitude..degrees_east.
       lat<-absorcion$Latitude..degrees_north.
       #unique(absorcion$CRUISE)
       #unique(absorcion$Provider)            
       
       nombres_lambda <- names(absorcion)[c(19:619)]
       numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
       absorcion[c(19:619)][absorcion[c(19:619)]==(-999)]<-NA
          absorcion[c(19:619)][absorcion[c(19:619)]<=0.0001]<-NA
          absorcion[c(19:619)][absorcion[c(19:619)]>10]<-NA
       nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(absorcion))
       
       for (i in c(1:length(modelo_lambda))){
         cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
         if (length(cuales[!is.na(cuales)])>1) {
           res<-rowMeans(absorcion[,18+cuales],na.rm=T)
           res[is.nan(res)]<-NA
           res[res==0]<-NA
           nueva_reflec[,i]<-res
         } else if (length(cuales[!is.na(cuales)])==1) {
           res<-absorcion[,18+cuales]
           res[is.nan(res)]<-NA
           res[res==0]<-NA
           nueva_reflec[,i]<-res
         } else {nueva_reflec[,i]<-rep(NA,length=nrow(absorcion))}
       }
       
       # absorcion[absorcion==0]<-NA
       colnames(nueva_reflec) <- paste("anap_",modelo_lambda, sep="")
       tabla_anap <- cbind(absorcion[,c(1:18)],nueva_reflec) 
       #colnames(tabla_anap)
       #unique(tabla_cp$Provider)
       #save(tabla_aph,modelo_lambda, paste(p_dir,"insitudb_rrs_modbands12.RData", sep=""))
       
#### ACDOM
       absorcion<-read.csv(paste(p_dir,"TOTAL_Casey_acdom.csv", sep=""))
       #names(absorcion)
       lon<-absorcion$Longitude..degrees_east.
       lat<-absorcion$Latitude..degrees_north.
       #unique(absorcion$CRUISE)
       #unique(absorcion$Provider)            
       
       nombres_lambda <- names(absorcion)[c(19:619)]
       numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
       absorcion[c(19:619)][absorcion[c(19:619)]==(-999)]<-NA
          absorcion[c(19:619)][absorcion[c(19:619)]<=0.0001]<-NA
          absorcion[c(19:619)][absorcion[c(19:619)]>10]<-NA
       nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(absorcion))
       
       for (i in c(1:length(modelo_lambda))){
         cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
         if (length(cuales[!is.na(cuales)])>1) {
           res<-rowMeans(absorcion[,18+cuales],na.rm=T)
           res[is.nan(res)]<-NA
           res[res==0]<-NA
           nueva_reflec[,i]<-res
         } else if (length(cuales[!is.na(cuales)])==1) {
           res<-absorcion[,18+cuales]
           res[is.nan(res)]<-NA
           res[res==0]<-NA
           nueva_reflec[,i]<-res
         } else {nueva_reflec[,i]<-rep(NA,length=nrow(absorcion))}
       }
       
       # absorcion[absorcion==0]<-NA
       colnames(nueva_reflec) <- paste("acdom_",modelo_lambda, sep="")
       tabla_acdom <- cbind(absorcion[,c(1:18)],nueva_reflec) 
       #colnames(tabla_acdom)
       #unique(tabla_cp$Provider)
       #save(tabla_aph,modelo_lambda, paste(p_dir,"insitudb_rrs_modbands12.RData", sep=""))       
       
#### BBP
       absorcion<-read.csv(paste(p_dir,"TOTAL_Casey_bbp.csv", sep=""))
       #names(absorcion)
       lon<-absorcion$Longitude..degrees_east.
       lat<-absorcion$Latitude..degrees_north.
       #unique(absorcion$CRUISE)
       #unique(absorcion$Provider)            
       
       nombres_lambda <- names(absorcion)[c(19:619)]
       numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
       absorcion[c(19:619)][absorcion[c(19:619)]==(-999)]<-NA
          absorcion[c(19:619)][absorcion[c(19:619)]<=0.0001]<-NA
          absorcion[c(19:619)][absorcion[c(19:619)]>10]<-NA
       nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(absorcion))
       
       for (i in c(1:length(modelo_lambda))){
         cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
         if (length(cuales[!is.na(cuales)])>1) {
           res<-rowMeans(absorcion[,18+cuales],na.rm=T)
           res[is.nan(res)]<-NA
           res[res==0]<-NA
           nueva_reflec[,i]<-res
         } else if (length(cuales[!is.na(cuales)])==1) {
           res<-absorcion[,18+cuales]
           res[is.nan(res)]<-NA
           res[res==0]<-NA
           nueva_reflec[,i]<-res
         } else {nueva_reflec[,i]<-rep(NA,length=nrow(absorcion))}
       }
       
       # absorcion[absorcion==0]<-NA
       colnames(nueva_reflec) <- paste("bbp_",modelo_lambda, sep="")
       tabla_bbp <- cbind(absorcion[,c(1:18)],nueva_reflec) 
       #colnames(tabla_bbp)
       #unique(tabla_cp$Provider)
       #save(tabla_aph,modelo_lambda, paste(p_dir,"insitudb_rrs_modbands12.RData", sep=""))       
       
########################            

#### Match all tables    
#####################      
      #dim(tabla_rss)
      #dim(tabla_ap) 
      #tabla_ap$ID%in%tabla_rss$ID
      #tabla_rss$ID%in%tabla_ap$ID    

      lista <- sort(unique(c(tabla_rss$ID,tabla_ap$ID,tabla_aph$ID,tabla_anap$ID,tabla_acdom$ID,tabla_bbp$ID)))
      #length(lista)

      tabla_rss_ordered <- tabla_rss[match(lista,tabla_rss$ID),]
      tabla_ap_ordered  <- tabla_ap[match(lista,tabla_ap$ID),] 
      tabla_aph_ordered  <- tabla_aph[match(lista,tabla_aph$ID),] 
      tabla_anap_ordered  <- tabla_anap[match(lista,tabla_anap$ID),] 
      tabla_acdom_ordered  <- tabla_acdom[match(lista,tabla_acdom$ID),] 
      tabla_bbp_ordered  <- tabla_bbp[match(lista,tabla_bbp$ID),] 
  
      ## Date.Time
      #colnames(tabla_rss_ordered)
      tabla_rss_ordered_date<-paste(tabla_rss_ordered[,10],"-",sprintf("%02d", tabla_rss_ordered[,11]),"-",sprintf("%02d", tabla_rss_ordered[,12]),sep="")
                                    #,"T",sprintf("%02d", tabla_rss_ordered[,13]),":",sprintf("%02d", tabla_rss_ordered[,14]),sep="")
      tabla_rss_ordered_date[is.na(tabla_rss_ordered[,10])]<-NA
      #colnames(tabla_ap_ordered)
      tabla_ap_ordered_date<-paste(tabla_ap_ordered[,10],"-",sprintf("%02d", tabla_ap_ordered[,11]),"-",sprintf("%02d", tabla_ap_ordered[,12]),sep="")
                                   #,"T",sprintf("%02d", tabla_ap_ordered[,13]),":",sprintf("%02d", tabla_ap_ordered[,14]),sep="")  
      tabla_ap_ordered_date[is.na(tabla_ap_ordered[,10])]<-NA
      #colnames(tabla_aph_ordered)
      tabla_aph_ordered_date<-paste(tabla_aph_ordered[,10],"-",sprintf("%02d", tabla_aph_ordered[,11]),"-",sprintf("%02d", tabla_aph_ordered[,12]),sep="")
                                    #,"T",sprintf("%02d", tabla_aph_ordered[,13]),":",sprintf("%02d", tabla_aph_ordered[,14]),sep="")  
      tabla_aph_ordered_date[is.na(tabla_aph_ordered[,10])]<-NA
      #colnames(tabla_anap_ordered)
      tabla_anap_ordered_date<-paste(tabla_anap_ordered[,10],"-",sprintf("%02d", tabla_anap_ordered[,11]),"-",sprintf("%02d", tabla_anap_ordered[,12]),sep="")
                                     #,"T",sprintf("%02d", tabla_anap_ordered[,13]),":",sprintf("%02d", tabla_anap_ordered[,14]),sep="") 
      tabla_anap_ordered_date[is.na(tabla_anap_ordered[,10])]<-NA
      #colnames(tabla_acdom_ordered)
      tabla_acdom_ordered_date<-paste(tabla_acdom_ordered[,10],"-",sprintf("%02d", tabla_acdom_ordered[,11]),"-",sprintf("%02d", tabla_acdom_ordered[,12]),sep="")
                                      #,"T",sprintf("%02d", tabla_acdom_ordered[,13]),":",sprintf("%02d", tabla_acdom_ordered[,14]),sep="")
      tabla_acdom_ordered_date[is.na(tabla_acdom_ordered[,10])]<-NA
      #colnames(tabla_bbp_ordered)
      tabla_bbp_ordered_date<-paste(tabla_bbp_ordered[,10],"-",sprintf("%02d", tabla_bbp_ordered[,11]),"-",sprintf("%02d", tabla_bbp_ordered[,12]),sep="")
                                    #,"T",sprintf("%02d", tabla_bbp_ordered[,13]),":",sprintf("%02d", tabla_bbp_ordered[,14]),sep="") 
      tabla_bbp_ordered_date[is.na(tabla_bbp_ordered[,10])]<-NA
      
      fechas <- data.frame(tabla_rss_ordered_date,   tabla_ap_ordered_date,     tabla_aph_ordered_date,
                           tabla_anap_ordered_date,  tabla_acdom_ordered_date,  tabla_bbp_ordered_date)
      fecha  <- apply(fechas, MARGIN=1, function(x)unique(x[!is.na(x)]))
      fecha3<-as.character(fecha)
      fecha3[fecha3=="character(0)"]<-NA  
      #length(fecha3)

      latitude <- apply(cbind(tabla_rss_ordered$Latitude,   tabla_ap_ordered$Latitude,
                              tabla_aph_ordered$Latitude,   tabla_anap_ordered$Latitude,
                              tabla_acdom_ordered$Latitude, tabla_bbp_ordered$Latitude),
                              MARGIN=1, FUN=mean, na.rm=TRUE)
      #sum(!is.na(latitude))
      #length(latitude)
      
      longitude <- apply(cbind(tabla_rss_ordered$Longitude, tabla_ap_ordered$Longitude,
                               tabla_aph_ordered$Longitude, tabla_anap_ordered$Longitude,
                               tabla_acdom_ordered$Longitude, tabla_bbp_ordered$Longitude),
                               MARGIN=1,FUN=mean, na.rm=TRUE)
      #sum(!is.na(longitude))
      #length(longitude)
      
      depth <- apply(cbind(tabla_rss_ordered$Depth,   tabla_ap_ordered$Depth,
                           tabla_aph_ordered$Depth,   tabla_anap_ordered$Depth,
                           tabla_acdom_ordered$Depth, tabla_bbp_ordered$Depth),
                           MARGIN=1,FUN=mean, na.rm=TRUE)
      depth[is.nan(depth)]<-NA
      #cbind(tabla_rss_ordered$Depth,   tabla_ap_ordered$Depth,
      #      tabla_aph_ordered$Depth,   tabla_anap_ordered$Depth,
      #      tabla_acdom_ordered$Depth, tabla_bbp_ordered$Depth, depth)      
      
      #sum(!is.na(depth))
      #sort(unique(depth))
      #length(depth)

      #plot(depth,tabla_rss_ordered[,22],xlim=c(0,5))
      
      #colnames(tabla_rss_ordered)
      metadata <- data.frame(lista, fecha3, latitude, longitude, depth)  
      tabla_total <- data.frame(metadata,
                                tabla_rss_ordered[,19:ncol(tabla_rss_ordered)],
                                tabla_ap_ordered[,19:ncol(tabla_ap_ordered)],
                                tabla_aph_ordered[,19:ncol(tabla_aph_ordered)],
                                tabla_anap_ordered[,19:ncol(tabla_anap_ordered)],
                                tabla_acdom_ordered[,19:ncol(tabla_acdom_ordered)],
                                tabla_bbp_ordered[,19:ncol(tabla_bbp_ordered)])
      #dim(tabla_total)
      #colnames(tabla_total)
      #head(tabla_total)
      tabla_total<-tabla_total[!is.nan(tabla_total$longitude),]
      #names(tabla_total)
      #tabla_total$depth
#####################      
save(tabla_total, file=paste(s_dir,"Casey_RssIop_modbands12.RData", sep=""))    
      


   
      
###############
###  CSIRO  ###
###############      
p_dir<-paste(global_path,"Dat_observations/optics/csvs/",sep="")      
      
#### AP
      reflectancia<-read.csv(paste(p_dir,"TOTAL_CSIRO_ap.csv", sep=""))   
      #names(reflectancia)
      #lon<-reflectancia$Longitude
      #lat<-reflectancia$Latitude
      #unique(reflectancia$CRUISE)
      #unique(reflectancia$Provider)
      nombres_lambda <- names(reflectancia)[c(17:37)]
      numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
      reflectancia[c(17:37)][reflectancia[c(17:37)]==(-999)]<-NA
          reflectancia[c(17:37)][reflectancia[c(17:37)]<=0.0001]<-NA
          reflectancia[c(17:37)][reflectancia[c(17:37)]>10]<-NA
      nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
      
      for (i in c(1:length(modelo_lambda))){
        cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
        if (length(cuales[!is.na(cuales)])>1) {
          res<-rowMeans(reflectancia[,16+cuales],na.rm=T)
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else if (length(cuales[!is.na(cuales)])==1) {
          res<-reflectancia[,16+cuales]
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
      }
      
     #reflectancia[reflectancia==0]<-NA
     #colnames(reflectancia) 
     colnames(nueva_reflec) <- paste("ap_",modelo_lambda, sep="")
     tabla_ap <- cbind(reflectancia[,c(1:16)],nueva_reflec) 
     #colnames(tabla_ap)

#### APH
      reflectancia<-read.csv(paste(p_dir,"TOTAL_CSIRO_aph.csv", sep=""))
      #names(reflectancia)
      nombres_lambda <- names(reflectancia)[c(17:38)]
      numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
      reflectancia[c(17:38)][reflectancia[c(17:38)]==(-999)]<-NA
          reflectancia[c(17:38)][reflectancia[c(17:38)]<=0.0001]<-NA
          reflectancia[c(17:38)][reflectancia[c(17:38)]>10]<-NA
      nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
      
      for (i in c(1:length(modelo_lambda))){
        cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
        if (length(cuales[!is.na(cuales)])>1) {
          res<-rowMeans(reflectancia[,16+cuales],na.rm=T)
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else if (length(cuales[!is.na(cuales)])==1) {
          res<-reflectancia[,16+cuales]
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
      }
      
      #reflectancia[reflectancia==0]<-NA
      colnames(nueva_reflec) <- paste("aph_",modelo_lambda, sep="")
      tabla_aph <- cbind(reflectancia[,c(1:16)],nueva_reflec) 
      #colnames(tabla_aph)
      
#### ACDOM
      reflectancia<-read.csv(paste(p_dir,"TOTAL_CSIRO_acdom.csv", sep=""))
      #names(reflectancia)
      nombres_lambda <- names(reflectancia)[c(17:38)]
      numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
      reflectancia[c(17:38)][reflectancia[c(17:38)]==(-999)]<-NA
          reflectancia[c(17:38)][reflectancia[c(17:38)]<=0.0001]<-NA
          reflectancia[c(17:38)][reflectancia[c(17:38)]>10]<-NA
      nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
      
      for (i in c(1:length(modelo_lambda))){
        cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
        if (length(cuales[!is.na(cuales)])>1) {
          res<-rowMeans(reflectancia[,16+cuales],na.rm=T)
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else if (length(cuales[!is.na(cuales)])==1) {
          res<-reflectancia[,16+cuales]
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
      }
      
      #reflectancia[reflectancia==0]<-NA
      colnames(nueva_reflec) <- paste("acdom_",modelo_lambda, sep="")
      tabla_acdom <- cbind(reflectancia[,c(1:16)],nueva_reflec) 
      #colnames(tabla_acdom)      
      
#### AD
      reflectancia<-read.csv(paste(p_dir,"TOTAL_CSIRO_ad.csv", sep=""))
      #names(reflectancia)
      nombres_lambda <- names(reflectancia)[c(17:38)]
      numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
      reflectancia[c(17:38)][reflectancia[c(17:38)]==(-999)]<-NA
          reflectancia[c(17:38)][reflectancia[c(17:38)]<=0.0001]<-NA
          reflectancia[c(17:38)][reflectancia[c(17:38)]>10]<-NA
      nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
      
      for (i in c(1:length(modelo_lambda))){
        cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
        if (length(cuales[!is.na(cuales)])>1) {
          res<-rowMeans(reflectancia[,16+cuales],na.rm=T)
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else if (length(cuales[!is.na(cuales)])==1) {
          res<-reflectancia[,16+cuales]
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
      }
      
      #reflectancia[reflectancia==0]<-NA
      colnames(nueva_reflec) <- paste("ad_",modelo_lambda, sep="")
      tabla_ad <- cbind(reflectancia[,c(1:16)],nueva_reflec) 
      #colnames(tabla_ad)  
##############      
      
#### Match all tables    
#####################      

      lista <- sort(unique(c(tabla_ap$ID,tabla_aph$ID,tabla_acdom$ID,tabla_ad$ID)))
      #length(lista)
      #sum(!is.na(match(lista,tabla_ap$ID)))
      #sum(!is.na(match(lista,tabla_aph$ID)))
      #sum(!is.na(match(lista,tabla_acdom$ID)))
      #sum(!is.na(match(lista,tabla_ad$ID)))
      
      tabla_ap_ordered     <- tabla_ap[match(lista,tabla_ap$ID),] 
      tabla_aph_ordered    <- tabla_aph[match(lista,tabla_aph$ID),] 
      tabla_acdom_ordered  <- tabla_acdom[match(lista,tabla_acdom$ID),] 
      tabla_ad_ordered     <- tabla_ad[match(lista,tabla_ad$ID),] 
      
      ## Date.Time
      #colnames(tabla_ap_ordered)
      tabla_ap_ordered_date<-paste(tabla_ap_ordered[,5],"-",sprintf("%02d", tabla_ap_ordered[,6]),"-",sprintf("%02d", tabla_ap_ordered[,7]),sep="")
      #,"T",sprintf("%02d", tabla_ap_ordered[,13]),":",sprintf("%02d", tabla_ap_ordered[,14]),sep="")  
      tabla_ap_ordered_date[is.na(tabla_ap_ordered[,5])]<-NA
      #colnames(tabla_aph_ordered)
      tabla_aph_ordered_date<-paste(tabla_aph_ordered[,5],"-",sprintf("%02d", tabla_aph_ordered[,6]),"-",sprintf("%02d", tabla_aph_ordered[,7]),sep="")
      #,"T",sprintf("%02d", tabla_aph_ordered[,13]),":",sprintf("%02d", tabla_aph_ordered[,14]),sep="")  
      tabla_aph_ordered_date[is.na(tabla_aph_ordered[,5])]<-NA
      #colnames(tabla_acdom_ordered)
      tabla_acdom_ordered_date<-paste(tabla_acdom_ordered[,5],"-",sprintf("%02d", tabla_acdom_ordered[,6]),"-",sprintf("%02d", tabla_acdom_ordered[,7]),sep="")
      #,"T",sprintf("%02d", tabla_acdom_ordered[,13]),":",sprintf("%02d", tabla_acdom_ordered[,14]),sep="")
      tabla_acdom_ordered_date[is.na(tabla_acdom_ordered[,5])]<-NA
      #colnames(tabla_bbp_ordered)
      tabla_ad_ordered_date<-paste(tabla_ad_ordered[,5],"-",sprintf("%02d", tabla_ad_ordered[,6]),"-",sprintf("%02d", tabla_ad_ordered[,7]),sep="")
      #,"T",sprintf("%02d", tabla_bbp_ordered[,13]),":",sprintf("%02d", tabla_bbp_ordered[,14]),sep="") 
      tabla_ad_ordered_date[is.na(tabla_ad_ordered[,5])]<-NA
      
      fechas <- data.frame(tabla_ap_ordered_date,     tabla_aph_ordered_date,
                           tabla_acdom_ordered_date,  tabla_ad_ordered_date)
      fecha  <- apply(fechas, MARGIN=1, function(x)unique(x[!is.na(x)]))
      fecha3<-as.character(fecha)
      fecha3[fecha3=="character(0)"]<-NA  
      #length(fecha3)
      
      latitude <- apply(cbind(tabla_ap_ordered$Latitude,    tabla_aph_ordered$Latitude,   
                              tabla_acdom_ordered$Latitude, tabla_ad_ordered$Latitude),
                              MARGIN=1, FUN=mean, na.rm=TRUE)
      #sum(!is.na(latitude))
      #length(latitude)
      
      longitude <- apply(cbind(tabla_ap_ordered$Longitude,     tabla_aph_ordered$Longitude, 
                               tabla_acdom_ordered$Longitude,  tabla_ad_ordered$Longitude),
                               MARGIN=1,FUN=mean, na.rm=TRUE)
      #sum(!is.na(longitude))
      #length(longitude)
      
      depth <- apply(cbind(tabla_ap_ordered$Depth,     tabla_aph_ordered$Depth, 
                           tabla_acdom_ordered$Depth,  tabla_ad_ordered$Depth),
                           MARGIN=1,FUN=mean, na.rm=TRUE)
      depth[is.nan(depth)]<-NA
      #sum(!is.na(depth))
      #sort(unique(depth))
      #length(depth)
      #cbind(tabla_ap_ordered$Depth,     tabla_aph_ordered$Depth, tabla_acdom_ordered$Depth,  tabla_ad_ordered$Depth, depth)      

      #colnames(tabla_aph_ordered)
      metadata <- data.frame(lista, fecha3, latitude, longitude, depth)  
      tabla_total <- data.frame(metadata,
                                tabla_ap_ordered[,17:ncol(tabla_ap_ordered)],
                                tabla_aph_ordered[,17:ncol(tabla_aph_ordered)],
                                tabla_acdom_ordered[,17:ncol(tabla_acdom_ordered)],
                                tabla_ad_ordered[,17:ncol(tabla_ad_ordered)])
      #dim(tabla_total)
      #colnames(tabla_total)
      #head(tabla_total)
      tabla_total<-tabla_total[!is.nan(tabla_total$longitude),]
      #dim(tabla_total)      
      
###############      
save(tabla_total, file=paste(s_dir,"CSIRO_Iop_modbands12.RData", sep=""))         
 
      
          
      
      
      
#####################
###  Astrid 2021  ###
#####################     
p_dir<-paste(global_path,"Dat_observations/optics/csvs/",sep="")     
      
#### AP 
      ## 2nm
      reflectancia<-read.csv(paste(p_dir,"TOTAL_Astrid4_ap_2nm.csv", sep=""))
      #names(reflectancia)
      nombres_lambda <- names(reflectancia)[c(20:295)]
      numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
      reflectancia[c(20:295)][reflectancia[c(20:295)]==(-999)]<-NA
          reflectancia[c(20:295)][reflectancia[c(20:295)]<=0.0001]<-NA
          reflectancia[c(20:295)][reflectancia[c(20:295)]>10]<-NA
      nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
      
      for (i in c(1:length(modelo_lambda))){
        cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
        if (length(cuales[!is.na(cuales)])>1) {
          res<-rowMeans(reflectancia[,19+cuales],na.rm=T)
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else if (length(cuales[!is.na(cuales)])==1) {
          res<-reflectancia[,19+cuales]
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
      }
      #reflectancia[reflectancia==0]<-NA
      colnames(nueva_reflec) <- paste("ap_",modelo_lambda, sep="")
      #names(reflectancia)[c(1:4,6,7,9:16,8,18:19)]
      tabla_ap2 <- cbind(reflectancia[,c(1:4,6,7,9:16,8,18:19)],nueva_reflec) 
      #dim(tabla_ap2) 
      
      ## 1nm
      # reflectancia<-read.csv(paste(p_dir,"TOTAL_Astrid4_ap_1nm.csv", sep="")) # corrected on March2022, 1 sample name error in MSM09
      reflectancia<-read.csv(paste(p_dir,"TOTAL_Astrid5_ap_1nm.csv", sep=""))
      #names(reflectancia)
      #unique(reflectancia$CRUISE)
      nombres_lambda <- names(reflectancia)[c(20:570)]
      numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
      reflectancia[c(20:570)][reflectancia[c(20:570)]==(-999)]<-NA
          reflectancia[c(20:570)][reflectancia[c(20:570)]<=0.0001]<-NA
          reflectancia[c(20:570)][reflectancia[c(20:570)]>10]<-NA
      nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
      
      for (i in c(1:length(modelo_lambda))){
        cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
        if (length(cuales[!is.na(cuales)])>1) {
          res<-rowMeans(reflectancia[,19+cuales],na.rm=T)
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else if (length(cuales[!is.na(cuales)])==1) {
          res<-reflectancia[,19+cuales]
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
      }
      #reflectancia[reflectancia==0]<-NA
      #colnames(tabla_ap2)
      #names(reflectancia)[c(1:4,6:16,18:19)]
      colnames(nueva_reflec) <- paste("ap_",modelo_lambda, sep="")
      tabla_ap1 <- cbind(reflectancia[,c(1:4,6:16,18:19)],nueva_reflec) 
      #dim(tabla_ap1)       
      #cbind(colnames(tabla_ap2),colnames(tabla_ap1))
      
      colnames(tabla_ap1)<-colnames(tabla_ap2)           
      tabla_ap<-rbind(tabla_ap2,tabla_ap1)
      #colnames(tabla_ap) 
      #dim(tabla_ap) 
      
      
#### APH
      # 2nm
      reflectancia<-read.csv(paste(p_dir,"TOTAL_Astrid4_aph_2nm.csv", sep=""))
      #names(reflectancia)
      nombres_lambda <- names(reflectancia)[c(20:245)]
      numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
      reflectancia[c(20:245)][reflectancia[c(20:245)]==(-999)]<-NA
          reflectancia[c(20:245)][reflectancia[c(20:245)]<=0.0001]<-NA
          reflectancia[c(20:245)][reflectancia[c(20:245)]>10]<-NA
      nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
      
      for (i in c(1:length(modelo_lambda))){
        cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
        if (length(cuales[!is.na(cuales)])>1) {
          res<-rowMeans(reflectancia[,19+cuales],na.rm=T)
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else if (length(cuales[!is.na(cuales)])==1) {
          res<-reflectancia[,19+cuales]
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
      }
      
      #reflectancia[reflectancia==0]<-NA
      colnames(nueva_reflec) <- paste("aph_",modelo_lambda, sep="")
      #names(reflectancia)[c(1:5,12,6:11,14:16,18:19)]
      tabla_aph2 <- cbind(reflectancia[,c(1:5,12,6:11,14:16,18:19)],nueva_reflec) 
      #colnames(tabla_aph2) 
      #dim(tabla_aph2) 
      
      # 1nm
      reflectancia<-read.csv(paste(p_dir,"TOTAL_Astrid4_aph_1nm.csv", sep=""))
      #names(reflectancia)
      #nique(reflectancia$CRUISE)
      nombres_lambda <- names(reflectancia)[c(18:568)]
      numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
      reflectancia[c(18:568)][reflectancia[c(18:568)]==(-999)]<-NA
          reflectancia[c(18:568)][reflectancia[c(18:568)]<=0.0001]<-NA
          reflectancia[c(18:568)][reflectancia[c(18:568)]>10]<-NA
      nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
      
      for (i in c(1:length(modelo_lambda))){
        cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
        if (length(cuales[!is.na(cuales)])>1) {
          res<-rowMeans(reflectancia[,17+cuales],na.rm=T)
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else if (length(cuales[!is.na(cuales)])==1) {
          res<-reflectancia[,17+cuales]
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
      }
      
      #reflectancia[reflectancia==0]<-NA
      colnames(nueva_reflec) <- paste("aph_",modelo_lambda, sep="")
      #names(reflectancia)[c(1:17)]
      tabla_aph1 <- cbind(reflectancia[,c(1:17)],nueva_reflec) 
      #colnames(tabla_aph1) 
      #dim(tabla_aph1) 
      #cbind(colnames(tabla_aph2),colnames(tabla_aph1))      
      
      #colnames(tabla_aph1)<-colnames(tabla_aph2)           
      tabla_aph<-rbind(tabla_aph2,tabla_aph1)      
      #colnames(tabla_aph) 
      #dim(tabla_aph) 
      
      
#### ANAP
      # 2nm
      reflectancia<-read.csv(paste(p_dir,"TOTAL_Astrid4_anap_2nm.csv", sep=""))
      #names(reflectancia)
      nombres_lambda <- names(reflectancia)[c(20:245)]
      numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
      reflectancia[c(20:245)][reflectancia[c(20:245)]==(-999)]<-NA
          reflectancia[c(20:245)][reflectancia[c(20:245)]<=0.0001]<-NA
          reflectancia[c(20:245)][reflectancia[c(20:245)]>10]<-NA
      nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
      
      for (i in c(1:length(modelo_lambda))){
        cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
        if (length(cuales[!is.na(cuales)])>1) {
          res<-rowMeans(reflectancia[,19+cuales],na.rm=T)
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else if (length(cuales[!is.na(cuales)])==1) {
          res<-reflectancia[,19+cuales]
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
      }
      
      #reflectancia[reflectancia==0]<-NA
      colnames(nueva_reflec) <- paste("anap_",modelo_lambda, sep="")
      #names(reflectancia)[c(1:5,12,6:11,14:16,18,19)]
      tabla_anap2 <- cbind(reflectancia[,c(1:5,12,6:11,14:16,18,19)],nueva_reflec) 
      #colnames(tabla_anap2)            
      
      # 1nm
      reflectancia<-read.csv(paste(p_dir,"TOTAL_Astrid4_anap_1nm.csv", sep=""))
      #names(reflectancia)
      nombres_lambda <- names(reflectancia)[c(17:417)]
      numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
      reflectancia[c(17:417)][reflectancia[c(17:417)]==(-999)]<-NA
          reflectancia[c(17:417)][reflectancia[c(17:417)]<=0.0001]<-NA
          reflectancia[c(17:417)][reflectancia[c(17:417)]>10]<-NA
      nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
      
      for (i in c(1:length(modelo_lambda))){
        cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
        if (length(cuales[!is.na(cuales)])>1) {
          res<-rowMeans(reflectancia[,16+cuales],na.rm=T)
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else if (length(cuales[!is.na(cuales)])==1) {
          res<-reflectancia[,16+cuales]
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
      }
      
      #reflectancia[reflectancia==0]<-NA
      colnames(nueva_reflec) <- paste("anap_",modelo_lambda, sep="")
      #names(reflectancia)[c(1:5,5,7:14,6,15:16)]
      tabla_anap1 <- cbind(reflectancia[,c(1:5,5,7:14,6,15:16)],nueva_reflec) 
      #colnames(tabla_anap2)
      #colnames(tabla_anap1)
      #cbind(colnames(tabla_anap2),colnames(tabla_anap1))      
      colnames(tabla_anap1)<-colnames(tabla_anap2)           
      tabla_anap<-rbind(tabla_anap2,tabla_anap1)      
      #colnames(tabla_anap)      
      
      
## ACDOM      
      # 2nm
      reflectancia<-read.csv(paste(p_dir,"TOTAL_Astrid4_acdom_2nm.csv", sep=""))
      #names(reflectancia)
      nombres_lambda <- names(reflectancia)[c(21:271)]
      numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
      reflectancia[c(21:271)][reflectancia[c(21:271)]==(-999)]<-NA
          reflectancia[c(21:271)][reflectancia[c(21:271)]<=0.0001]<-NA
          reflectancia[c(21:271)][reflectancia[c(21:271)]>10]<-NA
      nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
      
      for (i in c(1:length(modelo_lambda))){
        cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
        if (length(cuales[!is.na(cuales)])>1) {
          res<-rowMeans(reflectancia[,20+cuales],na.rm=T)
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else if (length(cuales[!is.na(cuales)])==1) {
          res<-reflectancia[,20+cuales]
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
      }
      
      #reflectancia[reflectancia==0]<-NA
      colnames(nueva_reflec) <- paste("acdom_",modelo_lambda, sep="")
      #names(reflectancia)[c(1:4,6:7,10:17,9,19,20)]
      tabla_acdom <- cbind(reflectancia[,c(1:4,6:7,10:17,9,19,20)],nueva_reflec)
      #cbind(colnames(tabla_anap2),colnames(tabla_acdom))      
      
    
## RRS
      # 1nm
      reflectancia<-read.csv(paste(p_dir,"TOTAL_Astrid4_rrs_1nm.csv", sep=""))
      #names(reflectancia)
      #reflectancia[,6]
      nombres_lambda <- names(reflectancia)[c(18:418)]
      numeros_lambda <- as.numeric(substring(nombres_lambda,2,4))          
      reflectancia[c(18:418)][reflectancia[c(18:418)]==(-999)]<-NA
        reflectancia[c(18:418)][reflectancia[c(18:418)]<=0]<-NA
        reflectancia[c(18:418)][reflectancia[c(18:418)]>0.15]<-NA
      nueva_reflec<-matrix(NA,ncol=length(modelo_lambda),nrow=nrow(reflectancia))
      
      for (i in c(1:length(modelo_lambda))){
        cuales<-which(numeros_lambda>=lambda_extremos[i] & numeros_lambda<lambda_extremos[i+1])
        if (length(cuales[!is.na(cuales)])>1) {
          res<-rowMeans(reflectancia[,17+cuales],na.rm=T)
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else if (length(cuales[!is.na(cuales)])==1) {
          res<-reflectancia[,17+cuales]
          res[is.nan(res)]<-NA
          res[res==0]<-NA
          nueva_reflec[,i]<-res
        } else {nueva_reflec[,i]<-rep(NA,length=nrow(reflectancia))}
      }
      
      #reflectancia[reflectancia==0]<-NA
      colnames(nueva_reflec) <- paste("rrs_",modelo_lambda, sep="")
      #colnames(tabla_anap2)
      #names(reflectancia)[c(1:17)]
      tabla_rrs <- cbind(reflectancia[,c(1:17)],nueva_reflec)
      #cbind(colnames(tabla_anap2),colnames(tabla_rrs))            
      
      #cbind(colnames(tabla_ap),colnames(tabla_aph),colnames(tabla_anap),colnames(tabla_acdom), colnames(tabla_rrs))
      
## HPLC
      #clorofila<-read.csv(paste(p_dir,"TOTAL_Astrid5_hplc.csv", sep="")) # corrected on March2022, 1 sample matched incorrectly to MSM09_ap_1nm
      clorofila<-read.csv(paste(p_dir,"TOTAL_Astrid6_hplc.csv", sep="")) 
      #names(clorofila)
      #dim(clorofila)
      clorofila[c(24)][clorofila[c(24)]<=0.001]<-NA
      clorofila[c(24)][clorofila[c(24)]>100]<-NA  
      
      
####################      
      
#### Match all tables    
#####################      

      lista <- sort(unique(c(tabla_ap$ID,tabla_aph$ID,tabla_anap$ID,clorofila$ID)))
      #length(lista)
      #sum(!is.na(match(lista,tabla_ap$ID)))
      #sum(!is.na(match(lista,tabla_aph$ID)))
      #sum(!is.na(match(lista,tabla_anap$ID)))
      #sum(!is.na(match(lista,tabla_acdom$ID)))
      #sum(!is.na(match(lista,tabla_rrs$ID)))
      #sum(!is.na(match(lista,clorofila$ID)))
      
      tabla_ap_ordered    <- tabla_ap[match(lista,tabla_ap$ID),] 
      tabla_aph_ordered   <- tabla_aph[match(lista,tabla_aph$ID),] 
      tabla_anap_ordered  <- tabla_anap[match(lista,tabla_anap$ID),] 
      tabla_acdom_ordered <- tabla_acdom[match(lista,tabla_acdom$ID),] 
      tabla_rrs_ordered   <- tabla_rrs[match(lista,tabla_rrs$ID),] 
      cloro_ordered       <- clorofila[match(lista,clorofila$ID),] 

      ## Date.Time
      #colnames(tabla_ap_ordered)
      tabla_ap_ordered_date<-paste(tabla_ap_ordered[,8],"-",sprintf("%02d", tabla_ap_ordered[,9]),"-",sprintf("%02d", tabla_ap_ordered[,10]),sep="")
      #,"T",sprintf("%02d", tabla_ap_ordered[,13]),":",sprintf("%02d", tabla_ap_ordered[,14]),sep="")  
      tabla_ap_ordered_date[is.na(tabla_ap_ordered[,8])]<-NA
      #colnames(tabla_aph_ordered)
      tabla_aph_ordered_date<-paste(tabla_aph_ordered[,8],"-",sprintf("%02d", tabla_aph_ordered[,9]),"-",sprintf("%02d", tabla_aph_ordered[,10]),sep="")
      #,"T",sprintf("%02d", tabla_aph_ordered[,13]),":",sprintf("%02d", tabla_aph_ordered[,14]),sep="")  
      tabla_aph_ordered_date[is.na(tabla_aph_ordered[,8])]<-NA
      #colnames(tabla_anap_ordered)
      tabla_anap_ordered_date<-paste(tabla_anap_ordered[,8],"-",sprintf("%02d", tabla_anap_ordered[,9]),"-",sprintf("%02d", tabla_anap_ordered[,10]),sep="")
      #,"T",sprintf("%02d", tabla_acdom_ordered[,13]),":",sprintf("%02d", tabla_acdom_ordered[,14]),sep="")
      tabla_anap_ordered_date[is.na(tabla_anap_ordered[,8])]<-NA
      
      #colnames(tabla_acdom_ordered)
      tabla_acdom_ordered_date<-paste(tabla_acdom_ordered[,8],"-",sprintf("%02d", tabla_acdom_ordered[,9]),"-",sprintf("%02d", tabla_acdom_ordered[,10]),sep="")
      #,"T",sprintf("%02d", tabla_acdom_ordered[,13]),":",sprintf("%02d", tabla_acdom_ordered[,14]),sep="")
      tabla_acdom_ordered_date[is.na(tabla_acdom_ordered[,8])]<-NA      
      #colnames(tabla_rrs_ordered)
      tabla_rrs_ordered_date<-paste(tabla_rrs_ordered[,8],"-",sprintf("%02d", tabla_rrs_ordered[,9]),"-",sprintf("%02d", tabla_rrs_ordered[,10]),sep="")
      #,"T",sprintf("%02d", tabla_acdom_ordered[,13]),":",sprintf("%02d", tabla_acdom_ordered[,14]),sep="")
      tabla_rrs_ordered_date[is.na(tabla_rrs_ordered[,8])]<-NA    
            
      #colnames(cloro_ordered)
      cloro_ordered_date<-paste(cloro_ordered[,11],"-",sprintf("%02d", cloro_ordered[,12]),"-",sprintf("%02d", cloro_ordered[,13]),sep="")
      #,"T",sprintf("%02d", tabla_bbp_ordered[,13]),":",sprintf("%02d", tabla_bbp_ordered[,14]),sep="") 
      cloro_ordered_date[is.na(cloro_ordered[,11])]<-NA
      
      fechas <- data.frame(tabla_ap_ordered_date,     tabla_aph_ordered_date,
                           tabla_anap_ordered_date,
                           tabla_acdom_ordered_date,  tabla_rrs_ordered_date,
                           cloro_ordered_date)
      #fecha  <- apply(fechas, MARGIN=1, function(x)unique(x[!is.na(x)]))
      fecha  <- apply(fechas, MARGIN=1, function(x)min(x[!is.na(x)]))
      fecha3<-as.character(fecha)
      fecha3[fecha3=="character(0)"]<-NA  
      #length(fecha3)
      #colnames(tabla_ap_ordered)
      #data.frame(tabla_ap_ordered$ID, fechas)
      #write.csv(resul,file=paste(s_dir,"check_dates.csv", sep=""))
      
      latitude <- apply(cbind(tabla_ap_ordered$Latitude,    tabla_aph_ordered$Latitude,   
                              tabla_anap_ordered$Latitude,
                              tabla_acdom_ordered$Latitude, tabla_rrs_ordered$Latitude,
                              cloro_ordered$Latitude),  MARGIN=1, FUN=mean, na.rm=TRUE)
      #sum(!is.na(latitude))
      #length(latitude)
      
      longitude <- apply(cbind(tabla_ap_ordered$Longitude,     tabla_aph_ordered$Longitude, 
                               tabla_anap_ordered$Longitude,
                               tabla_acdom_ordered$Longitude,  tabla_rrs_ordered$Longitude,                               
                               cloro_ordered$Longitude),     MARGIN=1,FUN=mean, na.rm=TRUE)
      #sum(!is.na(longitude))
      #length(longitude)
      
      depth <- apply(cbind(tabla_ap_ordered$Depth,     tabla_aph_ordered$Depth, 
                           tabla_anap_ordered$Depth,
                           tabla_acdom_ordered$Depth,  tabla_rrs_ordered$Depth,
                           cloro_ordered$Depth),  MARGIN=1,FUN=mean, na.rm=TRUE)
      depth[is.nan(depth)]<-NA
      #cbind(tabla_ap_ordered$Depth,     tabla_aph_ordered$Depth, 
      #      tabla_anap_ordered$Depth,
      #      tabla_acdom_ordered$Depth,  tabla_rrs_ordered$Depth,  cloro_ordered$Depth, depth)
      #sum(!is.na(depth))
      #sort(unique(depth))
      #length(depth)
      
      #colnames(tabla_anap_ordered)
      #colnames(cloro_ordered)
      metadata <- data.frame(lista, fecha3, latitude, longitude, depth)  
      tabla_total <- data.frame(metadata,
                                tabla_ap_ordered[,18:ncol(tabla_ap_ordered)],
                                tabla_aph_ordered[,18:ncol(tabla_aph_ordered)],
                                tabla_anap_ordered[,18:ncol(tabla_anap_ordered)],
                                tabla_acdom_ordered[,18:ncol(tabla_acdom_ordered)],
                                tabla_rrs_ordered[,18:ncol(tabla_rrs_ordered)],
                                cloro_ordered[,20:ncol(cloro_ordered)])
      
      #colnames(cloro_ordered)[20:ncol(cloro_ordered)]     
      #dim(tabla_total)
      #colnames(tabla_total)
      #head(tabla_total)
      tabla_total<-tabla_total[!is.nan(tabla_total$longitude),]
      #dim(tabla_total)
      
#####################      
#save(tabla_total, file=paste(s_dir,"Astrid5_IopHplc2nm_modbands12.RData", sep=""))       
      
save(tabla_total, file=paste(s_dir,"Astrid6_IopHplc2nm_modbands12.RData", sep=""))       
      
      
      
      
## Make grids:  takes time!
source(paste(global_path, "Scripts/05_Insitu_BO_Valente_Casey_Chase_CSIRO_Astrid_make_grids.R", sep=""))
      