source("00_Summary_EDITME.R")

##########################################
## SATELLITE ##  SEASONAL AND ANNUAL MEANS
##########################################
## ftp://oc-cci-data:ELaiWai8ae@ftp.rsg.pml.ac.uk/occci-v4.1/geographic/netcdf/monthly/





########################
####  SURFACE CHLA  #### 
########################

# Original archives from OCCCI: 
o_dir<-paste(global_path,"/Dat_observations/satellite/",sep="")
dir1 = "chla_MO_OCCCI/"

###################   
      archivo<-list.files(paste(o_dir,dir1, sep=""),pattern = ".nc")
      year1<-as.numeric(substring(archivo,57,60))  
      anito<-unique(year1)
      
      for (k in 1:length(anito)){   #length(anito)
        filename <- paste(o_dir,dir1,archivo[year1==anito[k]], sep="")
        #nchar(filename)
        year<-as.numeric(substring(filename,123,126))    
        ano <- anito[k]
        month<-as.numeric(substring(filename,127,128))
        day<-rep(15, length(filename))
        JD <- julian(as.Date(paste(year,month,day, sep="-")), origin = as.Date(paste(ano,"-01-01", sep=""))) #"%Y-%m-%d"       
        
        for (h in 1:length(filename)){      
          filenc <- open.nc(filename[h])
          tmp <- read.nc(filenc)
          z <- tmp$chlor_a
          lon <- tmp$lon
          lat <- tmp$lat
          colnames(z)=lat
          rownames(z)=lon
          z[which(is.nan(z))] <- NA 
          z[which(is.na(z))] <- NA
          rm(tmp)
          close.nc(filenc) 
          
          # Reducir resolucion del satelite para casar con recom
          z <- flipud(z)
          #dim(z)
          red=48    #8460/180  #4320/90   
          stepsi<-floor(nrow(z)/red)
          stepsj<-floor(ncol(z)/red)
          LAT<-rep(NA, length=stepsj)
          LON<-rep(NA, length=stepsi)
          ZETA<-matrix(NA, nrow=stepsi, ncol=stepsj)
          for (i in 1:stepsi){
            for (j in 1:stepsj){
              valores <- z[c((((i-1)*red)+1):(i*red)),c((((j-1)*red)+1):(j*red))]
              if (length(valores[!is.nan(valores)])<=(length(valores)/3)) {valores <- NA}
              ZETA[i,j]<-mean(valores, na.rm=TRUE)
              LAT[j]<-mean(lat[c((((j-1)*red)+1):(j*red))], na.rm=TRUE)} # end loop j
            LON[i]<-mean(lon[c((((i-1)*red)+1):(i*red))], na.rm=TRUE)}   # end loop i
          z<-ZETA 
          #dim(z)
          #image(z,col = jet(10), las=1)  #, zlim=c(-1.8,1.8)
          
          if (h==1){  m<-array(data = NA, dim=c(nrow(z),ncol(z),length(filename)), dimnames = NULL)
          dimnames(m) <- list(x=LON, y=LAT, z=seq(1,length(filename),by=1))
          m[,,h]<-z } else { m[,,h]<-z }  }  # end loop h
        res<-m
        dimnames(res) <- list(x=LON, y=LAT, z=JD)
        save(res,LAT,LON, file =paste(o_dir,"res_seasonal_OCCCI_chla/media_anual_chl_OCCCI_",anito[k],".RData",sep="") )
      } # end loop k 
###################


## CLIMATOLOGY LINEAL
####################   
   ## SEASONAL    
      o_dir <- paste(OS,"Datos/Dat_Satelite/",sep="")
      anito<-seq(2012,2018)
      load(paste(o_dir,"res_seasonal_OCCCI_chla/media_anual_chl_OCCCI_",anito[1],".RData",sep="") )
      m1<-res
      load(paste(o_dir,"res_seasonal_OCCCI_chla/media_anual_chl_OCCCI_",anito[2],".RData",sep="") )
      m2<-res
      load(paste(o_dir,"res_seasonal_OCCCI_chla/media_anual_chl_OCCCI_",anito[3],".RData",sep="") )
      m3<-res
      load(paste(o_dir,"res_seasonal_OCCCI_chla/media_anual_chl_OCCCI_",anito[4],".RData",sep="") )
      m4<-res
      load(paste(o_dir,"res_seasonal_OCCCI_chla/media_anual_chl_OCCCI_",anito[5],".RData",sep="") )
      m5<-res 
      load(paste(o_dir,"res_seasonal_OCCCI_chla/media_anual_chl_OCCCI_",anito[6],".RData",sep="") )
      m6<-res       
              JD <- as.numeric(unlist(dimnames(m1)[3]))    
              M <- array(NA, dim = c(dim(m1),length(anito)), dimnames=list(x=unlist(dimnames(m1)[1]), y=unlist(dimnames(m1)[2]), z=JD,
                                                                           w=anito))
              dim(M)
              M[,,,1] <- m1
              M[,,,2] <- m2
              M[,,,3] <- m3
              M[,,,4] <- m4
              M[,,,5] <- m5
              M[,,,6] <- m6
        dim(M)
        med <- apply(M , 1:3 , mean, na.rm=TRUE)
        rm(M)
        dim(med)
        dimnames(med) <- list(x=as.numeric(unlist(dimnames(m1)[1])),
                              y=as.numeric(unlist(dimnames(m1)[2])), z=JD)
        med[is.nan(med)] <- NA
        save(med, file =paste(o_dir,"climatologies_2012_2018/media_seasonal_chl_2deg_OCCCI_2012_2018.RData",sep="") )
      
   ## ANNUAL ##  
        source(paste(OS,"Programas/Funciones_varios/functions.r",sep=""))
        o_dir <- "/Users/ealvarez/Datos/Dat_Satelite/"
                   load(paste(o_dir,"climatologies_2012_2018/media_seasonal_chl_2deg_OCCCI_2012_2018.RData",sep="") )
        JD <- as.numeric(unlist(dimnames(med)[3])) 
        res <- apply(med, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)
        #res <- flipud(res)
        save(res, file =paste(o_dir,"climatologies_2012_2018/media_anual_chl_2deg_OCCCI_2012_2018.RData",sep="") )
####################        

## CLIMATOLOGY LOG
##################     
   ## SEASONAL    
   o_dir <- paste(OS,"Datos/Dat_Satelite/",sep="")
   anito<-c(2012:2018)
        load(paste(o_dir,"res_seasonal_OCCCI_chla/media_log_chl_OCCCI_",anito[1],".RData",sep="") )
        m1<-res
        load(paste(o_dir,"res_seasonal_OCCCI_chla/media_log_chl_OCCCI_",anito[2],".RData",sep="") )
        m2<-res
        load(paste(o_dir,"res_seasonal_OCCCI_chla/media_log_chl_OCCCI_",anito[3],".RData",sep="") )
        m3<-res
        load(paste(o_dir,"res_seasonal_OCCCI_chla/media_log_chl_OCCCI_",anito[4],".RData",sep="") )
        m4<-res
        load(paste(o_dir,"res_seasonal_OCCCI_chla/media_log_chl_OCCCI_",anito[5],".RData",sep="") )
        m5<-res 
        JD <- as.numeric(unlist(dimnames(m1)[3]))    
        M <- array(NA, dim = c(dim(m1),length(anito)), dimnames=list(x=unlist(dimnames(m1)[1]), y=unlist(dimnames(m1)[2]), z=JD, w=anito))
        M[,,,1] <- m1
        M[,,,2] <- m2
        M[,,,3] <- m3
        M[,,,4] <- m4
        M[,,,5] <- m5
        dim(M)
        med <- apply(M , 1:3 , mean, na.rm=TRUE)
        rm(M)
        dim(med)
        dimnames(med) <- list(x=as.numeric(unlist(dimnames(m1)[1])),
                              y=as.numeric(unlist(dimnames(m1)[2])), z=JD)
        med[is.nan(med)] <- NA
    save(med, file =paste(o_dir,"climatologies_2012_2018/media_seasonal_log_chl_OCCCI_2012_2018.RData",sep="") )
        
    ## ANNUAL ##  
    #source("/Users/ealvarez/Programas/Funciones_varios/functions.r")
    #o_dir <- "/Users/ealvarez/Datos/Dat_Satelite/"
    #load(paste(o_dir,"climatologies_2012_2018/media_seasonal_log_chl_OCCCI_2012_2018.RData",sep="") )
        JD <- as.numeric(unlist(dimnames(med)[3])) 
        res <- apply(med, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)
        #res <- flipud(res)
        save(res, file =paste(o_dir,"climatologies_2012_2018/media_anual_log_chl_OCCCI_2012_2018.RData",sep="") )
##################        




###############################
####  SURFACE REFLECTANCE  #### 
###############################

# Original archives from OCCCI: 
o_dir<-paste(global_path,"/Dat_observations/satellite/",sep="")
dir1 = "reflectance_MO_OCCCI/"

###################
    archivo<-list.files(paste(o_dir,dir1, sep=""),pattern = ".nc")
    year1<-as.numeric(substring(archivo,53,56))  
    anito<-unique(year1)
    
    custom_scale <- c("deepskyblue4","deepskyblue3","darkslategray1",
                      "mediumaquamarine","greenyellow",
                      "yellow", "gold1",
                      "orange", "tomato","orangered2","red3")
    landas_sat <- c("Rrs_412","Rrs_443","Rrs_490","Rrs_510","Rrs_555","Rrs_670")
    
    
    k=2
    
    
    
    for (k in c(1:2)){   # ano
    ano <- anito[k]
    
    for (w in 1:length(landas_sat)){   
        filename <- paste(o_dir,dir1,archivo[year1==anito[k]], sep="")
        #nchar(filename)
        year<-as.numeric(substring(filename,nchar(filename)-14,nchar(filename)-11))
        ano <- anito[k]
        month<-as.numeric(substring(filename,nchar(filename)-10,nchar(filename)-9))
        day<-rep(15, length(filename))
        JD <- julian(as.Date(paste(year,month,day, sep="-")), origin = as.Date(paste(ano,"-01-01", sep=""))) #"%Y-%m-%d"       
    
    #par(mfrow=c(6,2)) 
    #par(mar=c(2,2,1,1))
    #layout.show(n=12)
     
        for (h in 1:length(filename)){      
            filenc <- open.nc(filename[h])
            tmp <- read.nc(filenc)
            #names(tmp)
            z <- tmp[[which(names(tmp)==landas_sat[w])]]
            lon <- tmp$lon
            lat <- tmp$lat
            
            colnames(z)=lat
            rownames(z)=lon
            z[which(is.nan(z))] <- NA 
            z[which(is.na(z))] <- NA
            rm(tmp)
            close.nc(filenc) 
            
            # Reducir resolucion del satelite para casar con recom
            z <- flipud(z)
            #dim(z)
            red=48     # 48    #8460/180  #4320/90      #8460/360  #4320/180
            stepsi<-floor(nrow(z)/red)
            stepsj<-floor(ncol(z)/red)
            LAT<-rep(NA, length=stepsj)
            LON<-rep(NA, length=stepsi)
            ZETA<-matrix(NA, nrow=stepsi, ncol=stepsj)
            #dim(ZETA)
            for (i in 1:stepsi){
                for (j in 1:stepsj){
                    valores <- z[c((((i-1)*red)+1):(i*red)),c((((j-1)*red)+1):(j*red))]
                    if (length(valores[!is.nan(valores)])<=(length(valores)/3)) {valores <- NA}
                    ZETA[i,j]<-mean(valores, na.rm=TRUE)
                    LAT[j]<-mean(lat[c((((j-1)*red)+1):(j*red))], na.rm=TRUE)} # end loop j
                LON[i]<-mean(lon[c((((i-1)*red)+1):(i*red))], na.rm=TRUE)}   # end loop i
            z<-ZETA 
            
            #dim(z)
            #range(z, na.rm=TRUE)
            #image(z, col=custom_scale, zlim=c(0,0.05), las=1)  #, zlim=c(-1.8,1.8)
            
            if (h==1){  m<-array(data = NA, dim=c(nrow(z),ncol(z),length(filename)), dimnames = NULL)
            dimnames(m) <- list(x=LON, y=LAT, z=seq(1,length(filename),by=1))
            m[,,h]<-z } else { m[,,h]<-z }  }  # end loop h (month)
        res<-m
        #dim(res)
        dimnames(res) <- list(x=LON, y=-LAT, z=JD)
        save(res,LAT,LON, file=paste(o_dir,"res_seasonal_OCCCI_reflectance/media_anual_",
                     landas_sat[w],"_2deg_OCCCI_",ano,".RData",sep="") )
        } # end loop w (landa)
    
    } # end loop year
#####################    



## Unir landas en un solo archivo 
#################################    
ano=2013
landas_sat <- c("Rrs_412","Rrs_443","Rrs_490","Rrs_510","Rrs_555","Rrs_670")
year<-rep(ano, 12)
month<-c(1:12)
day<-rep(15, 12)
JD <- julian(as.Date(paste(year,month,day, sep="-")), origin = as.Date(paste(ano,"-01-01", sep=""))) #"%Y-%m-%d"       

for (w in 1:length(landas_sat)){
    load(paste(o_dir,"res_seasonal_OCCCI_reflectance/media_anual_",landas_sat[w],"_2deg_OCCCI_",ano,".RData",sep="") )
    
    if (w==1){  md<-array(data = NA, dim=c(dim(res),length(landas_sat)),
                          dimnames=list(x=LON, y=LAT, z=seq(1,12,by=1),l=landas_sat))
    md[,,,w]<-res } else { md[,,,w]<-res }       } # end loop w             

dim(md)
res <- md
save(res,LAT,LON,JD,landas_sat, file=paste(o_dir,"res_seasonal_OCCCI_reflectance/media_anual_Rrs_2deg_OCCCI_",ano,".RData",sep="") )
#################################



## CLIMATOLOGY
##############    
    ## SEASONAL    
    o_dir <- paste(OS,"Datos/Dat_Satelite/",sep="")
    anito<-c(2012:2018)
    landas_sat <- c("Rrs_412","Rrs_443","Rrs_490","Rrs_510","Rrs_555","Rrs_670")
    load(paste(o_dir, "res_seasonal_OCCCI_reflectance/media_anual_Rrs_2deg_OCCCI_",anito[1],".RData",sep="") )
    m1<-res
    dim(m1)
    load(paste(o_dir, "res_seasonal_OCCCI_reflectance/media_anual_Rrs_2deg_OCCCI_",anito[2],".RData",sep="") )
    m2<-res
    dim(m2)
    load(paste(o_dir, "res_seasonal_OCCCI_reflectance/media_anual_Rrs_2deg_OCCCI_",anito[3],".RData",sep="") )
    m3<-res
    dim(m3)
    load(paste(o_dir, "res_seasonal_OCCCI_reflectance/media_anual_Rrs_2deg_OCCCI_",anito[4],".RData",sep="") )
    m4<-res
    dim(m4)
    load(paste(o_dir, "res_seasonal_OCCCI_reflectance/media_anual_Rrs_2deg_OCCCI_",anito[5],".RData",sep="") )
    m5<-res
    dim(m5)
    load(paste(o_dir, "res_seasonal_OCCCI_reflectance/media_anual_Rrs_2deg_OCCCI_",anito[6],".RData",sep="") )
    m6<-res
    dim(m6)    
    
            JD <- as.numeric(unlist(dimnames(m1)[3]))    
            M <- array(NA, dim = c(dim(m1),length(anito)), dimnames=list(x=unlist(dimnames(m1)[1]), y=unlist(dimnames(m1)[2]), t=JD,
                                                                         l=landas_sat, w=anito))
            M[,,,,1] <- m1
            M[,,,,2] <- m2
            M[,,,,3] <- m3
            M[,,,,4] <- m4
            M[,,,,5] <- m5
            M[,,,,6] <- m6  
            
            dim(M)
            med <- apply(M , 1:4 , mean, na.rm=TRUE)
            rm(M)
            dim(med)
            dimnames(med) <- list(x=as.numeric(unlist(dimnames(m1)[1])),
                                  y=as.numeric(unlist(dimnames(m1)[2])),
                                  t=JD, l=landas_sat)
            med[is.nan(med)] <- NA
    save(med, file =paste(o_dir,"climatologies_2012_2018/media_seasonal_Rrs_2deg_OCCCI_2012_2018.RData",sep="") )
    

    ## ANNUAL ##  
    #source("/Users/ealvarez/Programas/Funciones_varios/functions.r")
    #o_dir <- "/Users/ealvarez/Datos/Dat_Satelite/"
    #load(paste(o_dir,"climatologies_2008_2012/media_seasonal_Rrs_2deg_OCCCI_2008_2012.RData",sep="") )
    #dim(med)
            JD <- as.numeric(unlist(dimnames(med))) 
            res <- apply(med, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)
            dim(res)
        save(res, file =paste(o_dir,"climatologies_2012_2018/media_anual_Rrs_2deg_OCCCI_2012_2018.RData",sep="") )
##############
    


            
            
            
                        
            
###########################
####   SURFACE IOP's   #### 
###########################
            aph_sat <- c("aph_412","aph_443","aph_490","aph_510","aph_555","aph_670")
            adg_sat <- c("adg_412","adg_443","adg_490","adg_510","adg_555","adg_670")
            atot_sat <- c("atot_412","atot_443","atot_490","atot_510","atot_555","atot_670")
            bbp_sat <- c("bbp_412","bbp_443","bbp_490","bbp_510","bbp_555","bbp_670")            
        
     # Original archives from OCCCI:                           
     o_dir<-paste(global_path,"/Dat_observations/satellite/",sep="")
     dir1 = "iops_MO_OCCCI/"

###########################        
       archivo<-list.files(paste(o_dir,dir1, sep=""),pattern = ".nc")
       year1<-as.numeric(substring(archivo,nchar(archivo)-14,nchar(archivo)-11))  
       anito<-unique(year1)
       
       k=2
       ano <- anito[k]       
       
       filename <- paste(o_dir,dir1,archivo[year1==ano], sep="")
       year<-as.numeric(substring(filename,nchar(filename)-14,nchar(filename)-11))
       month<-as.numeric(substring(filename,nchar(filename)-10,nchar(filename)-9))
       day<-rep(15, length(filename))
       JD <- julian(as.Date(paste(year,month,day, sep="-")), origin = as.Date(paste(ano,"-01-01", sep=""))) #"%Y-%m-%d"
        


  for (h in c(1:12)){           
        filenc <- open.nc(filename[h])
        tmp <- read.nc(filenc)
        #names(tmp)
        close.nc(filenc) 
            
            for (w in 1:length(aph_sat)){ 
            z <- tmp[[which(names(tmp)==aph_sat[w])]]
            lon <- tmp$lon
            lat <- tmp$lat
                colnames(z)=lat
                rownames(z)=lon
                z[which(is.nan(z))] <- NA 
                z[which(is.na(z))] <- NA
                # Reducir resolucion del satelite para casar con recom
                z <- flipud(z)
          #dim(z)
          #range(z, na.rm=TRUE)
          test <- hist(z, breaks=seq(min(z,na.rm=TRUE), max(z,na.rm=TRUE),length=25), plot=FALSE)
          z[z>min(test$mids[which(test$counts<10)], na.rm=TRUE)] <- NA
          #range(z, na.rm=TRUE)
          #hist(z)
          #points(test$mids,test$counts, type="b", ylim=c(0,30), col="green")
          
                red=48     # 48    #8460/180  #4320/90      #8460/360  #4320/180
                stepsi<-floor(nrow(z)/red)
                stepsj<-floor(ncol(z)/red)
                LAT<-rep(NA, length=stepsj)
                LON<-rep(NA, length=stepsi)
                ZETA<-matrix(NA, nrow=stepsi, ncol=stepsj)
                for (i in 1:stepsi){
                    for (j in 1:stepsj){
                        valores <- z[c((((i-1)*red)+1):(i*red)),c((((j-1)*red)+1):(j*red))]
                        if (length(valores[!is.nan(valores)])<=(length(valores)/3)) {valores <- NA}
                        ZETA[i,j]<-mean(valores, na.rm=TRUE)
                        LAT[j]<-mean(lat[c((((j-1)*red)+1):(j*red))], na.rm=TRUE)} # end loop j
                    LON[i]<-mean(lon[c((((i-1)*red)+1):(i*red))], na.rm=TRUE)}   # end loop i
                z<-ZETA 
                save(z,LAT,LON, file=paste(o_dir,"res_seasonal_OCCCI_IOPs/", h,"_",aph_sat[w],"_2deg_OCCCI_",ano,".RData",sep="") )
            } # end loop w (landa)
        
        
            for (w in 1:length(adg_sat)){             
            z <- tmp[[which(names(tmp)==adg_sat[w])]]
            lon <- tmp$lon
            lat <- tmp$lat
                colnames(z)=lat
                rownames(z)=lon
                z[which(is.nan(z))] <- NA 
                z[which(is.na(z))] <- NA
                    # Reducir resolucion del satelite para casar con recom
                    z <- flipud(z)
                    test <- hist(z, breaks=seq(min(z,na.rm=TRUE), max(z,na.rm=TRUE),length=25), plot=FALSE)
                    z[z>min(test$mids[which(test$counts<5)], na.rm=TRUE)] <- NA
                    red=48     # 48    #8460/180  #4320/90      #8460/360  #4320/180
                    stepsi<-floor(nrow(z)/red)
                    stepsj<-floor(ncol(z)/red)
                    LAT<-rep(NA, length=stepsj)
                    LON<-rep(NA, length=stepsi)
                    ZETA<-matrix(NA, nrow=stepsi, ncol=stepsj)
                    for (i in 1:stepsi){
                        for (j in 1:stepsj){
                            valores <- z[c((((i-1)*red)+1):(i*red)),c((((j-1)*red)+1):(j*red))]
                            if (length(valores[!is.nan(valores)])<=(length(valores)/3)) {valores <- NA}
                            ZETA[i,j]<-mean(valores, na.rm=TRUE)
                            LAT[j]<-mean(lat[c((((j-1)*red)+1):(j*red))], na.rm=TRUE)} # end loop j
                        LON[i]<-mean(lon[c((((i-1)*red)+1):(i*red))], na.rm=TRUE)}   # end loop i
                    z<-ZETA 
               save(z,LAT,LON, file=paste(o_dir,"res_seasonal_OCCCI_IOPs/",h,"_",adg_sat[w],"_2deg_OCCCI_",ano,".RData",sep="") )
            } # end loop w 
            
            for (w in 1:length(atot_sat)){             
                z <- tmp[[which(names(tmp)==atot_sat[w])]]
                lon <- tmp$lon
                lat <- tmp$lat
                colnames(z)=lat
                rownames(z)=lon
                z[which(is.nan(z))] <- NA 
                z[which(is.na(z))] <- NA
                # Reducir resolucion del satelite para casar con recom
                z <- flipud(z)
                test <- hist(z, breaks=seq(min(z,na.rm=TRUE), max(z,na.rm=TRUE),length=25), plot=FALSE)
                z[z>min(test$mids[which(test$counts<5)], na.rm=TRUE)] <- NA
                red=48     # 48    #8460/180  #4320/90      #8460/360  #4320/180
                stepsi<-floor(nrow(z)/red)
                stepsj<-floor(ncol(z)/red)
                LAT<-rep(NA, length=stepsj)
                LON<-rep(NA, length=stepsi)
                ZETA<-matrix(NA, nrow=stepsi, ncol=stepsj)
                for (i in 1:stepsi){
                    for (j in 1:stepsj){
                        valores <- z[c((((i-1)*red)+1):(i*red)),c((((j-1)*red)+1):(j*red))]
                        if (length(valores[!is.nan(valores)])<=(length(valores)/3)) {valores <- NA}
                        ZETA[i,j]<-mean(valores, na.rm=TRUE)
                        LAT[j]<-mean(lat[c((((j-1)*red)+1):(j*red))], na.rm=TRUE)} # end loop j
                    LON[i]<-mean(lon[c((((i-1)*red)+1):(i*red))], na.rm=TRUE)}   # end loop i
                z<-ZETA 
                save(z,LAT,LON, file=paste(o_dir,"res_seasonal_OCCCI_IOPs/",h,"_",atot_sat[w],"_2deg_OCCCI_",ano,".RData",sep="") )
            } # end loop w 
            
            for (w in 1:length(bbp_sat)){             
                z <- tmp[[which(names(tmp)==bbp_sat[w])]]
                lon <- tmp$lon
                lat <- tmp$lat
                colnames(z)=lat
                rownames(z)=lon
                z[which(is.nan(z))] <- NA 
                z[which(is.na(z))] <- NA
                # Reducir resolucion del satelite para casar con recom
                z <- flipud(z)
                test <- hist(z, breaks=seq(min(z,na.rm=TRUE), max(z,na.rm=TRUE),length=25), plot=FALSE)
                z[z>min(test$mids[which(test$counts<5)], na.rm=TRUE)] <- NA
                red=48     # 48    #8460/180  #4320/90      #8460/360  #4320/180
                stepsi<-floor(nrow(z)/red)
                stepsj<-floor(ncol(z)/red)
                LAT<-rep(NA, length=stepsj)
                LON<-rep(NA, length=stepsi)
                ZETA<-matrix(NA, nrow=stepsi, ncol=stepsj)
                for (i in 1:stepsi){
                    for (j in 1:stepsj){
                        valores <- z[c((((i-1)*red)+1):(i*red)),c((((j-1)*red)+1):(j*red))]
                        if (length(valores[!is.nan(valores)])<=(length(valores)/3)) {valores <- NA}
                        ZETA[i,j]<-mean(valores, na.rm=TRUE)
                        LAT[j]<-mean(lat[c((((j-1)*red)+1):(j*red))], na.rm=TRUE)} # end loop j
                    LON[i]<-mean(lon[c((((i-1)*red)+1):(i*red))], na.rm=TRUE)}   # end loop i
                z<-ZETA 
                save(z,LAT,LON, file=paste(o_dir,"res_seasonal_OCCCI_IOPs/",h,"_",bbp_sat[w],"_2deg_OCCCI_",ano,".RData",sep="") )
            } # end loop w  (lambda)     
        rm(z)
        rm(tmp)
     } # end loop h (month)       
     
###########################     
     
       
       
       
       
     aph_sat <- c("aph_412", "aph_443", "aph_490", "aph_510", "aph_555", "aph_670")
     adg_sat <- c("adg_412", "adg_443", "adg_490", "adg_510", "adg_555", "adg_670")
     atot_sat<- c("atot_412","atot_443","atot_490","atot_510","atot_555","atot_670")
     bbp_sat <- c("bbp_412", "bbp_443", "bbp_490", "bbp_510", "bbp_555", "bbp_670")            
     
     ano <- 2013 # 2015
     o_dir <- "E:/Datos/Dat_Satelite/"
     #w=6
     
            for (w in 1:length(atot_sat)){  #1
                year<-rep(ano, 12)
                month<-c(1:12)
                day<-rep(15, 12)
                JD <- julian(as.Date(paste(year,month,day, sep="-")), origin = as.Date(paste(ano,"-01-01", sep=""))) #"%Y-%m-%d"       
                
                for (h in c(1:12)){
                load(paste(o_dir,"res_seasonal_OCCCI_IOPs/", h,"_",
                           atot_sat[w],"_2deg_OCCCI_",ano,".RData",sep="") ) #2
                    #dim(z)
                    #range(z, na.rm=TRUE)
                    #image(z, col=custom_scale, zlim=c(0,0.05), las=1)  #, zlim=c(-1.8,1.8)
                    if (h==1){  m<-array(data = NA, dim=c(nrow(z),ncol(z),12), dimnames = NULL)
                    dimnames(m) <- list(x=LON, y=LAT, z=seq(1,12,by=1))
                    m[,,h]<-z } else { m[,,h]<-z }  }  # end loop h
                res<-m
                #dim(res)
                dimnames(res) <- list(x=LON, y=-LAT, z=JD)
                
                if (w==1){  md<-array(data = NA, dim=c(dim(res),
                                length(atot_sat)), #3
                                dimnames=list(x=LON, y=LAT, z=seq(1,12,by=1),
                                l=atot_sat)) #4
                #dimnames(md) <- list(x=LON, y=LAT, z=seq(1,12,by=1),l=aph_sat)
                md[,,,w]<-res } else { md[,,,w]<-res }  
                } # end loop w             
        
#dim(md)
res <- md
save(res,LAT,LON,JD, atot_sat, file=paste(o_dir, #5
"res_seasonal_OCCCI_IOPs/media_anual_Atot_2deg_OCCCI_", #6
ano,".RData",sep="") )

           




## CLIMATOLOGY Aph
###################    
    ## SEASONAL    
    o_dir <- paste(OS,"Datos/Dat_Satelite/",sep="")
    anito<-c(2012:2018)
    landas_sat <- c("Rrs_412","Rrs_443","Rrs_490","Rrs_510","Rrs_555","Rrs_670")
    
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Aph_2deg_OCCCI_",anito[1],".RData",sep="") )
    m1<-res
    dim(m1)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Aph_2deg_OCCCI_",anito[2],".RData",sep="") )
    m2<-res
    dim(m2)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Aph_2deg_OCCCI_",anito[3],".RData",sep="") )
    m3<-res
    dim(m3)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Aph_2deg_OCCCI_",anito[4],".RData",sep="") )
    m4<-res
    dim(m4)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Aph_2deg_OCCCI_",anito[5],".RData",sep="") )
    m5<-res
    dim(m5)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Aph_2deg_OCCCI_",anito[6],".RData",sep="") )
    m6<-res
    dim(m6)
    
    JD <- as.numeric(unlist(dimnames(m1)[3]))    
    M <- array(NA, dim = c(dim(m1),length(anito)), dimnames=list(x=unlist(dimnames(m1)[1]), y=unlist(dimnames(m1)[2]), t=JD,
                                                                 l=landas_sat, w=anito))
    M[,,,,1] <- m1
    M[,,,,2] <- m2
    M[,,,,3] <- m3
    M[,,,,4] <- m4
    M[,,,,5] <- m5
    M[,,,,6] <- m6   
    
    dim(M)
    med <- apply(M , 1:4 , mean, na.rm=TRUE)
    rm(M)
    dim(med)
    dimnames(med) <- list(x=as.numeric(unlist(dimnames(m1)[1])),
                          y=as.numeric(unlist(dimnames(m1)[2])),
                          t=JD, l=landas_sat)
    med[is.nan(med)] <- NA
    save(med, file =paste(o_dir,"climatologies_2012_2018/media_seasonal_Aph_2deg_OCCCI_2012_2018.RData",sep="") )
    


    ## ANNUAL ##  
    #source("/Users/ealvarez/Programas/Funciones_varios/functions.r")
    #o_dir <- "/Users/ealvarez/Datos/Dat_Satelite/"
    #load(paste(o_dir,"climatologies_2012_2018/media_seasonal_Aph_2deg_OCCCI_2012_2018.RData",sep="") )
    #dim(med)
    JD <- as.numeric(unlist(dimnames(med))) 
    res <- apply(med, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)
    dim(res)
    save(res, file =paste(o_dir,"climatologies_2012_2018/media_anual_Aph_2deg_OCCCI_2012_2018.RData",sep="") )
    
###################

            
            
            
## CLIMATOLOGY Adg
###################    
    ## SEASONAL    
    o_dir <- paste(OS,"Datos/Dat_Satelite/",sep="")
    anito<-c(2012:2017)
    landas_sat <- c("Rrs_412","Rrs_443","Rrs_490","Rrs_510","Rrs_555","Rrs_670")
    
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Adg_2deg_OCCCI_",anito[1],".RData",sep="") )
    m1<-res
    dim(m1)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Adg_2deg_OCCCI_",anito[2],".RData",sep="") )
    m2<-res
    dim(m2)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Adg_2deg_OCCCI_",anito[3],".RData",sep="") )
    m3<-res
    dim(m3)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Adg_2deg_OCCCI_",anito[4],".RData",sep="") )
    m4<-res
    dim(m4)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Adg_2deg_OCCCI_",anito[5],".RData",sep="") )
    m5<-res
    dim(m5)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Adg_2deg_OCCCI_",anito[6],".RData",sep="") )
    m6<-res
    dim(m6)
    
    JD <- as.numeric(unlist(dimnames(m1)[3]))    
    M <- array(NA, dim = c(dim(m1),length(anito)), dimnames=list(x=unlist(dimnames(m1)[1]), y=unlist(dimnames(m1)[2]), t=JD,
                                                                 l=landas_sat, w=anito))
    M[,,,,1] <- m1
    M[,,,,2] <- m2
    M[,,,,3] <- m3
    M[,,,,4] <- m4
    M[,,,,5] <- m5
    M[,,,,6] <- m6   
    
    dim(M)
    med <- apply(M , 1:4 , mean, na.rm=TRUE)
    rm(M)
    dimnames(med) <- list(x=as.numeric(unlist(dimnames(m1)[1])),
                          y=as.numeric(unlist(dimnames(m1)[2])),
                          t=JD, l=landas_sat)
    med[is.nan(med)] <- NA
    save(med, file =paste(o_dir,"climatologies_2012_2018/media_seasonal_Adg_2deg_OCCCI_2012_2018.RData",sep="") )
    
    ## ANNUAL ##  
    JD <- as.numeric(unlist(dimnames(med))) 
    res <- apply(med, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)
    dim(res)
    save(res, file =paste(o_dir,"climatologies_2012_2018/media_anual_Adg_2deg_OCCCI_2012_2018.RData",sep="") )
    
###################    
   
    
    
    
## CLIMATOLOGY Atot
###################    
    ## SEASONAL    
    o_dir <- paste(OS,"Datos/Dat_Satelite/",sep="")
    anito<-c(2012:2017)
    landas_sat <- c("Rrs_412","Rrs_443","Rrs_490","Rrs_510","Rrs_555","Rrs_670")
    
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Atot_2deg_OCCCI_",anito[1],".RData",sep="") )
    m1<-res
    dim(m1)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Atot_2deg_OCCCI_",anito[2],".RData",sep="") )
    m2<-res
    dim(m2)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Atot_2deg_OCCCI_",anito[3],".RData",sep="") )
    m3<-res
    dim(m3)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Atot_2deg_OCCCI_",anito[4],".RData",sep="") )
    m4<-res
    dim(m4)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Atot_2deg_OCCCI_",anito[5],".RData",sep="") )
    m5<-res
    dim(m5)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Atot_2deg_OCCCI_",anito[6],".RData",sep="") )
    m6<-res
    dim(m6)
    
    JD <- as.numeric(unlist(dimnames(m1)[3]))    
    M <- array(NA, dim = c(dim(m1),length(anito)), dimnames=list(x=unlist(dimnames(m1)[1]), y=unlist(dimnames(m1)[2]), t=JD,
                                                                 l=landas_sat, w=anito))
    M[,,,,1] <- m1
    M[,,,,2] <- m2
    M[,,,,3] <- m3
    M[,,,,4] <- m4
    M[,,,,5] <- m5
    M[,,,,6] <- m6   
    
    dim(M)
    med <- apply(M , 1:4 , mean, na.rm=TRUE)
    rm(M)
    dimnames(med) <- list(x=as.numeric(unlist(dimnames(m1)[1])),
                          y=as.numeric(unlist(dimnames(m1)[2])),
                          t=JD, l=landas_sat)
    med[is.nan(med)] <- NA
    save(med, file =paste(o_dir,"climatologies_2012_2018/media_seasonal_Atot_2deg_OCCCI_2012_2018.RData",sep="") )
    
    ## ANNUAL ##  
    JD <- as.numeric(unlist(dimnames(med))) 
    res <- apply(med, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)
    dim(res)
    save(res, file =paste(o_dir,"climatologies_2012_2018/media_anual_Atot_2deg_OCCCI_2012_2018.RData",sep="") )
    
###################     
    
    
    
    
## CLIMATOLOGY Bbp
###################    
    ## SEASONAL    
    o_dir <- paste(OS,"Datos/Dat_Satelite/",sep="")
    anito<-c(2012:2017)
    landas_sat <- c("Rrs_412","Rrs_443","Rrs_490","Rrs_510","Rrs_555","Rrs_670")
    
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Bbp_2deg_OCCCI_",anito[1],".RData",sep="") )
    m1<-res
    dim(m1)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Bbp_2deg_OCCCI_",anito[2],".RData",sep="") )
    m2<-res
    dim(m2)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Bbp_2deg_OCCCI_",anito[3],".RData",sep="") )
    m3<-res
    dim(m3)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Bbp_2deg_OCCCI_",anito[4],".RData",sep="") )
    m4<-res
    dim(m4)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Bbp_2deg_OCCCI_",anito[5],".RData",sep="") )
    m5<-res
    dim(m5)
    load(paste(o_dir,"res_seasonal_OCCCI_IOPs/media_anual_Bbp_2deg_OCCCI_",anito[6],".RData",sep="") )
    m6<-res
    dim(m6)
    
    JD <- as.numeric(unlist(dimnames(m1)[3]))    
    M <- array(NA, dim = c(dim(m1),length(anito)), dimnames=list(x=unlist(dimnames(m1)[1]), y=unlist(dimnames(m1)[2]), t=JD,
                                                                 l=landas_sat, w=anito))
    M[,,,,1] <- m1
    M[,,,,2] <- m2
    M[,,,,3] <- m3
    M[,,,,4] <- m4
    M[,,,,5] <- m5
    M[,,,,6] <- m6   
    
    dim(M)
    med <- apply(M , 1:4 , mean, na.rm=TRUE)
    rm(M)
    dimnames(med) <- list(x=as.numeric(unlist(dimnames(m1)[1])),
                          y=as.numeric(unlist(dimnames(m1)[2])),
                          t=JD, l=landas_sat)
    med[is.nan(med)] <- NA
    save(med, file =paste(o_dir,"climatologies_2012_2018/media_seasonal_Bbp_2deg_OCCCI_2012_2018.RData",sep="") )
    
    ## ANNUAL ##  
    JD <- as.numeric(unlist(dimnames(med))) 
    res <- apply(med, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)
    dim(res)
    save(res, file =paste(o_dir,"climatologies_2012_2018/media_anual_Bbp_2deg_OCCCI_2012_2018.RData",sep="") )
    
###################     
    