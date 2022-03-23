if (Sys.info()['sysname']=="Windows") {OS<-"C:/"} else {OS<-paste("/Users/",Sys.info()['user'],"/", sep="")}
#source(paste(OS,"Programas/Funciones_varios/functions.r",sep=""))
library(RNetCDF)
library(akima)
library(abind)
library(plot3D)

########################################################
### Interpolate SURFACE variables from model output  ###
########################################################

# Reflectance (x,y,l)
# Biomass: Chla (log and lineal), PPC:TChla, CDOM and Detr (x,y)
# IOP's surface (x,y,l)

    o_dir<-paste(OS,"Datos/Ser_Stan/global_14_aphyt/Res_model/",sep="")        
    p_dir<-paste(OS,"Datos/Ser_Stan/global_14_aphyt/Res_model/interpolated/",sep="")
# Input: 
# Output:     
    
    #path_figures<-paste(OS,"Datos/Res_C20_radtrans/radtrans_v5_marshall/",sep="")
    
    # List of simulations
    experimentos<-read.csv(paste(OS,"Datos/Ser_Stan/global_14_aphyt/run_log_marshall_PPC_PS_SA_G100.csv",sep=""), sep=",")
    modelnames<-experimentos$name[c(9:29)]
    ruta<-o_dir
    
    ##  experimentos$name[c(30,31)] are identical to experimentos$name[c(21,9)], just run in another machine (G100)
    
    
#################
## REFLECTANCE ##
#################    
    for (w in 1:length(modelnames)){
        modelname<-modelnames[w]
        o_dir<-ruta[w]     

        # Reflectance
        archivo <-"/RADiags3DFlux.nc"
        filename <- paste(o_dir,modelname,archivo, sep="")
        filenc <- open.nc(filename)
        filerc <- read.nc(filenc)            
              #names(filerc)        
              lamda_new <- seq(400,700,by=25)
              m<-filerc$reflecta
                    # convert from irradiance reflectance to remote sensing reflectance
                    # using a bidirectional function Q (Q=4 in Dutkiewicz 2015, Q=3 in Dutkiewicz 2018)
                    # convert to above-surface remotely sensed reflectance using the formula of Lee et al. (2002)
                    m <- m/3
                    m<-(0.52*m)/(1-1.7*m) 
                    m[m>0.5]<-NA
                    m[m<0.0]<-NA
                    
              m[m==0]<-NA
              #m <- m[,,,2:dim(m)[4]]
              #dim(m)
              dimnames(m) <- list(x=filerc$X, y=filerc$Y, l=lamda_new, t=(((filerc$T/60)/60)/24))
              fechas_recom_real  <-(((filerc$T/60)/60)/24)
              fechas_recom       <-as.numeric(fechas_recom_real)-((4*360))
              meses_recom        <-ceiling(fechas_recom/30)
    
        close.nc(filenc)
        
    # 2x2 degree horizontal resolution
    md<-array(data = NA, dim=c(180,90,dim(m)[3],dim(m)[4]), dimnames = NULL)   
    for (k in 1:dim(m)[4]) {         # date
        z<-m[,,,k]                   #  180 126 13
        ZETA<-abind(z[91:180,,],z[1:90,,], along=1)   #  180 126 10
        X<-c(filerc$X[91:180],filerc$X[1:90])   
        interpoladico<-array(data = NA, dim=c(180,90,dim(ZETA)[3]), dimnames = NULL) 
        for (j in 1:dim(ZETA)[3]) {     # lambda 
            # Interpolate in 2D: remove distortion in the SHemis
            interm<-ZETA[,,j]           #  180 126
            nuevo<-seq(-79,79,by=2)
            viejas<-as.numeric(rownames(interm))
            ZETA_NEW<-matrix(NA,ncol=length(nuevo), nrow=nrow(interm))
            rownames(ZETA_NEW)<-viejas
            colnames(ZETA_NEW)<-nuevo
            indices<-as.numeric(colnames(ZETA))
            
            for (i in 1:(length(nuevo)-1)){
                cuales <- which(indices>=nuevo[i] & indices<nuevo[i+1])
                mat<-interm[,cuales]
                mat <- matrix(mat,nrow=nrow(ZETA_NEW))
                ZETA_NEW[,i]<-apply(mat, MARGIN=1, mean, na.rm=TRUE)
                }  # end loop i
            
            # Add NAs in the Arctic (>80)    
            cajasur<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
            colnames(cajasur)<-seq(-89,-81,by=2)
            cajanorte<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
            colnames(cajanorte)<-seq(81,89,by=2)
            res<-cbind(cajasur, ZETA_NEW, cajanorte)   
            rownames(res)<-rownames(ZETA_NEW)          
            interpoladico[,,j] <- res  }  # end loop j lambda  #  180 90 13
        #range(res,na.rm=TRUE)
        md[,,,k] <- interpoladico    }  # end loop k date      #  180 90 13 12
    
    dimnames(md) <- list(x=X, y=colnames(res), l=lamda_new, t=fechas_recom)
    
    # Monthly resolution
    res <-array(data = NA, dim=c(dim(md)[1:3],length(unique(meses_recom))),
                dimnames = list(x=dimnames(md)$x, y=dimnames(md)$y,
                                l=dimnames(md)$l, t=unique(meses_recom)))  
    for (q in unique(meses_recom)){
        res[,,,q] <- apply(md[,,,which(meses_recom==q)], MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)}
    
    # SEASONAL
    save(res, file =paste(p_dir,"2reflectance/",modelname,".RData",sep="") )
    
    # ANNUAL 
    res <- apply(md, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)
    save(res, file =paste(p_dir,"1reflectance/",modelname,".RData",sep="") )
    
    } # end w modelrun
#################    
    

    
  
    
    
                
#######################        
## CLORO SURFACE LOG ##
#######################    
            
        for (w in 1:length(modelnames)){
        modelname<-modelnames[w]
        o_dir<-ruta[w]
                    # depth
                    archivo<-"/grid.nc"
                    filename <- paste(o_dir,modelname,archivo, sep="")
                    filen <- open.nc(filename)
                    filedepth <- read.nc(filen)       
                    profun<-filedepth$Z
                    profunU<-filedepth$Zu        
                    
                    archivo <-"/recomDiags3D10daily.nc"
                    filename <- paste(o_dir,modelname,archivo, sep="")
                    filenc <- open.nc(filename)
                    filerc <- read.nc(filenc)
                    
                    # Chla todo
                    m<-(filerc$TRAC06+filerc$TRAC15)
                    m[m==0]<-NA
                    dimnames(m) <- list(x=filerc$X, y=filerc$Y, z=profun, t=(((filerc$T/60)/60)/24))
                    fechas_recom_real  <-(((filerc$T/60)/60)/24)
                    fechas_recom       <-as.numeric(fechas_recom_real)-((4*360))
                    meses_recom        <-ceiling(fechas_recom/30)
                    
                    close.nc(filenc)
                    close.nc(filen)
                    rm(filedepth)
                    
         md<-array(data = NA, dim=c(180,90,dim(m)[4]), dimnames = NULL)   
         for (k in 1:dim(m)[4]) {  # una k por fecha
              z<-m[,,profun>=-(10),k]  #dim(z)
              z[which(z==0)] <- NA
              z[z<1.000000e-11] <- 0
              ZETA<-abind(z[91:180,],z[1:90,], along=1)   #  180 126 30
              X<-c(filerc$X[91:180],filerc$X[1:90])     
                        
                  # Interpolate in 2D
                  nuevo<-seq(-79,79,by=2)
                  viejas<-as.numeric(rownames(ZETA))
                  ZETA_NEW<-matrix(NA,ncol=length(nuevo), nrow=nrow(ZETA))
                  rownames(ZETA_NEW)<-viejas
                  colnames(ZETA_NEW)<-nuevo
                  indices<-as.numeric(colnames(ZETA))
                        
                        for (i in 1:(length(nuevo)-1)){
                            cuales <- which(indices>=nuevo[i] & indices<nuevo[i+1])
                            mat<-ZETA[,cuales]
                            mat <- matrix(mat,nrow=nrow(ZETA_NEW))
                            ZETA_NEW[,i]<-apply(mat, MARGIN=1, mean, na.rm=TRUE)
                        }  # end loop i
                        
                        # Add NAs Arctic 
                        cajasur<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                        colnames(cajasur)<-seq(-89,-81,by=2)
                        cajanorte<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                        colnames(cajanorte)<-seq(81,89,by=2)
                        res<-cbind(cajasur, ZETA_NEW, cajanorte)   
                        rownames(res)<-rownames(ZETA_NEW)          
                        md[,,k] <- res  }  # end loop k
          dimnames(md) <- list(x=X, y=colnames(res), t=fechas_recom)
          md <- log10(md)
          
          # Mothly resolution
          res <-array(data = NA, dim=c(180,90,length(unique(meses_recom))),
                                dimnames = list(x=X, y=colnames(res), t=unique(meses_recom)))  
          for (q in unique(meses_recom)){
               res[,,q] <- apply(md[,,which(meses_recom==q)], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)}
                    
          # SEASONAL
          save(res, file =paste(p_dir,"2cloro_log/",modelname,".RData",sep="") )
                    
          # ANNUAL 
          res <- apply(md, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)
          save(res, file =paste(p_dir,"1cloro_log/",modelname,".RData",sep="") )
                    
          }  # end loop w   
#######################                
                
               
     
              

##########################        
## CLORO SURFACE LINEAL ##
##########################   
   
         for (w in 1:length(modelnames)){
             modelname<-modelnames[w]
             o_dir<-ruta[w]
                    # depth
                    archivo<-"/grid.nc"
                    filename <- paste(o_dir,modelname,archivo, sep="")
                    filen <- open.nc(filename)
                    filedepth <- read.nc(filen)       
                    profun<-filedepth$Z
                    profunU<-filedepth$Zu        
                    
               archivo <-"/recomDiags3D10daily.nc"
               filename <- paste(o_dir,modelname,archivo, sep="")
               filenc <- open.nc(filename)
               filerc <- read.nc(filenc)
               #names(filerc)
                   
               # Chla todo
       m<-(filerc$TRAC06+filerc$TRAC15)
       #m<-filerc$TRAC15  # diatoms
       #m<-filerc$TRAC06  # phyto
               m[m==0]<-NA
               #dim(m)
               #image(m[,,1,3])
               # range(m,na.rm=TRUE)
               dimnames(m) <- list(x=filerc$X, y=filerc$Y, z=profun, t=(((filerc$T/60)/60)/24))
               fechas_recom_real  <-(((filerc$T/60)/60)/24)
               fechas_recom       <-as.numeric(fechas_recom_real)-((4*360))
               meses_recom        <-ceiling(fechas_recom/30)
                    
                    close.nc(filenc)
                    close.nc(filen)
                    rm(filedepth)
                    
        md<-array(data = NA, dim=c(180,90,dim(m)[4]), dimnames = NULL)   
        for (k in 1:dim(m)[4]) {
             z<-m[,,profun>=-(10),k]
             z[which(z==0)] <- NA
             z[z<1.000000e-11] <- 0
             ZETA<-abind(z[91:180,],z[1:90,], along=1)   #  180 126 30
             X<-c(filerc$X[91:180],filerc$X[1:90])     
                        
                        # Interpolate in 2D
                        nuevo<-seq(-79,79,by=2)
                        viejas<-as.numeric(rownames(ZETA))
                        ZETA_NEW<-matrix(NA,ncol=length(nuevo), nrow=nrow(ZETA))
                        rownames(ZETA_NEW)<-viejas
                        colnames(ZETA_NEW)<-nuevo
                        indices<-as.numeric(colnames(ZETA))
                        
                        for (i in 1:(length(nuevo)-1)){
                            cuales <- which(indices>=nuevo[i] & indices<nuevo[i+1])
                            mat<-ZETA[,cuales]
                            mat <- matrix(mat,nrow=nrow(ZETA_NEW))
                            ZETA_NEW[,i]<-apply(mat, MARGIN=1, mean, na.rm=TRUE)
                        }  # end loop i
                        
                        # Add NAs Arctic     
                        cajasur<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                        colnames(cajasur)<-seq(-89,-81,by=2)
                        cajanorte<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                        colnames(cajanorte)<-seq(81,89,by=2)
                        res<-cbind(cajasur, ZETA_NEW, cajanorte)   
                        rownames(res)<-rownames(ZETA_NEW)          
                        md[,,k] <- res  }  # end loop k
          dimnames(md) <- list(x=X, y=colnames(res), t=fechas_recom)
        
          # Monthly res
          res <-array(data = NA, dim=c(180,90,length(unique(meses_recom))),
                      dimnames = list(x=X, y=colnames(res), t=unique(meses_recom)))  
          for (q in unique(meses_recom)){
              res[,,q] <- apply(md[,,which(meses_recom==q)], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)}
                    
           # SEASONAL
        save(res, file =paste(p_dir,"2cloro_lin/",modelname,".RData",sep="") )
        #save(res, file =paste(p_dir,"2cloro_dia/",modelname,".RData",sep="") )
        #save(res, file =paste(p_dir,"2cloro_phy/",modelname,".RData",sep="") )
           rm(filerc)
           # ANNUAL 
        #res <- apply(md, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)
        #save(res, file =paste(p_dir,"1cloro_lin/",modelname,".RData",sep="") )
                    
           }  # end loop w   
##########################                
                
         
    
    

###################################        
## CDOM and DETR. SURFACE LINEAL ##
###################################   

      for (w in 1:length(modelnames)){
          modelname<-modelnames[w]
          o_dir<-ruta[w]             
             # depth
             archivo<-"/grid.nc"
             filename <- paste(o_dir,modelname,archivo, sep="")
             filen <- open.nc(filename)
             filedepth <- read.nc(filen)       
             profun<-filedepth$Z
             profunU<-filedepth$Zu        
             
             archivo <-"/recomDiags3D10daily.nc"
             filename <- paste(o_dir,modelname,archivo, sep="")
             filenc <- open.nc(filename)
             filerc <- read.nc(filenc)
             #names(filerc)
             
             # CDOM
        #m<-filerc$TRAC22
             # Detritus
        m<-filerc$TRAC08
             m[m==0]<-NA
             dimnames(m) <- list(x=filerc$X, y=filerc$Y, z=profun, t=(((filerc$T/60)/60)/24))
             fechas_recom_real  <-(((filerc$T/60)/60)/24)
             fechas_recom       <-as.numeric(fechas_recom_real)-((4*360))
             meses_recom        <-ceiling(fechas_recom/30)
             
             close.nc(filenc)
             close.nc(filen)
             rm(filedepth)
             
             md<-array(data = NA, dim=c(180,90,dim(m)[4]), dimnames = NULL)   
             for (k in 1:dim(m)[4]) {
                 z<-m[,,profun>=-(10),k]
                 z[which(z==0)] <- NA
                 z[z<1.000000e-11] <- 0
                 ZETA<-abind(z[91:180,],z[1:90,], along=1)   #  180 126 30
                 X<-c(filerc$X[91:180],filerc$X[1:90])     
                 
                 # Interpolate in 2D
                 nuevo<-seq(-79,79,by=2)
                 viejas<-as.numeric(rownames(ZETA))
                 ZETA_NEW<-matrix(NA,ncol=length(nuevo), nrow=nrow(ZETA))
                 rownames(ZETA_NEW)<-viejas
                 colnames(ZETA_NEW)<-nuevo
                 indices<-as.numeric(colnames(ZETA))
                 
                 for (i in 1:(length(nuevo)-1)){
                     cuales <- which(indices>=nuevo[i] & indices<nuevo[i+1])
                     mat<-ZETA[,cuales]
                     mat <- matrix(mat,nrow=nrow(ZETA_NEW))
                     ZETA_NEW[,i]<-apply(mat, MARGIN=1, mean, na.rm=TRUE)
                 }  # end loop i
                 
                 # Add NAs Arctic
                 cajasur<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                 colnames(cajasur)<-seq(-89,-81,by=2)
                 cajanorte<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                 colnames(cajanorte)<-seq(81,89,by=2)
                 res<-cbind(cajasur, ZETA_NEW, cajanorte)   
                 rownames(res)<-rownames(ZETA_NEW)          
                 md[,,k] <- res  }  # end loop k
             dimnames(md) <- list(x=X, y=colnames(res), t=fechas_recom)

             # Monthly res
             res <-array(data = NA, dim=c(180,90,length(unique(meses_recom))),
                         dimnames = list(x=X, y=colnames(res), t=unique(meses_recom)))  
             for (q in unique(meses_recom)){
                 res[,,q] <- apply(md[,,which(meses_recom==q)], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)}
             
             # SEASONAL
        #save(res, file =paste(p_dir,"2cdom_surface/",modelname,".RData",sep="") )
        save(res, file =paste(p_dir,"2detr_surface/",modelname,".RData",sep="") )

             # ANNUAL 
        #res <- apply(md, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)
        #save(res, file =paste(p_dir,"1cloro_lin/",modelname,".RData",sep="") )
             
             }  # end loop w   
###################################             
             
             
             
             

#################        
### IOP's 4D  ###   spectral IOPs in (x,y): it needs to export "surface (sur)" diagnostics (data.diagnostics)
#################      

    for (g in 1:length(modelnames)){
          modelname<-modelnames[g]
          o_dir<-ruta[g]
             
             # Depth
             archivo<-"/grid.nc"
             filename <- paste(o_dir,modelname,archivo, sep="")
             filen <- open.nc(filename)
             filedepth <- read.nc(filen)       
             profun<-filedepth$Z
             
             # Spectral IOP's in surface # m-1
             archivo <-"/RADiags3DFlux.nc"
             filename <- paste(o_dir,modelname,archivo, sep="")
             filenc <- open.nc(filename)
             filerc <- read.nc(filenc)
             #names(filerc)
             longitude<-filerc$X
             latitude<-filerc$Y
             lamda_new <- seq(400,700,by=25)
             fechas_recom_real  <-(((filerc$T/60)/60)/24)
             fechas_recom       <-as.numeric(fechas_recom_real)-((4*360))
             meses_recom        <-ceiling(fechas_recom/30)            

             
        variables <-c("abtotsur","acdomsur","abparsur","abphysur",
                      #"bttotsur","btparsur","btphysur",
                      "bbtotsur")
        nombres <-c("atot", "acd", "apt", "aph",
                    #"btot", "bpt", "bph",
                    "bbp")        
        #cbind(variables, nombres)
        
        for (w in 1:length(variables)){
            
          if (file.exists(paste(p_dir,"2iop_",nombres[w],"/", sep=""))==FALSE){
          dir.create(paste(p_dir,"2iop_",nombres[w],"/", sep=""), recursive=TRUE)}
         
                 m <- filerc[[which(names(filerc)==variables[w])]]   # m-1
                 m[m==0]<-NA
        dimnames(m) <- list(x=filerc$X, y=filerc$Y, l=lamda_new, t=fechas_recom_real)
        md<-array(data = NA, dim=c(180,90,dim(m)[3],dim(m)[4]), dimnames = NULL)   
                 for (k in 1:dim(m)[4]) {     # k date
                     z<-m[,,,k]               #  180 126 13
                     ZETA<-abind(z[91:180,,],z[1:90,,], along=1)
                     X<-c(filerc$X[91:180],filerc$X[1:90])   
                     interpoladico<-array(data = NA, dim=c(180,90,dim(ZETA)[3]), dimnames = NULL) 
                     for (j in 1:dim(ZETA)[3]) {    # j lambda 
                         # Interpolate in 2D
                         interm<-ZETA[,,j]          #  180 126
                         nuevo<-seq(-79,79,by=2)
                         viejas<-as.numeric(rownames(interm))
                         ZETA_NEW<-matrix(NA,ncol=length(nuevo), nrow=nrow(interm))
                         rownames(ZETA_NEW)<-viejas
                         colnames(ZETA_NEW)<-nuevo
                         indices<-as.numeric(colnames(ZETA))
                         
                         for (i in 1:(length(nuevo)-1)){
                             cuales <- which(indices>=nuevo[i] & indices<nuevo[i+1])
                             mat<-interm[,cuales]
                             mat <- matrix(mat,nrow=nrow(ZETA_NEW))
                             ZETA_NEW[,i]<-apply(mat, MARGIN=1, mean, na.rm=TRUE)
                          }  # end loop i
                         
                         # Add NAs Arctic    
                         cajasur<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                         colnames(cajasur)<-seq(-89,-81,by=2)
                         cajanorte<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                         colnames(cajanorte)<-seq(81,89,by=2)
                         res<-cbind(cajasur, ZETA_NEW, cajanorte)   
                         rownames(res)<-rownames(ZETA_NEW)          
                         interpoladico[,,j] <- res  }  # end loop j    #  180 90 13
                     md[,,,k] <- interpoladico    }  # end loop k      #  180 90 13 12
                 
                 dimnames(md) <- list(x=X, y=colnames(res), l=lamda_new, t=fechas_recom)
                 
                 # Monthly res
                 res <-array(data = NA, dim=c(dim(md)[1:3],length(unique(meses_recom))),
                             dimnames = list(x=dimnames(md)$x, y=dimnames(md)$y,
                                             l=dimnames(md)$l, t=unique(meses_recom)))  
                 for (q in unique(meses_recom)){
                     res[,,,q] <- apply(md[,,,which(meses_recom==q)], MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)}
                 
        # SEASONAL
        save(res, file=paste(p_dir,"2iop_",nombres[w],"/",modelname,".RData",sep="") )
        }  # end loop w    (variable)       
             

             
##### Adg  (sum detr + cdom) #############
   
        if (file.exists(paste(p_dir,"2iop_adg/", sep=""))==FALSE){
          dir.create(paste(p_dir,"2iop_adg/", sep=""), recursive=TRUE)}
        
            #names(filerc)
            m <- filerc$acdomsur+filerc$abparsur
            m[m==0]<-NA
            dimnames(m) <- list(x=filerc$X, y=filerc$Y, l=lamda_new, t=fechas_recom_real)
            md<-array(data = NA, dim=c(180,90,dim(m)[3],dim(m)[4]), dimnames = NULL)   
            for (k in 1:dim(m)[4]) {      # k date
                z<-m[,,,k]                # 180 126 13
                ZETA<-abind(z[91:180,,],z[1:90,,], along=1)
                X<-c(filerc$X[91:180],filerc$X[1:90])   
                interpoladico<-array(data = NA, dim=c(180,90,dim(ZETA)[3]), dimnames = NULL) 
                for (j in 1:dim(ZETA)[3]) {     # j lambda 
                    # Interpolate in 2D
                    interm<-ZETA[,,j]           # 180 126
                    nuevo<-seq(-79,79,by=2)
                    viejas<-as.numeric(rownames(interm))
                    ZETA_NEW<-matrix(NA,ncol=length(nuevo), nrow=nrow(interm))
                    rownames(ZETA_NEW)<-viejas
                    colnames(ZETA_NEW)<-nuevo
                    indices<-as.numeric(colnames(ZETA))
                    
                    for (i in 1:(length(nuevo)-1)){
                        cuales <- which(indices>=nuevo[i] & indices<nuevo[i+1])
                        mat<-interm[,cuales]
                        mat <- matrix(mat,nrow=nrow(ZETA_NEW))
                        ZETA_NEW[,i]<-apply(mat, MARGIN=1, mean, na.rm=TRUE)
                      }  # end loop i

                    # Add NAs Arctic    
                    cajasur<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                    colnames(cajasur)<-seq(-89,-81,by=2)
                    cajanorte<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                    colnames(cajanorte)<-seq(81,89,by=2)
                    res<-cbind(cajasur, ZETA_NEW, cajanorte)   
                    rownames(res)<-rownames(ZETA_NEW)          
                    interpoladico[,,j] <- res  }  # end loop j     #  180 90 13
                md[,,,k] <- interpoladico    }  # end loop k       #  180 90 13 12
            
            dimnames(md) <- list(x=X, y=colnames(res), l=lamda_new, t=fechas_recom)
            
            # Monthly res
            res <-array(data = NA, dim=c(dim(md)[1:3],length(unique(meses_recom))),
                        dimnames = list(x=dimnames(md)$x, y=dimnames(md)$y,l=dimnames(md)$l, t=unique(meses_recom)))  
            for (q in unique(meses_recom)){
                res[,,,q] <- apply(md[,,,which(meses_recom==q)], MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)}
 
       save(res, file=paste(p_dir,"2iop_adg/",modelname,".RData",sep="") )
       # end Adg      
           
       
       
##### App  (sum phyto + detr) #############

       if (file.exists(paste(p_dir,"2iop_app/", sep=""))==FALSE){
         dir.create(paste(p_dir,"2iop_app/", sep=""), recursive=TRUE)}
       
       #names(filerc)
       m <- filerc$abphysur+filerc$abparsur
       m[m==0]<-NA
       dimnames(m) <- list(x=filerc$X, y=filerc$Y, l=lamda_new, t=fechas_recom_real)
       md<-array(data = NA, dim=c(180,90,dim(m)[3],dim(m)[4]), dimnames = NULL)   
       for (k in 1:dim(m)[4]) {      # k date
         z<-m[,,,k]                  # 180 126 13
         ZETA<-abind(z[91:180,,],z[1:90,,], along=1)
         X<-c(filerc$X[91:180],filerc$X[1:90])   
         interpoladico<-array(data = NA, dim=c(180,90,dim(ZETA)[3]), dimnames = NULL) 
         for (j in 1:dim(ZETA)[3]) {   # j lambda 
           # Interpolate in 2D
           interm<-ZETA[,,j]           #  180 126
           nuevo<-seq(-79,79,by=2)
           viejas<-as.numeric(rownames(interm))
           ZETA_NEW<-matrix(NA,ncol=length(nuevo), nrow=nrow(interm))
           rownames(ZETA_NEW)<-viejas
           colnames(ZETA_NEW)<-nuevo
           indices<-as.numeric(colnames(ZETA))
           
           for (i in 1:(length(nuevo)-1)){
             cuales <- which(indices>=nuevo[i] & indices<nuevo[i+1])
             mat<-interm[,cuales]
             mat <- matrix(mat,nrow=nrow(ZETA_NEW))
             ZETA_NEW[,i]<-apply(mat, MARGIN=1, mean, na.rm=TRUE)
            }  # end loop i
           
           # Add NAs Arctic     
           cajasur<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
           colnames(cajasur)<-seq(-89,-81,by=2)
           cajanorte<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
           colnames(cajanorte)<-seq(81,89,by=2)
           res<-cbind(cajasur, ZETA_NEW, cajanorte)   
           rownames(res)<-rownames(ZETA_NEW)          
           interpoladico[,,j] <- res  }  # end loop j    #  180 90 13
         md[,,,k] <- interpoladico    }  # end loop k    #  180 90 13 12
       
       dimnames(md) <- list(x=X, y=colnames(res), l=lamda_new, t=fechas_recom)
       
       # Monthly res
       res <-array(data = NA, dim=c(dim(md)[1:3],length(unique(meses_recom))),
                   dimnames = list(x=dimnames(md)$x, y=dimnames(md)$y,l=dimnames(md)$l, t=unique(meses_recom)))  
       for (q in unique(meses_recom)){
         res[,,,q] <- apply(md[,,,which(meses_recom==q)], MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)}
       
       save(res, file=paste(p_dir,"2iop_app/",modelname,".RData",sep="") )
       # end App     
       
 
    }  # end loop g   (modelname)                
#################                         
             
        
        
        
        
########        
## PPC        
########           
        
        for (w in 1:length(modelnames)){
          modelname<-modelnames[w]
          o_dir<-ruta[w]
          # depth
          archivo<-"/grid.nc"
          filename <- paste(o_dir,modelname,archivo, sep="")
          filen <- open.nc(filename)
          filedepth <- read.nc(filen)       
          profun<-filedepth$Z
          profunU<-filedepth$Zu        
          
                    # png(file=paste(o_dir,modelname,"/D1_alpha_PPC.png",sep=""),width = 700, height = 700, units = "px", pointsize = 16, bg = "white")  
                    # par(mfrow=c(3,2))
                    # par(mar=c(2,2,1,1))
                    # layout.show(6)
                    
              # Biomass
              archivo <-"/recomDiags3D10daily.nc"
              filename <- paste(o_dir,modelname,archivo, sep="")
              filenc <- open.nc(filename)
              filerc <- read.nc(filenc)
              carbon_phy<-filerc$TRAC05
              carbon_dia<-filerc$TRAC14          
              
                    ### D1          
                    # m<-filerc$TRAC23
                    # m[m==0]<-NA
                    # dim(m)
                    # range(m,na.rm=T)
                    # image2D(m[,,1,1],  col=jet(9), las=1, cex.axis=0.8, main="D1 small phyto", colkey=FALSE,xaxt="n", yaxt="n",zlim=c(0.0,1))          
                    # m<-filerc$TRAC24
                    # m[m==0]<-NA
                    # dim(m)
                    # range(m,na.rm=T)
                    # image2D(m[,,1,1],  col=jet(9), las=1, cex.axis=0.8, main="D1 diatoms", colkey=FALSE,xaxt="n", yaxt="n",zlim=c(0.0,1))

              archivo <-"/MAREDiags3D10daily.nc"
              filename <- paste(o_dir,modelname,archivo, sep="")
              filenc <- open.nc(filename)
              filerc <- read.nc(filenc)
              #names(filerc)
              
                    ### alpha
                    # m<-filerc$alphaphy
                    # m[m==0]<-NA
                    # dim(m)
                    # range(m,na.rm=T)
                    # image2D(m[,,1,1],  col=jet(9), las=1, cex.axis=0.8, main="alpha small phyto", colkey=FALSE,xaxt="n", yaxt="n",zlim=c(0.0,0.14))          
                    # m<-filerc$alphadia
                    # m[m==0]<-NA
                    # dim(m)
                    # range(m,na.rm=T)
                    # image2D(m[,,1,1],  col=jet(9), las=1, cex.axis=0.8, main="alpha diatoms", colkey=FALSE,xaxt="n", yaxt="n",zlim=c(0.0,0.19))

                    ### PPC
                    # m<-filerc$ppcphy
                    # m[m==0]<-NA
                    # dim(m)
                    # range(m,na.rm=T)
                    # image2D(m[,,1,1],  col=jet(9), las=1, cex.axis=0.8, main="PPC small phyto", colkey=FALSE,xaxt="n", yaxt="n",zlim=c(0.0,1))          
                    # m<-filerc$ppcdia
                    # m[m==0]<-NA
                    # dim(m)
                    # range(m,na.rm=T)
                    # image2D(m[,,1,1],  col=jet(9), las=1, cex.axis=0.8, main="PPC diatoms", colkey=FALSE,xaxt="n", yaxt="n",zlim=c(0.0,1))

              #dev.off()
          
          # PPC total
          m_m<-((filerc$alphaphy*carbon_phy)+(filerc$alphadia*carbon_dia))/(carbon_phy+carbon_dia)
          m_g<-((0.14*carbon_phy)+(0.19*carbon_dia))/(carbon_phy+carbon_dia)
          m<-1-(m_m/m_g)
          m[m==0]<-NA
          dimnames(m) <- list(x=filerc$X, y=filerc$Y, z=profun, t=(((filerc$T/60)/60)/24))
          fechas_recom_real  <-(((filerc$T/60)/60)/24)
          fechas_recom       <-as.numeric(fechas_recom_real)-((4*360))
          meses_recom        <-ceiling(fechas_recom/30)
              close.nc(filenc)
              close.nc(filen)
              rm(filedepth)
              
          md<-array(data = NA, dim=c(180,90,dim(m)[4]), dimnames = NULL)   
              for (k in 1:dim(m)[4]) {
                z<-m[,,profun>=-(10),k]
                z[which(z==0)] <- NA
                z[z<1.000000e-11] <- 0
                ZETA<-abind(z[91:180,],z[1:90,], along=1)     #  180 126 30
                X<-c(filerc$X[91:180],filerc$X[1:90])     
                
                # Interpolate in 2D
                nuevo<-seq(-79,79,by=2)
                viejas<-as.numeric(rownames(ZETA))
                ZETA_NEW<-matrix(NA,ncol=length(nuevo), nrow=nrow(ZETA))
                rownames(ZETA_NEW)<-viejas
                colnames(ZETA_NEW)<-nuevo
                indices<-as.numeric(colnames(ZETA))
                
                for (i in 1:(length(nuevo)-1)){
                  cuales <- which(indices>=nuevo[i] & indices<nuevo[i+1])
                  mat<-ZETA[,cuales]
                  mat <- matrix(mat,nrow=nrow(ZETA_NEW))
                  ZETA_NEW[,i]<-apply(mat, MARGIN=1, mean, na.rm=TRUE)
                }  # end loop i
                
                # Add NAs Arctic    
                cajasur<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                colnames(cajasur)<-seq(-89,-81,by=2)
                cajanorte<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                colnames(cajanorte)<-seq(81,89,by=2)
                res<-cbind(cajasur, ZETA_NEW, cajanorte)   
                rownames(res)<-rownames(ZETA_NEW)          
                md[,,k] <- res  }  # end loop k
          dimnames(md) <- list(x=X, y=colnames(res), t=fechas_recom)

          # Monthly res
          res <-array(data = NA, dim=c(180,90,length(unique(meses_recom))),
                      dimnames = list(x=X, y=colnames(res), t=unique(meses_recom)))  
          for (q in unique(meses_recom)){
            res[,,q] <- apply(md[,,which(meses_recom==q)], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)}
          
          rm(filerc)
          # SEASONAL
          save(res, file =paste(p_dir,"2ppc_total/",modelname,".RData",sep="") )
          # ANNUAL 
          res <- apply(md, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)
          save(res, file =paste(p_dir,"1ppc_total/",modelname,".RData",sep="") )
          
        }  # end loop w   
########                
        
        
  