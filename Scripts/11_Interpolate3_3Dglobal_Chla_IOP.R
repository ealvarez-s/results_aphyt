if (Sys.info()['sysname']=="Windows") {OS<-"C:/"} else {OS<-paste("/Users/",Sys.info()['user'],"/", sep="")}
source(paste(OS,"Programas/Funciones_varios/functions.r",sep=""))
library(RNetCDF)
library(akima)
library(abind)
library(plot3D)

######################################################
## Interpolate GLOBAL 3D fields from model output  ###
######################################################

# Biomass: Chla 3D, Chl:C and PPC:TChla.
# IOP's 3D: only one waveband (450nm)

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

        
    
    
###################
##  CHLA global  ##
###################   

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
        
        ## Clorofila ##
        archivo <-"/recomDiags3D10daily.nc"
        filename <- paste(o_dir,modelname,archivo, sep="")
        filenc <- open.nc(filename)
        filerc <- read.nc(filenc)
        #names(filerc)
        m <- filerc$TRAC06 + filerc$TRAC15   # mg Chla
        dimnames(m)      <- list(x=filerc$X, y=filerc$Y, z=profun, t=(((filerc$T/60)/60)/24))
        
        fechas_recom_real  <-(((filerc$T/60)/60)/24)
        fechas_recom       <-as.numeric(fechas_recom_real)-((4*360))
        meses_recom        <-ceiling(fechas_recom/30)
        
        md<-array(data = NA, dim=c(180,90,30,dim(m)[4]), dimnames = NULL)   
        for (k in 1:dim(m)[4]) {   # date
            z<-m[,,,k]             # 180 126 30
            range(z,na.rm=TRUE)
            ZETA<-abind(z[91:180,,],z[1:90,,], along=1)   #  180 126 30
            X<-c(filerc$X[91:180],filerc$X[1:90])   
            range(ZETA,na.rm=TRUE)
            interpoladico<-array(data = NA, dim=c(180,90,dim(ZETA)[3]), dimnames = NULL) 
            for (j in 1:dim(ZETA)[3]) {       #  depth 
                # Interpolate in 2D
                interm<-ZETA[,,j]             #  180 126
                range(interm,na.rm=TRUE)
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
                range(ZETA_NEW,na.rm=TRUE)
                
                # Add NAs Arctic  
                cajasur<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                colnames(cajasur)<-seq(-89,-81,by=2)
                cajanorte<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                colnames(cajanorte)<-seq(81,89,by=2)
                res<-cbind(cajasur, ZETA_NEW, cajanorte)   
                rownames(res)<-rownames(ZETA_NEW)          
                interpoladico[,,j] <- res  }  # end loop j    #  180 90 30
            
            md[,,,k] <- interpoladico  }  # end loop k        #  180 90 30 36

        dimnames(md) <- list(x=X, y=colnames(res), z=profun, t=fechas_recom)
        # 10daily 
        #res <- md
        #save(res, file =paste(p_dir,"3cloro_global/",modelname,".RData",sep="") )
        
        # Monthly res
        res <-array(data = NA, dim=c(180,90,30,length(unique(meses_recom))),
                    dimnames = list(x=X, y=colnames(res), z=profun, t=unique(meses_recom)))  
        for (q in unique(meses_recom)){
            res[,,,q] <- apply(md[,,,which(meses_recom==q)], MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)}
        
        # SEASONAL
        save(res, file =paste(p_dir,"2cloro_global/",modelname,".RData",sep="") )
        
        # ANNUAL 
        res <- apply(md, MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)
        save(res, file =paste(p_dir,"1cloro_global/",modelname,".RData",sep="") )
        
    }  # end loop w        
###################    

    

###########################
##  Chla:C (g:g) global  ##    
###########################   

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
        
        ## Clorofila:Carbon ##
        archivo <-"/recomDiags3D10daily.nc"
        filename <- paste(o_dir,modelname,archivo, sep="")
        filenc <- open.nc(filename)
        filerc <- read.nc(filenc)
        m <- (filerc$TRAC06+filerc$TRAC15)/((filerc$TRAC05+filerc$TRAC14)*12.0107000000057)
        dimnames(m)      <- list(x=filerc$X, y=filerc$Y, z=profun, t=(((filerc$T/60)/60)/24))
        m[m==0]<-NA
        m[m>1] <- NA
        
        fechas_recom_real  <-(((filerc$T/60)/60)/24)
        fechas_recom       <-as.numeric(fechas_recom_real)-((4*360))
        meses_recom        <-ceiling(fechas_recom/30)
        
        md<-array(data = NA, dim=c(180,90,10,dim(m)[4]), dimnames = NULL)   
        for (k in 1:dim(m)[4]) {          # date
            z<-m[,,profun>=-(220),k]      # 180 126 10
            range(z,na.rm=TRUE)
            #z[which(z<=1.0e-11)] <- 0
            #z[which(z==0)] <- NA
            ZETA<-abind(z[91:180,,],z[1:90,,], along=1)   #  180 126 10
            X<-c(filerc$X[91:180],filerc$X[1:90])   
            range(ZETA,na.rm=TRUE)
            interpoladico<-array(data = NA, dim=c(180,90,dim(ZETA)[3]), dimnames = NULL) 
            for (j in 1:dim(ZETA)[3]) {    #  depth
                # Interpolate in 2D
                interm<-ZETA[,,j]          #  180 126
                range(interm,na.rm=TRUE)
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
                range(ZETA_NEW,na.rm=TRUE)
                
                # Anadir bloque NAS artico     
                cajasur<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                colnames(cajasur)<-seq(-89,-81,by=2)
                cajanorte<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                colnames(cajanorte)<-seq(81,89,by=2)
                res<-cbind(cajasur, ZETA_NEW, cajanorte)   
                rownames(res)<-rownames(ZETA_NEW)          
                interpoladico[,,j] <- res  }  # end loop j    #  180 90 10
            
            md[,,,k] <- interpoladico  }  # end loop k        #  180 90 10 36

        dimnames(md) <- list(x=X, y=colnames(res), z=profun[profun>=-(220)], t=fechas_recom)
        
        # 10daily 
        #res <- md
        #save(res, file =paste(p_dir,"3ratio_global/",modelname,".RData",sep="") )
        
        # Monthly res
        res <-array(data = NA, dim=c(180,90,10,length(unique(meses_recom))),
                    dimnames = list(x=X, y=colnames(res), z=profun[profun>=-(220)], t=unique(meses_recom)))  
        for (q in unique(meses_recom)){
            res[,,,q] <- apply(md[,,,which(meses_recom==q)], MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)}
        
        # SEASONAL
        save(res, file =paste(p_dir,"2ratio_global/",modelname,".RData",sep="") )
        
        # ANNUAL 
        res <- apply(md, MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)
        save(res, file =paste(p_dir,"1ratio_global/",modelname,".RData",sep="") )
        
        }  # end loop w        
###########################    


    
##############################
##  PPC:TChla (g:g) global  ##    
##############################   
    
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
      
      ## PPC:Clorofila ##
      archivo <-"/recomDiags3D10daily.nc"
      filename <- paste(o_dir,modelname,archivo, sep="")
      filenc <- open.nc(filename)
      filerc <- read.nc(filenc)
      carbon_phy<-filerc$TRAC05
      carbon_dia<-filerc$TRAC14          
      
      # PPC total
      archivo <-"/MAREDiags3D10daily.nc"
      filename <- paste(o_dir,modelname,archivo, sep="")
      filenc <- open.nc(filename)
      filerc <- read.nc(filenc)
      #names(filerc)
          m_m<-((filerc$alphaphy*carbon_phy)+(filerc$alphadia*carbon_dia))/(carbon_phy+carbon_dia)
          m_g<-((0.14*carbon_phy)+(0.19*carbon_dia))/(carbon_phy+carbon_dia)
          m<-1-(m_m/m_g)
          m[m==0]<-NA      
          dim(m)
          #range(m_g, na.rm=TRUE)
          fechas_recom_real  <-(((filerc$T/60)/60)/24)
          fechas_recom       <-as.numeric(fechas_recom_real)-((4*360))
          meses_recom        <-ceiling(fechas_recom/30)
          #if (dim(m)[4]!=36) { m<-m[,,,2:37]}
          dimnames(m)<-list(x=c(filerc$X),x=c(filerc$Y),z=profun, t=fechas_recom)
          
      md<-array(data = NA, dim=c(180,90,10,dim(m)[4]), dimnames = NULL)   
      for (k in 1:dim(m)[4]) {        # date
        z<-m[,,profun>=-(220),k]      # 180 126 10
        range(z,na.rm=TRUE)
        ZETA<-abind(z[91:180,,],z[1:90,,], along=1)    #  180 126 10
        X<-c(filerc$X[91:180],filerc$X[1:90])   
        range(ZETA,na.rm=TRUE)
        interpoladico<-array(data = NA, dim=c(180,90,dim(ZETA)[3]), dimnames = NULL) 
        for (j in 1:dim(ZETA)[3]) {      # depth
          # Interpolate in 2D
          interm<-ZETA[,,j]              #  180 126
          range(interm,na.rm=TRUE)
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
          interpoladico[,,j] <- res  }  # end loop j      #  180 90 10
        
        md[,,,k] <- interpoladico  }  # end loop k        #  180 90 10 36

      dimnames(md) <- list(x=X, y=colnames(res), z=profun[profun>=-(220)], t=fechas_recom)
      
      # 10daily 
      #res <- md
      #save(res, file =paste(p_dir,"3ppc_global/",modelname,".RData",sep="") )
      
      # Monthly res
      res <-array(data = NA, dim=c(180,90,10,length(unique(meses_recom))),
                  dimnames = list(x=X, y=colnames(res), z=profun[profun>=-(220)], t=unique(meses_recom)))  
      for (q in unique(meses_recom)){
        res[,,,q] <- apply(md[,,,which(meses_recom==q)], MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)}
      
      # SEASONAL
      save(res, file =paste(p_dir,"2ppc_global/",modelname,".RData",sep="") )
      
      # ANNUAL 
      res <- apply(md, MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)
      save(res, file =paste(p_dir,"1ppc_global/",modelname,".RData",sep="") )
      
    }  # end loop w        
##############################    
    
    
    
#######################                
###  IOP's GLOBAL   ###     only one waveband (450 nm)                         
#######################     
    
    modelnames<-experimentos$name[c(21,9)]
    ruta<-o_dir  
    
    for (g in 1:length(modelnames)){
        modelname<-modelnames[g]
        o_dir<-ruta[g]
        # depth
        archivo<-"/grid.nc"
        filename <- paste(o_dir,modelname,archivo, sep="")
        filen <- open.nc(filename)
        filedepth <- read.nc(filen)       
        profun<-filedepth$Z
        profunU<-filedepth$Zu        
        
        ## IOP's ##
        archivo <-"/RADiags3D10daily.nc"
        filename <- paste(o_dir,modelname,archivo, sep="")
        filenc <- open.nc(filename)
        filerc <- read.nc(filenc)
        fechas_recom_real  <-(((filerc$T/60)/60)/24)
        fechas_recom       <-as.numeric(fechas_recom_real)-((4*360))
        meses_recom        <-ceiling(fechas_recom/30) 
        
        variables <-c("abtot_wb", "acdom_wb", "abpar_wb", "abphy_wb",
                      "bttot_wb", "btpar_wb", "btphy_wb", "bbtot_wb")
        
 
            for (w in 1:length(variables)){
                if (file.exists(paste(p_dir,"2iop_",variables[w],"/", sep=""))==FALSE){
                dir.create(paste(p_dir,"2iop_",variables[w],"/", sep=""), recursive=TRUE)}
                
                m <- filerc[[which(names(filerc)==variables[w])]]   # m-1
                dimnames(m)  <- list(x=filerc$X, y=filerc$Y, z=profun, t=fechas_recom)
                
                md<-array(data = NA, dim=c(180,90,dim(m)[3],dim(m)[4]), dimnames = NULL)   
                for (k in 1:dim(m)[4]) {    # one k per date
                    z<-m[,,,k]              # 180 126 30
                    ZETA<-abind(z[91:180,,],z[1:90,,], along=1)   #  180 126 30
                    X<-c(filerc$X[91:180],filerc$X[1:90])   
                    interpoladico<-array(data = NA, dim=c(180,90,dim(ZETA)[3]), dimnames = NULL) 
                    for (j in 1:dim(ZETA)[3]) {    # depth
                        # Interpolate in 2D
                        interm<-ZETA[,,j]          # 180 126
                        range(interm,na.rm=TRUE)
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
                        interpoladico[,,j] <- res  }  # end loop j      #  180 90 7
                    md[,,,k] <- interpoladico    }  # end loop k        #  180 90 7 12
                
                dimnames(md) <- list(x=X, y=colnames(res), z=profun, t=fechas_recom)
                
                # Monthly res
                res <-array(data = NA, dim=c(180,90,30,length(unique(meses_recom))),
                            dimnames = list(x=X, y=colnames(res), z=profun, t=unique(meses_recom)))  
                for (q in unique(meses_recom)){
                    res[,,,q] <- apply(md[,,,which(meses_recom==q)], MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)}
                
                # SEASONAL
                save(res, file=paste(p_dir,"2iop_",variables[w],"/",modelname,"_450.RData",sep="") )
                }  # end loop w variable  
                
                
      #  Adg
        if (file.exists(paste(p_dir,"2iop_abdg_wb", sep=""))==FALSE){
        dir.create(paste(p_dir,"2iop_abdg_wb", sep=""))}
        m <- filerc$acdom_wb + filerc$abpar_wb
        dimnames(m)  <- list(x=filerc$X, y=filerc$Y, z=profun, t=fechas_recom)
        
            md<-array(data = NA, dim=c(180,90,dim(m)[3],dim(m)[4]), dimnames = NULL)   
            for (k in 1:dim(m)[4]) {     # date
                z<-m[,,,k]               # 180 126 30
                ZETA<-abind(z[91:180,,],z[1:90,,], along=1)   #  180 126 30
                X<-c(filerc$X[91:180],filerc$X[1:90])   
                interpoladico<-array(data = NA, dim=c(180,90,dim(ZETA)[3]), dimnames = NULL) 
                for (j in 1:dim(ZETA)[3]) {    # depth
                    # Interpolate in 2D
                    interm<-ZETA[,,j]          # 180 126
                    range(interm,na.rm=TRUE)
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
                    
                    # Anadir bloque NAS artico     
                    cajasur<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                    colnames(cajasur)<-seq(-89,-81,by=2)
                    cajanorte<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
                    colnames(cajanorte)<-seq(81,89,by=2)
                    res<-cbind(cajasur, ZETA_NEW, cajanorte)   
                    rownames(res)<-rownames(ZETA_NEW)          
                    interpoladico[,,j] <- res  }  # end loop j    #  180 90 7
                md[,,,k] <- interpoladico    }  # end loop k    #  180 90 7 12
            
            dimnames(md) <- list(x=X, y=colnames(res), z=profun, t=fechas_recom)
            
            # Monthly res
            res <-array(data = NA, dim=c(180,90,30,length(unique(meses_recom))),
                        dimnames = list(x=X, y=colnames(res), z=profun, t=unique(meses_recom)))  
            for (q in unique(meses_recom)){
                res[,,,q] <- apply(md[,,,which(meses_recom==q)], MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)}
            
            # SEASONAL
            save(res, file=paste(p_dir,"2iop_abdg_wb/",modelname,".RData",sep="") )
      
    #  App
            if (file.exists(paste(p_dir,"2iop_abpp_wb", sep=""))==FALSE){
            dir.create(paste(p_dir,"2iop_abpp_wb", sep=""))}
            m <- filerc$abphy_wb + filerc$abpar_wb
            dimnames(m)  <- list(x=filerc$X, y=filerc$Y, z=profun, t=fechas_recom)
            
            md<-array(data = NA, dim=c(180,90,dim(m)[3],dim(m)[4]), dimnames = NULL)   
            for (k in 1:dim(m)[4]) {     # date
              z<-m[,,,k]                 # 180 126 30
              ZETA<-abind(z[91:180,,],z[1:90,,], along=1)   #  180 126 30
              X<-c(filerc$X[91:180],filerc$X[1:90])   
              interpoladico<-array(data = NA, dim=c(180,90,dim(ZETA)[3]), dimnames = NULL) 
              for (j in 1:dim(ZETA)[3]) {  # depth
                # Interpolate in 2D
                interm<-ZETA[,,j]          # 180 126
                range(interm,na.rm=TRUE)
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
                interpoladico[,,j] <- res  }  # end loop j    #  180 90 7
              md[,,,k] <- interpoladico    }  # end loop k    #  180 90 7 12

            dimnames(md) <- list(x=X, y=colnames(res), z=profun, t=fechas_recom)
            
            # Monthly res
            res <-array(data = NA, dim=c(180,90,30,length(unique(meses_recom))),
                        dimnames = list(x=X, y=colnames(res), z=profun, t=unique(meses_recom)))  
            for (q in unique(meses_recom)){
              res[,,,q] <- apply(md[,,,which(meses_recom==q)], MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)}
            
            # SEASONAL
            save(res, file=paste(p_dir,"2iop_abpp_wb/",modelname,".RData",sep="") )
            
    }  # end loop g modelname
#######################    
