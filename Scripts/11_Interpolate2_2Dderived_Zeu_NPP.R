if (Sys.info()['sysname']=="Windows") {OS<-"C:/"} else {OS<-paste("/Users/",Sys.info()['user'],"/", sep="")}
source(paste(OS,"Programas/Funciones_varios/functions.r",sep=""))
library(RNetCDF)
library(akima)
library(abind)
library(plot3D)
library(Hmisc)
library(plotrix)

#########################################
## 2D derived variables from model output
#########################################
# Zeu from broadband PAR (for spectral Zeu, it is necessary one run per waveband)
# NPP and EXP (total annual)

o_dir<-paste(OS,"Datos/Ser_Stan/global_14_aphyt/Res_model/",sep="")        
p_dir<-paste(OS,"Datos/Ser_Stan/global_14_aphyt/Res_model/interpolated/",sep="")
# Input: 
# Output:     

#path_figures<-paste(OS,"Datos/Res_C20_radtrans/radtrans_v5_marshall/",sep="")

    # List of simulations
    experimentos<-read.csv(paste(OS,"Datos/Ser_Stan/global_14_aphyt/run_log_marshall_PPC_PS_SA_G100.csv",sep=""), sep=",")
    modelnames<-experimentos$name[c(21,9)]
    ruta<-o_dir
    
    ##  experimentos$name[c(30,31)] are identical to experimentos$name[c(21,9)], just run in another machine (G100)

      
      
###########################
## NON-SPECTRAL PAR and Zeu
###########################  

      for (w in 1:length(modelnames)){  
          modelname<-modelnames[w]
          o_dir <- ruta[w]
          
          # Depth
          archivo<-"/grid.nc"
          filename <- paste(o_dir,modelname,archivo, sep="")
          filen <- open.nc(filename)
          filedepth <- read.nc(filen)       
          profun<-filedepth$Z
          
          # Total Par surface wm2
          archivo <-"/diags2D.nc"
          filename <- paste(o_dir,modelname,archivo, sep="")
          filenc <- open.nc(filename)
          filerc <- read.nc(filenc)
          longitude<-filerc$X
          latitude<-filerc$Y
          fechas_recom_real  <-(((filerc$T/60)/60)/24)
          fechas_recom       <-as.numeric(fechas_recom_real)-((4*360))
          meses_recom        <-ceiling(fechas_recom/30)    
          #names(filerc)
          m<-filerc$PARSURF
          m[m<=1e-8] <- NA
          dimnames(m) <- list(x=longitude,y=latitude,t=fechas_recom)
          # Monthly res
          res <-array(data=NA, dim=c(dim(m)[1:2],length(unique(meses_recom))))
          for (q in unique(meses_recom)){
              res[,,q] <- apply(m[,,which(meses_recom==q)], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)}
          TPARS <- res
          dimnames(TPARS) <- list(x=longitude,y=latitude, t=unique(meses_recom))            
          
          # Total PAR 3D   # W m-2
          archivo <-"/diags3D10daily.nc"
          #archivo <-"/MAREDiags3D10daily.nc"
          filename <- paste(o_dir,modelname,archivo, sep="")
          filenc <- open.nc(filename)
          filerc <- read.nc(filenc)
          #names(filerc)
          m<-filerc$par3d 
          dimnames(m) <- list(x=longitude,y=latitude,z=profun,t=fechas_recom)
          # Monthly res
          res <-array(data=NA, dim=c(dim(m)[1:3],length(unique(meses_recom))))
          for (q in unique(meses_recom)){
              res[,,,q] <- apply(m[,,,which(meses_recom==q)], MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)}
          TPARD <- res
          dimnames(TPARD) <- list(x=longitude,y=latitude, z=profun, t=unique(meses_recom))
          
          save(TPARS,TPARD,file=paste(o_dir,modelname,"/PAR_4D10daily.RData",sep=""))     
              
          ## Euphotic depth TOTAL (m)     
          input <- abind(TPARS, TPARD, along=3)    
          resul <- apply(input, MARGIN=c(1,2,4), FUN=find_eudepth, prof=profun)     
          dimnames(resul) <- list(x=longitude, y=latitude, t=unique(meses_recom))
          save(resul, file=paste(o_dir,modelname,"/Zeu_3D10daily.RData",sep="") )
      
      }  # end modelnames
#########################              



        ## Set of Zeu (broadband) per model
        source("/Users/ealvarez/Programas/Funciones_varios/functions.r")
        library(akima)
        library(abind)
        o_dir<-"/Users/ealvarez/Datos/Ser_Stan/global_14_aphyt/"        
        p_dir<-"/Users/ealvarez/Datos/Ser_Stan/global_14_aphyt/interpolated/"

        experimentos<-read.csv("/Users/ealvarez/Datos/Ser_Stan/global_14_aphyt/run_log_marshall_PPC_PS_SA_G100.csv", sep=",")
        modelnames<-experimentos$name[c(21,9,30,31)]
        ruta<-experimentos$path[c(21,9,30,31)]
            
        nombre <- substring(modelnames,1,100)
        ZEU <- array(NA,dim=c(180,126,12,length(modelnames)))
        for (w in 1:length(modelnames)){  
            modelname<-modelnames[w]
            load(paste(ruta[w],modelname,"/Zeu_3D10daily.RData",sep=""))
                longitude <- dimnames(resul)$x
                latitude <- dimnames(resul)$y
                fechas <- dimnames(resul)$t
            zeu_recom <- resul
            ZEU[,,,w] <- zeu_recom }
        dimnames(ZEU) <- list(x=longitude, y=latitude, t=fechas, w=modelnames)
        save(modelnames,longitude,latitude,ZEU, file=paste(o_dir,"Zeu_3D10daily_final4.RData", sep=""))
        
        

        
        
               
################################
##  NPP and EXP total annual  ##   PgC year-1    
################################    

    o_dir<-"/Users/ealvarez/Datos/Ser_Stan/global_14_aphyt/"        
    p_dir<-"/Users/ealvarez/Datos/Ser_Stan/global_14_aphyt/interpolated/"
    experimentos<-read.csv(paste(OS,"Datos/Ser_Stan/global_14_aphyt/run_log_marshall_PPC_PS_SA_G100.csv",sep=""), sep=",")
    modelnames<-experimentos$name[c(21,9,30,31)]     
    ruta<-experimentos$path[c(21,9,30,31)]
    
    # RESULTS
    NPPtotal <- rep(NA,length=length(modelnames))        # NPP & Export production (Pg year-1)
    EXPtotal  <- rep(NA,length=length(modelnames))

    
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
                areas <- (filedepth$rA)
                
                ### NPP & EXP annual total ###
                archivo <-"/recomDiags2D.nc"
                filename <- paste(o_dir,modelname,archivo, sep="")
                filenc <- open.nc(filename)
                filerc <- read.nc(filenc)
                fechas_recom_real  <-(((filerc$T/60)/60)/24)
                fechas_recom       <-as.numeric(fechas_recom_real)-((4*360)-1)
                meses_recom        <-ceiling(fechas_recom/30) 
                tiempito <- (((filerc$T/60)/60)/24)
                
                # NPP   mmolC m-2 d-1
                NPP<-(filerc$NETPPVIS+filerc$NETPPVID)
                NPP[NPP==0]<-NA
                dimnames(NPP) <- list(x=filerc$X, y=filerc$Y, t=(((filerc$T/60)/60)/24))

                # EXPORT  mmolC m-2 d-1
                if (sum(match(names(filerc),"EXPORTC"),na.rm=T)==1){
                EXP<-(filerc$EXPORTC)
                EXP[EXP==0]<-NA
                dimnames(EXP) <- list(x=filerc$X, y=filerc$Y, t=(((filerc$T/60)/60)/24)) }
                
                AREAS <- array(rep(areas,times=12),dim=dim(NPP),dimnames=list(x=filerc$X, y=filerc$Y, t=tiempito) )  
                NPPtotal[w] <- sum(NPP*AREAS*(360/12), na.rm=TRUE)*1e-18*12     #  mmolC m-2 d-1  to  Pg C year-1
                if (sum(match(names(filerc),"EXPORTC"),na.rm=T)==1){
                EXPtotal[w] <- sum(EXP*AREAS*(360/12), na.rm=TRUE)*1e-18*12 }   #  mmolC m-2 d-1  to  Pg C year-1
                
                close.nc(filenc)
                rm(filerc)             
            } # end loop w modelnames

    data.frame(modelnames,NPPtotal,EXPtotal) 
    save(modelnames,NPPtotal,EXPtotal,CNPPtotal,file=paste(o_dir,"Tannual_NPP_final4.RData",sep=""))     

########################