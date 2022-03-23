source("00_Summary_EDITME.R")

# OASIM Input Light Below Ocean Surface
# Input: original NetCDF files obtained from /gmaoftp/NOBM/monthly/eds
# User-defined: initial/final wavebands to extract, and initial/final years to average
# Output: one binary file with one climatological year for each waveband
## The format can be checked in the example binary files (loc1_oasim) in MITgcm_contrib/darwin2/verification/radtrans_1d
## file.info(paste(pdir,"loc1_oasim_edp_below.bin", sep=""))
## scan(paste(pdir,"loc1_oasim_edp_below.bin", sep=""), raw())


##########    
### OASIM files per lamda in one degree boxes
##########

      rdir <-paste(global_path,"/Dat_WB_OASIM/originals/",sep="")
      nx <- as.double(400) # first lambda
      ny <- as.double(700) # last lambda
      ano1 <- as.character(as.double(2006)) # first year
      ano2 <- as.character(as.double(2008)) # last year 
      sdir <-paste(global_path,"/Dat_WB_OASIM/binary/",sep="")

  # YEARS
  anos<-seq(as.numeric(ano1), as.numeric(ano2), by=1)
  
  # MONTHS
  meses <- sprintf("%02d.nc", seq(1,12)) 
  #archivos <- sprintf(paste("edsEd",ano1,"%02d.nc", sep=""), seq(1,12))
  longitud_old <- seq(0.5,359.5, by=1)   #360
  latitud_old  <- seq(-89.5,89.5,by=1)   #180
  
  # LAMBDAS
  lamda_old <- c(250,325,350,375,400,425,450,475,500,525,550,575,600,625,650,675,700,
                 725,775,850,950,1050,1150,1250,1350,1450,1550,1650,1750,1900,2200,2900,3700)
  lamda_new <-lamda_old[lamda_old>=nx & lamda_old<=ny]
  o <- match(lamda_new,lamda_old)
  savearchivoEd <-paste("Ed_",lamda_new, "_360x180x12", sep="")
  savearchivoEs <-paste("Es_",lamda_new, "_360x180x12", sep="")
  #data.frame(archivos, savearchivo)    


      # loop landa
      for (i in c(1:length(savearchivoEd))){    
        lambda <- o[i]
        finalEd <- array(NA,dim=c(360,180,12))
        finalEs <- array(NA,dim=c(360,180,12))
        filenameSaveEd <- paste(sdir,savearchivoEd[i], sep="")   
        filenameSaveEs <- paste(sdir,savearchivoEs[i], sep="")
        
        # loop mes
        for (j in c(1:length(meses))){
          archivosEd<-intersect(list.files(rdir,pattern="edsEd"),list.files(rdir,pattern=meses[j]))
          archivosEs<-intersect(list.files(rdir,pattern="edsEs"),list.files(rdir,pattern=meses[j]))
          
          # loop year
          for (h in c(1:length(anos))){        
            archivoEd <- archivosEd[h]
            filename <- paste(rdir,archivoEd, sep="")
            filenc <- open.nc(filename)
            filerc <- read.nc(filenc)
            tmp <- filerc$eds_Ed
            tmp[tmp>9e9] <- NA
            tmp <- tmp[,,lambda]
            ZETA<-abind(tmp[181:360,],tmp[1:180,], along=1)   #  180 126 30
            if (h==1){ med<-array(data = NA, dim=c(dim(ZETA),length(anos)),
                                  dimnames=list(x=longitud_old,y=latitud_old,t=anos))
            med[,,h]<-ZETA } else { med[,,h]<-ZETA }
            
            archivoEs <- archivosEs[h]
            filename <- paste(rdir,archivoEs, sep="")
            filenc <- open.nc(filename)
            filerc <- read.nc(filenc)
            tmp <- filerc$eds_Es
            tmp[tmp>9e9] <- NA
            tmp <- tmp[,,lambda]
            ZETA<-abind(tmp[181:360,],tmp[1:180,], along=1)   #  180 126 30
            if (h==1){ mes<-array(data = NA, dim=c(dim(ZETA),length(anos)),
                                  dimnames=list(x=longitud_old,y=latitud_old,t=anos))
            mes[,,h]<-ZETA } else { mes[,,h]<-ZETA }
          } # end loop h anos
          
          ZEd<-apply(med, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)
          finalEd[,,j] <- ZEd
          ZEs<-apply(mes, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)
          finalEs[,,j] <- ZEs
        } # end loop j month
        
        dimnames(finalEd) <- list(x=longitud_old, y=latitud_old, t=c(1:12))    
        dimnames(finalEs) <- list(x=longitud_old, y=latitud_old, t=c(1:12))
        
        ### Save as binary Ed
        finalEd[is.na(finalEd)] <- 0
        finalEd[is.nan(finalEd)] <- 0
        filename <- file(paste(sdir,savearchivoEd[i],"_32b.bin", sep=""), "wb")
        objeto <- c(finalEd)
        # size==4 # 32bit, single precision
        writeBin(object=objeto, con=filename, size = 4, endian = "big") #, useBytes=TRUE)  
        close(filename)       
        
        ### Save as binary Es
        finalEs[is.na(finalEs)] <- 0
        finalEs[is.nan(finalEs)] <- 0
        filename <- file(paste(sdir,savearchivoEs[i],"_32b.bin", sep=""), "wb")
        objeto <- c(finalEs)
        # size==4        # 32bit, single precision   # 
        writeBin(object=objeto, con=filename, size = 4, endian = "big") #, useBytes=TRUE)  
        close(filename)
        
      } # end loop i lambda
      