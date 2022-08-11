source("00_Summary_EDITME.R")
source("functions_misc.R")

################################
###  Chla and PPC concentrations
################################   
# Input: .csv tables of hplc samples, quality controlled (FLAG3 and FLAG5)
# Output:  grids with 2x2 degrees and model depth layers

p_dir<-paste(global_path,"Dat_observations/HPLC/",sep="")
# one for MAREDAT: "pigments_MAREDAT_and_Claustre_quality_controled.csv"
# one for CSIRO: "TOTAL_CSIRO_hplc.csv" 
# other for Bracher : "TOTAL_Astrid6_hplc.csv"
# other for Chase: "TOTAL_Chase_hplc.csv"
# combine them in ALL 

        # define grid in depth : use the model
        experimentos<-read.csv(paste(global_path,"/Res_model/run_log_marshall_PPC_PS_SA_G100.csv",sep=""), sep=",")
        modelnames<-experimentos$name[c(21,9)]
        ruta<-experimentos$path[c(21,9)]
        load(paste(global_path,"/Res_model/interpolated/","1cloro_global/",modelnames[1],".RData",sep="")) 
        prof_field <-(-as.numeric(unlist(dimnames(res)$z)))      
        Zu <- prof_field-(c(10,diff(prof_field))/2)
        Zl <- prof_field+(c(diff(prof_field),500)/2)
        #cbind(prof_field, ((Zu+Zl)/2), Zu, Zl)
        latitud  <- cbind(seq(-90,88,by=2),seq(-89,89,by=2),seq(-88,90,by=2))     
        longitud <- cbind(seq(-180,178,by=2),seq(-179,179,by=2),seq(-178,180,by=2))
        clases_lat  <-seq(-90,90,by=2)
        clases_log  <-seq(-180,180,by=2)
        clases_pro  <-c(0,Zl)
        clases_time <-c(0:12)
        month<-c(1:12)   
    
        
###########################################
###  Grid 2x2 with Astrid, CSIRO and Chase
########################################### 

# Tabla datos Astrid
astrid<-read.csv(paste(p_dir,"TOTAL_Astrid6_hplc.csv", sep=""))    
#names(astrid)
nombres_astrid<-names(astrid)[c(3,9,16,11:13,18,19,24,28,31,37,51,54:56,53,58)]
astrid2<-astrid[,c(3,9,16,11:13,18,19,24,28,31,37,51,54:56,53,58)] 
#dim(astrid2)
colnames(astrid2)
sum(complete.cases(astrid2[,c(9,13)])) 

# Tabla datos CSIRO
csiro<-read.csv(paste(p_dir,"TOTAL_CSIRO_hplc.csv", sep=""))    
#names(csiro)   
nombres_csiro<-names(csiro)[c(3,9:13,15,16,21,25,28,34,48,51:53,50,55)]    
csiro2<-csiro[,c(3,9:13,15,16,21,25,28,34,48,51:53,50,55)]  
#dim(csiro2)
colnames(csiro2)
sum(complete.cases(csiro2[,c(9,13)])) 

# Tabla datos Chase
chase<-read.csv(paste(p_dir,"TOTAL_Chase_hplc.csv", sep=""))
#names(chase)       
nombres_chase<-c(names(chase)[c(1,14,4:7,12,13,16:18,44,43,46:48)],NA,NA)  
chase2<-cbind(chase[,c(1,14,4:7,12,13,16:18,44,43,46:48)],matrix(NA,ncol=2,nrow=nrow(chase)))
colnames(chase2)<-nombres_astrid
#dim(chase2)

#cbind(nombres_astrid,nombres_csiro,nombres_chase)
experimentos<-rbind(astrid2, csiro2)    
#colnames(experimentos)    
    
    ## The three tables together    
    cruise <- experimentos$DATABASE
    lati<-experimentos$lat   # plot(long, lati)
    long<-experimentos$lon
    depth<-as.numeric(as.character(experimentos$Depth.m.)) # CUIDADO: lo lee como factor
    #range(depth, na.rm=TRUE)

    YEAR  <- experimentos$Year    #range(YEAR, na.rm=TRUE)
    times <- experimentos$Month   #sort(unique(times))
    DAY <- experimentos$Day       #sort(unique(DAY))
    datetime <- as.Date(paste(DAY,times,YEAR,sep="/"), format="%d/%m/%Y")  #range(datetime, na.rm=TRUE)    
    JULIAN <- julian(datetime, origin = as.Date("2000-01-01"))             #range(JULIAN, na.rm=TRUE)

    #range(datetime[!is.na(experimentos$MTChla)],na.rm=T)
    
    ## Name aggregated pigments (original) and invividual pigments (review1)
    namePigments <- colnames(experimentos)[9:16]
    posPigments <- c(9:16)     
    unidades    <- rep(1e+3,8)
    #cbind(colnames(experimentos)[posPigments],namePigments,unidades) 

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

    # Fill grid cells for each group of pigments 
    # one per pigment, all in ng L-1
    for (w in c(1:length(namePigments))){
      NAME=namePigments[w] 
      pigment <- experimentos[,posPigments[w]]   # ug L-1  
      pigment[experimentos$MYFLAG3==3] <- NA
      if (NAME=="MPPC"){pigment[experimentos$MYFLAG5==5]<-NA}
         #range(datetime[!is.na(pigment)])
      value   <- c(pigment)*unidades[w]      # change from ug L-1 to ng L-1, to merge later with MAREDAT
      experiment <- data.frame(value, factor1, factor2, factor3, factor4)

                jpeg(file=paste(p_dir,"quality_check/distr_",NAME,"_Bracher_CSIRO_Chase.jpg", sep=""))                     
                par(mfrow=c(2,1))
                hist(log10(value), las=1, breaks=seq(-3,6,by=0.1), main=NAME, xaxt="n")
                axis(1, at=seq(-4,6,by=1), labels=10^seq(-4,6,by=1), las=1)
                res <- aggregate(value~factor1+factor2+factor3+factor4,experiment,
                                 mean, na.rm=TRUE, drop=TRUE, simplify=FALSE,
                                 na.action = na.omit)
                chla_diat <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(niveles4)),
                                   dimnames=list(x=niveles1,y=niveles2,z=niveles3,t=niveles4))
                valor <- unlist(res$value)
    
                for (i in 1:nrow(res)){
                  chla_diat[which(niveles1==res[i,1]),
                            which(niveles2==res[i,2]),
                            which(niveles3==res[i,3]),
                            which(niveles4==res[i,4])] <- valor[i]} # end loop i
                chla_diat[chla_diat==0] <-NA
                hist(log10(chla_diat), las=1, breaks=seq(-3,6,by=0.1), main=NAME, xaxt="n")
                axis(1, at=seq(-4,6,by=1), labels=10^seq(-4,6,by=1), las=1)
                dev.off()       
                
      longitude=seq(-179,179,length=180)
      latitude=seq(-89,89,by=2)
      dimnames(chla_diat) <- list(x=seq(-179,179,length=180),y=seq(-89,89,by=2),z=prof_field,t=c(1:12))
      res <- chla_diat        

      # SEASONAL 
      save(longitude, latitude, prof_field,  month,  res,
           file=paste(p_dir, "Bracher_CSIRO_Chase/", "media_seasonal_", NAME,".Rdata", sep=""))
      
      # ANNUAL
      res[res==0] <- NA
      res <- apply(res, MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)
      save(longitude, latitude, prof_field, res,
           file=paste(p_dir,"Bracher_CSIRO_Chase/", "media_annual_", NAME,".Rdata", sep=""))
      
      # jpeg(file=paste(p_dir,"Bracher_CSIRO_Chase/quality_maps/map_",NAME,".jpg", sep=""), width = 400, height = 300, pointsize = 10, quality=100)                     
      # par(mfrow=c(1,1))
      # latit <- latitude
      # longi <- longitude
      # res2 <- apply(res, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)  
      # image(x=longi,y=latit,res2, breaks=c(0,0,10000),col=c("transparent","red3"),ylab="", xlab="", las=1, main=NAME)
      # dev.off()         
      
    } # end loop w (pigment)
###########






################################
### Grid 2x2 pigments MAREDAT  #
################################    


# Tabla datos MAREDAT
experimentos<-read.csv(paste(p_dir,"pigments_MAREDAT_and_Claustre_quality_controled.csv", sep=""))
#names(experimentos)
        #dim(experimentos)
        colnames(experimentos)
        sum(complete.cases(experimentos[,c(19,62)]))         
        
        cruise <- experimentos$Cruise_des
        lati<-experimentos$lat   # plot(long, lati)
        long<-experimentos$lon
        depth<-as.numeric(as.character(experimentos$Depth.m.)) # CUIDADO: lo lee como factor
        YEAR  <- experimentos$Year       # range(YEAR, na.rm=TRUE)
        times <- experimentos$Month
        DAY <- experimentos$Day        
        datetime <- as.Date(paste(DAY,times,YEAR, sep="/"), format="%d/%m/%Y")         
        JULIAN <- julian(datetime, origin = as.Date("1988-01-01")) # range(JULIAN, na.rm=TRUE)  
        #range(datetime[!is.na(experimentos$MTChla)],na.rm=T)
        
      ## Name aggregated pigments
      namePigments <-  c("MTChla","MTChlb","MTChlc","MPSC","MPPC")         
      namePigments2 <- c("MTChla","MTChlb","MTChlc","MPSC","MPPC")
      posPigments <- c(19,24,33,43,62)         
      #cbind(colnames(experimentos)[posPigments],namePigments2,namePigments) 

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


      # Fill grid cells for each pigment group     
      for (w in c(1:length(namePigments))){
        NAME=namePigments[w]      
        pigment <- experimentos[,posPigments[w]]   # ng L-1  
        pigment[experimentos$MYFLAG3==3] <- NA
        if (NAME=="MPPC"){pigment[experimentos$MYFLAG5==5]<-NA}
        #range(datetime[!is.na(pigment)])
        
        value   <- c(pigment)             # already in ng L-1  !!!!!!
        experiment <- data.frame(value, factor1, factor2, factor3, factor4)

                #jpeg(file=paste(p_dir,"quality_check/distr_",NAME,"_MAREDAT.jpg", sep="")) 
                #par(mfrow=c(2,1))
                #hist(log10(value), las=1, breaks=seq(-4,6,by=0.1), main=NAME, xaxt="n")
                #axis(1, at=seq(-4,6,by=1), labels=10^seq(-4,6,by=1), las=1)
                
                res <- aggregate(value~factor1+factor2+factor3+factor4,experiment,
                                 mean, na.rm=TRUE, drop=TRUE, simplify=FALSE,
                                 na.action = na.omit)

                chla_diat <- array(NA, dim=c(length(niveles1),length(niveles2),length(niveles3),length(niveles4)),
                                   dimnames=list(x=niveles1,y=niveles2,z=niveles3,t=niveles4))
                valor <- unlist(res$value)

                for (i in 1:nrow(res)){
                  chla_diat[which(niveles1==res[i,1]),
                            which(niveles2==res[i,2]),
                            which(niveles3==res[i,3]),
                            which(niveles4==res[i,4])] <- valor[i]} # end loop i
                chla_diat[chla_diat==0] <-NA
                #hist(log10(chla_diat), las=1, breaks=seq(-4,6,by=0.1), main=NAME, xaxt="n")  
                #axis(1, at=seq(-4,6,by=1), labels=10^seq(-4,6,by=1), las=1)
                #dev.off()       
                
        longitude=seq(-179,179,length=180)
        latitude=seq(-89,89,by=2)
        dimnames(chla_diat) <- list(x=seq(-179,179,length=180),y=seq(-89,89,by=2),z=prof_field, t=c(1:12))
        res <- chla_diat        
    
        # SEASONAL 
        save(longitude, latitude, prof_field, month, res,file=paste(p_dir,"MAREDAT/","media_seasonal_", NAME,".Rdata", sep=""))
        
        # ANNUAL
        res[res==0] <- NA
        res <- apply(res, MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)
        save(longitude, latitude, prof_field, res,file=paste(p_dir,"MAREDAT/","media_annual_", NAME,".Rdata", sep=""))
  
                      # jpeg(file=paste(p_dir,"MAREDAT/quality_maps/map_",NAME,".jpg", sep=""), width = 400, height = 300, pointsize = 10, quality=100)                     
                      # par(mfrow=c(1,1))
                      # latit <- latitude
                      # longi <- longitude
                      # res2 <- apply(res, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)  
                      # image(x=longi,y=latit,res2, breaks=c(0,0,10000),col=c("transparent","red3"), ylab="", xlab="", las=1, main=NAME)
                      # dev.off()         
  
      } # end loop w (pigment)
###########


        
        
     
        
        
#######################################################
###  Grid 2x2 with Astrid_CSIRO_Chase + MAREDAT : ALL #
####################################################### 

## Name pigments (aggregated revised + individual)
namePigments <- c("MTChla","MTChlb","MTChlc","MPSC","MPPC")

for (k in 1:length(namePigments)){
  NAME <-namePigments[k] 
  load(paste(p_dir,"MAREDAT/media_seasonal_", NAME,".RData", sep=""))
  #dim(res)
  maredat <- res
  load(paste(p_dir,"Bracher_CSIRO_Chase/media_seasonal_", NAME,".Rdata", sep=""))
  #dim(res)  
  astrid <- res
  
  res2 <- abind(maredat,astrid, along=0.5)
  #dim(res2)
  dimnames(res2) <- list(h=c(1:2),x=longitude, y=latitude, z=prof_field, t=month)
  res <- apply(res2, MARGIN=c("x","y","z","t"), FUN=mean, na.rm=TRUE)
  
  # SEASONAL 
  save(longitude, latitude, prof_field, month, res,file=paste(p_dir,"media_seasonal_", NAME,"_ALL.Rdata", sep=""))
  
  # ANNUAL
  res[res==0] <- NA
  res <- apply(res, MARGIN=c("x","y","z"), FUN=mean, na.rm=TRUE)
  save(longitude, latitude, prof_field, res,file=paste(p_dir,"media_annual_",  NAME,"_ALL.Rdata", sep=""))
  
  
            # jpeg(file=paste(p_dir,"quality_check/map_",NAME,".jpg", sep=""), width = 400, height = 300, pointsize = 10, quality=100)                     
            # par(mfrow=c(1,1))
            # latit <- latitude
            # longi <- longitude
            # res2 <- apply(res, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)  
            # image(x=longi,y=latit,res2, breaks=c(0,0,10000),col=c("transparent","red3"),ylab="", xlab="", las=1, main=NAME)
            # 
            # maredat2 <- apply(maredat, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)
            # image(x=longi,y=latit,maredat2, breaks=c(0,0,10000),col=c("transparent","blue"),ylab="", xlab="", las=1, add=TRUE)        
            # 
            # astrid2 <- apply(astrid, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)
            # image(x=longi,y=latit,astrid2, breaks=c(0,0,10000),col=c("transparent","green"),ylab="", xlab="", las=1, add=TRUE)  
            # 
            # legend(x="topright", legend=c("MAREDAT","Astrid"), pch=22, col=c("blue","green"), pt.bg=c("blue","green"))
            # dev.off()         
            
      } # end loop k
###############    





        
        
        
###############################
### Grids PPC:TChla in situ ###
###############################   

## Name pigments
namePigments <- c("MTChla", "MTChlb", "MTChlc", "MPSC", "MPPC")


### In situ MAREDAT PPC
for (k in 1:length(namePigments)){
  fileguardar1 <- paste(p_dir,"MAREDAT/media_annual_",namePigments[k],".RData", sep="")        
  load(fileguardar1)        
  res[is.na(res)] <- 0
  assign(namePigments[k],res)}

          clorofila <- MTChla
          pigmentos <- MPPC
          pigmentos[clorofila==0] <- NA 
          clorofila[clorofila==0] <- NA 
          PPCinsitu_MAREDAT <- pigmentos/clorofila
          #dim(PPCinsitu_MAREDAT)

for (k in 1:length(namePigments)){
  fileguardar1 <- paste(p_dir, "MAREDAT/media_seasonal_",namePigments[k],".RData", sep="")        
  load(fileguardar1)        
  res[is.na(res)] <- 0
  assign(namePigments[k],res)}

          clorofila <- MTChla
          pigmentos <- MPPC
          pigmentos[clorofila==0] <- NA 
          clorofila[clorofila==0] <- NA 
          PPCinsitu_MAREDAT_s <- pigmentos/clorofila
          #dim(PPCinsitu_MAREDAT_s)                    

### In situ NEW            
for (k in 1:length(namePigments)){
  fileguardar1 <- paste(p_dir,"Bracher_CSIRO_Chase/media_annual_", namePigments[k],".RData", sep="")        
  load(fileguardar1)        
  res[is.na(res)] <- 0
  assign(namePigments[k],res)}

          clorofila <- MTChla
          pigmentos <- MPPC
          pigmentos[clorofila==0] <- NA 
          clorofila[clorofila==0] <- NA 
          PPCinsitu_NEW <- pigmentos/clorofila
          #dim(PPCinsitu_NEW)

for (k in 1:length(namePigments)){
  fileguardar1 <- paste(p_dir,"Bracher_CSIRO_Chase/media_seasonal_", namePigments[k],".RData", sep="")        
  load(fileguardar1)        
  res[is.na(res)] <- 0
  assign(namePigments[k],res)}

          clorofila <- MTChla
          pigmentos <- MPPC
          pigmentos[clorofila==0] <- NA 
          clorofila[clorofila==0] <- NA 
          PPCinsitu_NEW_s <- pigmentos/clorofila
          #dim(PPCinsitu_NEW_s)    
          

### In situ ALL           
for (k in 1:length(namePigments)){
  fileguardar1 <- paste(p_dir, "media_annual_", namePigments[k],"_ALL.RData", sep="")        
  load(fileguardar1)        
  res[is.na(res)] <- 0
  assign(namePigments[k],res)}

          clorofila <- MTChla
          pigmentos <- MPPC
          pigmentos[clorofila==0] <- NA 
          clorofila[clorofila==0] <- NA 
          PPCinsitu_ALL <- pigmentos/clorofila
          #dim(PPCinsitu_ALL)
          

for (k in 1:length(namePigments)){
  fileguardar1 <- paste(p_dir,"media_seasonal_",namePigments[k],"_ALL.RData", sep="")        
  load(fileguardar1)        
  res[is.na(res)] <- 0
  assign(namePigments[k],res)}

          clorofila <- MTChla
          pigmentos <- MPPC
          pigmentos[clorofila==0] <- NA 
          clorofila[clorofila==0] <- NA 
          PPCinsitu_ALL_s <- pigmentos/clorofila
          #dim(PPCinsitu_ALL_s)               

save(PPCinsitu_MAREDAT,    PPCinsitu_MAREDAT_s,
     PPCinsitu_NEW,        PPCinsitu_NEW_s,
     PPCinsitu_ALL,        PPCinsitu_ALL_s,
     file=paste(p_dir,"grids_pigments/","PPCinsitu_new3.RData", sep=""))                          

###############################    









###########################
### Grids TChla in situ ### 
###########################    
    ## Name pigments
    namePigments <- c("MTChla")
   
##### In situ MAREDAT TChla
    # Annual
    for (k in 1:length(namePigments)){
      fileguardar1 <- paste(p_dir,"MAREDAT/media_annual_",namePigments[k],".RData", sep="")        
      load(fileguardar1)        
      res[is.na(res)] <- 0
      assign(namePigments[k],res)}
            clorofila <- MTChla
            clorofila[clorofila==0] <- NA 
            TChla_MAREDAT <- clorofila
            #dim(TChla_MAREDAT)
            
    # Seasonal
    for (k in 1:length(namePigments)){
      fileguardar1 <- paste(p_dir, "MAREDAT/media_seasonal_",namePigments[k],".RData", sep="")        
      load(fileguardar1)        
      res[is.na(res)] <- 0
      assign(namePigments[k],res)}
            clorofila <- MTChla
            clorofila[clorofila==0] <- NA 
            TChla_MAREDAT_s <- clorofila
            #dim(TChla_MAREDAT_s)                    
    
    
##### In situ NEW
    # Annual
    for (k in 1:length(namePigments)){
      fileguardar1 <- paste(p_dir,"Bracher_CSIRO_Chase/media_annual_", namePigments[k],".RData", sep="")        
      load(fileguardar1)        
      res[is.na(res)] <- 0
      assign(namePigments[k],res)}
            clorofila <- MTChla
            clorofila[clorofila==0] <- NA 
            TChla_NEW <- clorofila
            #dim(TChla_NEW)

    # Seasonal
    for (k in 1:length(namePigments)){
      fileguardar1 <- paste(p_dir,"Bracher_CSIRO_Chase/media_seasonal_", namePigments[k],".RData", sep="")        
      load(fileguardar1)        
      res[is.na(res)] <- 0
      assign(namePigments[k],res)}
            clorofila <- MTChla
            clorofila[clorofila==0] <- NA 
            TChla_NEW_s <- clorofila
            #dim(TChla_NEW_s)    
        
 
##### In situ ALL
    # Annual
    for (k in 1:length(namePigments)){
      fileguardar1 <- paste(p_dir, "media_annual_", namePigments[k],"_ALL.RData", sep="")        
      load(fileguardar1)        
      res[is.na(res)] <- 0
      assign(namePigments[k],res)}
            clorofila <- MTChla
            clorofila[clorofila==0] <- NA 
            TChla_ALL <- clorofila
            #dim(TChla_ALL)
                     
    # Seasonal
    for (k in 1:length(namePigments)){
      fileguardar1 <- paste(p_dir,"media_seasonal_",namePigments[k],"_ALL.RData", sep="")        
      load(fileguardar1)        
      res[is.na(res)] <- 0
      assign(namePigments[k],res)}
            clorofila <- MTChla
            clorofila[clorofila==0] <- NA 
            TChla_ALL_s <- clorofila
            #dim(TChla_ALL_s)               

    save(TChla_MAREDAT,    TChla_MAREDAT_s,
         TChla_NEW,        TChla_NEW_s,
         TChla_ALL,        TChla_ALL_s,
         file=paste(p_dir,"grids_pigments/","TChla_insitu_new2.RData",sep=""))                           
    
###########################    
    
