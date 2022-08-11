source("Scripts/00_Summary_EDITME.R")
source("Scripts/functions_misc.R")

###############################################################
### COMPARE SPECTRAL SHAPES Rrs MODEL vs SATELLITE/OBSERVATIONS
###############################################################

  o_dir<-paste(global_path,"Res_model/",sep="")        
  p_dir<-paste(global_path,"Res_model/interpolated/",sep="")
  path_figures<-paste(global_path,"Figures/",sep="")
  # List of simulations
  experimentos<-read.csv(paste(o_dir,"run_log_marshall_PPC_PS_SA_G100.csv",sep=""), sep=",")
  s_dir <- paste(global_path,"Dat_observations/satellite/",sep="")
  is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
    modelnames<-experimentos$name[c(21,9)]
    ruta<-o_dir
    nombre <- c(expression(paste({a^{"*"}}[PH],"(", lambda, ") constant",sep="")),
                expression(paste({a^{"*"}}[PH],"(", lambda, ") variable",sep="")))

      # Colors
      mycol1<-usecol(pal_bordeaux,n=4)
      mycol2<-usecol(pal_petrol,n=4)
      colores_puntos<-c(mycol1[3],mycol2[3])
      colores<-c(mycol1[1],mycol2[1])
      colores_obs<-c("grey50","grey30")

       
      
##############################
### Spectral Shape REFLECTANCE
##############################

        # MODEL
        modelname<-modelnames[1]        
        load(paste(p_dir,"1reflectance/",modelname,".RData",sep=""))
        lamda_new <- as.numeric(dimnames(res)$l)
        model_mask <- res[,,6]
        model_mask[model_mask>1.0e-9]  <- 1
        model_mask[model_mask<=1.0e-9] <- 0
        
        # SAT
        landas_sat <- c("Rrs_412","Rrs_443","Rrs_490","Rrs_510","Rrs_555","Rrs_670")
        landa_sat <- c(412,443,490,510,555,670)
        load(paste(s_dir,"climatologies_2012_2018/media_seasonal_Rrs_2deg_OCCCI_2012_2018.RData",sep=""))
        res <- med
            # ANNUAL 
            res <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)
            for (w in 1:length(landa_sat)){res[,,w][model_mask==0] <- NA}
            sat <- apply(res, MARGIN="l", FUN=mean, na.rm=TRUE)
            #sat_SO <- apply(res[,1:25,], MARGIN="l", FUN=mean, na.rm=TRUE)
            #sat_AO <- apply(res[,66:90,], MARGIN="l", FUN=mean, na.rm=TRUE)
            #sat_TR <- apply(res[,26:65,], MARGIN="l", FUN=mean, na.rm=TRUE)

        # In situ
        load(paste(is_dir,"/media_annual_Rrs.Rdata",sep=""))
        lambda_obs <- lambda
            # ANNUAL
            for (w in 1:length(lambda_obs)){matriz_annual[,,w][model_mask==0] <- NA}
            obs <- apply(matriz_annual, MARGIN=c("l"), FUN=mean, na.rm=TRUE)
            #obs_SO <- apply(matriz_annual[,1:25,],  MARGIN=c("l"), FUN=mean, na.rm=TRUE)
            #obs_AO <- apply(matriz_annual[,66:90,], MARGIN=c("l"), FUN=mean, na.rm=TRUE)
            #obs_TR <- apply(matriz_annual[,26:65,], MARGIN=c("l"), FUN=mean, na.rm=TRUE)

            

############    
## FIGURE 10 RRS(lambda) by TChla concentration
############    
            
png(file=paste(path_figures,"Figure10_Reflectance_spectra_lambda.png",sep=""), width = 1200, height = 360, units = "px", pointsize = 18,bg = "white")
        par(mfrow=c(1,4))
        par(mar=c(4,3,1,1))
        par(oma=c(2,3,2,0))
        layout.show(n=4)      
        labtext=1.4
        lettersize=2
        textinner=1.2
        cexaxis=1.2        
        
#### All pixels 
###############        
       #  In situ             
       plot(landa_sat, sat, las=1, type="b", col=colores_obs[1], bg=colores_obs[1], pch=21, cex=1.2, lwd=2,
             ylim=c(0,0.012), xlim=c(400,700), ylab=expression(paste(R[RS]," (",sr^-1,")")), xlab=expression(lambda),
             mgp=c(3,1,0), cex.lab=labtext, yaxs="i", cex.axis=cexaxis)
       points(lambda_obs, obs, las=1, type="b", col=colores_obs[2], bg=colores_obs[2], pch=21, lwd=2, cex=1.2)
              
       legend(x=540, y=0.010, legend=c("Satellite","Obs. in situ"),col = c(colores_obs[1],colores_obs[2]), lwd=2, pch=21,
                       pt.bg = c(colores_obs[1],colores_obs[2]), pt.cex=c(1.2,1.2),bty="n")
       legend(x=520, y=0.008, legend=nombre, lwd=2,col=colores, bty="n")        
       text(x=400,y=0.0111,labels="Global means", cex=textinner, pos=4)
       mtext(3, at=390, adj=0, line=0.5, text="a)", cex=lettersize, font=1)
              
       # MODEL
       for (h in c(1:length(modelnames))){
               modelname<-modelnames[h]       
               load(paste(p_dir,"1reflectance/",modelname,".RData",sep=""))
               lamda_new <- as.numeric(dimnames(res)$l)
               for (w in 1:length(lamda_new)){res[,,w][model_mask==0] <- NA}
               mod <- apply(res, MARGIN=c("l"), FUN=mean, na.rm=TRUE)
               points(lamda_new, mod, las=1, type="l", col=colores[h], bg=colores[h], pch=21, lwd=2,ylim=c(0,0.01), xlim=c(400,700),lty=1)
               }
############          
        
#### Only pixels with 1 +- 0.1 mg Chla m-3       
#############        
        # In situ
        load(paste(is_dir,"/media_annual_Rrs.Rdata",sep=""))
        lambda_obs <- lambda
          for (w in 1:length(lambda_obs)){matriz_annual[,,w][model_mask==0] <- NA}
          for (w in 1:length(lambda_obs)){matriz_annual[,,w][model_mask==0] <- NA}
          load(paste(is_dir,"/media_annual_Chla.Rdata",sep=""))
          cloro <- res[,,1]
          for (w in 1:length(lambda_obs)){matriz_annual[,,w][cloro<=0.9 | cloro>=1.1] <- NA}
          obs <- apply(matriz_annual, MARGIN=c("l"), FUN=mean, na.rm=TRUE)
        plot(lambda_obs, obs, las=1, type="b", col=colores_obs[2], bg=colores_obs[2], pch=21, lwd=2, cex=1.2,
             ylim=c(0,0.012), xlim=c(400,700), ylab=expression(paste(R[RS]," (",sr^-1,")")), xlab=expression(lambda),
             mgp=c(3,1,0), cex.lab=labtext, yaxs="i", cex.axis=cexaxis)
        text(x=400,y=0.0111, labels=expression(paste("Pixels with 1±0.1 mg Chla ",m^-3)), cex=textinner, pos=4)
        mtext(3, at=390, adj=0, line=0.5, text="b)", cex=lettersize, font=1)
        
        # MODEL
        for (h in c(1:length(modelnames))){
            modelname<-modelnames[h]       
            load(paste(p_dir,"2cloro_lin/",modelname,".RData",sep=""))
            cloro <- apply(res, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)
            load(paste(p_dir,"1reflectance/",modelname,".RData",sep=""))
            lamda_new <- as.numeric(dimnames(res)$l)
            for (w in 1:length(lamda_new)){res[,,w][model_mask==0] <- NA}         
            for (w in 1:length(lamda_new)){res[,,w][cloro<=0.9 | cloro>=1.1] <- NA} 
            mod <- apply(res, MARGIN=c("l"), FUN=mean, na.rm=TRUE)
            points(lamda_new, mod, las=1, type="l", col=colores[h], pch=21, lwd=2)
        }
#############
        
#### Only pixels with 0.5 +- 0.05 mg Chla m-3       
#############        
        # In situ
        load(paste(is_dir,"/media_annual_Rrs.Rdata",sep=""))
        lambda_obs <- lambda
          for (w in 1:length(lambda_obs)){matriz_annual[,,w][model_mask==0] <- NA}
          for (w in 1:length(lambda_obs)){matriz_annual[,,w][model_mask==0] <- NA}
          load(paste(is_dir,"/media_annual_Chla.Rdata",sep=""))
          cloro <- res[,,1]
          for (w in 1:length(lambda_obs)){matriz_annual[,,w][cloro<=0.45 | cloro>=0.55] <- NA}
          obs <- apply(matriz_annual, MARGIN=c("l"), FUN=mean, na.rm=TRUE)
        plot(lambda_obs, obs, las=1, type="b", col=colores_obs[2], bg=colores_obs[2], pch=21, lwd=2,cex=1.2,
             ylim=c(0,0.012), xlim=c(400,700),ylab=expression(paste(R[RS]," (",sr^-1,")")), xlab=expression(lambda),
             mgp=c(3,1,0), cex.lab=labtext, yaxs="i", cex.axis=cexaxis)
        text(x=400,y=0.0111,labels=expression(paste("Pixels with 0.5±0.05 mg Chla ",m^-3)), cex=textinner, pos=4)
        mtext(3, at=390, adj=0, line=0.5, text="c)", cex=lettersize, font=1)
        
        # MODEL
        for (h in c(1:length(modelnames))){
          modelname<-modelnames[h]       
          load(paste(p_dir,"2cloro_lin/",modelname,".RData",sep=""))
          cloro <- apply(res, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)
          load(paste(p_dir,"1reflectance/",modelname,".RData",sep=""))
          lamda_new <- as.numeric(dimnames(res)$l)
          for (w in 1:length(lamda_new)){res[,,w][model_mask==0] <- NA}         
          for (w in 1:length(lamda_new)){res[,,w][cloro<=0.45 | cloro>=0.55] <- NA} 
          mod <- apply(res, MARGIN=c("l"), FUN=mean, na.rm=TRUE)
          points(lamda_new, mod, las=1, type="l", col=colores[h], pch=21, lwd=2)
        }        
#############
        
#### Only pixels with 0.1 +- 0.01 mg Chla m-3       
#############        
        # In situ
        load(paste(is_dir,"/media_annual_Rrs.Rdata",sep=""))
        lambda_obs <- lambda
          for (w in 1:length(lambda_obs)){matriz_annual[,,w][model_mask==0] <- NA}
          for (w in 1:length(lambda_obs)){matriz_annual[,,w][model_mask==0] <- NA}
          load(paste(is_dir,"/media_annual_Chla.Rdata",sep=""))
          cloro <- res[,,1]
          for (w in 1:length(lambda_obs)){matriz_annual[,,w][cloro<=0.09 | cloro>=0.11] <- NA}
          obs <- apply(matriz_annual, MARGIN=c("l"), FUN=mean, na.rm=TRUE)
        plot(lambda_obs, obs, las=1, type="b", col=colores_obs[2], bg=colores_obs[2], pch=21, lwd=2,cex=1.2,
             ylim=c(0,0.012), xlim=c(400,700),ylab=expression(paste(R[RS]," (",sr^-1,")")), xlab=expression(lambda),
             mgp=c(3,1,0), cex.lab=labtext, yaxs="i", cex.axis=cexaxis)
        text(x=400,y=0.0111,labels=expression(paste("Pixels with 0.1±0.01 mg Chla ",m^-3)), cex=textinner, pos=4)
        mtext(3, at=390, adj=0, line=0.5, text="d)", cex=lettersize, font=1)
        
        # MODEL
        for (h in c(1:length(modelnames))){
          modelname<-modelnames[h]       
          load(paste(p_dir,"2cloro_lin/",modelname,".RData",sep=""))
          cloro <- apply(res, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)
          load(paste(p_dir,"1reflectance/",modelname,".RData",sep=""))
          lamda_new <- as.numeric(dimnames(res)$l)
          for (w in 1:length(lamda_new)){res[,,w][model_mask==0] <- NA}         
          for (w in 1:length(lamda_new)){res[,,w][cloro<=0.09 | cloro>=0.11] <- NA} 
          mod <- apply(res, MARGIN=c("l"), FUN=mean, na.rm=TRUE)
          points(lamda_new, mod, las=1, type="l", col=colores[h], pch=21, lwd=2)
        }          
#############
        
  #mtext(side=1,at=0.55,line=0.5,outer=T,text=expression(lambda)) 
  mtext(side=2,at=0.55,line=0.6,outer=T,text=expression(paste(R[RS]," (",sr^-1,")")))   
  
        dev.off()
#########        




  
  
  
#############  
### FIGURE 11: Spectral shape by Biomes for both mod and sat/obs
#############  
 
png(file=paste(path_figures,"Figure11_Reflectance_spectral_shape_biomes.png",sep=""), width=1000, height=660, units="px", pointsize=18,bg="white")    

      layout(matrix(c(1,5,9,12,15,
                      2,6,10,13,16,
                      3,7,11,14,17,
                      4,8,18,18,18),ncol=5,byrow=TRUE))
      par(mar=c(2,4,1,1))
      par(oma=c(4,4,3,3))
      layout.show(n=18)     
      lettersize=2.0
      textinner=1.4
      cexaxis=1.2   
      cexlab=1.2
      cextitle=1.2
        modelname<-modelnames[1]       
        load(paste(p_dir,"1reflectance/",modelname,".RData",sep="")) 
        model_mask <- res[,,6]
        model_mask[model_mask>1.0e-9]  <- 1
        model_mask[model_mask<=1.0e-9] <- 0            
        
        
### Find Biome in each grid point: both sat and model are in 2 degrees
###################
        # Time_Varying_Biomes.nc from (Fay & McKinley 2014)
        load(paste(global_path,"Misc/Biomes_global.RData",sep=""))
        long<-c(seq(0.5,179.5,by=1),seq(-179.5,-0.5,by=1))
        matrix_lon<-matrix(rep(seq(-179,179,by=2),times=90),ncol=90,byrow=F)
        matrix_lat<-matrix(rep(seq(-89,89,by=2),times=180), ncol=90,byrow=T)
        res <- abind(matrix_lon, matrix_lat, along=0.5)
        dimnames(res) <- list(h=c(1:2),x=c(1:180), y=c(1:90))
        find_biome <- function(v){
          xn <- v[1]
          yn <- v[2]
          if (xn==0 & yn==0) {res <- NA }else{ 
            minimos <- sort(abs(long-xn),na.last=NA)[1]
            cuales <- match(minimos,abs(long-xn))
            intervalo_lon <- long[cuales]
            cuales_lon <- cuales
            minimos <- sort(abs(lati-yn),na.last=NA)[1]
            cuales <- match(minimos, abs(lati-yn))
            intervalo_lat <- lati[cuales]
            cuales_lat <- cuales
            res <- b[cuales_lon, cuales_lat]}
          res <- round(res,0)
          return(res)}
        res5 <- apply(res, MARGIN=c("x","y"), FUN=find_biome)
        biomes<-res5
###################         
        
        
        
## Spectral shape by Biome
########################## 
        for (q in c(4,6,5,7,11,12,13,14,2,3,15,9,10,16,1,8,17)){
          if (q==4|q==6|q==5|q==8) par(mar=c(0,2,0,0))
          if (q==2|q==3|q==15) par(mar=c(0,2,0,0))
          if (q==11|q==12|q==13|q==14) par(mar=c(0,0,0,2))
          if (q==9|q==10|q==16) par(mar=c(0,0,0,2))
          if (q==1|q==8|q==17) par(mar=c(0,2,0,0))
          
          # SAT
          landas_sat <- c("Rrs_412","Rrs_443","Rrs_490","Rrs_510","Rrs_555","Rrs_670")
          landa_sat <- c(412,443,490,510,555,670)
          load(paste(s_dir,"climatologies_2012_2018/media_seasonal_Rrs_2deg_OCCCI_2012_2018.RData",sep=""))
          res <- med
              # ANNUAL 
              res <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)
              for (w in 1:length(landa_sat)){res[,,w][model_mask==0] <- NA}
              for (w in 1:length(landa_sat)){res[,,w][is.na(biomes)] <- NA}
              for (w in 1:length(landa_sat)){res[,,w][biomes!=q] <- NA}
              sat <- apply(res, MARGIN="l", FUN=mean, na.rm=TRUE)
              
          # In situ
          load(paste(is_dir,"/media_annual_Rrs.Rdata",sep=""))
          lambda_obs <- lambda
              for (w in 1:length(lambda_obs)){matriz_annual[,,w][model_mask==0] <- NA}
              for (w in 1:length(lambda_obs)){matriz_annual[,,w][is.na(biomes)] <- NA}
              for (w in 1:length(lambda_obs)){matriz_annual[,,w][biomes!=q] <- NA}
              obs <- apply(matriz_annual, MARGIN=c("l"), FUN=mean, na.rm=TRUE)
              
          plot(landa_sat, sat, las=1, type="b", col=colores_obs[1], pch=19, lwd=2,ylim=c(0,0.015), xlim=c(400,700),
               main="", col.main="black", ylab="",xlab="", xaxt="n", yaxt="n", cex.axis=cexaxis, cex.lab=cexlab)
            legend(x="topright", legend=biome_labels[q],bty="n", cex=textinner)
            points(lambda_obs, obs, las=1, type="b", col=colores_obs[2], pch=19, lwd=2)

          if (q==7|q==14) {         axis(1,seq(450,650,by=50), cex.axis=cexaxis)
                                    mtext(side=1,at=550,line=3,outer=F,text=expression(lambda))}    #else{axis(1,seq(400,700,by=50), labels=F)}
          if (q==4|q==6|q==5|q==7) {axis(2,seq(0.002,0.014,length=4), las=1, cex.axis=cexaxis)} #else{axis(2,seq(0.002,0.012,length=6), labels=F)}
          if (q==15|q==16) {        axis(1,seq(450,650,by=50), cex.axis=cexaxis)
                                    mtext(side=1,at=550,line=3,outer=F,text=expression(lambda))}
          if (q==2|q==3|q==15) {    axis(2,seq(0.002,0.014,length=4), las=1, cex.axis=cexaxis)} #else{axis(2,seq(0.002,0.012,length=6), labels=F)}
          if (q==17) {              axis(1,seq(450,650,by=50), cex.axis=cexaxis)
                                    mtext(side=1,at=550,line=3,outer=F,text=expression(lambda))}
          if (q==1|q==8|q==17) {    axis(2,seq(0.002,0.014,length=4), las=1, cex.axis=cexaxis)} #else{axis(2,seq(0.002,0.012,length=6), labels=F)}
          
          # MODEL
          for (h in c(1:length(modelnames))){    
            modelname<-modelnames[h]       
            load(paste(p_dir,"1reflectance/",modelname,".RData",sep=""))
            lamda_new <- as.numeric(dimnames(res)$l)
                for (w in 1:length(lamda_new)){res[,,w][model_mask==0] <- NA}
                for (w in 1:length(lamda_new)){res[,,w][is.na(biomes)] <- NA}
                for (w in 1:length(lamda_new)){res[,,w][biomes!=q] <- NA}
                mod <- apply(res,  MARGIN=c("l"), FUN=mean, na.rm=TRUE)
             points(lamda_new, mod, las=1, type="l", col=colores[h], pch=19, lwd=2,ylim=c(0,0.01), xlim=c(400,700),lty=1)
          } # end loop h modelname
        } # end loop q biome 
##########################
  # Labels      
  mtext(3,adj=0,line=0.5,text=c("a)","b)","c)"),at=c(0.00,0.40,0.8035),cex=lettersize, font=1, outer=T)    
  mtext(3,adj=0,line=0.5,text=c("Permanently-stratified biomes ", "Seasonally-stratified biomes", "Ice biomes"),
        at=c(0.05,0.45,0.85), cex=cextitle, font=1, outer=T)    
  mtext(side=2,at=0.5,line=1.8,outer=T,text=expression(paste(R[RS],"(",sr^-1,")")))         

### LEGEND
    par(mar=c(0,6,3,4))        
    plot(1,1,type="n",xaxt="n",yaxt="n",ylab="", xlab="", bty="n")        
    legend(x="bottomleft", legend=c("Satellite","Observations in situ"),col = colores_obs, lwd=2, pch=21,
                   pt.bg = colores_obs, pt.cex=c(1.2,1.2),bty="n", cex=textinner)
    legend(x="bottomright", legend=nombre,
           col=colores, lwd=2, bty="n", cex=textinner)
    
    dev.off()
    
####################################################