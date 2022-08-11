source("Scripts/00_Summary_EDITME.R")
source("Scripts/functions_misc.R")

#######################         
#### FIGURES in SURFACE constant vs variable
#######################         
# Chla, Reflectance, IOP's, PPC
# Summary metrics

o_dir<-paste(global_path,"Res_model/",sep="")        
p_dir<-paste(global_path,"Res_model/interpolated/",sep="")
path_figures<-paste(global_path,"Figures/",sep="")
# List of simulations
experimentos<-read.csv(paste(o_dir,"run_log_marshall_PPC_PS_SA_G100.csv",sep=""), sep=",")
s_dir <- paste(global_path,"Dat_observations/satellite/",sep="")
is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")

         modelnames<-experimentos$name[c(21,9)]
         ruta<-experimentos$path[c(21,9)]        
         nombre <- c(expression(paste("EXP-noPPC: ",{{"a*"}}[PH],"(", lambda, ") constant",sep="")),
                     expression(paste("EXP-PPC: ",{{"a*"}}[PH],"(", lambda, ") variable",sep="")))

         ## Color palette
         mycol1<-usecol(pal_bordeaux,n=4)
         mycol2<-usecol(pal_petrol,n=4)
         colores_puntos<-c(mycol1[3],mycol2[3])
         colores<-c(mycol1[1],mycol2[1])
         #custom_scale <- wes_palette(length(custom_scale),name = "Zissou1", type = "continuous")
         custom_scale <- c("deepskyblue4","deepskyblue3","darkslategray1","mediumaquamarine","greenyellow","yellow",
                           "gold1","orange", "tomato","orangered2")         
         
                         # ## Chlorophyll and PPC
                         # mat <- matrix(c(1:200), ncol=200)
                         # colores<- viridis(n=200, alpha = 1)
                         # image(mat, col=colores, las=1, xaxt="n", yaxt="n", main="", cex.main=1.2)
                         # 
                         # ## Optics
                         # mat <- matrix(c(1:length(custom_scale)), ncol=length(custom_scale))
                         # image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", cex.main=0.8)
                         # 
                         # ## Constant vs. variable
                         # mat <- matrix(c(1:4), ncol=4)
                         # image(mat, col=mycol1, las=1, xaxt="n", yaxt="n", main="", cex.main=1.2)
                         # image(mat, col=mycol2, las=1, xaxt="n", yaxt="n", main="", cex.main=1.2)
                         # ## Constituents
                         # mat <- matrix(c(1:4), ncol=4)
                         # paleta<-wes_palette(4, name = "Zissou1", type = "continuous")
                         # image(mat, col=paleta, las=1, xaxt="n", yaxt="n", main="", cex.main=1.2) 
                         # axis(4,at=seq(0+0.04,1-0.04, length=4),labels=c("water","phyto","cdom","nap"), cex.axis=1.2,las=1, font=2)        
                   
         
                 
                           
########################################
#### FIGURE S4: log(TChla) and PPC:TChla
########################################         

    png(file=paste(path_figures,"FigureS4_Chla_log_PPC_corr.png",sep=""),width = 1200, height = 400, units = "px", pointsize = 17, bg = "white") 
        layout(matrix(c(1:6,10,7:9), ncol=5, byrow=TRUE), widths = c(1,1,1,1,0.3))
        par(family="")
        par(mar=c(1,1,1,0))
        par(oma=c(3,3,2,1))                 
        #layout.show(n=10)
        lettersize=1.2
############################         
         
  ## CHLA log
  ##############    
  ##### In situ MAREDAT and others
         is_dir <- paste(global_path,"Dat_observations/HPLC/grids_pigments",sep="") 
         load(paste(is_dir,"/TChla_insitu_new2.RData",sep=""))
         longitude<-as.numeric(dimnames(TChla_ALL)$x)
         latitude<-as.numeric(dimnames(TChla_ALL)$y)        
         matriz1  <- (apply(TChla_ALL[,,1:2], MARGIN=c("x","y"), mean, na.rm=TRUE))/1e+3 # ng L-1 to ug L-1
  ##### In situ Valente and others: No usaria esta!!
         is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
         load(paste(is_dir,"/media_annual_Chla.Rdata",sep=""))
         matriz2  <- apply(res[,,1:2], MARGIN=c("x","y"), mean, na.rm=TRUE) # ug L-1
         nueva<-abind(matriz1, matriz2, along=0.5)
         matriz<-apply(nueva,MARGIN=c(2,3), mean, na.rm=TRUE)

         image2D(x=longitude, y=latitude, log10(matriz), col=viridis(100),
                 zlim=c(-2,1.2),las=1,bty="n", main="TChla in situ 1988-2019", xaxt="n", yaxt="n",
                 cex.axis=1.0, cex.lab=1.0, cex=1.2,xlab="", ylab="", colkey=F)
             # Map
             load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
             contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
             axis(2, at=seq(-90,90,length=10), labels=T,   las=1)
             axis(1, at=seq(-180,180,length=19), labels=FALSE,  las=1)
             mtext(3,text="a)",at=-175,cex=lettersize, font=2)
             text(x=80,y=-76,labels="n = 2255", font=2, cex=1.1)

  #### SAT
         modelname<-modelnames[1]
         o_dir<-ruta[1]              
         # MODEL
         load(paste(p_dir,"1cloro_log/",modelname,".RData",sep=""))
         mod <- res
         load(paste(s_dir,"climatologies_2012_2018/media_anual_log_chl_OCCCI_2012_2018.RData", sep=""))
         media <- apply(res, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
         media[is.na(mod)] <- NA
         longitude<-as.numeric(dimnames(res)$x)
         latitude<-(-as.numeric(dimnames(res)$y))
         image2D(x=longitude, y=latitude, media, col=viridis(100), las=1, cex.axis=0.8, colkey=FALSE,main="Chla OC-CCI 2012-2018",xaxt="n", yaxt="n",zlim=c(-2,1.2))
             # Map
             load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
             contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
             axis(2, at=seq(-90,90,length=10), labels=FALSE,   las=1)
             axis(1, at=seq(-180,180,length=19), labels=FALSE,  las=1)
             mtext(3,text="b)",at=-175,cex=lettersize, font=2)
             text(x=80,y=-76,labels="n = 9466", font=2, cex=1.1)
             
R_chla_ins<-c(NA,NA)      
AE_chla_ins<-c(NA,NA)  
n_chla_ins<-c(NA,NA) 
R_chla_sat<-c(NA,NA)     
AE_chla_sat<-c(NA,NA)  
n_chla_sat<-c(NA,NA)  
  
         for (w in 1:length(modelnames)){ 
           modelname<-modelnames[w]
           o_dir<-ruta[w]              
           # MODEL
           load(paste(p_dir,"1cloro_log/",modelname,".RData",sep=""))
           mod <- res
           
         ## Correlation in situ
           x<-cbind(c(mod), c(log10(matriz)))
           x[x=="-Inf"] <- NA
           x<-x[complete.cases(x)==TRUE,]
           R_chla_ins<- rcorr(x, type="pearson")$r[1,2]
           AE_chla_ins<- (sum(x[,1]-x[,2]))/nrow(x)  
           n_chla_ins<- nrow(x)
           
          ## Correlation sat
           x<-cbind(c(mod), c(media))
           x[x=="-Inf"] <- NA
           x<-x[complete.cases(x)==TRUE,]
           R_chla_sat<- rcorr(x, type="pearson")$r[1,2]
           AE_chla_sat<- (sum(x[,1]-x[,2]))/nrow(x)
           n_chla_sat<- nrow(x)
           image2D(x=longitude, y=latitude, mod,  col=viridis(100), las=1, cex.axis=0.8, main=nombre[w], colkey=FALSE,xaxt="n", yaxt="n",zlim=c(-2,1.2))
           
           # Map
           load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
           contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50", xaxt="n", yaxt="n", bty="n",add=T)
           axis(2, at=seq(-90,90,length=10), labels=FALSE,   las=1)
           axis(1, at=seq(-180,180,length=19), labels=FALSE,  las=1)

           # Metrics           
           text(x=-110, y=-70, labels=paste("R = ",round(R_chla_ins,3)), font=2)
           text(x=-110, y=-80, labels=paste("Bias = ",round(AE_chla_ins,3)), font=2)
           text(x=75, y=-70, labels=paste("R = ",round(R_chla_sat,3)), font=2)
           text(x=75, y=-80, labels=paste("Bias = ",round(AE_chla_sat,3)), font=2)

           if (w==1){mtext(3,text="c)",at=-175,cex=lettersize, font=2)}
           if (w==2){mtext(3,text="d)",at=-175,cex=lettersize, font=2)}
           }  # end loop model
         
         
         ###ESCALA
         par(mar=c(1,1,1,3))
         mat <- matrix(c(1:100), ncol=100)
         image(mat, col=viridis(100), las=1, xaxt="n", yaxt="n", main=expression(paste("mg ",m^-3,sep="")), cex.main=0.8) 
         axis(4,at=seq(0-0.005,1+0.005, length=7)[seq(1,7,by=2)],labels=10^seq(-2,1.2,by=0.5)[seq(1,7,by=2)], cex.axis=1.2, las=1)
  ##############         
         
  ## PPC:TChla      
  ############      
         # In situ
         is_dir <- paste(global_path,"Dat_observations/HPLC/grids_pigments",sep="") 
         load(paste(is_dir,"/PPCinsitu_new3.RData",sep=""))         
         matriz1  <- apply(PPCinsitu_ALL[,,1:2], MARGIN=c("x","y"), mean, na.rm=TRUE)
         matriz11 <- apply(PPCinsitu_NEW[,,1:2], MARGIN=c("x","y"), mean, na.rm=TRUE)
         matriz111 <- apply(PPCinsitu_MAREDAT[,,1:2], MARGIN=c("x","y"), mean, na.rm=TRUE)

         par(mar=c(1,1,1,0))
         image2D(x=longitude, y=latitude, matriz1, col=viridis(100), las=1,bty="n", main="PPC:TChla in situ 1992-2019", xaxt="n", yaxt="n", zlim=c(0.0,1),cex.axis=1.0, cex.lab=1.0, cex=1.2,xlab="", ylab="", mgp=c(1.5,1,0), colkey=FALSE)
             # Map
             load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
             contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
             axis(2, at=seq(-90,90,length=10), labels=T,   las=1)
             axis(1, at=seq(-180,180,length=19), labels=T,  las=1)
             mtext(3,text="e)",at=-175,cex=lettersize, font=2)
             text(x=80, y=-76,labels="n = 1386", font=2, cex=1.1)
             

         #Mod
R_ppc<-c(NA,NA)      
AE_ppc<-c(NA,NA) 
n_ppc<-c(NA,NA)  
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
           R_ppc<- rcorr(x, type="pearson")$r[1,2]
           AE_ppc<- (sum(x[,1]-x[,2]))/nrow(x)
           n_ppc<- nrow(x)
           image2D(x=longitude, y=latitude, mod,  col=viridis(100), las=1, cex.axis=0.8, main=nombre[w], colkey=FALSE,xaxt="n", yaxt="n",zlim=c(0.0,1))
               # Map
               load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
               contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50", xaxt="n", yaxt="n", bty="n",add=T)
               axis(2, at=seq(-90,90,length=10), labels=FALSE,   las=1)
               axis(1, at=seq(-180,180,length=19), labels=T,  las=1)
               
           # Metrics
           text(x=-110, y=-70, labels=paste("R = ",round(R_ppc,3)), font=2)
           text(x=-110, y=-80, labels=paste("Bias = ",round(AE_ppc,3)), font=2)

           if (w==1){mtext(3,text="f)",at=-175,cex=lettersize, font=2)}
           if (w==2){mtext(3,text="g)",at=-175,cex=lettersize, font=2)}           
         }  # end loop model
         
         

         ###ESCALA
          par(mar=c(1,1,1,3))
          mat <- matrix(c(1:100), ncol=100)
          image(mat, col=viridis(100), las=1, xaxt="n", yaxt="n", main="g:g", cex.main=0.8) 
          axis(4,at=seq(0-0.005,1+0.005, length=6),labels=seq(0,1,length=6), cex.axis=1.2, las=1)
  ############      
        
############################ 
          
          ## Labels
          mtext(1,at=seq(0.12,0.815,length=4)[-2], line=1.5,outer=T,text=expression(paste(degree, "E")))
          mtext(2,at=c(0.25,0.75),                 line=1.2,outer=T,text=expression(paste(degree, "N")), las=1)          
          
     dev.off()
        
         
         
          
          
          
        
         
###########################    
#### FIGURE S5: REFLECTANCE
########################### 
          
     png(file=paste(path_figures,"FigureS5_Rrs_corr.png",sep=""),width = 1200, height = 1200, units = "px", pointsize =20, bg = "white")  
         layout(cbind(matrix(c(1:24),ncol=4,byrow=F),c(25:30)),widths=c(rep(1,4),0.3))
         #par(mfcol=c(6,4))
         par(mar=c(1,0.5,1,0))
         par(oma=c(2,2,2,1))
         layout.show(n=30)
         lettersize=1.2
############################          
         
         
##### In situ Rrs
#################         
         is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
         load(paste(is_dir,"/media_annual_Rrs.Rdata",sep=""))           
         
         #400
         image2D(x=longitude, y=latitude,matriz_annual[,,1], col=custom_scale, zlim=c(0,0.020), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                 main=expression(paste(R[RS], " (400nm) ",italic("in situ"), sep="")))
                 # Map
                 load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                 contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                 axis(2, at=seq(-90,90,length=5), labels=T,   las=1)
                 axis(1, at=seq(-180,180,length=10), labels=FALSE,  las=1)
                 mtext(3,text="a)",at=-175,cex=lettersize, font=2) 
                 text(x=75,y=-76,labels="n = 659", font=2, cex=1.1)
                 
         #450 
         image2D(x=longitude, y=latitude,matriz_annual[,,3], col=custom_scale, zlim=c(0,0.020), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                 main=expression(paste(R[RS], " (450nm) ",italic("in situ"), sep="")))
                 # Map
                 load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                 contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                 axis(2, at=seq(-90,90,length=5), labels=T,   las=1)
                 axis(1, at=seq(-180,180,length=10), labels=FALSE,  las=1)
                 mtext(3,text="b)",at=-175,cex=lettersize, font=2) 
                 text(x=75,y=-76,labels="n = 747", font=2, cex=1.1)
                 
         #475
         image2D(x=longitude, y=latitude, matriz_annual[,,4], col=custom_scale, zlim=c(0,0.015), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                 main=expression(paste(R[RS], " (475nm) ",italic("in situ"), sep="")))
                 # Map
                 load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                 contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                 axis(2, at=seq(-90,90,length=5), labels=T,   las=1)
                 axis(1, at=seq(-180,180,length=10), labels=FALSE,  las=1)
                 mtext(3,text="c)",at=-175,cex=lettersize, font=2)
                 text(x=75,y=-76,labels="n = 219", font=2, cex=1.1)

         #500
         image2D(x=longitude, y=latitude, matriz_annual[,,5], col=custom_scale, zlim=c(0,0.010), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                 main=expression(paste(R[RS], " (500nm) ",italic("in situ"), sep="")))
                 # Map
                 load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                 contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                 axis(2, at=seq(-90,90,length=5), labels=T,   las=1)
                 axis(1, at=seq(-180,180,length=10), labels=FALSE,  las=1)
                 mtext(3,text="d)",at=-175,cex=lettersize, font=2)  
                 text(x=75,y=-76,labels="n = 754", font=2, cex=1.1)
                 
         #550
         image2D(x=longitude, y=latitude, matriz_annual[,,7], col=custom_scale, zlim=c(0,0.005), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                 main=expression(paste(R[RS], " (550nm) ",italic("in situ"), sep="")))
                 # Map
                 load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                 contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                 axis(2, at=seq(-90,90,length=5), labels=T,   las=1)
                 axis(1, at=seq(-180,180,length=10), labels=FALSE,  las=1)
                 mtext(3,text="e)",at=-175,cex=lettersize, font=2) 
                 text(x=75,y=-76,labels="n = 722", font=2, cex=1.1)
         #675
         image2D(x=longitude, y=latitude, matriz_annual[,,12], col=custom_scale, zlim=c(0,0.001), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                 main=expression(paste(R[RS], " (675nm) ",italic("in situ"), sep="")))
                 # Map
                 load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                 contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                 axis(2, at=seq(-90,90,length=5), labels=T,   las=1)
                 axis(1, at=seq(-180,180,length=10),
                 labels=c("",seq(-180,180,length=10)[2],"",seq(-180,180,length=10)[4],"",seq(-180,180,length=10)[6],"",seq(-180,180,length=10)[8],"",seq(-180,180,length=10)[10]),las=1)
                 mtext(3,text="f)",at=-175,cex=lettersize, font=2) 
                 text(x=75,y=-76,labels="n = 699", font=2, cex=1.1)
############         

####### SAT
###########         
         modelname<-modelnames[1]
         load(paste(p_dir,"1reflectance/",modelname,".RData",sep=""))
         mod <- res
         z <- mod[,,12]
           landas_sat <- c("Rrs_412","Rrs_443","Rrs_490","Rrs_510","Rrs_555","Rrs_670")
           load(paste(s_dir,"climatologies_2012_2018/media_seasonal_Rrs_2deg_OCCCI_2012_2018.RData",sep=""))
           res <- med
           longitude<-as.numeric(dimnames(res)$x)
           latitude<-(-as.numeric(dimnames(res)$y))
           
           media412 <- apply(res[,,,1], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
           media412[is.na(z)] <- NA
           image2D(x=longitude, y=latitude,media412, col=custom_scale, zlim=c(0,0.020), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                   main=expression(paste(R[RS], " (412nm) ",italic("OC-CCI"), sep="")))
                   # Map
                   load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                   contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                   axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                   axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                   mtext(3,text="g)",at=-175,cex=lettersize, font=2) 
                   text(x=75,y=-76,labels="n = 9465", font=2, cex=1.1)
                   
           media443 <- apply(res[,,,2], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
           media443[is.na(z)] <- NA
           image2D(x=longitude, y=latitude,media443, col=custom_scale, zlim=c(0,0.020), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                   main=expression(paste(R[RS], " (443nm) ",italic("OC-CCI"), sep="")))
                   # Map
                   load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                   contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                   axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                   axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                   mtext(3,text="h)",at=-175,cex=lettersize, font=2) 
                   text(x=75,y=-76,labels="n = 9465", font=2, cex=1.1)
                   
           media490 <- apply(res[,,,3], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
           media490[is.na(z)] <- NA
           image2D(x=longitude, y=latitude,media490, col=custom_scale, zlim=c(0,0.015), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                   main=expression(paste(R[RS], " (490nm) ",italic("OC-CCI"), sep="")))
                   # Map
                   load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                   contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                   axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                   axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                   mtext(3,text="i)",at=-175,cex=lettersize, font=2) 
                   text(x=75,y=-76,labels="n = 9465", font=2, cex=1.1)
                   
           media510 <- apply(res[,,,4], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
           media510[is.na(z)] <- NA
           image2D(x=longitude, y=latitude,media510, col=custom_scale, zlim=c(0,0.010), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                   main=expression(paste(R[RS], " (510nm) ",italic("OC-CCI"), sep="")))
                   # Map
                   load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                   contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                   axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                   axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                   mtext(3,text="j)",at=-175,cex=lettersize, font=2) 
                   text(x=75,y=-76,labels="n = 9465", font=2, cex=1.1)
                   
           media555 <- apply(res[,,,5], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
           media555[is.na(z)] <- NA
           image2D(x=longitude, y=latitude,media555, col=custom_scale, zlim=c(0,0.005), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                   main=expression(paste(R[RS], " (555nm) ",italic("OC-CCI"), sep="")))
                   # Map
                   load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                   contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                   axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                   axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                   mtext(3,text="k)",at=-175,cex=lettersize, font=2) 
                   text(x=75,y=-76,labels="n = 9465", font=2, cex=1.1)
                   
           media670 <- apply(res[,,,6], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
           media670[is.na(z)] <- NA
           image2D(x=longitude, y=latitude,media670, col=custom_scale, zlim=c(0,0.001), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                   main=expression(paste(R[RS], " (670nm) ",italic("OC-CCI"), sep="")))
                   # Map
                   load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                   contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                   axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                   axis(1, at=seq(-180,180,length=10), labels=c("",seq(-180,180,length=10)[2],"",seq(-180,180,length=10)[4],"",seq(-180,180,length=10)[6],"",seq(-180,180,length=10)[8],"",seq(-180,180,length=10)[10]),  las=1)
                   mtext(3,text="l)",at=-175,cex=lettersize, font=2) 
                   text(x=75,y=-76,labels="n = 9465", font=2, cex=1.1)
               
#######################         

                   
#### MODEL
##########          
         for (w in 1:length(modelnames)){  
           modelname<-modelnames[w]
           #o_dir<-ruta[w]    
           # MODEL
           load(paste(p_dir,"1reflectance/",modelname,".RData",sep=""))
           mod <- res
           z <- mod[,,12]
           if (w==1) {image2D(x=longitude, y=latitude, mod[,,1],  col=custom_scale, zlim=c(0,0.020), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                   main=expression(paste(R[RS], " (400nm) ",italic("EXP-noPPC"), sep="")))}else{
                     image2D(x=longitude, y=latitude, mod[,,1],  col=custom_scale, zlim=c(0,0.020), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                             main=expression(paste(R[RS], " (400nm) ",italic("EXP-PPC"), sep=""))) }
                   # Correlation Insitu 400 vs 400
                   x<-cbind(c(mod[,,1]), c(matriz_annual[,,1]))
                   x[x=="-Inf"] <- NA
                   x<-x[complete.cases(x)==TRUE,]
                   R_412<- rcorr(x, type="pearson")$r[1,2  ]
                   nrow(x)
                   text(x=-75,y=-76,labels=paste("R = ",round(R_412,3)), font=2, cex=1.1) 
                   
                   # Correlation Sat 412 vs 400
                   x<-cbind(c(mod[,,1]), c(media412))
                   x[x=="-Inf"] <- NA
                   x<-x[complete.cases(x)==TRUE,]
                   R_412<- rcorr(x, type="pearson")$r[1,2  ]
                   nrow(x)
                   text(x=75,y=-76,labels=paste("R = ",round(R_412,3)), font=2, cex=1.1)
                        # Map
                        load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                        contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                        axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                        axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                        if (w==1){mtext(3,text="m)",at=-175,cex=lettersize, font=2)} 
                        if (w==2){mtext(3,text="s)",at=-175,cex=lettersize, font=2)}  
           
        if (w==1) {image2D(x=longitude, y=latitude,mod[,,3],  col=custom_scale, zlim=c(0,0.020), las=1, xaxt="n", yaxt="n",colkey=FALSE,
                   main=expression(paste(R[RS], " (450nm) ",italic("EXP-noPPC"), sep="")))}else{
                   image2D(x=longitude, y=latitude,mod[,,3],  col=custom_scale, zlim=c(0,0.020), las=1, xaxt="n", yaxt="n",colkey=FALSE,
                   main=expression(paste(R[RS], " (450nm) ",italic("EXP-PPC"), sep="")))}
                   # Correlation Insitu 450 vs 450
                   x<-cbind(c(mod[,,3]), c(matriz_annual[,,3]))
                   x[x=="-Inf"] <- NA
                   x<-x[complete.cases(x)==TRUE,]
                   R_443<- rcorr(x, type="pearson")$r[1,2  ]
                   nrow(x)
                   text(x=-75,y=-76,labels=paste("R = ",round(R_443,3)), font=2, cex=1.1)           
                   # Correlation 443 vs 450
                   x<-cbind(c(mod[,,3]), c(media443))
                   x[x=="-Inf"] <- NA
                   x<-x[complete.cases(x)==TRUE,]
                   R_443<- rcorr(x, type="pearson")$r[1,2  ]
                   nrow(x)
                   text(x=75,y=-76,labels=paste("R = ",round(R_443,3)), font=2, cex=1.1)
                         # Map
                         load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                         contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                         axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                         axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                         if (w==1){mtext(3,text="n)",at=-175,cex=lettersize, font=2)} 
                         if (w==2){mtext(3,text="t)",at=-175,cex=lettersize, font=2)}                          
                 
         if (w==1) {image2D(x=longitude, y=latitude,mod[,,4],  col=custom_scale, zlim=c(0,0.015), las=1, xaxt="n", yaxt="n",colkey=FALSE,
                   main=expression(paste(R[RS], " (475nm) ",italic("EXP-noPPC"), sep=""))) }else{ 
                   image2D(x=longitude, y=latitude,mod[,,4],  col=custom_scale, zlim=c(0,0.015), las=1, xaxt="n", yaxt="n",colkey=FALSE,
                   main=expression(paste(R[RS], " (475nm) ",italic("EXP-PPC"), sep=""))) }
                   # Correlation Insitu 475 vs 475
                   x<-cbind(c(mod[,,4]), c(matriz_annual[,,4]))
                   x[x=="-Inf"] <- NA
                   x<-x[complete.cases(x)==TRUE,]
                   R_490<- rcorr(x, type="pearson")$r[1,2  ]
                   text(x=-75,y=-76,labels=paste("R = ",round(R_490,3)), font=2, cex=1.1)           
                   # Correlation 490 vs 475
                   x<-cbind(c(mod[,,4]), c(media490))
                   x[x=="-Inf"] <- NA
                   x<-x[complete.cases(x)==TRUE,]
                   R_490<- rcorr(x, type="pearson")$r[1,2  ]
                   text(x=75,y=-76,labels=paste("R = ",round(R_490,3)), font=2, cex=1.1)
                         # Map
                         load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                         contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                         axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                         axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                         if (w==1){mtext(3,text="o)",at=-175,cex=lettersize, font=2)} 
                         if (w==2){mtext(3,text="u)",at=-175,cex=lettersize, font=2)}                   
                 
        if (w==1) {image2D(x=longitude, y=latitude,mod[,,5],  col=custom_scale, zlim=c(0,0.010), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                   main=expression(paste(R[RS], " (500nm) ",italic("EXP-noPPC"), sep="")))}else{
                   image2D(x=longitude, y=latitude,mod[,,5],  col=custom_scale, zlim=c(0,0.010), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                   main=expression(paste(R[RS], " (500nm) ",italic("EXP-PPC"), sep="")))}
                   # Correlation Insitu 500 vs 500
                   x<-cbind(c(mod[,,5]), c(matriz_annual[,,5]))
                   x[x=="-Inf"] <- NA
                   x<-x[complete.cases(x)==TRUE,]
                   R_510<- rcorr(x, type="pearson")$r[1,2  ]
                   text(x=-75,y=-76,labels=paste("R = ",round(R_510,3)), font=2, cex=1.1)              
                   # Correlation 510 vs 500
                   x<-cbind(c(mod[,,5]), c(media510))
                   x[x=="-Inf"] <- NA
                   x<-x[complete.cases(x)==TRUE,]
                   R_510<- rcorr(x, type="pearson")$r[1,2  ]
                   text(x=75,y=-76,labels=paste("R = ",round(R_510,3)), font=2, cex=1.1)           
                         # Map
                         load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                         contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                         axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                         axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                         if (w==1){mtext(3,text="p)",at=-175,cex=lettersize, font=2)} 
                         if (w==2){mtext(3,text="v)",at=-175,cex=lettersize, font=2)}                      
                 
           
        if (w==1) {image2D(x=longitude, y=latitude,mod[,,7],  col=custom_scale, zlim=c(0,0.005), las=1, xaxt="n", yaxt="n",colkey=FALSE,
                   main=expression(paste(R[RS], " (550nm) ",italic("EXP-noPPC"), sep="")))}else{
                   image2D(x=longitude, y=latitude,mod[,,7],  col=custom_scale, zlim=c(0,0.005), las=1, xaxt="n", yaxt="n",colkey=FALSE,
                   main=expression(paste(R[RS], " (550nm) ",italic("EXP-PPC"), sep="")))}
                   # Correlation Insitu 550 vs 550
                   x<-cbind(c(mod[,,7]), c(matriz_annual[,,7]))
                   x[x=="-Inf"] <- NA
                   x<-x[complete.cases(x)==TRUE,]
                   R_555<- rcorr(x, type="pearson")$r[1,2  ]
                   text(x=-75,y=-76,labels=paste("R = ",round(R_555,3)), font=2, cex=1.1)             
                   # Correlation 555 vs 550
                   x<-cbind(c(mod[,,7]), c(media555))
                   x[x=="-Inf"] <- NA
                   x<-x[complete.cases(x)==TRUE,]
                   R_555<- rcorr(x, type="pearson")$r[1,2]
                   text(x=75,y=-76,labels=paste("R = ",round(R_555,3)), font=2, cex=1.1) 
                           # Map
                           load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep="")) 
                           contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                           axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                           axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                           if (w==1){mtext(3,text="q)",at=-175,cex=lettersize, font=2)} 
                           if (w==2){mtext(3,text="w)",at=-175,cex=lettersize, font=2)}                     
                   
        if (w==1) {image2D(x=longitude, y=latitude,mod[,,12],  col=custom_scale, zlim=c(0,0.001), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                   main=expression(paste(R[RS], " (675nm) ",italic("EXP-noPPC"), sep="")))}else{
                   image2D(x=longitude, y=latitude,mod[,,12],  col=custom_scale, zlim=c(0,0.001), las=1, xaxt="n", yaxt="n", colkey=FALSE,
                   main=expression(paste(R[RS], " (675nm) ",italic("EXP-PPC"), sep="")))}
                   # Correlation Insitu 675 vs 675
                   x<-cbind(c(mod[,,12]), c(matriz_annual[,,12]))
                   x[x=="-Inf"] <- NA
                   x<-x[complete.cases(x)==TRUE,]
                   R_670<- rcorr(x, type="pearson")$r[1,2  ]
                   text(x=-75,y=-76,labels=paste("R = ",round(R_670,3)), font=2, cex=1.1)            
                   # Correlation 670 vs 675
                   x<-cbind(c(mod[,,12]), c(media670))
                   x[x=="-Inf"] <- NA
                   x<-x[complete.cases(x)==TRUE,]
                   R_670<- rcorr(x, type="pearson")$r[1,2]
                   text(x=75,y=-76,labels=paste("R = ",round(R_670,3)), font=2, cex=1.1)
                           # Map
                           load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                           contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                           axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                           axis(1, at=seq(-180,180,length=10), labels=c("",seq(-180,180,length=10)[2],"",seq(-180,180,length=10)[4],"",seq(-180,180,length=10)[6],"",seq(-180,180,length=10)[8],"",seq(-180,180,length=10)[10]),  las=1)
                           if (w==1){mtext(3,text="r)",at=-175,cex=lettersize, font=2)} 
                           if (w==2){mtext(3,text="x)",at=-175,cex=lettersize, font=2)}                      
             }
#######################
           
           ### SCALE
           par(mar=c(1,0.5,1,3))
           mat <- matrix(c(1:length(custom_scale)), ncol=length(custom_scale))
     image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(sr^-1), cex.main=0.8) 
     axis(4,at=seq(0-0.05,1+0.05, length=6),labels=round(seq(0,0.02,length=6),3), cex.axis=1.2, las=1)
           
     image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(sr^-1), cex.main=0.8) 
     axis(4,at=seq(0-0.05,1+0.05, length=6),labels=round(seq(0,0.02,length=6),3), cex.axis=1.2, las=1)

     image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(sr^-1), cex.main=0.8) 
     axis(4,at=seq(0-0.05,1+0.05, length=6),labels=round(seq(0,0.015,length=6),3), cex.axis=1.2, las=1)

     image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(sr^-1), cex.main=0.8) 
     axis(4,at=seq(0-0.05,1+0.05, length=6),labels=round(seq(0,0.01,length=6),3), cex.axis=1.2, las=1)

     image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(sr^-1), cex.main=0.8) 
     axis(4,at=seq(0-0.05,1+0.05, length=6),labels=round(seq(0,0.005,length=6),3), cex.axis=1.2, las=1)

     image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(sr^-1), cex.main=0.8) 
     axis(4,at=seq(0-0.05,1+0.05, length=6),labels=round(seq(0,0.001,length=6),4), cex.axis=1.2, las=1)
     
     dev.off()
     
         
      
     
     
        
         
#############################
#### FIGURE S6: IOP's surface
############################# 
         
     png(file=paste(path_figures,"FigureS6_IOPs_corr.png",sep=""),width = 1200, height = 1200, units = "px", pointsize = 20, bg = "white")  
           layout(matrix(c(1:30),ncol=5,byrow=T),widths=c(rep(1,4),0.3))
           par(mar=c(1,0.5,1,0))
           par(oma=c(2,2,2,1))
           layout.show(n=30)
           textinner=1.2
           lettersize=1.2
############################         
         
#### Atot
###########################           
           ##### In situ
           plot(1,1,type="n",yaxt="n",xaxt="n",bty="n")
           
           #### Sat
           modelname<-modelnames[1]
           load(paste(p_dir,"2iop_atot/",modelname,".RData",sep=""))
           z <- apply(res,MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
           z[z<0.01] <- NA
           MASCARA <- z
           load(paste(s_dir,"climatologies_2012_2018/media_seasonal_Atot_2deg_OCCCI_2012_2018.RData",sep=""))
           res <- med           
           landa_sat <- dimnames(res)$l
           longitude<-as.numeric(dimnames(res)$x)
           latitude<-as.numeric(dimnames(res)$y)
           sat <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)        
           media <- sat[,,2]
           media[is.na(z)] <- NA
           image2D(x=longitude, y=latitude,media, col=custom_scale, zlim=c(0,0.2), las=1, xaxt="n", yaxt="n",colkey=F,
                   main=expression(paste(a[TOT], " (443nm) ",italic("OC-CCI"), sep="")))
                   # Map
                   load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                   contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                   axis(2, at=seq(-90,90,length=5), labels=T,   las=1)
                   axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                   mtext(3,text="a)",at=-175,cex=lettersize, font=2)  
                   text(x=75,y=-76,labels="n = 9466", font=2, cex=1.1)
                           
           # Mod  
           for (w in 1:length(modelnames)){ 
             modelname<-modelnames[w]
             load(paste(p_dir,"2iop_atot/",modelname,".RData",sep=""))
             z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
             z[z<0.01] <- NA
           if (w==1){ image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.2), las=1, xaxt="n", yaxt="n",colkey=F,
                      main=expression(paste(a[TOT], " (450nm) ",italic("EXP-noPPC"), sep="")))}else{
                      image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.2), las=1, xaxt="n", yaxt="n",colkey=F,
                      main=expression(paste(a[TOT], " (450nm) ",italic("EXP-PPC"), sep="")))}
                   # Correlation
                   x<-cbind(c(z), c(media))
                   x[x=="-Inf"] <- NA
                   x<-x[complete.cases(x)==TRUE,]
                   tot<- rcorr(x, type="pearson")$r[1,2]
                   text(x=75,y=-77,labels=paste("R = ",round(tot,3)), font=2, cex=textinner)
                         # Map
                         load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                         contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                         axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                         axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
             if (w==1){mtext(3,text="b)",at=-175,cex=lettersize, font=2)   }
             if (w==2){mtext(3,text="c)",at=-175,cex=lettersize, font=2)   } 
           }
###########################    
           ### ESCALA
           par(mar=c(1,0.5,1,3))
           mat <- matrix(c(1:length(custom_scale)), ncol=length(custom_scale))
           image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(m^-1), cex.main=0.8) 
           axis(4,at=seq(0-0.05,1+0.05, length=6),labels=round(seq(0,0.2,length=6),2), cex.axis=1.2, las=1)
           
          
#### Aph
           par(mar=c(1,0.5,1,0))
############################           
           ##### In situ
           is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
           load(paste(is_dir,"/media_annual_Aph.Rdata",sep=""))           
           image2D(x=longitude, y=latitude,matriz_annual[,,1,3], col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n", colkey=F,
                  main=expression(paste(a[PH], " (450nm) ",italic("in situ"), sep="")))
                  # Map
                  load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                  contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                  axis(2, at=seq(-90,90,length=5), labels=T,   las=1)
                  axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                  mtext(3,text="d)",at=-175,cex=lettersize, font=2)  
                  text(x=75,y=-76,labels="n = 457", font=2, cex=1.1)

           # Sat
           modelname<-modelnames[1]
           load(paste(p_dir,"2iop_aph/",modelname,".RData",sep=""))
           z <- apply(res,MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
           z[z<1e-5] <- NA
           MASCARA<- z
           load(paste(s_dir,"climatologies_2012_2018/media_seasonal_Aph_2deg_OCCCI_2012_2018.RData",sep=""))
           res <- med            
           landa_sat <- dimnames(res)$l
           sat <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)        
           media <- sat[,,2]
           media[is.na(z)] <- NA
           image2D(x=longitude, y=latitude,media, col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n",colkey=F,
                   main=expression(paste(a[PH], " (443nm) ",italic("OC-CCI"), sep="")))
                   # Map
                   load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                   contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                   axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                   axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                   mtext(3,text="e)",at=-175,cex=lettersize, font=2)
                   text(x=75,y=-76,labels="n = 9466", font=2, cex=1.1)
                   
           # Mod  
           for (w in 1:length(modelnames)){ 
             modelname<-modelnames[w]
             load(paste(p_dir,"2iop_aph/",modelname,".RData",sep=""))
             z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
             z[z<1e-5] <- NA
           if (w==1){ image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n",colkey=F,
                      main=expression(paste(a[PH], " (450nm) ",italic("EXP-noPPC"), sep="")))}else{
                      image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n",colkey=F,
                      main=expression(paste(a[PH], " (450nm) ",italic("EXP-PPC"), sep="")))}            
                     # Correlation Insitu 450 vs 450
                     x<-cbind(c(z), c(matriz_annual[,,1,3]))
                     x[x=="-Inf"] <- NA
                     x<-x[complete.cases(x)==TRUE,]
                     aph<- rcorr(x, type="pearson")$r[1,2  ]
                     text(x=-110,y=-77,labels=paste("R = ",round(aph,3)), font=2, cex=textinner)               
                     # Correlation
                     x<-cbind(c(z), c(media))
                     x[x=="-Inf"] <- NA
                     x<-x[complete.cases(x)==TRUE,]
                     aph<- rcorr(x, type="pearson")$r[1,2]
                     text(x=75,y=-77,labels=paste("R = ",round(aph,3)), font=2, cex=textinner)
                           # Map
                           load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                           contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                           axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                           axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                           if (w==1){mtext(3,text="f)",at=-175,cex=lettersize, font=2)   }
                           if (w==2){mtext(3,text="g)",at=-175,cex=lettersize, font=2)   }              
              }        
#############################           
           ### ESCALA
           par(mar=c(1,0.5,1,3))
           mat <- matrix(c(1:length(custom_scale)), ncol=length(custom_scale))
           image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(m^-1), cex.main=0.8) 
           axis(4,at=seq(0-0.05,1+0.05, length=6),labels=round(seq(0,0.1,length=6),2), cex.axis=1.2, las=1)
           
           
#### Anap
           par(mar=c(1,0.5,1,0))
############################           
           ##### In situ
           is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
           load(paste(is_dir,"/media_annual_Anap.Rdata",sep=""))
           image2D(x=longitude, y=latitude,matriz_annual[,,1,3], col=custom_scale, zlim=c(0,0.03), las=1, xaxt="n", yaxt="n", colkey=F,
                    main=expression(paste(a[NAP], " (450nm) ",italic("in situ"), sep="")))
                    # Map
                    load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                    contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                    axis(2, at=seq(-90,90,length=5), labels=T,   las=1)
                    axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                    mtext(3,text="h)",at=-175,cex=lettersize, font=2)  
                    text(x=75,y=-76,labels="n = 502", font=2, cex=1.1)
                    
           # Sat
           plot(1,1,type="n",yaxt="n",xaxt="n",bty="n")
           
           # Mod  
           for (w in 1:length(modelnames)){ 
             modelname<-modelnames[w]
             load(paste(p_dir,"2iop_apt/",modelname,".RData",sep=""))
             z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
             z[z<1e-5] <- NA
             if (w==1){image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.03), las=1, xaxt="n", yaxt="n",colkey=F,
                       main=expression(paste(a[NAP], " (450nm) ",italic("EXP-noPPC"), sep="")))}else{
                       image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.03), las=1, xaxt="n", yaxt="n",colkey=F,
                       main=expression(paste(a[NAP], " (450nm) ",italic("EXP-PPC"), sep="")))}        
                     # Correlation Insitu 450 vs 450
                     x<-cbind(c(z), c(matriz_annual[,,1,3]))
                     x[x=="-Inf"] <- NA
                     x<-x[complete.cases(x)==TRUE,]
                     aph<- rcorr(x, type="pearson")$r[1,2  ]
                     text(x=-110,y=-77,labels=paste("R = ",round(aph,3)), font=2, cex=textinner)
                           # Map
                           load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                           contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                           axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                           axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                           if (w==1){mtext(3,text="i)",at=-175,cex=lettersize, font=2)   }
                           if (w==2){mtext(3,text="j)",at=-175,cex=lettersize, font=2)   }                
           }        
############################           
           ### ESCALA
           par(mar=c(1,0.5,1,3))
           mat <- matrix(c(1:length(custom_scale)), ncol=length(custom_scale))
           image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(m^-1), cex.main=0.8) 
           axis(4,at=seq(0-0.05,1+0.05, length=6),labels=round(seq(0,0.03,length=6),3), cex.axis=1.2, las=1)
           
           
#### Ap
           par(mar=c(1,0.5,1,0))
############################           
           ##### In situ
           is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
           load(paste(is_dir,"/media_annual_Ap.Rdata",sep=""))           
           image2D(x=longitude, y=latitude,matriz_annual[,,1,3], col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n", colkey=F,
                    main=expression(paste(a[P], " (450nm) ",italic("in situ"), sep="")))
                    # Map
                    load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                    contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                    axis(2, at=seq(-90,90,length=5), labels=T,   las=1)
                    axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                    mtext(3,text="k)",at=-175,cex=lettersize, font=2)    
                    text(x=75,y=-76,labels="n = 622", font=2, cex=1.1)
           # Sat
           plot(1,1,type="n",yaxt="n",xaxt="n",bty="n")
           
           # Mod  
           for (w in 1:length(modelnames)){ 
             modelname<-modelnames[w]
             load(paste(p_dir,"2iop_app/",modelname,".RData",sep=""))
             z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
             z[z<1e-5] <- NA
           if (w==1){image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n",colkey=F,
                     main=expression(paste(a[P], " (450nm) ",italic("EXP-noPPC"), sep="")))}else{
                     image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n",colkey=F,
                     main=expression(paste(a[P], " (450nm) ",italic("EXP-PPC"), sep="")))}        
                     # Correlation Insitu 450 vs 450
                     x<-cbind(c(z), c(matriz_annual[,,1,3]))
                     x[x=="-Inf"] <- NA
                     x<-x[complete.cases(x)==TRUE,]
                     aph<- rcorr(x, type="pearson")$r[1,2  ]
                     text(x=-110,y=-77,labels=paste("R = ",round(aph,3)), font=2, cex=textinner)
                           # Map
                           load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                           contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                           axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                           axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                           if (w==1){mtext(3,text="l)",at=-175,cex=lettersize, font=2)   }
                           if (w==2){mtext(3,text="m)",at=-175,cex=lettersize, font=2)   }   
           }        
############################           
           ### ESCALA
           par(mar=c(1,0.5,1,3))
           mat <- matrix(c(1:length(custom_scale)), ncol=length(custom_scale))
           image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(m^-1), cex.main=0.8) 
           axis(4,at=seq(0-0.05,1+0.05, length=6),labels=round(seq(0,0.1,length=6),3), cex.axis=1.2, las=1)
           
                      
#### Acdom
           par(mar=c(1,0.5,1,0))
#############################           
           ##### In situ
           is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
           load(paste(is_dir,"/media_annual_Acdom.Rdata",sep=""))           
           image2D(x=longitude, y=latitude,matriz_annual[,,1,3], col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n", colkey=F,
                      main=expression(paste(a[CDOM], " (450nm) ",italic("in situ"), sep="")))
                      # Map
                      load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                      contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                      axis(2, at=seq(-90,90,length=5), labels=T,   las=1)
                      axis(1, at=seq(-180,180,length=10), labels=c("",seq(-180,180,length=10)[2],"",seq(-180,180,length=10)[4],"",seq(-180,180,length=10)[6],"",seq(-180,180,length=10)[8],"",seq(-180,180,length=10)[10]),  las=1)
                      mtext(3,text="n)",at=-175,cex=lettersize, font=2) 
                      text(x=75,y=-76,labels="n = 279", font=2, cex=1.1)
           # Sat
           plot(1,1,type="n",yaxt="n",xaxt="n",bty="n")
           
           # Mod  
           for (w in 1:length(modelnames)){ 
             modelname<-modelnames[w]
             load(paste(p_dir,"2iop_acd/",modelname,".RData",sep=""))
             z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
             z[z<1e-5] <- NA
        if (w==1){   image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n",colkey=F,
                     main=expression(paste(a[CDOM], " (450nm) ",italic("EXP-noPPC"), sep="")))}else{
                     image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n",colkey=F,
                     main=expression(paste(a[CDOM], " (450nm) ",italic("EXP-PPC"), sep="")))}                
                     # Correlation Insitu 450 vs 450
                     x<-cbind(c(z), c(matriz_annual[,,1,3]))
                     x[x=="-Inf"] <- NA
                     x<-x[complete.cases(x)==TRUE,]
                     aph<- rcorr(x, type="pearson")$r[1,2  ]
                     text(x=-110,y=-77,labels=paste("R = ",round(aph,3)), font=2, cex=textinner)
                             # Map
                             load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                             contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                             axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                             axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
                             if (w==1){mtext(3,text="o)",at=-175,cex=lettersize, font=2)   }
                             if (w==2){mtext(3,text="p)",at=-175,cex=lettersize, font=2)   }  
           }        
############################           
           ### ESCALA
           par(mar=c(1,0.5,1,3))
           mat <- matrix(c(1:length(custom_scale)), ncol=length(custom_scale))
           image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(m^-1), cex.main=0.8) 
           axis(4,at=seq(0-0.05,1+0.05, length=6),labels=round(seq(0,0.1,length=6),2), cex.axis=1.2, las=1)

           
           
#### Adg
           par(mar=c(1,0.5,1,0))
############################           
           ##### In situ
           plot(1,1,type="n",yaxt="n",xaxt="n",bty="n")
           
           # SAT
           modelname<-modelnames[1]
           load(paste(p_dir,"2iop_adg/",modelname,".RData",sep=""))
           z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
           z[z<1e-5] <- NA
           load(paste(s_dir,"climatologies_2012_2018/media_seasonal_Adg_2deg_OCCCI_2012_2018.RData",sep=""))
           res <- med 
           landa_sat <- dimnames(res)$l
           sat <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)        
           media <- sat[,,2]
           media[is.na(z)] <- NA
           image2D(x=longitude, y=latitude,media, col=custom_scale, zlim=c(0,0.1), las=1,xaxt="n", yaxt="n",colkey=F,
                   main=expression(paste(a[DG], " (443nm) ",italic("OC-CCI"), sep="")))
                   # Map
                   load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                   contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                   axis(2, at=seq(-90,90,length=5), labels=T,   las=1)
                   axis(1, at=seq(-180,180,length=10), labels=c("",seq(-180,180,length=10)[2],"",seq(-180,180,length=10)[4],"",seq(-180,180,length=10)[6],"",seq(-180,180,length=10)[8],"",seq(-180,180,length=10)[10]),  las=1)
                   mtext(3,text="q)",at=-175,cex=lettersize, font=2) 
                   text(x=75,y=-76,labels="n = 9466", font=2, cex=1.1)
                   
             # Mod
             for (w in 1:length(modelnames)){ 
               modelname<-modelnames[w]
               load(paste(p_dir,"2iop_adg/",modelname,".RData",sep=""))
               z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
               z[z<1e-5] <- NA
          if (w==1){   image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n",colkey=F,
                       main=expression(paste(a[DG], " (450nm) ",italic("EXP-noPPC"), sep="")))}else{
                       image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n",colkey=F,
                       main=expression(paste(a[DG], " (450nm) ",italic("EXP-PPC"), sep="")))}   
                       # Correlation
                       x<-cbind(c(z), c(media))
                       x[x=="-Inf"] <- NA
                       x<-x[complete.cases(x)==TRUE,]
                       adg<- rcorr(x, type="pearson")$r[1,2]
                       text(x=75,y=-77,labels=paste("R = ",round(adg,3)), font=2, cex=textinner)
                             # Map
                             load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
                             contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
                             axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
                             axis(1, at=seq(-180,180,length=10), labels=c("",seq(-180,180,length=10)[2],"",seq(-180,180,length=10)[4],"",seq(-180,180,length=10)[6],"",seq(-180,180,length=10)[8],"",seq(-180,180,length=10)[10]),  las=1)
                             if (w==1){mtext(3,text="r)",at=-175,cex=lettersize, font=2)   }
                             if (w==2){mtext(3,text="s)",at=-175,cex=lettersize, font=2)   }  
               }        
############################ 
           ### ESCALA
           par(mar=c(1,0.5,1,3))
           mat <- matrix(c(1:length(custom_scale)), ncol=length(custom_scale))
           image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(m^-1), cex.main=0.8) 
           axis(4,at=seq(0-0.05,1+0.05, length=6),labels=round(seq(0,0.1,length=6),3), cex.axis=1.2, las=1)
           
      
     dev.off()
      
         
         
         
     
     
     
     
     
     
     modelnames<-experimentos$name[c(30,32)]
     ruta<-experimentos$path[c(30,32)]        
     nombre <- c(expression(paste("EXP-noPPC: ",{{"a*"}}[NAP],"(", lambda, ") [", m^-2," ",mmolC^-1,"]",sep="")),
                 expression(paste("EXP-D: ",{{"a*"}}[NAP],"(", lambda, ") [", m^-2," ",particle^-1,"]",sep="")))     
     
     
#############################
#### FIGURE S8: aNAP surface ( review R1 (2022-06-30) )
############################# 
     
     png(file=paste(path_figures,"FigureS8_Anap_corr.png",sep=""),width = 1200, height = 900, units = "px", pointsize = 20, bg = "white")  
     layout(matrix(c(1:16),ncol=4,byrow=T),widths=c(rep(1,3),0.25))
     par(mar=c(1,0.5,1,0))
     par(oma=c(2,2,2,1))
     layout.show(n=16)
     textinner=1.2
     lettersize=1.2
     ############################         
     
     #### Atot
     ###########################           

     #### Sat
     modelname<-modelnames[1]
     load(paste(p_dir,"2iop_atot/",modelname,".RData",sep=""))
     z <- apply(res,MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
     z[z<0.01] <- NA
     MASCARA <- z
     load(paste(s_dir,"climatologies_2012_2018/media_seasonal_Atot_2deg_OCCCI_2012_2018.RData",sep=""))
     res <- med           
     landa_sat <- dimnames(res)$l
     longitude<-as.numeric(dimnames(res)$x)
     latitude<-as.numeric(dimnames(res)$y)
     sat <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)        
     media <- sat[,,2]
     media[is.na(z)] <- NA
     image2D(x=longitude, y=latitude,media, col=custom_scale, zlim=c(0,0.2), las=1, xaxt="n", yaxt="n",colkey=F,
             main=expression(paste(a[TOT], " (443nm) ",italic("OC-CCI"), sep="")))
     # Map
     load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
     contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
     axis(2, at=seq(-90,90,length=5), labels=T,   las=1)
     axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
     mtext(3,text="a)",at=-175,cex=lettersize, font=2)  
     text(x=75,y=-76,labels="n = 9466", font=2, cex=1.1)
     
     # Mod  
     for (w in 1:length(modelnames)){ 
       modelname<-modelnames[w]
       load(paste(p_dir,"2iop_atot/",modelname,".RData",sep=""))
       z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
       z[z<0.01] <- NA
       if (w==1){ image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.2), las=1, xaxt="n", yaxt="n",colkey=F,
                          main=expression(paste(a[TOT], " (450nm) ",italic("EXP-noPPC"), sep="")))}else{
                            image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.2), las=1, xaxt="n", yaxt="n",colkey=F,
                                    main=expression(paste(a[TOT], " (450nm) ",italic("EXP-D"), sep="")))}
       # Correlation
       x<-cbind(c(z), c(media))
       x[x=="-Inf"] <- NA
       x<-x[complete.cases(x)==TRUE,]
       tot<- rcorr(x, type="pearson")$r[1,2]
       text(x=75,y=-77,labels=paste("R = ",round(tot,3)), font=2, cex=textinner)
       # Map
       load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
       contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
       axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
       axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
       if (w==1){mtext(3,text="b)",at=-175,cex=lettersize, font=2)   }
       if (w==2){mtext(3,text="c)",at=-175,cex=lettersize, font=2)   } 
     }
     ###########################    
     ### ESCALA
     par(mar=c(1,0.5,1,3))
     mat <- matrix(c(1:length(custom_scale)), ncol=length(custom_scale))
     image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(m^-1), cex.main=0.8) 
     axis(4,at=seq(0-0.05,1+0.05, length=6),labels=round(seq(0,0.2,length=6),2), cex.axis=1.2, las=1)
     
     
     #### Anap
     par(mar=c(1,0.5,1,0))
     ############################           
     ##### In situ
     is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
     load(paste(is_dir,"/media_annual_Anap.Rdata",sep=""))
     image2D(x=longitude, y=latitude,matriz_annual[,,1,3], col=custom_scale, zlim=c(0,0.03), las=1, xaxt="n", yaxt="n", colkey=F,
             main=expression(paste(a[NAP], " (450nm) ",italic("in situ"), sep="")))
     # Map
     load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
     contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
     axis(2, at=seq(-90,90,length=5), labels=T,   las=1)
     axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
     mtext(3,text="d)",at=-175,cex=lettersize, font=2)  
     text(x=75,y=-76,labels="n = 502", font=2, cex=1.1)
     

     # Mod  
     for (w in 1:length(modelnames)){ 
       modelname<-modelnames[w]
       load(paste(p_dir,"2iop_apt/",modelname,".RData",sep=""))
       z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
       #z[z<1e-5] <- NA
       z[z<0] <- 0
       if (w==1){image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.03), las=1, xaxt="n", yaxt="n",colkey=F,
                         main=expression(paste(a[NAP], " (450nm) ",italic("EXP-noPPC"), sep="")))}else{
                           image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.03), las=1, xaxt="n", yaxt="n",colkey=F,
                                   main=expression(paste(a[NAP], " (450nm) ",italic("EXP-D"), sep="")))}        
       # Correlation Insitu 450 vs 450
       x<-cbind(c(z), c(matriz_annual[,,1,3]))
       x[x=="-Inf"] <- NA
       x<-x[complete.cases(x)==TRUE,]
       aph<- rcorr(x, type="pearson")$r[1,2  ]
       text(x=-110,y=-77,labels=paste("R = ",round(aph,3)), font=2, cex=textinner)
       # Map
       load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
       contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
       axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
       axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
       if (w==1){mtext(3,text="e)",at=-175,cex=lettersize, font=2)   }
       if (w==2){mtext(3,text="f)",at=-175,cex=lettersize, font=2)   }                
     }        
     ############################           
     ### ESCALA
     par(mar=c(1,0.5,1,3))
     mat <- matrix(c(1:length(custom_scale)), ncol=length(custom_scale))
     image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(m^-1), cex.main=0.8) 
     axis(4,at=seq(0-0.05,1+0.05, length=6),labels=round(seq(0,0.03,length=6),3), cex.axis=1.2, las=1)
     
     
     #### Ap
     par(mar=c(1,0.5,1,0))
     ############################           
     ##### In situ
     is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
     load(paste(is_dir,"/media_annual_Ap.Rdata",sep=""))           
     image2D(x=longitude, y=latitude,matriz_annual[,,1,3], col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n", colkey=F,
             main=expression(paste(a[P], " (450nm) ",italic("in situ"), sep="")))
     # Map
     load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
     contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
     axis(2, at=seq(-90,90,length=5), labels=T,   las=1)
     axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
     mtext(3,text="g)",at=-175,cex=lettersize, font=2)    
     text(x=75,y=-76,labels="n = 622", font=2, cex=1.1)

     # Mod  
     for (w in 1:length(modelnames)){ 
       modelname<-modelnames[w]
       load(paste(p_dir,"2iop_app/",modelname,".RData",sep=""))
       z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
       z[z<1e-5] <- NA
       if (w==1){image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n",colkey=F,
                         main=expression(paste(a[P], " (450nm) ",italic("EXP-noPPC"), sep="")))}else{
                           image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n",colkey=F,
                                   main=expression(paste(a[P], " (450nm) ",italic("EXP-D"), sep="")))}        
       # Correlation Insitu 450 vs 450
       x<-cbind(c(z), c(matriz_annual[,,1,3]))
       x[x=="-Inf"] <- NA
       x<-x[complete.cases(x)==TRUE,]
       aph<- rcorr(x, type="pearson")$r[1,2  ]
       text(x=-110,y=-77,labels=paste("R = ",round(aph,3)), font=2, cex=textinner)
       # Map
       load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
       contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
       axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
       axis(1, at=seq(-180,180,length=10), labels=F,  las=1)
       if (w==1){mtext(3,text="h)",at=-175,cex=lettersize, font=2)   }
       if (w==2){mtext(3,text="i)",at=-175,cex=lettersize, font=2)   }   
     }        
     ############################           
     ### ESCALA
     par(mar=c(1,0.5,1,3))
     mat <- matrix(c(1:length(custom_scale)), ncol=length(custom_scale))
     image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(m^-1), cex.main=0.8) 
     axis(4,at=seq(0-0.05,1+0.05, length=6),labels=round(seq(0,0.1,length=6),3), cex.axis=1.2, las=1)
     

     #### Adg
     par(mar=c(1,0.5,1,0))
     ############################           

     # SAT
     modelname<-modelnames[1]
     load(paste(p_dir,"2iop_adg/",modelname,".RData",sep=""))
     z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
     z[z<1e-5] <- NA
     load(paste(s_dir,"climatologies_2012_2018/media_seasonal_Adg_2deg_OCCCI_2012_2018.RData",sep=""))
     res <- med 
     landa_sat <- dimnames(res)$l
     sat <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)        
     media <- sat[,,2]
     media[is.na(z)] <- NA
     image2D(x=longitude, y=latitude,media, col=custom_scale, zlim=c(0,0.1), las=1,xaxt="n", yaxt="n",colkey=F,
             main=expression(paste(a[DG], " (443nm) ",italic("OC-CCI"), sep="")))
     # Map
     load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
     contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
     axis(2, at=seq(-90,90,length=5), labels=T,   las=1)
     axis(1, at=seq(-180,180,length=10), labels=c("",seq(-180,180,length=10)[2],"",seq(-180,180,length=10)[4],"",seq(-180,180,length=10)[6],"",seq(-180,180,length=10)[8],"",seq(-180,180,length=10)[10]),  las=1)
     mtext(3,text="j)",at=-175,cex=lettersize, font=2) 
     text(x=75,y=-76,labels="n = 9466", font=2, cex=1.1)
     
     # Mod
     for (w in 1:length(modelnames)){ 
       modelname<-modelnames[w]
       load(paste(p_dir,"2iop_adg/",modelname,".RData",sep=""))
       z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
       z[z<1e-5] <- NA
       if (w==1){   image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n",colkey=F,
                            main=expression(paste(a[DG], " (450nm) ",italic("EXP-noPPC"), sep="")))}else{
                              image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.1), las=1, xaxt="n", yaxt="n",colkey=F,
                                      main=expression(paste(a[DG], " (450nm) ",italic("EXP-D"), sep="")))}   
       # Correlation
       x<-cbind(c(z), c(media))
       x[x=="-Inf"] <- NA
       x<-x[complete.cases(x)==TRUE,]
       adg<- rcorr(x, type="pearson")$r[1,2]
       text(x=75,y=-77,labels=paste("R = ",round(adg,3)), font=2, cex=textinner)
       # Map
       load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
       contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
       axis(2, at=seq(-90,90,length=5), labels=F,   las=1)
       axis(1, at=seq(-180,180,length=10), labels=c("",seq(-180,180,length=10)[2],"",seq(-180,180,length=10)[4],"",seq(-180,180,length=10)[6],"",seq(-180,180,length=10)[8],"",seq(-180,180,length=10)[10]),  las=1)
       if (w==1){mtext(3,text="k)",at=-175,cex=lettersize, font=2)   }
       if (w==2){mtext(3,text="l)",at=-175,cex=lettersize, font=2)   }  
     }        
     ############################ 
     ### ESCALA
     par(mar=c(1,0.5,1,3))
     mat <- matrix(c(1:length(custom_scale)), ncol=length(custom_scale))
     image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(m^-1), cex.main=0.8) 
     axis(4,at=seq(0-0.05,1+0.05, length=6),labels=round(seq(0,0.1,length=6),3), cex.axis=1.2, las=1)
     
     
     dev.off()
     
     
     
    
     
     
     modelnames<-experimentos$name[c(21,9)]
     ruta<-experimentos$path[c(21,9)]        
     nombre <- c(expression(paste("EXP-noPPC: ",{{"a*"}}[PH],"(", lambda, ") constant",sep="")),
                 expression(paste("EXP-PPC: ",{{"a*"}}[PH],"(", lambda, ") variable",sep="")))
     
 
           
#########################################################################
#### Taylor diagram with all output EXP-noPPC and EXP-PPC: constant vs variable
#########################################################################
    
## Compute positions first (Taylor 2001), then plot diagram 
     
  ## MODEL vs FIELD          
  ############################             
  posX<-matrix(NA, ncol=2, nrow=12)
  posY<-matrix(NA, ncol=2, nrow=12)  

           for (w in c(1:length(modelnames))){ 
             modelname<-modelnames[w]
             
      #### CHLA log
      ##############             
                   load(paste(p_dir,"1cloro_log/",modelname,".RData",sep=""))
                   mod <- res
                   longitude<-as.numeric(dimnames(res)$x)
                   latitude<-as.numeric(dimnames(res)$y)
                   
                   ##### In situ MAREDAT and others
                   is_dir <- paste(global_path,"Dat_observations/HPLC/grids_pigments",sep="") 
                   load(paste(is_dir,"/TChla_insitu_new2.RData",sep=""))                           
                   longitude<-as.numeric(dimnames(TChla_ALL)$x)
                   latitude<-as.numeric(dimnames(TChla_ALL)$y)        
                   matriz1  <- (apply(TChla_ALL[,,1:2], MARGIN=c("x","y"), mean, na.rm=TRUE))/1e+3 # ng L-1 to ug L-1
      
                   ##### In situ Valente and others
                   is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
                   load(paste(is_dir,"/media_annual_Chla.Rdata",sep=""))
                   dim(res)
                   matriz2  <- apply(res[,,1:2], MARGIN=c("x","y"), mean, na.rm=TRUE) # ug L-1
                   
                   nueva<-abind(matriz1, matriz2, along=0.5)
                   matriz<-apply(nueva,MARGIN=c(2,3), mean, na.rm=TRUE)
      
                   # Positions for Taylor Diagram
                   tabla=cbind(c(log10(matriz)),c(mod))
                   tabla[is.nan(tabla)]<-NA
                   tabla=tabla[complete.cases(tabla),]
                   ref=tabla[,1]
                   model=tabla[,2]
                   sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                   R <- cor(ref,model, method= c('pearson'))
                   posX[1,w] <- sd.f * R
                   posY[1,w] <- sd.f * sin(acos(R))

      ##############
                   
      #### PPC:TChla
      ##############             
                   # In situ
                   is_dir <- paste(global_path,"Dat_observations/HPLC/grids_pigments",sep="") 
                   load(paste(is_dir,"/PPCinsitu_new3.RData",sep=""))  
                   matriz1  <- apply(PPCinsitu_ALL[,,1:2], MARGIN=c("x","y"), mean, na.rm=TRUE)
                   matriz1[matriz1=="Inf" | matriz1=="-Inf"] <- NA
                   
                   # MODEL
                   load(paste(p_dir,"1ppc_total/",modelname,".RData",sep=""))
                   mod <- res
                   
                   # Positions for Taylor Diagram
                           tabla=cbind(c(matriz1),c(mod))
                           tabla[is.nan(tabla)]<-NA
                           tabla=tabla[complete.cases(tabla),]
                           ref=tabla[,1]
                           model=tabla[,2]
                           sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                           R <- cor(ref,model, method= c('pearson'))
                           posX[2,w] <- sd.f * R
                           posY[2,w] <- sd.f * sin(acos(R)) 
      ###############             
                   
      #### Reflectance
      ################           
                   # MODEL
                   load(paste(p_dir,"1reflectance/",modelname,".RData",sep=""))
                   mod <- res
                   z <- mod[,,12]

                   ##### In situ Rrs
                   is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
                   load(paste(is_dir,"/media_annual_Rrs.Rdata",sep=""))           
      
                   # Positions for Taylor Diagram
                                 tabla=cbind(c(matriz_annual[,,1]),c(mod[,,1]))
                                 tabla[is.nan(tabla)]<-NA
                                 tabla=tabla[complete.cases(tabla),]
                                 ref=tabla[,1]
                                 model=tabla[,2]
                                 sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                                 R <- cor(ref,model, method= c('pearson'))
                                 posX[3,w] <- sd.f * R
                                 posY[3,w] <- sd.f * sin(acos(R)) 
                                 
                   # Positions for Taylor Diagram
                                 tabla=cbind(c(matriz_annual[,,3]),c(mod[,,3]))
                                 tabla[is.nan(tabla)]<-NA
                                 tabla=tabla[complete.cases(tabla),]
                                 ref=tabla[,1]
                                 model=tabla[,2]
                                 sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                                 R <- cor(ref,model, method= c('pearson'))
                                 posX[4,w] <- sd.f * R
                                 posY[4,w] <- sd.f * sin(acos(R))  
                                 
                   # Positions for Taylor Diagram
                                 tabla=cbind(c(matriz_annual[,,4]),c(mod[,,4]))
                                 tabla[is.nan(tabla)]<-NA
                                 tabla=tabla[complete.cases(tabla),]
                                 ref=tabla[,1]
                                 model=tabla[,2]
                                 sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                                 R <- cor(ref,model, method= c('pearson'))
                                 posX[5,w] <- sd.f * R
                                 posY[5,w] <- sd.f * sin(acos(R))  

                   # Positions for Taylor Diagram
                                 tabla=cbind(c(matriz_annual[,,5]),c(mod[,,5]))
                                 tabla[is.nan(tabla)]<-NA
                                 tabla=tabla[complete.cases(tabla),]
                                 ref=tabla[,1]
                                 model=tabla[,2]
                                 sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                                 R <- cor(ref,model, method= c('pearson'))
                                 posX[6,w] <- sd.f * R
                                 posY[6,w] <- sd.f * sin(acos(R))  

                   # Positions for Taylor Diagram
                                 tabla=cbind(c(matriz_annual[,,7]),c(mod[,,7]))
                                 tabla[is.nan(tabla)]<-NA
                                 tabla=tabla[complete.cases(tabla),]
                                 ref=tabla[,1]
                                 model=tabla[,2]
                                 sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                                 R <- cor(ref,model, method= c('pearson'))
                                 posX[7,w] <- sd.f * R
                                 posY[7,w] <- sd.f * sin(acos(R))  

                   # Positions for Taylor Diagram
                                 tabla=cbind(c(matriz_annual[,,12]),c(mod[,,12]))
                                 tabla[is.nan(tabla)]<-NA
                                 tabla=tabla[complete.cases(tabla),]
                                 ref=tabla[,1]
                                 model=tabla[,2]
                                 sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                                 R <- cor(ref,model, method= c('pearson'))
                                 posX[8,w] <- sd.f * R
                                 posY[8,w] <- sd.f * sin(acos(R))

      ################             
                   
      #### IOPs
      #############             
                   
               #### Aph
                   load(paste(p_dir,"2iop_aph/",modelname,".RData",sep=""))
                   z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
                   z[z<1e-5] <- NA
                   ##### In situ
                   is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
                   load(paste(is_dir,"/media_annual_Aph.Rdata",sep=""))           

                   # Positions for Taylor Diagram
                                 tabla=cbind(c(matriz_annual[,,1,3]),c(z))
                                 tabla[is.nan(tabla)]<-NA
                                 tabla=tabla[complete.cases(tabla),]
                                 ref=tabla[,1]
                                 model=tabla[,2]
                                 sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                                 R <- cor(ref,model, method= c('pearson'))
                                 posX[9,w] <- sd.f * R
                                 posY[9,w] <- sd.f * sin(acos(R))  
               #### Anap
                   load(paste(p_dir,"2iop_apt/",modelname,".RData",sep=""))
                   z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
                   z[z<1e-5] <- NA
                   ##### In situ
                   is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
                   load(paste(is_dir,"/media_annual_Anap.Rdata",sep=""))           
                   
                   # Positions for Taylor Diagram
                                 tabla=cbind(c(matriz_annual[,,1,3]),c(z))
                                 tabla[is.nan(tabla)]<-NA
                                 tabla=tabla[complete.cases(tabla),]
                                 ref=tabla[,1]
                                 model=tabla[,2]
                                 sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                                 R <- cor(ref,model, method= c('pearson'))
                                 posX[10,w] <- sd.f * R
                                 posY[10,w] <- sd.f * sin(acos(R))              
               #### Ap
                   load(paste(p_dir,"2iop_app/",modelname,".RData",sep=""))
                   z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
                   z[z<1e-5] <- NA
                   ##### In situ
                   is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
                   load(paste(is_dir,"/media_annual_Ap.Rdata",sep=""))           
                   
                   # Positions for Taylor Diagram
                                 tabla=cbind(c(matriz_annual[,,1,3]),c(z))
                                 tabla[is.nan(tabla)]<-NA
                                 tabla=tabla[complete.cases(tabla),]
                                 ref=tabla[,1]
                                 model=tabla[,2]
                                 sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                                 R <- cor(ref,model, method= c('pearson'))
                                 posX[11,w] <- sd.f * R
                                 posY[11,w] <- sd.f * sin(acos(R))              
              #### Acdom
                   load(paste(p_dir,"2iop_acd/",modelname,".RData",sep=""))
                   z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
                   z[z<1e-5] <- NA
                   ##### In situ
                   is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
                   load(paste(is_dir,"/media_annual_Acdom.Rdata",sep=""))           
                   
                   # Positions for Taylor Diagram
                                 tabla=cbind(c(matriz_annual[,,1,3]),c(z))
                                 tabla[is.nan(tabla)]<-NA
                                 tabla=tabla[complete.cases(tabla),]
                                 ref=tabla[,1]
                                 model=tabla[,2]
                                 sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                                 R <- cor(ref,model, method= c('pearson'))
                                 posX[12,w] <- sd.f * R
                                 posY[12,w] <- sd.f * sin(acos(R))
      ################             
                 
                   } # end loop modelnames              
                 
         
  ## MODEL vs SAT
  ############################ 
  possX<-matrix(NA, ncol=2, nrow=10)
  possY<-matrix(NA, ncol=2, nrow=10) 

           for (w in c(1:length(modelnames))){ 

      #### CHLA log
      #################  
                   modelname<-modelnames[w]
                   load(paste(p_dir,"1cloro_log/",modelname,".RData",sep=""))
                   mod <- res
      
                   # SAT
                   load(paste(s_dir,"climatologies_2012_2018/media_anual_log_chl_OCCCI_2012_2018.RData", sep=""))
                   media <- apply(res, MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
                   media[is.na(mod)] <- NA
                   longitude<-as.numeric(dimnames(res)$x)
                   latitude<-(-as.numeric(dimnames(res)$y))
                   
                  # Positions for Taylor Diagram
                   tabla=cbind(c(media),c(mod))
                   tabla[is.nan(tabla)]<-NA
                   tabla=tabla[complete.cases(tabla),]
                   ref=tabla[,1]
                   model=tabla[,2]
                   sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                   R <- cor(ref,model, method= c('pearson'))
                   possX[1,w] <- sd.f * R
                   possY[1,w] <- sd.f * sin(acos(R))
      #################               
                     
      
      #### Reflectance
      ################             
               # MODEL
               load(paste(p_dir,"1reflectance/",modelname,".RData",sep=""))
               mod <- res
               z <- mod[,,12]

               # SAT
               landas_sat <- c("Rrs_412","Rrs_443","Rrs_490","Rrs_510","Rrs_555","Rrs_670")
               load(paste(s_dir,"climatologies_2012_2018/media_seasonal_Rrs_2deg_OCCCI_2012_2018.RData",sep=""))
               res <- med
               
               media412 <- apply(res[,,,1], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
               media412[is.na(z)] <- NA
               # Positions for Taylor Diagram
                     tabla=cbind(c(media412),c(mod[,,1]))
                     tabla[is.nan(tabla)]<-NA
                     tabla=tabla[complete.cases(tabla),]
                     ref=tabla[,1]
                     model=tabla[,2]
                     sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                     R <- cor(ref,model, method= c('pearson'))
                     possX[2,w] <- sd.f * R
                     possY[2,w] <- sd.f * sin(acos(R))                 
                     
                 media443 <- apply(res[,,,2], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
                 media443[is.na(z)] <- NA
                 # Positions for Taylor Diagram
                             tabla=cbind(c(media443),c(mod[,,3]))
                             tabla[is.nan(tabla)]<-NA
                             tabla=tabla[complete.cases(tabla),]
                             ref=tabla[,1]
                             model=tabla[,2]
                             sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                             R <- cor(ref,model, method= c('pearson'))
                             possX[3,w] <- sd.f * R
                             possY[3,w] <- sd.f * sin(acos(R))
                 
                   media490 <- apply(res[,,,3], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
                   media490[is.na(z)] <- NA
                   # Positions for Taylor Diagram
                               tabla=cbind(c(media490),c(mod[,,4]))
                               tabla[is.nan(tabla)]<-NA
                               tabla=tabla[complete.cases(tabla),]
                               ref=tabla[,1]
                               model=tabla[,2]
                               sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                               R <- cor(ref,model, method= c('pearson'))
                               possX[4,w] <- sd.f * R
                               possY[4,w] <- sd.f * sin(acos(R))                 
                               
                       media510 <- apply(res[,,,4], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
                       media510[is.na(z)] <- NA
                       # Positions for Taylor Diagram
                                     tabla=cbind(c(media510),c(mod[,,5]))
                                     tabla[is.nan(tabla)]<-NA
                                     tabla=tabla[complete.cases(tabla),]
                                     ref=tabla[,1]
                                     model=tabla[,2]
                                     sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                                     R <- cor(ref,model, method= c('pearson'))
                                     possX[5,w] <- sd.f * R
                                     possY[5,w] <- sd.f * sin(acos(R))  
                                     
                          media555 <- apply(res[,,,5], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
                          media555[is.na(z)] <- NA
                          # Positions for Taylor Diagram
                                    tabla=cbind(c(media555),c(mod[,,7]))
                                    tabla[is.nan(tabla)]<-NA
                                    tabla=tabla[complete.cases(tabla),]
                                    ref=tabla[,1]
                                    model=tabla[,2]
                                    sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                                    R <- cor(ref,model, method= c('pearson'))
                                    possX[6,w] <- sd.f * R
                                    possY[6,w] <- sd.f * sin(acos(R)) 
                                    
                            media670 <- apply(res[,,,6], MARGIN=c("x","y"), FUN=mean, na.rm=TRUE)        
                            media670[is.na(z)] <- NA
                            # Positions for Taylor Diagram
                                        tabla=cbind(c(media670),c(mod[,,12]))
                                        tabla[is.nan(tabla)]<-NA
                                        tabla=tabla[complete.cases(tabla),]
                                        ref=tabla[,1]
                                        model=tabla[,2]
                                        sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                                        R <- cor(ref,model, method= c('pearson'))
                                        possX[7,w] <- sd.f * R
                                        possY[7,w] <- sd.f * sin(acos(R))       
      ################                     
            
      #### IOPs
      ##############                     
           #### Atot
               load(paste(p_dir,"2iop_atot/",modelname,".RData",sep=""))
               z<- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
               z[z<0.01] <- NA
               MASCARA<-z
                    # SAT
                    load(paste(s_dir,"climatologies_2012_2018/media_seasonal_Atot_2deg_OCCCI_2012_2018.RData",sep=""))
                    res <- med
                    landa_sat <- dimnames(res)$l
                    longitude<-as.numeric(dimnames(res)$x)
                    latitude<-(-as.numeric(dimnames(res)$y))
                    sat <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)        
                    media <- sat[,,2]
                    media[is.na(z)] <- NA
                    # Positions for Taylor Diagram
                              tabla=cbind(c(media),c(z))
                              tabla[is.nan(tabla)]<-NA
                              tabla=tabla[complete.cases(tabla),]
                              ref=tabla[,1]
                              model=tabla[,2]
                              sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                              R <- cor(ref,model, method= c('pearson'))
                              possX[8,w] <- sd.f * R
                              possY[8,w] <- sd.f * sin(acos(R))                            
            #### Aph
                load(paste(p_dir,"2iop_aph/",modelname,".RData",sep=""))
                z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
                z[z<1e-5] <- NA
                      # SAT
                      load(paste(s_dir,"climatologies_2012_2018/media_seasonal_Aph_2deg_OCCCI_2012_2018.RData",sep=""))
                      res <- med
                      landa_sat <- dimnames(res)$l
                      sat <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)        
                      media <- sat[,,2]
                      media[is.na(z)] <- NA
                      # Positions for Taylor Diagram
                               tabla=cbind(c(media),c(z))
                               tabla[is.nan(tabla)]<-NA
                               tabla=tabla[complete.cases(tabla),]
                               ref=tabla[,1]
                               model=tabla[,2]
                               sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                               R <- cor(ref,model, method= c('pearson'))
                               possX[9,w] <- sd.f * R
                               possY[9,w] <- sd.f * sin(acos(R))   
                               
            #### Adg
                 load(paste(p_dir,"2iop_adg/",modelname,".RData",sep=""))
                 z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
                 z[z<1e-5] <- NA
                       # SAT
                       load(paste(s_dir,"climatologies_2012_2018/media_seasonal_Adg_2deg_OCCCI_2012_2018.RData",sep=""))
                       res <- med
                       landa_sat <- dimnames(res)$l
                       sat <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)        
                       media <- sat[,,2]
                       media[is.na(z)] <- NA
                       # Positions for Taylor Diagram
                                 tabla=cbind(c(media),c(z))
                                 tabla[is.nan(tabla)]<-NA
                                 tabla=tabla[complete.cases(tabla),]
                                 ref=tabla[,1]
                                 model=tabla[,2]
                                 sd.f <- (sd(model)/(mean(model)))/(sd(ref)/mean(ref)) 
                                 R <- cor(ref,model, method= c('pearson'))
                                 possX[10,w] <- sd.f * R
                                 possY[10,w] <- sd.f * sin(acos(R))       
      ##############            
                            
                 } # end loop modelnames              

  
  
   
###########################################################
## FIGURE 12: Taylor diagrams Model vs In situ and Satellite (before review was Figure 4)
###########################################################

  png(file=paste(path_figures,"Figure12_Taylor_diagrams_split.png",sep=""),width = 1210, height = 730, units = "px", pointsize = 24, bg = "white") 
    par(mfcol=c(2,3)) 
    par(family="")
    par(mar=c(5,5,2,2))  
    par(oma=c(1,2,2,1))                 
    layout.show(n=6)
    lettersize=3
    lettersizeS=1.9
    labsize=0.8
    legsize=1.0
    cexpoints=1.4
  
    # BGC products In situ
    ############################         
    taylor.diagram.mine(ref=1, model=1, ref.sd=TRUE,
                        normalize=T,show.gamma=T, gamma.col="limegreen", gamma.cex=0.7,
                        pch=19, col=colores[w], main="", las=1,
                        mar=c(3,3,3,3), pcex=0.01,
                        sd.arcs=c(0.4,0.6,0.8,1.2), font=2, maxsd=1.5)
    points(1,0, pch=21, cex=cexpoints, font=2)
    text(0.98,0.05,pos=4, labels="Ref")
              points(posX[1,1],posY[1,1], pch="a", col=colores[1], cex=cexpoints, font=2)  
              points(posX[1,2],posY[1,2], pch="a", col=colores[2], cex=cexpoints, font=2)   
                     #segments(x0=posX[1,1],x1=posX[1,2],y0=posY[1,1],y1=posY[1,2], lwd=2)
              points(posX[2,1],posY[2,1], pch="b", col=colores[1], cex=cexpoints, font=2)  
              points(posX[2,2],posY[2,2], pch="b", col=colores[2], cex=cexpoints, font=2)    
                     #segments(x0=posX[2,1],x1=posX[2,2],y0=posY[2,1],y1=posY[2,2], lwd=2)
              ## Legends
              legend(x=1.0, y=1.8,  xjust=0, 
                     legend=c(expression(paste(log[10],"(TChla)", sep="")),"PPC:TChla"), y.intersp=1.1,
                     pch=c("a","b"), bty="n", cex=legsize, pt.cex =1.0, pt.lwd = 2)
              text(x=-0.15,y=1.60, labels="a)",cex=lettersizeS, pos=4, font=1)        
    ############################          
              
                        
    # BGC products Satellite
    ############################          
    taylor.diagram.mine(ref=1, model=1, ref.sd=TRUE,
                        normalize=T,show.gamma=T, gamma.col="limegreen", gamma.cex=0.7,
                        pch=19, col=colores[w], main="", las=1,
                        mar=c(3,3,3,3), pcex=0.01,
                        sd.arcs=c(0.4,0.6,0.8,1.2), font=2, maxsd=1.5)
    points(1,0, pch=21, cex=cexpoints, font=2)
    text(0.98,0.05,pos=4, labels="Ref")
              points(possX[1,1],possY[1,1], pch="a", col=colores[1], cex=cexpoints, font=2)  
              points(possX[1,2],possY[1,2], pch="a", col=colores[2], cex=cexpoints, font=2)     
                #segments(x0=possX[1,1],x1=possX[1,2],y0=possY[1,1],y1=possY[1,2], lwd=2)
              ## Legends 
              #legend(x="topright", legend=c("TChla","Rrs412","Rrs443","Rrs490","Rrs510","Rrs555","Rrs670","Atot","Aph","Adg"),pch=simbolos[c(1,3:10,14)], bty="n", cex=0.75)
              legend(x=1.0, y=1.8,  xjust=0,  #0.66,y=1.5,y.intersp=1.1,
                     legend=expression(paste(log[10],"(TChla)", sep="")),
                     pch="a", bty="n", cex=legsize, pt.cex =1.0, pt.lwd = 2)
              text(x=-0.15,y=1.60, labels="b)",cex=lettersizeS, pos=4, font=1)        
    ############################          
      
    # Rrs in situ
    ############################            
    taylor.diagram.mine(ref=1, model=1, ref.sd=TRUE,
                        normalize=T,show.gamma=T, gamma.col="limegreen", gamma.cex=0.7,
                        pch=19, col=colores[w], main="", las=1,
                        mar=c(3,3,3,3), pcex=0.01,
                        sd.arcs=c(0.4,0.6,0.8,1.2), font=2, maxsd=1.2)  
    points(1,0, pch=21, cex=cexpoints, font=2)
    text(0.98,0.05,pos=4, labels="Ref")
                        points(posX[3,1],posY[3,1], pch="a", col=colores[1], cex=cexpoints, font=2)  
                        points(posX[3,2],posY[3,2], pch="a", col=colores[2], cex=cexpoints, font=2)     
                          #segments(x0=posX[3,1],x1=posX[3,2],y0=posY[3,1],y1=posY[3,2], lwd=2)
                        points(posX[4,1],posY[4,1], pch="b", col=colores[1], cex=cexpoints, font=2)  
                        points(posX[4,2],posY[4,2], pch="b", col=colores[2], cex=cexpoints, font=2)     
                          #segments(x0=posX[4,1],x1=posX[4,2],y0=posY[4,1],y1=posY[4,2], lwd=2) 
                        points(posX[5,1],posY[5,1], pch="c", col=colores[1], cex=cexpoints, font=2)  
                        points(posX[5,2],posY[5,2], pch="c", col=colores[2], cex=cexpoints, font=2)     
                          #segments(x0=posX[5,1],x1=posX[5,2],y0=posY[5,1],y1=posY[5,2], lwd=2)
                        points(posX[6,1],posY[6,1], pch="d", col=colores[1], cex=cexpoints, font=2)  
                        points(posX[6,2],posY[6,2], pch="d", col=colores[2], cex=cexpoints, font=2)  
                          #segments(x0=posX[6,1],x1=posX[6,2],y0=posY[6,1],y1=posY[6,2], lwd=2)
                        points(posX[7,1],posY[7,1], pch="e", col=colores[1], cex=cexpoints, font=2)  
                        points(posX[7,2],posY[7,2], pch="e", col=colores[2], cex=cexpoints, font=2)     
                          #segments(x0=posX[7,1],x1=posX[7,2],y0=posY[7,1],y1=posY[7,2], lwd=2)
                        points(posX[8,1],posY[8,1], pch="f", col=colores[1], cex=cexpoints, font=2)  
                        points(posX[8,2],posY[8,2], pch="f", col=colores[2], cex=cexpoints, font=2)        
                          #segments(x0=posX[8,1],x1=posX[8,2],y0=posY[8,1],y1=posY[8,2], lwd=2)
                        legend(x=0.9, y=1.5,  xjust=0,  # 1.27,y=1.5,y.intersp=1.1,
                               legend=c(expression(R[RS]("400nm")),expression(R[RS]("450nm")),expression(R[RS]("475nm")),
                                        expression(R[RS]("500nm")),expression(R[RS]("550nm")),expression(R[RS]("675nm"))),
                               pch=c("a","b","c","d","e","f"), bty="n", cex=legsize, pt.cex =1.0, pt.lwd = 2)
                        text(x=-0.12,y=1.30, labels="c)",cex=lettersizeS, pos=4, font=1)        
    ############################                    
                                                
    # Rrs satellite
    ############################            
    taylor.diagram.mine(ref=1, model=1, ref.sd=TRUE,
                        normalize=T,show.gamma=T, gamma.col="limegreen", gamma.cex=0.7,
                        pch=19, col=colores[w], main="", las=1,
                        mar=c(3,3,3,3), pcex=0.01,
                        sd.arcs=c(0.4,0.6,0.8,1.2), font=2, maxsd=1.2)   
    points(1,0, pch=21, cex=cexpoints, font=2)
    text(0.98,0.05,pos=4, labels="Ref")
                        points(possX[2,1],possY[2,1], pch="a", col=colores[1], cex=cexpoints, font=2)  
                        points(possX[2,2],possY[2,2], pch="a", col=colores[2], cex=cexpoints, font=2)  
                          #segments(x0=possX[2,1],x1=possX[2,2],y0=possY[2,1],y1=possY[2,2], lwd=2)
                        points(possX[3,1],possY[3,1], pch="b", col=colores[1], cex=cexpoints, font=2)  
                        points(possX[3,2],possY[3,2], pch="b", col=colores[2], cex=cexpoints, font=2)      
                          #segments(x0=possX[3,1],x1=possX[3,2],y0=possY[3,1],y1=possY[3,2], lwd=2)
                        points(possX[4,1],possY[4,1], pch="c", col=colores[1], cex=cexpoints, font=2)  
                        points(possX[4,2],possY[4,2], pch="c", col=colores[2], cex=cexpoints, font=2)     
                          #segments(x0=possX[4,1],x1=possX[4,2],y0=possY[4,1],y1=possY[4,2], lwd=2)
                        points(possX[5,1],possY[5,1], pch="d", col=colores[1], cex=cexpoints, font=2)  
                        points(possX[5,2],possY[5,2], pch="d", col=colores[2], cex=cexpoints, font=2)     
                          #segments(x0=possX[5,1],x1=possX[5,2],y0=possY[5,1],y1=possY[5,2], lwd=2)
                        points(possX[6,1],possY[6,1], pch="e", col=colores[1], cex=cexpoints, font=2)  
                        points(possX[6,2],possY[6,2], pch="e", col=colores[2], cex=cexpoints, font=2)  
                          #segments(x0=possX[6,1],x1=possX[6,2],y0=possY[6,1],y1=possY[6,2], lwd=2)
                        points(possX[7,1],possY[7,1], pch="f", col=colores[1], cex=cexpoints, font=2)  
                        points(possX[7,2],possY[7,2], pch="f", col=colores[2], cex=cexpoints, font=2)     
                          #segments(x0=possX[7,1],x1=possX[7,2],y0=possY[7,1],y1=possY[7,2], lwd=2)
                        
                        legend(x=0.9, y=1.5,  xjust=0, #=1.25,y=1.5,y.intersp=1.1,
                               #       legend=c("Rrs412","Rrs443","Rrs490","Rrs510","Rrs555","Rrs670"),
                               legend=c(expression(R[RS]("412nm")),expression(R[RS]("443nm")),
                                        expression(R[RS]("490nm")),expression(R[RS]("510nm")),
                                        expression(R[RS]("555nm")),expression(R[RS]("670nm"))),
                               pch=c("a","b","c","d","e","f"), bty="n", cex=legsize, pt.cex =1.0, pt.lwd = 2)
                        text(x=-0.12,y=1.30, labels="d)",cex=lettersizeS, pos=4, font=1)        
    ############################
                        
    # IOPs in situ
    ############################            
    taylor.diagram.mine(ref=1, model=1, ref.sd=TRUE,
                        normalize=T,show.gamma=T, gamma.col="limegreen", gamma.cex=0.7,
                        pch=19, col=colores[w], main="", las=1,
                        mar=c(3,3,3,3), pcex=0.01,
                        sd.arcs=c(0.4,0.6,0.8,1.2), font=2, maxsd=1.2)   
    points(1,0, pch=21, cex=cexpoints, font=2)
    text(0.98,0.05,pos=4, labels="Ref")
                    points(posX[9,1],posY[9,1], pch="a", col=colores[1], cex=cexpoints, font=2)  
                    points(posX[9,2],posY[9,2], pch="a", col=colores[2], cex=cexpoints, font=2)     
                     #segments(x0=posX[9,1],x1=posX[9,2],y0=posY[9,1],y1=posY[9,2], lwd=2)
                    
                    points(posX[10,1],posY[10,1], pch="c", col=colores[1], cex=cexpoints, font=2)  
                    points(posX[10,2],posY[10,2], pch="c", col=colores[2], cex=cexpoints, font=2)     
                     #segments(x0=posX[10,1],x1=posX[10,2],y0=posY[10,1],y1=posY[10,2], lwd=2)
                    
                    points(posX[11,1],posY[11,1], pch="b", col=colores[1], cex=cexpoints, font=2)  
                    points(posX[11,2],posY[11,2], pch="b", col=colores[2], cex=cexpoints, font=2)     
                     #segments(x0=posX[11,1],x1=posX[11,2],y0=posY[11,1],y1=posY[11,2], lwd=2)
                    
                    points(posX[12,1],posY[12,1], pch="d", col=colores[1], cex=cexpoints, font=2)  
                    points(posX[12,2],posY[12,2], pch="d", col=colores[2], cex=cexpoints, font=2)  
                     #segments(x0=posX[12,1],x1=posX[12,2],y0=posY[12,1],y1=posY[12,2], lwd=2)
                    
                    legend(x=0.85, y=1.5,  xjust=0, # 0.89,y=1.5,y.intersp=1.1,
                           legend=c(expression(a[PH]("450nm")),
                                    expression(a[P]("450nm")),
                                    expression(a[NAP]("450nm")),
                                    expression(a[CDOM]("450nm"))),
                           pch=c("a","b","c","d"), bty="n", cex=legsize, pt.cex =1.0, pt.lwd = 2)
                    text(x=-0.12,y=1.30, labels="e)",cex=lettersizeS, pos=4, font=1)        
    ############################                
                    
    # IOPs satellite
    ############################           
    taylor.diagram.mine(ref=1, model=1, ref.sd=TRUE,
                        normalize=T,show.gamma=T, gamma.col="limegreen", gamma.cex=0.7,
                        pch=19, col=colores[w], main="", las=1,
                        mar=c(3,3,3,3), pcex=0.01,
                        sd.arcs=c(0.4,0.6,0.8,1.2), font=2, maxsd=1.2)    
    points(1,0, pch=21, cex=cexpoints, font=2)
    text(0.98,0.05,pos=4, labels="Ref")
                      points(possX[8,1],possY[8,1], pch="f", col=colores[1], cex=cexpoints, font=2)  
                      points(possX[8,2],possY[8,2], pch="f", col=colores[2], cex=cexpoints, font=2)      
                        #segments(x0=possX[8,1],x1=possX[8,2],y0=possY[8,1],y1=possY[8,2], lwd=2)
                      
                      points(possX[9,1],possY[9,1], pch="a", col=colores[1], cex=cexpoints, font=2)  
                      points(possX[9,2],possY[9,2], pch="a", col=colores[2], cex=cexpoints, font=2)     
                        #segments(x0=possX[9,1],x1=possX[9,2],y0=possY[9,1],y1=possY[9,2], lwd=2)
                      
                      points(possX[10,1],possY[10,1], pch="e", col=colores[1], cex=cexpoints, font=2)  
                      points(possX[10,2],possY[10,2], pch="e", col=colores[2], cex=cexpoints, font=2)     
                        #segments(x0=possX[10,1],x1=possX[10,2],y0=possY[10,1],y1=possY[10,2], lwd=2)
                      
                      legend(x=0.85, y=1.5,  xjust=0, #0.89,y=1.5, y.intersp=1.1,
                             legend=c(expression(a[PH]("443nm")),
                                      expression(a[DG]("443nm")),
                                      expression(a[TOT]("443nm"))),
                             pch=c("a","e","f"), bty="n", cex=legsize, pt.cex =1.0, pt.lwd = 2)
                      text(x=-0.12,y=1.30, labels="f)",cex=lettersizeS, pos=4, font=1)        
    ############################
                      
## Labels
mtext(1,at=seq(0.16,0.83,length=3),line=-0.5,outer=T,text="normalized st. dev.", cex=labsize)
mtext(2,at=c(0.25,0.75),  line=-0.5,outer=T,text="normalized st. dev.", cex=labsize)
mtext(4,at=c(0.2,0.7),    line=-1.2,outer=T,text="Correlation", cex=labsize)

############################                  
dev.off()           
           
              






#########################
#### FIGURE 4 Aph surface
######################### 

png(file=paste(path_figures,"Figure4_Aph_surface.png",sep=""), width = 1210, height = 730, units = "px", pointsize = 20, bg = "white")  
        layout(matrix(c(1,2,5,3,4,5),ncol=3,byrow=T),widths=c(1,1,0.18))
        par(mar=c(2,2,3,1))
        par(oma=c(3,3,1,1))
        layout.show(n=5)
        textinner=1.2
        lettersize=2
        cexaxis=1.2
        labsize=1.4
#############################         

#### Sat
        modelname<-modelnames[1]
        load(paste(p_dir,"2iop_aph/",modelname,".RData",sep=""))
        z <- apply(res,MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
        z[z<1e-5] <- NA
        MASCARA<- z
        load(paste(s_dir,"climatologies_2012_2018/media_seasonal_Aph_2deg_OCCCI_2012_2018.RData",sep=""))
        res <- med            
        landa_sat <- dimnames(res)$l
        sat <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)        
        media <- sat[,,2]
        media[is.na(z)] <- NA
        image2D(x=longitude, y=latitude,  media, col=custom_scale, zlim=c(0,0.1),
                las=1,xaxt="n", yaxt="n",colkey=F, cex.axis=cexaxis, cex.main=labsize, 
                main=expression(paste(a[PH], " (443nm) OC-CCI 2012-2018",sep="")))  
              # Map
              load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
              contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
              axis(2, at=seq(-90,90,length=10),   cex.axis=cexaxis,   las=1)
              axis(1, at=seq(-160,160,length=9),  cex.axis=cexaxis,   las=1)     
              mtext(3,at=-170,line=0.5, text="a)",cex=lettersize, outer=F, font=1) 
              text(x=82,y=-78, labels="n = 9466", font=2, cex=1.1)

##### In situ
        is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
        load(paste(is_dir,"/media_annual_Aph.Rdata",sep=""))  
        latitude<-seq(-89,89,by=2)
        image2D(x=longitude, y=latitude,matriz_annual[,,1,3], col=custom_scale, cex.main=labsize,
                zlim=c(0,0.1), las=1, xaxt="n", yaxt="n", colkey=F, main=expression(paste(a[PH], " (450nm) in situ 1997-2019",sep="")), cex.axis=cexaxis)
              # Map
              load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
              contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
              axis(2, at=seq(-90,90,length=10),   cex.axis=cexaxis,   las=1)
              axis(1, at=seq(-160,160,length=9), cex.axis=cexaxis,   las=1)
              mtext(3,at=-170,line=0.5, text="b)",cex=lettersize, outer=F, font=1) 
              text(x=82,y=-78,labels="n = 457", font=2, cex=1.1)
              
#### Mod  
        for (w in 1:length(modelnames)){ 
          modelname<-modelnames[w]
          load(paste(p_dir,"2iop_aph/",modelname,".RData",sep=""))
          z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
          z[z<1e-5] <- NA
          image2D(x=longitude, y=latitude,z, col=custom_scale, zlim=c(0,0.1), cex.main=labsize,
                  las=1, xaxt="n", yaxt="n", main=nombre[w], colkey=F, cex.axis=cexaxis)  #, zlim=c(-1.8,1.8)
          # Map
          load(paste(global_path,"Misc/batimetria_world_degree6.Rdata", sep=""))
          contour(x=lon, y=lat, z=prof,levels=c(0), lwd=0.4, ylab="", xlab="", drawlabels=FALSE, las=1, xaxs="i", yaxs="i", col="grey50",xaxt="n", yaxt="n", bty="n", add=T)
          axis(2, at=seq(-90,90,length=10),   cex.axis=cexaxis,   las=1)
          axis(1, at=seq(-160,160,length=9), cex.axis=cexaxis,   las=1)
          
          # Correlation Insitu 450 vs 450
          x<-cbind(c(z), c(matriz_annual[,,1,3]))
          x[x=="-Inf"] <- NA
          x<-x[complete.cases(x)==TRUE,]
          aph<- rcorr(x, type="pearson")$r[1,2]
          ae <- (sum(x[,1]-x[,2]))/nrow(x)   
          legend(x="topright",legend=c(paste("R = ",round(aph,3)), paste("Bias = ",round(ae,4))),cex=textinner,
                 bty="n",text.font=2)
          
          # Correlation
          x<-cbind(c(z), c(media))
          x[x=="-Inf"] <- NA
          x<-x[complete.cases(x)==TRUE,]
          aph<- rcorr(x, type="pearson")$r[1,2]
          ae <- (sum(x[,1]-x[,2]))/nrow(x)
          legend(x="topleft",legend=c(paste("R = ",round(aph,3)), paste("Bias = ",round(ae,4))),cex=textinner,
                 bty="n",text.font=2)
          if (w==1){ mtext(3,at=-170,line=0.5, text="c)",cex=lettersize, outer=F, font=1) }
          if (w==2){ mtext(3,at=-170,line=0.5, text="d)",cex=lettersize, outer=F, font=1) }
        }        
#############################           
          ###ESCALA
          par(mar=c(2,0.5,3,4))
          mat <- matrix(c(1:length(custom_scale)), ncol=length(custom_scale))
          image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(m^-1), cex.main=1.2) 
          axis(4,at=seq(0-0.055,1+0.055, length=11),labels=round(seq(0,0.1,length=11),3), cex.axis=1.2, las=1)

          ## Labels
          mtext(1,at=c(0.23,0.69),line=1,outer=T,text=expression(paste(degree, "E")))
          mtext(2,at=c(0.25,0.75),line=0.8,outer=T,text=expression(paste(degree, "N")), las=1)

dev.off()










###################################
### FIGURE 6  Gradient aPH latitude
###################################


png(file=paste(path_figures,"Figure6_Aph_latitude2.png",sep=""), width = 1100, height = 720, units = "px", pointsize = 18, bg = "white") 
    par(mfrow=c(1,2))
    par(mar=c(2,3,1,2))
    par(oma=c(3,2,2,1))
    lettersize=2.0
    textinner=1.2
    layout.show(2)

        ##### In situ
        is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
        load(paste(is_dir,"/media_seasonal_Aph.Rdata",sep=""))   
        #dim(matriz_seasonal)                      
        situ<-matriz_seasonal[,,1,,3]
        insitu <- apply(situ,MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
        localizacion<-rep(latitud[,2],each=length(longitud[,2]))

##### A: EXP-1
##################          
        modelname<-modelnames[1]
        load(paste(p_dir,"2iop_aph/",modelname,".RData",sep=""))
        z <- apply(res,MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
        z[z<1e-5] <- NA
        MASCARA<- z
        load(paste(s_dir,"climatologies_2012_2018/media_seasonal_Aph_2deg_OCCCI_2012_2018.RData",sep=""))
        res <- med            
        landa_sat <- dimnames(res)$l
        longitude<-as.numeric(dimnames(res)$x)
        latitude<-as.numeric(dimnames(res)$y)
        sat <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)        
        media <- sat[,,2]
        media[is.na(z)] <- NA
        midpoint<-apply(media,MARGIN=c(2),FUN=mean, na.rm=T)
            
            plot(y=latitude,x=rep(0.1,length=length(latitude)),las=1,xlim=c(0,0.08),ylim=c(-80,90), xaxs="i", yaxs="i",type="n",bty="l",yaxt="n", xlab="", ylab="")
            axis(2,at=seq(-80,90,by=10),labels=F, las=1)
            axis(2,at=seq(-80,80,by=20), las=1)
                # In situ
                localizacion2<-localizacion[!is.nan(c(insitu))]
                valores<-c(insitu)[!is.nan(c(insitu))]
                    grupos<-cut(localizacion2,breaks=seq(-90,90,by=2))
                    grupitos<-split(valores,grupos)
                    intervalos<-seq(-89,89,by=2)
                    longitud<-sapply(grupitos,length)
                na<-sapply(grupitos,quantile,probs=0.5, na.rm=T)
                na[longitud<=2]<-NA
                            segments(x0=na,x1=na,y0=intervalos-1,y1=intervalos+1, lwd=2, col="grey40")   
                                        na3<-sapply(grupitos,quantile,probs=0.90, na.rm=T)
                                        na3[longitud<=2]<-NA
                                        na1<-sapply(grupitos,quantile,probs=0.10, na.rm=T)
                                        na1[longitud<=2]<-NA
                            #segments(x0=na1,x1=na3,y0=intervalos,y1=intervalos, col="grey40", lwd=1)             
                                        na3<-sapply(grupitos,quantile,probs=0.75, na.rm=T)
                                        na3[longitud<=2]<-NA
                                        na1<-sapply(grupitos,quantile,probs=0.25, na.rm=T)
                                        na1[longitud<=2]<-NA
                            segments(x0=na1,x1=na3,y0=intervalos,y1=intervalos, col="grey40", lwd=2)
                                
                # Sat
                mycol<-usecol(pal_grau,n=4,alpha=0.2)
                upper<-apply(media,MARGIN=c(2),FUN=quantile, probs=c(0.6,0.7,0.8,0.9), na.rm=T)
                lower<-apply(media,MARGIN=c(2),FUN=quantile, probs=c(0.4,0.3,0.2,0.1), na.rm=T)
                equis<-c(latitude,latitude[length(latitude):1])
                ies<-c(upper[4,],lower[4,][length(lower[4,]):1])
                polygon(y=equis[!is.na(ies)], x=ies[!is.na(ies)], col=mycol[4], border=NA)
                ies<-c(upper[3,],lower[3,][length(lower[3,]):1])
                polygon(y=equis[!is.na(ies)], x=ies[!is.na(ies)], col=mycol[3], border=NA)
                ies<-c(upper[2,],lower[2,][length(lower[2,]):1])
                polygon(y=equis[!is.na(ies)], x=ies[!is.na(ies)], col=mycol[2], border=NA)
                ies<-c(upper[1,],lower[1,][length(lower[1,]):1])
                polygon(y=equis[!is.na(ies)], x=ies[!is.na(ies)], col=mycol[1], border=NA)
                mtext(3,at=0.005, line=1,text="a) EXP-noPPC",cex=lettersize, outer=F, font=1) 
                
                # Constant
                w=1
                modelname<-modelnames[w]
                      load(paste(p_dir,"2iop_aph/",modelname,".RData",sep=""))
                      dim(res)
                      z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
                      z[z<1e-5] <- NA
                      dim(z)
                      midpoint<-apply(z,MARGIN=c(2),FUN=mean, na.rm=T)
                      #seecol("grad_all")
                      mycol<-usecol(pal_bordeaux,n=4,alpha=0.2)
                      upper<-apply(z,MARGIN=c(2),FUN=quantile, probs=c(0.6,0.7,0.8,0.9), na.rm=T)
                      lower<-apply(z,MARGIN=c(2),FUN=quantile, probs=c(0.4,0.3,0.2,0.1), na.rm=T)
                      equis<-c(latitude,latitude[length(latitude):1])
                      ies<-c(upper[4,],lower[4,][length(lower[4,]):1])
                      polygon(y=equis[!is.na(ies)], x=ies[!is.na(ies)], col=mycol[4], border=NA)
                      ies<-c(upper[3,],lower[3,][length(lower[3,]):1])
                      polygon(y=equis[!is.na(ies)], x=ies[!is.na(ies)], col=mycol[3], border=NA)
                      ies<-c(upper[2,],lower[2,][length(lower[2,]):1])
                      polygon(y=equis[!is.na(ies)], x=ies[!is.na(ies)], col=mycol[2], border=NA)
                      ies<-c(upper[1,],lower[1,][length(lower[1,]):1])
                      polygon(y=equis[!is.na(ies)], x=ies[!is.na(ies)], col=mycol[1], border=NA)
                  
                          
                ### SCALE
                par(new=T)
                mycol<-usecol(pal_grau,n=4,alpha=0.8)
                par(mar=c(10,14,14,9))
                mat <- matrix(c(1:4), ncol=4)
                image(mat, col=mycol, las=1, xaxt="n", yaxt="n", main="", cex.main=1.2, bty="n") 
                axis(4,at=seq(0.00,1.00,length=4),labels=F, cex.axis=1.2, las=1, col=mycol[1]) 
                mtext(3,at=0.1, line=0.02,text="Satellite", cex=0.88)
                
                par(new=T)
                mycol<-usecol(pal_bordeaux,n=4,alpha=0.8)
                par(mar=c(10,18,14,5))
                mat <- matrix(c(1:4), ncol=4)
                image(mat, col=mycol, las=1, xaxt="n", yaxt="n", main="", cex.main=1.2, bty="n") 
                axis(4,at=seq(0.00,1.00,length=4),labels=c("[40-60]%","[30-70]%","[20-80]%","[10-90]%"), cex.axis=1.0, las=1,col=mycol[1])    
                mtext(3,at=0.1, line=0,text=expression(paste(a[PH],"*(",lambda,")"," cnt",sep="")), cex=0.88)
##################                
                                
                       
##### B: EXP-2
##################          
          par(mar=c(2,2,1,1))
          modelname<-modelnames[2]
          load(paste(p_dir,"2iop_aph/",modelname,".RData",sep=""))
          z <- apply(res,MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
          z[z<1e-5] <- NA
          MASCARA<- z
          load(paste(s_dir,"climatologies_2012_2018/media_seasonal_Aph_2deg_OCCCI_2012_2018.RData",sep=""))
          res <- med            
          landa_sat <- dimnames(res)$l
          longitude<-as.numeric(dimnames(res)$x)
          latitude<-as.numeric(dimnames(res)$y)
          sat <- apply(res, MARGIN=c("x","y","l"), FUN=mean, na.rm=TRUE)        
          media <- sat[,,2]
          media[is.na(z)] <- NA
          midpoint<-apply(media,MARGIN=c(2),FUN=mean, na.rm=T)

            plot(y=latitude,x=rep(0.1,length=length(latitude)),las=1,xlim=c(0,0.08),ylim=c(-80,90), xaxs="i", yaxs="i",type="n",
                 bty="l",yaxt="n", xlab="", ylab="")
            axis(2,at=seq(-80,90,by=10),labels=F, las=1)
            axis(2,at=seq(-80,80,by=20), las=1)
                # In situ
                localizacion2<-localizacion[!is.nan(c(insitu))]
                valores<-c(insitu)[!is.nan(c(insitu))]
                grupos<-cut(localizacion2,breaks=seq(-90,90,by=2))
                grupitos<-split(valores,grupos)
                intervalos<-seq(-89,89,by=2)
                longitud<-sapply(grupitos,length)
                na<-sapply(grupitos,quantile,probs=0.5, na.rm=T)
                na[longitud<=2]<-NA
                    segments(x0=na,x1=na,y0=intervalos-1,y1=intervalos+1, lwd=2, col="grey40")   
                              na3<-sapply(grupitos,quantile,probs=0.90, na.rm=T)
                              na3[longitud<=2]<-NA
                              na1<-sapply(grupitos,quantile,probs=0.10, na.rm=T)
                              na1[longitud<=2]<-NA
                    #segments(x0=na1,x1=na3,y0=intervalos,y1=intervalos, col="grey40", lwd=1)             
                                  na3<-sapply(grupitos,quantile,probs=0.75, na.rm=T)
                                  na3[longitud<=2]<-NA
                                  na1<-sapply(grupitos,quantile,probs=0.25, na.rm=T)
                                  na1[longitud<=2]<-NA
                    segments(x0=na1,x1=na3,y0=intervalos,y1=intervalos, col="grey40", lwd=2)
                      
                ### Sat            
                mycol<-usecol(pal_grau,n=4,alpha=0.2)
                upper<-apply(media,MARGIN=c(2),FUN=quantile, probs=c(0.6,0.7,0.8,0.9), na.rm=T)
                lower<-apply(media,MARGIN=c(2),FUN=quantile, probs=c(0.4,0.3,0.2,0.1), na.rm=T)
                equis<-c(latitude,latitude[length(latitude):1])
                ies<-c(upper[4,],lower[4,][length(lower[4,]):1])
                polygon(y=equis[!is.na(ies)], x=ies[!is.na(ies)], col=mycol[4], border=NA)
                ies<-c(upper[3,],lower[3,][length(lower[3,]):1])
                polygon(y=equis[!is.na(ies)], x=ies[!is.na(ies)], col=mycol[3], border=NA)
                ies<-c(upper[2,],lower[2,][length(lower[2,]):1])
                polygon(y=equis[!is.na(ies)], x=ies[!is.na(ies)], col=mycol[2], border=NA)
                ies<-c(upper[1,],lower[1,][length(lower[1,]):1])
                polygon(y=equis[!is.na(ies)], x=ies[!is.na(ies)], col=mycol[1], border=NA)
                mtext(3,at=0.005, line=1,text="b) EXP-PPC",cex=lettersize, outer=F, font=1) 
                
                # Variable
                w=2
                modelname<-modelnames[w]
                    load(paste(p_dir,"2iop_aph/",modelname,".RData",sep=""))
                    z <- apply(res[,,3,],MARGIN=c("x","y"),FUN=mean, na.rm=TRUE)
                    z[z<1e-5] <- NA
                    midpoint<-apply(z,MARGIN=c(2),FUN=mean, na.rm=T)
                    #seecol("grad_all")
                    mycol<-usecol(pal_petrol,n=4,alpha=0.2)
                    upper<-apply(z,MARGIN=c(2),FUN=quantile, probs=c(0.6,0.7,0.8,0.9), na.rm=T)
                    lower<-apply(z,MARGIN=c(2),FUN=quantile, probs=c(0.4,0.3,0.2,0.1), na.rm=T)
                    equis<-c(latitude,latitude[length(latitude):1])
                    ies<-c(upper[4,],lower[4,][length(lower[4,]):1])
                    polygon(y=equis[!is.na(ies)], x=ies[!is.na(ies)], col=mycol[4], border=NA)
                    ies<-c(upper[3,],lower[3,][length(lower[3,]):1])
                    polygon(y=equis[!is.na(ies)], x=ies[!is.na(ies)], col=mycol[3], border=NA)
                    ies<-c(upper[2,],lower[2,][length(lower[2,]):1])
                    polygon(y=equis[!is.na(ies)], x=ies[!is.na(ies)], col=mycol[2], border=NA)
                    ies<-c(upper[1,],lower[1,][length(lower[1,]):1])
                    polygon(y=equis[!is.na(ies)], x=ies[!is.na(ies)], col=mycol[1], border=NA)
                    

                    ### SCALE
                    par(new=T)
                    mycol<-usecol(pal_grau,n=4,alpha=0.8)
                    par(mar=c(10,14,14,9))
                    mat <- matrix(c(1:4), ncol=4)
                    image(mat, col=mycol, las=1, xaxt="n", yaxt="n", main="", cex.main=1.2, bty="n") 
                    axis(4,at=seq(0.00,1.00,length=4),labels=F, cex.axis=1.2, las=1,col=mycol[1]) 
                    mtext(3,at=0.1, line=0.02,text="Satellite", cex=0.88)
                    par(new=T)
                    mycol<-usecol(pal_petrol,n=4,alpha=0.8)
                    par(mar=c(10,18,14,5))
                    mat <- matrix(c(1:4), ncol=4)
                    image(mat, col=mycol, las=1, xaxt="n", yaxt="n", main="", cex.main=1.2, bty="n") 
                    axis(4,at=seq(0.00,1.00,length=4),labels=c("[40-60]%","[30-70]%","[20-80]%","[10-90]%"), cex.axis=1.0, las=1,col=mycol[1])    
                    mtext(3,at=0.1, line=0,text=expression(paste(a[PH],"*(",lambda,")"," var",sep="")), cex=0.88)
##################                    

      ## Labels
      mtext(2,at=0.50,line=-0.5,outer=T,text=expression(paste(degree, "N")),las=1, cex=1.2)
      mtext(1,at=c(0.25,0.75),line=0.75,outer=T,text=expression(paste(a[PH]," (450nm) [", m^-1,"] in situ & model")))
      mtext(1,at=c(0.25,0.75),line=1.75,outer=T,text=expression(paste(a[PH]," (443nm) [", m^-1,"] satellite")))

      dev.off()        
                








#################################################
#### FIGURE 7: metrics aPH observed vs. simulated
#################################################         

png(file=paste(path_figures,"Figure7_aph_metrics.png",sep=""), width = 900, height = 1100, units = "px", pointsize = 18, bg = "white") 
     
      layout(matrix(c(2,1,9,4,3,10,6,5,11,8,7,12),ncol=3, byrow=TRUE), widths = c(0.4,1,1))
      par(family="")
      par(mar=c(4,4,0,0))  
      par(oma=c(1,3,3,2))                 
      layout.show(n=12)
          mycol1<-usecol(pal_bordeaux,n=4)
          mycol2<-usecol(pal_petrol,n=4)
          colores_puntos<-c(mycol1[3],mycol2[3])
          colores<-c(mycol1[1],mycol2[1])
          lettersize=2.0

## MODEL vs FIELD/SAT         
############################             
          simbolos<-letters[seq( from = 1, to = 14)]
          lambda[c(1,3,4,5,7,12)]
          landa_correspondiente<-c(1,1,2,3,4,4,5,5,5,5,5,6,6)
          landa_correspondiente[c(1,3,4,5,7,12)]
          letra<-c("a)","b)","c)","d)","e)")


        for (w in c(1:length(modelnames))){ 
          modelname<-modelnames[w]  
          ### Mod             
          load(paste(p_dir,"2iop_aph/",modelname,".RData",sep=""))
          #dim(res)
          ### In situ
          is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")
          load(paste(is_dir,"/media_seasonal_Aph.Rdata",sep=""))   
          dim(matriz_seasonal[,,1,,])
          ### Sat
          load(paste(s_dir,"climatologies_2012_2018/media_seasonal_Aph_2deg_OCCCI_2012_2018.RData",sep=""))
          #dim(med) 
          
          for (l in c(1,3,4,5)){    
                #### Aph mod
                z <- apply(res[,,l,],MARGIN=c("x","y","t"),FUN=mean, na.rm=TRUE)
                z[z<1e-5] <- NA
                dim(z)
                ##### In situ
                situ<-matriz_seasonal[,,1,,l]
                dim(situ)
                ### SAT
                landa_sat <- dimnames(med)$l
                landa_correspondiente[l]
                media <- med[,,,landa_correspondiente[l]]
                dim(media)
                media[is.na(z)] <- NA
          
         plot(log10(c(z)), log10(c(media)), xlim=c(-4,0), ylim=c(-4,0), pch=19,cex=0.2,
              mgp=c(2,1,0), xaxt="n", yaxt="n",  col="grey30", xaxs="i",yaxs="i",
              xlab=substitute(paste(a[PH],"(",wl,"nm) [",m^-1,"] model"), list(wl=lambda[l])), ylab="")
              abline(a=0,b=1, col="black")
         points(log10(c(z)), log10(c(situ)), xlim=c(-4,0), ylim=c(-4,0), pch=19, cex=0.3, col="grey70")
                    axis(1, at=seq(-4,1), labels=10^seq(-4,1), las=1)
                    axis(2, at=seq(-4,1), labels=10^seq(-4,1), las=1)
                    regre1<-lm(log10(c(media))~log10(c(z)))
                    regre2<-lm(log10(c(situ))~log10(c(z)))

         histogram<-hist(log10(z),breaks=c(-6,seq(-4,0,length=98),2),plot=F)
         par(new=T)
         barplot(histogram$density[-c(1,length(histogram$density))], xlab="", ylab="", ylim=c(0,4), xaxt="n", yaxt="n", col="grey50",border="grey50")
                    
                    
         # Correlation
         x<-cbind(log10(z), log10(media))
                    x<-x[complete.cases(x)==TRUE,]
                    ene <- nrow(x)
                    R_chla_s     <- rcorr(x, type="pearson")$r[1,2]
                    RI_chla_s    <- exp(sqrt((sum((log(x[,2]/x[,1]))^2, na.rm=TRUE))/nrow(x)))                  
                    MEF_chla_s   <- (sum((x[,2]-mean(x[,2]))^2) - sum((x[,1]-x[,2])^2)) / sum((x[,2]-mean(x[,2]))^2)             
                    # igual a: 1-(RMSE_chla^2)/var(x[,2])
                    RMSE_chla_s  <- sqrt((sum((x[,1]-x[,2])^2))/nrow(x))
                    AE_chla_s    <- (sum(x[,1]-x[,2]))/nrow(x)
                    AAE_chla_s   <- (sum(abs(x[,1]-x[,2])))/nrow(x)                 
                    x<-cbind(c(z), c(media))
                    x<-x[complete.cases(x)==TRUE,]
                    AE_chla_s    <- (sum(x[,1]-x[,2]))/nrow(x)
         legend(x="bottomright",legend=paste(c("R=","RMSE=","Bias="),c(round(R_chla_s,3),round(RMSE_chla_s,3),round(AE_chla_s,3)), sep=""),bty="n", text.col="black")            
                    
         # Correlation
         x<-cbind(log10(z), log10(situ))
                    x<-x[complete.cases(x)==TRUE,]
                    ene <- nrow(x)
                    R_chla_s     <- rcorr(x, type="pearson")$r[1,2]
                    RI_chla_s    <- exp(sqrt((sum((log(x[,2]/x[,1]))^2, na.rm=TRUE))/nrow(x)))                  
                    MEF_chla_s   <- (sum((x[,2]-mean(x[,2]))^2) - sum((x[,1]-x[,2])^2)) / sum((x[,2]-mean(x[,2]))^2)             
                    # igual a: 1-(RMSE_chla^2)/var(x[,2])
                    RMSE_chla_s  <- sqrt((sum((x[,1]-x[,2])^2))/nrow(x))
                    AE_chla_s    <- (sum(x[,1]-x[,2]))/nrow(x)
                    AAE_chla_s   <- (sum(abs(x[,1]-x[,2])))/nrow(x)                 
                    x<-cbind(c(z), c(situ))
                    x<-x[complete.cases(x)==TRUE,]
                    AE_chla_s    <- (sum(x[,1]-x[,2]))/nrow(x)
         legend(x="topleft",legend=paste(c("R=","RMSE=","Bias="),c(round(R_chla_s,3),round(RMSE_chla_s,3),round(AE_chla_s,3)), sep=""),bty="n", text.col="grey60")
              
                    
                    if (w==1){
                 histogram1<-hist(log10(media),breaks=c(-6,seq(-4,0,length=98),2),plot=F)
                 escalado<-histogram1$counts/max(histogram1$counts,na.rm=T)
                 barplot(escalado[-c(1,length(histogram1$density))],
                         xlab="", xlim=c(0,1), xaxt="n", yaxt="n", yaxs="i", horiz=T, col="grey30",border="grey30", bty="o",
                         ylab=substitute(paste(a[PH],"(",wl,"nm) [",m^-1,"] satellite"), list(wl=substring(landa_sat[landa_correspondiente[l]],5,7))))

                 histogram2<-hist(log10(situ),breaks=c(-6,seq(-4,0,length=98),2),plot=F)
                 escalado<-histogram2$counts/max(histogram2$counts,na.rm=T)
                 barplot(escalado[-c(1,length(histogram2$density))]/4,
                         xlab="", ylab="", xlim=c(0,2), xaxt="n", yaxt="n", horiz=T, add=T, col="grey70", border="grey70")
                         axis(2, at=seq(1,116,length=5), labels=10^seq(-4,0), las=1)     

                 mtext(2, line=4, at=60, cex=0.70, text=substitute(paste(a[PH],"(",wl,"nm) [",m^-1,"] in situ"), list(wl=lambda[l])), col="grey60")
                 mtext(3,at=0.5, line=-2,text=letra[landa_correspondiente[l]],cex=lettersize, outer=F, font=1)
                 mtext(3,at=1.1, line=-2.2,text=paste(lambda[l],"nm", sep=""),cex=0.9, outer=F, font=3)
                 }                    
                    
          } # end lambda
          
          } # end loop modelnames              
        
############################  

    ## Labels
    mtext(3,at=c(0.42,0.82),line=1,outer=T,text=nombre)
          
dev.off()           
