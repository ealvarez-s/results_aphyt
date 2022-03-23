if (Sys.info()['sysname']=="Windows") {OS<-"C:/"} else {OS<-paste("/Users/",Sys.info()['user'],"/", sep="")}
source(paste(OS,"Programas/Funciones_varios/functions.r",sep=""))
library(Hmisc)
library(plotrix)
library(plot3D)
library(RNetCDF)
library(akima)
library(abind)
library(unikn)      
library(viridisLite)
library(wesanderson)


#######################         
#### FIGURES in DEPTH constant vs variable
#######################  
# Chla, IOPs, Dominance Abs/Scat

      o_dir<-paste(OS,"Datos/Ser_Stan/global_14_aphyt/",sep="")        
      p_dir<-paste(OS,"Datos/Ser_Stan/global_14_aphyt/interpolated/",sep="")
      #path_figures<-paste(OS,"Documentos/5_Trabajos/20_Radtrans_aph/Figuras/Figuras_2021_review/",sep="")
      path_figures<-paste(OS,"Documentos/5_Trabajos/20_Radtrans_aph/reviews_coauthors2/para_enviar_JAMES/Figuras/",sep="")
      experimentos<-read.csv(paste(OS,"Datos/Ser_Stan/global_14_aphyt/run_log_marshall_PPC_PS_SA_G100.csv",sep=""), sep=",")
      s_dir <- paste(OS,"Datos/Dat_Satelite/",sep="")

      modelnames<-experimentos$name[c(21,9)]
      ruta<-experimentos$path[c(21,9)]        
      nombre <- c(expression(paste("EXP-1: ",{a^{"*"}}[PH],"(", lambda, ") constant",sep="")),
                  expression(paste("EXP-2: ",{a^{"*"}}[PH],"(", lambda, ") variable",sep="")))
      
      
            # NPP
            load(paste(OS,"Datos/Ser_G100/global_14_aphyt/Tannual_NPP_final4.RData",sep=""))     
            NPP<-NPPtotal
            EXP<-EXPtotal
            # Zeu
            load(paste(o_dir,"Zeu_3D10daily_final4.RData",sep=""))    
    

            ## Color palette
            mycol1<-usecol(pal_bordeaux,n=4)
            mycol2<-usecol(pal_petrol,n=4)
            colores_puntos<-c(mycol1[3],mycol2[3])
            colores<-c(mycol1[1],mycol2[1])
            custom_scale<-wes_palette(length(custom_scale),name = "Zissou1", type = "continuous")
            custom_scale <- c("deepskyblue4","deepskyblue3","darkslategray1","mediumaquamarine","greenyellow","yellow",
                              "gold1","orange", "tomato","orangered2")
    
    
    
####################
## FIGURE 6 Absorption of constituents in depth
####################
    
    myfunc <- function(x){
      if (sum(is.na(x))!=0){ return(NA)}
      else {return(as.numeric(unlist(which(x==max(x, na.rm=TRUE)))))}}
  
  
png(file=paste(path_figures,"Figure6_IOPs_depth_dominance.png",sep=""), width = 1200, height = 800, pointsize = 18,units = "px",bg = "white") 
    layout(matrix(c(1,4:7,2,8:11,3,12:15),ncol=3, byrow=FALSE),widths=c(1,1,0.18))
    par(mar=c(2,2,1,1))
    par(oma=c(2,3,2,2))
    lettersize=1.8
    textinner=1.2
    layout.show(n=15)

   
## Chla    
####################   
    for (w in c(1:length(modelnames))){
      zeu_mean <- apply(ZEU[,,,w], MARGIN=c("y"), FUN=mean,na.rm=TRUE)
      load(paste(p_dir,"1cloro_global/",modelnames[w],".RData",sep="")) 
      cloro <- res    
      
      ## DCM CHLOROPHYLL
      matrizC <- apply(cloro[,,], MARGIN=c("y","z"), mean, na.rm=TRUE)
      depthC   <- as.numeric(colnames(matrizC))
      latitudC <- as.numeric(rownames(matrizC))   
      matrizC[matrizC>1] <- 1
      image2D(x=latitudC, y=depthC, z=matrizC, col=viridis(100), ylim=c(-200.0,0), zlim=c(0.0,1.0),
              main="", las=1, xlim=c(-82,80), colkey=FALSE, xaxt="n", cex.main=2.0)
      axis(1, at=seq(-80,80,by=20), labels=seq(-80,80,by=20), las=1)
      points(x=latitude, y=zeu_mean, col="black", type="l", lty=1, lwd=2)
      text(x=28,y=-150,  labels=paste("NPP =",round(NPP[w],2)), col="white", font=1, pos=4, cex=1.2)
      text(x=28,y=-180,  labels=paste("EXP =",-round(EXP[w],2)), col="white", font=1, pos=4, cex=1.2)
      text(x=79,y=-170,  labels=expression(paste("PgC ",y^-1, sep="")), col="white", font=2, pos=2, cex=1.2)
   
      mtext(3, at=2, line=0.5, text=nombre[w], cex=1.3, outer=F, font=1)
      
      if(w==1){mtext(3,at=-72,text="a)",cex=lettersize, outer=F, font=1)}      
      if(w==2){mtext(3,at=-72,text="b)",cex=lettersize, outer=F, font=1) }  
      }
    
    # SCALE
    par(mar=c(2,0.5,1,3))
    escala <- matrix(c(1:100), ncol=100, byrow=TRUE)
    image(escala, las=1,col=viridis(100), xaxt="n", yaxt="n", bty="n", main=expression(paste("mg ",m^-3, sep="")), cex.main=1.0)  
    axis(4,at=seq(0-0.001,1+0.001,length=6),labels=c(seq(0,1.0,length=6)[-6],">1"), las=1,cex.axis=1.2)      
####################
    
    
## Constituents
##################    
    par(mar=c(2,2,1,1))
    for (w in 1:length(modelnames)){    
      modelname<-modelnames[w]
      zeu <- apply(ZEU[,,,w], MARGIN=c(2), FUN=mean,na.rm=TRUE)   
      latitude_zeu <- latitude    
      
      ## ABSORPTION
      load(paste(o_dir,"interpolated/2iop_abtot_wb/",modelname,"_450.RData",sep=""))
      abtotave <- res
      load(paste(o_dir,"interpolated/2iop_acdom_wb/",modelname,"_450.RData",sep=""))
      acdomave <- res
      load(paste(o_dir,"interpolated/2iop_abpar_wb/",modelname,"_450.RData",sep=""))
      abparave <- res
      load(paste(o_dir,"interpolated/2iop_abphy_wb/",modelname,"_450.RData",sep=""))
      abphyave <- res
      

      matrizT <- apply(abtotave, MARGIN=c("y","z"), mean, na.rm=TRUE)
      depthC   <- as.numeric(colnames(matrizT))
      latitudC <- as.numeric(rownames(matrizT))
      matrizA <- apply(abphyave, MARGIN=c("y","z"), mean, na.rm=TRUE)
      matrizB <- apply(acdomave, MARGIN=c("y","z"), mean, na.rm=TRUE)
      matrizC <- apply(abparave, MARGIN=c("y","z"), mean, na.rm=TRUE)
      matrizW <- matrizT-(matrizA+matrizB+matrizC)
      

      # aPH      
      image2D(x=latitudC, y=depthC, z=matrizA, col=custom_scale, ylim=c(-200.0,0), xlim=c(-82,80),main="", las=1, zlim=c(0,0.05), colkey=F,xaxt="n")
          axis(1, at=seq(-80,80,by=20), labels=seq(-80,80,by=20), las=1)
          points(x=latitude_zeu, y=zeu, col="black", type="l", lty=1, lwd=2)
          text(x=75,y=-165,labels=expression(a[PH]("450nm")), font=2, cex=1.4, col="white", pos=2)
          if(w==1){mtext(3,at=-72,text="c)",cex=lettersize, outer=F, font=1) }      
          if(w==2){mtext(3,at=-72,text="d)",cex=lettersize, outer=F, font=1) }  
          
      # aCDOM    
      image2D(x=latitudC, y=depthC, z=matrizB, col=custom_scale, ylim=c(-200.0,0), xlim=c(-82,80),main="", las=1, zlim=c(0,0.05), colkey=F, xaxt="n")
          axis(1, at=seq(-80,80,by=20), labels=seq(-80,80,by=20), las=1)
          points(x=latitude_zeu, y=zeu, col="black", type="l", lty=1, lwd=2)
          text(x=75,y=-165,labels=expression(a[CDOM]("450nm")), font=2, cex=1.4, pos=2)
          if(w==1){mtext(3,at=-72,text="e)",cex=lettersize, outer=F, font=1) }      
          if(w==2){mtext(3,at=-72,text="f)",cex=lettersize, outer=F, font=1) }
            
      # aNAP
      image2D(x=latitudC, y=depthC, z=matrizC, col=custom_scale, ylim=c(-200.0,0), xlim=c(-82,80), main="", las=1, zlim=c(0,0.01), colkey=F,xaxt="n")
          axis(1, at=seq(-80,80,by=20), labels=seq(-80,80,by=20), las=1)
          points(x=latitude_zeu, y=zeu, col="black", type="l", lty=1, lwd=2)
          text(x=75,y=-165,labels=expression(a[NAP]("450nm")), font=2, cex=1.4, col="black", pos=2)
          if(w==1){mtext(3,at=-72,text="g)",cex=lettersize, outer=F, font=1) }      
          if(w==2){mtext(3,at=-72,text="h)",cex=lettersize, outer=F, font=1) }
            
      # Dominance Absorption 
      matrizW <- matrizW/matrizT
      matrizA <- matrizA/matrizT
      matrizB <- matrizB/matrizT
      matrizC <- matrizC/matrizT 
      res <- abind(matrizW,matrizB,matrizA,matrizC,along=0.5)
      dimnames(res) <- list(h=c(1:4), y=dimnames(matrizT)$y, z=dimnames(matrizT)$z)
      m <- apply(res, MARGIN=c("y","z"), FUN=myfunc)
      #dim(m)
      #range(m, na.rm=TRUE)    
      #unique(m)
      depthC   <- as.numeric(colnames(matrizT))
      latitudC <- as.numeric(rownames(matrizT))
      paleta<-wes_palette(4, name = "Zissou1", type = "continuous")
      #c("lightsteelblue3","goldenrod3","palegreen3","indianred3")
      image2D(x=latitudC, y=depthC, z=m, breaks=seq(0.5,4.5,by=1),
              col=paleta[c(1,3,2,4)],ylim=c(-200.0,0), main="", las=1,colkey =F, xlim=c(-82,80),xaxt="n")
              axis(1, at=seq(-80,80,by=20), labels=seq(-80,80,by=20), las=1)
              points(x=latitude_zeu, y=zeu, col="black", type="l", lty=1, lwd=2)
              if(w==1){mtext(3,at=-72,text="i)",cex=lettersize, outer=F, font=1) }      
              if(w==2){mtext(3,at=-72,text="j)",cex=lettersize, outer=F, font=1) } 
    } # end loop w
##################       
    
    ### SCALE
    par(mar=c(2,0.5,1,3))
    mat <- matrix(c(1:length(custom_scale)), ncol=length(custom_scale))
    image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(m^-1), cex.main=1.2) 
      axis(4,at=seq(0-0.05,1+0.05, length=11)[seq(1,11,by=2)],labels=round(seq(0,0.01,length=11),3)[seq(1,11,by=2)], cex.axis=1.2, las=1)    
    image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(m^-1), cex.main=1.2) 
      axis(4,at=seq(0-0.05,1+0.05, length=11)[seq(1,11,by=2)],labels=round(seq(0,0.05,length=11),3)[seq(1,11,by=2)], cex.axis=1.2, las=1)    
    image(mat, col=custom_scale, las=1, xaxt="n", yaxt="n", main=expression(m^-1), cex.main=1.2) 
      axis(4,at=seq(0-0.05,1+0.05, length=11)[seq(1,11,by=2)],labels=round(seq(0,0.05,length=11),3)[seq(1,11,by=2)], cex.axis=1.2, las=1)
    mat <- matrix(c(1:4), ncol=4)
    paleta<-wes_palette(4, name = "Zissou1", type = "continuous")
    image(mat, col=paleta, las=1, xaxt="n", yaxt="n", main="", cex.main=1.2) 
      axis(4,at=seq(0,1, length=4),labels=c("water","phyto","cdom","nap"), cex.axis=1.2,las=1, font=2)        

##################    
    ## Labels
    mtext(1,at=c(0.24,0.70),line=1,outer=T,text=expression(paste(degree,"N",sep="")))
    mtext(2,at=0.5,line=1,outer=T,text="depth (m)")  

dev.off()      
