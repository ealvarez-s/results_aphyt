if (Sys.info()['sysname']=="Windows") {OS<-"C:/"} else {OS<-paste("/Users/",Sys.info()['user'],"/", sep="")}
library(unikn) 
library(wesanderson)

###################################################
### INITIAL OPTICAL PROPERTIES for Constituents ###                                                  
###################################################

# Average aPH curves from literature (Alvarez et al. 2022 PiO under review)
# Weight-specific absorption coefficients for pigments (Bidigare et al. 1990)
# Curve decomposition for aPS (Bidigare 1990, Babin et al 1996, Hickman et al 2010)
# Mass-specific b and bb_to_b (Dutkiewicz, 2018)


paleta<-wes_palette(3, name = "Chevalier1", type = "continuous")
path_figures<-paste(OS,"Documentos/5_Trabajos/20_Radtrans_aph/reviews_coauthors2/para_enviar_JAMES/Figuras/",sep="") 


      nombres <-c("Prochlorococcus", "Synechococcus", "Small Eukaryotes", "Green algae", "Cocolithophorids",    "Brown algae", "Diatoms")  
      clusters<-c("Prochlorococcus", "Synechococcus", "SmallEuk",         "Chloro",      "Emiliania",           "Brown",       "Diatoms")      
      length(nombres)
      length(clusters)

      oscuros<-c(paleta[1],paleta[3], paleta[2], paleta[1],  paleta[2], paleta[2],  paleta[2])
      #oscuros<-c("green2","blue2", "orange2", "green2", "orange2", "orange2", "orange2")
      claros <-c("green", "blue",  "orange",  "green",  "orange",  "orange",  "orange")
      colores<-oscuros
      colores2<-claros
      length(colores)
      
      
################################################################################
### Average Chla-specific aPH (m2 mgChla-1) and reconstruct aPS (m2 mmolC-1) ###
################################################################################
      
      ### ABSORPTION
      dir<-"/Users/ealvarez/Datos/Dat_Aph/"
      datos<-read.csv(paste(dir,"csvs/absorption_spectra.csv", sep=""), sep=";")
      datos<-datos[datos$FLAG_1==0,]
      datos<-datos[datos$FLAG_2==0,]
      datos<-datos[!is.na(datos[,52]),]
      lambdas_in<-as.numeric(substring(colnames(datos)[22:102],2,4))
      lambdas_out<-seq(400,700)
      resultado_aph <- matrix(NA, nrow=length(nombres), ncol=length(lambdas_out))
      resultado_aps <- matrix(NA, nrow=length(nombres), ncol=length(lambdas_out))
      
      ### SPECTRA RECONSTRUCTION   
      ## Weight-specific absorption coefficients for pigments downloaded from:
      ## http://www.oceanopticsbook.info/view/optical_constituents_of_the_ocean/_phytoplankton
      espectros<-read.csv(paste(dir,"Bidigare_et_al_1990.csv",sep=""))
      lambda_bidigare<-espectros$lAMBDA
      spectra_reconstruction<-function(x,P1="CHLA",P2="PEB",P3="PSC",P4="PPC"){
        fila<-unlist(x)
        fila_new<-approx(x=lambdas_in, y=fila, xout=lambda_bidigare)
        aph<-fila_new$y
        # weight-specific abs. coeff. : Bidigare 1990
        final <- data.frame(espectros[,c(which(colnames(espectros)==P1),
                                         which(colnames(espectros)==P2),
                                         which(colnames(espectros)==P3),
                                         which(colnames(espectros)==P4))],aph, espectros$lAMBDA)     
        colnames(final)<-c("P1","P2","P3","P4","aph","lambda")
        final <- final[complete.cases(final),]
        
        # fit linear model: Hickman 2010
        fit <- lm(aph ~ P1 + P2 + P3 + P4 -1, data=final)
        #summary(fit,correlation =TRUE)
        coeffs <- round(coefficients(fit),3)
        
        # pigment ratio : Babin et al 1996
        cnp<-(coeffs[4]*final[,4])/((coeffs[1]*final[,1])+(coeffs[2]*final[,2])+(coeffs[3]*final[,3])+(coeffs[4]*final[,4]))
        aps<-final[,5]*(1-cnp)
        aps_new<-approx(x=final[,6], y=aps, xout=lambdas_in)
        aps<-aps_new$y
        return(aps)}    
   

      
##############################################################         
### Aph spectra averaged per cluster and aps reconstructed ###
##############################################################     
      
      for (j in seq(1,length(clusters))){
        #datos$Genera[datos$Cluster==j]
        genero<-datos[datos$aph_class==clusters[j] & !is.na(datos$aph_class),22:102]
        genero[,lambdas_in<400]<-NA
        resul<-approx(x=lambdas_in, y=colMeans(genero,na.rm=T), xout=lambdas_out)
        resultado_aph[j,]<-resul$y        
      
        # Pigments of each group:       
        if (nombres[j]=="Prochlorococcus"){
          genero2<-apply(genero,MARGIN=1,FUN=spectra_reconstruction,
                         P1="CHLA",P2="CHLB",  P3="PSC",P4="PPC")
                         etiqueta<-c("CHLA","CHLB","PSC","PPC")}
        if (nombres[j]=="Synechococcus"){
          genero2<-apply(genero,MARGIN=1,FUN=spectra_reconstruction,
                         P1="CHLA",P2="PEB",   P3="PSC",P4="PPC")
                         etiqueta<-c("CHLA","PEB","PSC","PPC")}
        if (nombres[j]=="Small Eukaryotes"){
          genero2<-apply(genero,MARGIN=1,FUN=spectra_reconstruction,
                         P1="CHLA",P2="CHLC",   P3="PSC",P4="PPC")
                         etiqueta<-c("CHLA","CHLC","PSC","PPC")}
        if (nombres[j]=="Green algae"){
          genero2<-apply(genero,MARGIN=1,FUN=spectra_reconstruction,
                         P1="CHLA",P2="CHLB",  P3="PSC",P4="PPC")
                         etiqueta<-c("CHLA","CHLB","PSC","PPC")}      
        if (nombres[j]=="Brown algae"){
          genero2<-apply(genero,MARGIN=1,FUN=spectra_reconstruction,
                         P1="CHLA",P2="CHLC",   P3="PSC",P4="PPC")
                         etiqueta<-c("CHLA","CHLC","PSC","PPC")}       
        if (nombres[j]=="Cocolithophorids"){
          genero2<-apply(genero,MARGIN=1,FUN=spectra_reconstruction,
                         P1="CHLA",P2="CHLC",   P3="PSC",P4="PPC")
                         etiqueta<-c("CHLA","CHLC","PSC","PPC")}   
        
        if (nombres[j]=="Diatoms"){
          genero2<-apply(genero,MARGIN=1,FUN=spectra_reconstruction,
                         P1="CHLA",P2="CHLC",   P3="PSC",P4="PPC")
                         etiqueta<-c("CHLA","CHLC","PSC","PPC")}      
    
        genero2<-t(genero2)
        resul<-approx(x=lambdas_in, y=colMeans(genero2,na.rm=T), xout=lambdas_out)
        resultado_aps[j,]<-resul$y
      }  
      
      rownames(resultado_aph)<-nombres
      colnames(resultado_aph)<-lambdas_out
      rownames(resultado_aps)<-nombres
      colnames(resultado_aps)<-lambdas_out

    
##########################################################
### Average mass-specific b (m2 mgC-1) to (m2 mmolC-1) ###
##########################################################
     
     dir<-"/Users/ealvarez/Datos/Dat_Aph/"
     spectra<-read.csv(paste(dir,"csvs/scatter_mass_Dut2015.csv",sep=""), sep=";")
     #names(spectra)
      # Agregar por grupo
      resultado_bph  <- rbind(spectra$Small_B,  spectra$Diatom_B)
      resultado_bb   <- rbind(spectra$Small_BB, spectra$Diatom_BB)
      lambdas_small<-seq(400,700,by=25)
      rownames(resultado_bph)<-c("SmallPhyto","Diatoms")
      colnames(resultado_bph)<-lambdas_small        
      
    
      
#############################################
### Create initial optics_phyto_recom.dat ###
#############################################   
  
  modelo_lambda <- c(400.0, 425.0, 450.0, 475.0,500.0, 525.0,
                     550.0, 575.0, 600.0, 625.0, 650.0, 675.0,700.0)
  lambda_extremos<-sort(unique(c(modelo_lambda[-length(modelo_lambda)]+(diff(modelo_lambda)/2),400,700)))
  
      # Datos
      lambdas_out
      datos11<-colMeans(resultado_aph[1:6,],na.rm=T)
      datos21<-resultado_aph[7,]
      datos12<-colMeans(resultado_aps[1:6,],na.rm=T)
      datos22<-resultado_aps[7,]      
      
      datos13<-resultado_bph[1,]
      datos23<-resultado_bph[2,]
      datos14<-resultado_bb[1,]
      datos24<-resultado_bb[2,]
      
      datos11b<-rep(NA,length(modelo_lambda))
      datos12b<-rep(NA,length(modelo_lambda))
      datos21b<-rep(NA,length(modelo_lambda))
      datos22b<-rep(NA,length(modelo_lambda))
      datos13b<-rep(NA,length(modelo_lambda))
      datos23b<-rep(NA,length(modelo_lambda))
      
      for (i in c(1:length(modelo_lambda))){
        cuales<-which(lambdas_out>=lambda_extremos[i] & lambdas_out<lambda_extremos[i+1])
        datos11b[i]<-mean(datos11[cuales],na.rm=T)
        datos12b[i]<-mean(datos12[cuales],na.rm=T)
        datos21b[i]<-mean(datos21[cuales],na.rm=T)
        datos22b[i]<-mean(datos22[cuales],na.rm=T)
        datos13b[i]<-mean(datos13[cuales],na.rm=T)
        datos23b[i]<-mean(datos23[cuales],na.rm=T)
      }      

      datosothers  <- cbind(paste("",modelo_lambda, sep=" "),
                            format(round(datos11b,4), nsmall = 4),
                            format(round(datos12b,4), nsmall = 4),
                            format(round(datos13,4), nsmall = 4),
                            format(round(datos14,9), nsmall = 9))
      datosdiatoms <- cbind(paste("",modelo_lambda, sep=" "),
                            format(round(datos22b,4), nsmall = 4),
                            format(round(datos22b,4), nsmall = 4),
                            format(round(datos23,4), nsmall = 4),
                            format(round(datos24,9), nsmall = 9))    

      # SAVE phyto_optics.dat file
      pdir<-"/Users/ealvarez/Datos/Res_C20_radtrans/"
      sink(paste(pdir,"phyto_optics/optics_phyto_recom_carbon.dat",sep=""))
      # HEADER
      cat("Computed in 10_Optics1_phyto_REcoM_carbon.R\n")
      cat("# aPH compiled from literature (absorption_spectra.csv) (m2 mgChla-1)\n")    
      cat("# aPS reconstructed based on Bidigare/Babin/Hickman (m2 mgChla-1)\n")
      cat("# b converted to mol-specific from Dutkiewicz2015 (m2 mmolC-1)\n")
      cat("# non-spectral bb converted to mol-specific from Dutkiewicz2015 (m2 molC-1)\n")
      cat("Format I4,3F10.4,F20.14\n")
      
      cat("*** Others ***\n")    
      write.table(datosothers, file=paste(pdir,"phyto_optics/optics_phyto_recom_carbon.dat",sep=""),
                  sep = "    ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
      #sep = "\t"
      cat("*** Diatom ***\n",file=paste(pdir,"phyto_optics/optics_phyto_recom_carbon.dat",sep=""), append = TRUE)
      write.table(datosdiatoms, file=paste(pdir,"phyto_optics/optics_phyto_recom_carbon.dat",sep=""),
                  sep = "    ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
      sink()    
      
#############################################      
      
      write.csv(datosothers,  file="/Users/ealvarez/Datos/Res_C20_radtrans/phyto_optics/datosothers_carbon.csv")
      write.csv(datosdiatoms, file="/Users/ealvarez/Datos/Res_C20_radtrans/phyto_optics/datosdiatoms_carbon.csv")
      
      
      
      
      
################        
### FIGURE 2 ###
################      
    
png(file=paste(path_figures, "Figure2_1_optics_PFTs_recom2_carbon.png", sep=""),width = 1210, height = 490, pointsize=22)                  
      
      layout(matrix(c(1,3,4,2,3,4),ncol=3,byrow=T))
      par(mar=c(4,4.5,1,1))
      par(oma=c(1,1,1,0))
      layout.show(n=4)
      lettersize=1.8
      lettersizeS=1.0
      
##### A  Small Phyto 
      plot(1,1,type="n", col=colores[j], xlab="",ylab=expression(paste(m^2," mg ",Chla^-1,sep="")),
           las=1, xlim=c(400,700),ylim=c(0,0.17), main="", cex.main=1.0,
           xaxt="n",yaxt="n", cex=0.5, pch=19, lwd=0.6, las=1)
      axis(1,at=seq(400,800,by=50), labels=T)
      axis(2,at=seq(0,0.18,by=0.03), labels=T, las=1)

      for (j in seq(1,6,by=1)){
        #datos$Genera[datos$Cluster==j]
        genero<-datos[(datos$bph_class==clusters[j]|datos$bph_class==clusters[j+1]|datos$bph_class==clusters[j+2]) & !is.na(datos$bph_class),22:102]
        genero[,lambdas_in<400]<-NA
        matpoints(cbind(lambdas_in),t(genero),type="b", col=colores[j],
                  las=1, xlim=c(400,700),ylim=c(0,0.18), main=nombres[j], cex.main=1.0,
                  xaxt="n",yaxt="n", cex=0.5, pch=19, lwd=0.6, las=1)
      }        
      points(lambdas_out,colMeans(resultado_aph[1:6,],na.rm=T),type="l",lwd=3, col="black", cex=0.5, pch=19)
      points(lambdas_out,colMeans(resultado_aps[1:6,],na.rm=T),type="l",lwd=3, lty=3, col="black", cex=0.5, pch=19)
      
      legend(x="topright",legend=c("types with CHLB","types with PEB","types with CHLC"), bty="n", pch=19, col=colores[c(1,2,3)]) 
      mtext(3,at=380,line=0.4,text=c("a) Small phytoplankton"),font=1, cex=lettersizeS, adj=0) 
      
    
##### B Diatoms
      j=7
      genero<-datos[datos$aph_class==clusters[j] & !is.na(datos$aph_class),22:102]
      genero[,lambdas_in<400]<-NA
      matplot(cbind(lambdas_in),t(genero),type="b", col=colores[j], ylab=expression(paste(m^2," mg ",Chla^-1,sep="")), xlab=expression(lambda),
              las=1, xlim=c(400,700),ylim=c(0,0.17), main="", cex.main=1.0,
              xaxt="n",yaxt="n", cex=0.5, pch=19, lwd=0.6, las=1)
      axis(1,at=seq(400,800,by=50), labels=T)
      axis(2,at=seq(0,0.18,by=0.03), labels=T, las=1)
      
      points(lambdas_out,resultado_aph[7,],type="l",lwd=3, col="black", cex=0.5, pch=19)
      points(lambdas_out,resultado_aps[7,],type="l",lwd=3, lty=3, col="black", cex=0.5, pch=19)
      
      legend(x="topright",legend=c("types with CHLC"), bty="n", pch=19, col=colores[c(3)]) 
      mtext(3,at=380,line=0.4,text=c("b) Diatoms"),font=1,cex=lettersizeS, adj=0) 
      
      
### Average mass-specific b (m2 mgC-1) to (m2 mmolC-1) ###
      #dir<-"/Users/ealvarez/Datos/Dat_Aph/"
      spectra<-read.csv(paste(dir,"csvs/scatter_mass_Dut2015.csv",sep=""), sep=";")
      #names(spectra)
      # Agregar por grupo
      resultado_bph  <- rbind(spectra$Small_B,  spectra$Diatom_B)
      resultado_bb   <- rbind(spectra$Small_BB, spectra$Diatom_BB)
      lambdas_small<-seq(400,700,by=25)
      rownames(resultado_bph)<-c("SmallPhyto","Diatoms")
      colnames(resultado_bph)<-lambdas_small        
      
### Average aph and aps ###
###########################      
      modelo_lambda <- c(400.0, 425.0, 450.0, 475.0,500.0, 525.0,550.0, 575.0, 600.0, 625.0, 650.0, 675.0,700.0)
      lambda_extremos<-sort(unique(c(modelo_lambda[-length(modelo_lambda)]+(diff(modelo_lambda)/2),400,700)))
      
      # Datos
      datos11<-colMeans(resultado_aph[1:6,],na.rm=T)
      datos21<-resultado_aph[7,]
      datos12<-colMeans(resultado_aps[1:6,],na.rm=T)
      datos22<-resultado_aps[7,]      
      
      datos13<-resultado_bph[1,]
      datos23<-resultado_bph[2,]
      datos14<-resultado_bb[1,]
      datos24<-resultado_bb[2,]
      
      datos11b<-rep(NA,length(modelo_lambda))
      datos12b<-rep(NA,length(modelo_lambda))
      datos21b<-rep(NA,length(modelo_lambda))
      datos22b<-rep(NA,length(modelo_lambda))
      datos13b<-rep(NA,length(modelo_lambda))
      datos23b<-rep(NA,length(modelo_lambda))
      
      for (i in c(1:length(modelo_lambda))){
        cuales<-which(lambdas_out>=lambda_extremos[i] & lambdas_out<lambda_extremos[i+1])
        datos11b[i]<-mean(datos11[cuales],na.rm=T)
        datos12b[i]<-mean(datos12[cuales],na.rm=T)
        datos21b[i]<-mean(datos21[cuales],na.rm=T)
        datos22b[i]<-mean(datos22[cuales],na.rm=T)
        datos13b[i]<-mean(datos13[cuales],na.rm=T)
        datos23b[i]<-mean(datos23[cuales],na.rm=T)
      }      
      
     
##### C Constituents Absorption 
###############################      
      # ## Constituents
      #mat <- matrix(c(1:6), ncol=6)
      paleta<-wes_palette(6, name = "Zissou1", type = "continuous")
      #image(mat, col=paleta, las=1, xaxt="n", yaxt="n", main="", cex.main=1.2) 
      # axis(4,at=seq(0+0.04,1-0.04, length=4),labels=c("water","phyto","cdom","nap"), cex.axis=1.2,las=1, font=2)        
      
      # Water
      pdir <- paste(OS,"Datos/Res_C20_radtrans/phyto_optics/",sep="")
      datos <- read.csv(paste(pdir,"absorb_Files_darwin3.csv", sep=""))
      colores<-c("lightblue","darkblue","goldenrod1","darkslategray4","seagreen3","seagreen1","salmon3","salmon1")
      colores<-c(paleta[1],paleta[1],paleta[4],paleta[6],paleta[3],paleta[3],paleta[5],paleta[5])
      plot(datos$lambda, datos$sw_ab.m.1., las=1, type="l",
           ylim=c(0,0.5),pch=19,lty=1, col=colores[1], ylab="", xlab=expression(lambda), lwd=2)
      abline(v=lambda_extremos, col="grey70")
      points(datos$lambda, datos$sw_ab.m.1., type="p",pch=19,col=colores[1])
      
      # ACDOM: Dut2015
      acdom= 0.18*exp(-0.021*(lambdas_out-450)) 
      points(lambdas_out, acdom, las=1, type="l", ylim=c(0,0.06),
             pch=19,lty=1, col=colores[3], ylab=expression(paste(m^2," mg ",Chla^-1,sep="")), xlab=expression(lambda),lwd=2)
      acdom= 0.18*exp(-0.021*(modelo_lambda-450))
      abline(v=lambda_extremos, col="grey70")
      points(modelo_lambda, acdom, las=1, type="p",pch=19,lty=1,col=colores[3],lwd=1)
      
      # APart
      apart= 0.016*exp(-0.013*(lambdas_out-440)) 
      points(lambdas_out, apart, las=1, type="l",pch=19,lty=1, col=colores[4], lwd=2)
      apart= 0.016*exp(-0.013*(modelo_lambda-440)) 
      points(modelo_lambda, apart, las=1, type="p",pch=19,lty=1, col=colores[4], lwd=1)
      
      # APhyto
      points(lambdas_out,datos11, type="l", las=1, col=colores[5], lwd=2)
      #points(lambdas_out,datos12, type="l", las=1, col=colores[6], lwd=2)
      points(lambdas_out,datos21, type="l", las=1, col=colores[7], lwd=2)
      #points(lambdas_out,datos22, type="l", las=1, col=colores[8], lwd=2)
      points(modelo_lambda, datos11b, las=1, col=colores[5], pch=19)
      #points(modelo_lambda, datos12b, las=1, bg=colores[6], pch=21)
      points(modelo_lambda, datos21b, las=1, col=colores[7], pch=19)
      #points(modelo_lambda, datos22b, las=1, bg=colores[8], pch=21)      
      
      mtext(side=3,at=380, line=0.4,text="c) Absorption",cex=lettersizeS, font=1, adj=0)
      #axis(4,at=seq(0,0.1,0.02),labels=T, las=1)
      #mtext(4,at=0.03,outer=F,text="m2 mmolC-1")         
      
      
      legend(x=425,y=0.5,pch=19,col=colores[c(1,3,4,5,7)],
             legend=c(expression(paste("water (",m^-1,")", sep="")),
                      expression(paste("CDOM (",m^2," ",mmolC^-1,")", sep="")),
                      expression(paste("NAP (",m^2," ",mmolC^-1,")", sep="")),
                      expression(paste("small phyt. (",m^2," ",mgChl^-1,")", sep="")),
                      expression(paste("diatoms (",m^2," ",mgChl^-1,")", sep=""))))
      ############################      
      
      
      
##### D Constituents Scatter
############################      
      par(mar=c(4,3.5,1,2))
      
      # Water
      pdir <- paste(OS, "Datos/Res_C20_radtrans/phyto_optics/", sep="")
      datos <- read.csv(paste(pdir,"absorb_Files_darwin3.csv", sep=""))
      plot(datos$lambda, datos$sw_st.m.1., las=1, type="l",
           ylim=c(0,0.4),pch=19,lty=1, col=colores[1], ylab="", xlab=expression(lambda), lwd=2)
      abline(v=lambda_extremos, col="grey70")
      points(datos$lambda, datos$sw_st.m.1., type="p",pch=19,col=colores[1])
      
      # BPART    
      apart= 0.345*(550/lambdas_out)^0.5
      points(lambdas_out, apart, las=1, type="l",pch=19,lty=1, col=colores[4], ylim=c(0,0.1),ylab="", xlab=expression(lambda), lwd=2, yaxt="n")
      apart= 0.345*(550/modelo_lambda)^0.5
      points(modelo_lambda, apart, las=1, type="p",pch=19,lty=1, col=colores[4])
      
      # BPhyto
      points(modelo_lambda,datos13, type="l", las=1, col=colores[5],lwd=2, ylim=c(0,0.01))
      points(modelo_lambda,datos23, type="l", las=1, col=colores[7],lwd=2)
      abline(v=lambda_extremos, col="grey70")      
      points(modelo_lambda, datos13, las=1, col=colores[5], pch=19)
      points(modelo_lambda, datos23, las=1, col=colores[7], pch=19)          
      
      mtext(side=3,at=380,line=0.4,text="d) Scattering",cex=lettersizeS, font=1, adj=0)
      #axis(4,at=seq(0,0.1,0.02),labels=T, las=1)
      #mtext(4,at=0.03,outer=F,text="m2 mmolC-1")         
      
      legend(x=450,y=0.3,pch=19,col=colores[c(1,4,5,7)],
             legend=c(expression(paste("water (",m^-1,")", sep="")),
                      expression(paste("NAP (",m^2," ",mmolC^-1,")", sep="")),
                      expression(paste("small phyt. (",m^2," ",mmolC^-1,")", sep="")),
                      expression(paste("diatoms (",m^2," ",mmolC^-1,")", sep="")) ))
      #########################      
      dev.off()      
      ################      
      
    