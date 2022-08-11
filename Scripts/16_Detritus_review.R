source("Scripts/00_Summary_EDITME.R")
source("Scripts/functions_misc.R")

p_dir<-paste(global_path,"Phyto_optics/",sep="")
path_figures<-paste(global_path,"Figures/",sep="")
datos <- read.csv(paste(p_dir,"absorb_Files_darwin3.csv", sep=""))

### DETRITUS          
  png(file=paste(path_figures,"FigureS7_Absorption_Scattering_detritus.png",sep=""),
      width=1150, height =400, pointsize=17, family="Helvetica") 
         par(family="")
         par(mfcol=c(1,2))
         par(mar=c(4,4.2,1,1))
         cexletter=1.8
         
    # Absorption         
    #RECOM_CALC_APART=F, Dut2015  
    plot(datos$lambda, datos[,7]/1.00e-15, las=1, type="b",
            pch=19,lty=1, col="black", ylim=c(0,0.04),
            ylab=expression(paste(m^2," ",mmolC^-1,sep="")),
            xlab=expression(lambda), mgp=c(2.8,1,0), lwd=2, yaxs="i")
            text(x=700,y=0.035,pos=2, labels=expression(paste("a*"[NAP],(lambda)," = ",m^2, particle^-1, "/ 1·",10^-15," ",mmolC," ", particle^-1, sep="")), cex=1.1, col="black")
    # RECOM_CALC_APART=T, Gre2017 
    apart= 0.016*exp(-0.013*(datos$lambda-440)) 
    points(datos$lambda, apart, las=1, type="b",pch=19,lty=1, col="red4", lwd=2)
    text(x=700,y=0.028,pos=2,labels=expression(paste("a*"[NAP],(lambda)," = 0.016 · ",e^-0.013(lambda-440))), cex=1.1, col="red4")
    text(x=410,y=0.0, pos=3, labels="a)", cex=cexletter, font=2)
    
    # Scattering
    #RECOM_CALC_APART=F, Dut2015  
    plot(datos$lambda, datos[,8]/1.00e-15, las=1, type="b",
         pch=19,lty=1, col="black", ylim=c(0,0.7),
         ylab=expression(paste(m^2," ",mmolC^-1,sep="")),
         xlab=expression(lambda), mgp=c(2.8,1,0), lwd=2, yaxs="i")
         #media <- sum(25*(datos[,8]/1.06e-13))/(700-400)
         text(x=700,y=0.6,pos=2,labels=expression(paste("b*"[NAP],(lambda)," = ",m^2, particle^-1, "/ 1·",10^-15," ",mmolC," ", particle^-1)), cex=1.1, col="black")
    # RECOM_CALC_APART=T, Gre2017    
    apart= 0.345*(550/datos$lambda)^0.5
    points(datos$lambda, apart, las=1, type="b",pch=19,lty=1, col="red4", lwd=2)
    text(x=700,y=0.5,pos=2,labels=expression(paste("b*"[NAP],(lambda)," = 0.345 · (550/", lambda,"",")"^0.5)), cex=1.1, col="red4")
    text(x=410,y=0.0, pos=3, labels="b)", cex=cexletter, font=2)
    
    dev.off()
    