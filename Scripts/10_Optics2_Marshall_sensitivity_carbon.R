
###########################################################
### SCALE OPTICAL PROPERTIES to mean(aph) in literature ###                                                  
###########################################################

## Starting from the initial: optics_phyto_recom_carbon.dat
if (Sys.info()['sysname']=="Windows") {OS<-"C:/"} else {OS<-paste("/Users/",Sys.info()['user'],"/", sep="")}
library(unikn)
library(plot3D)

dir<-paste(OS, "Datos/Res_C20_radtrans/phyto_optics/",sep="")

    
###################################################################    
###    Combinations aph and aps for small phyto and diatoms    ####
### all of them keep alpha in 0.14 (small phyto) and 0.19 (diatoms)
################################################################### 
    
  mycol<-usecol(pal = c(rev(pal_karpfenblau), pal_peach),n=20, alpha = 0.75)
  par(mfrow=c(1,3))
  par(mar=c(2,4,2,2))
  par(oma=c(1,1,1,1))
  layout.show(n=3)
        
        # ALL combinations alpha
        ab <- seq(0.0030,0.041, length=100)    
        qy <- seq(2.1e-5, 4.0e-4, length=100)  
        res <- matrix(c(rep(ab, times=length(qy))*rep(qy,each=length(ab))*86400),ncol=length(ab), byrow=TRUE)
        colnames(res) <- ab
        rownames(res) <- qy
        #range(res,na.rm=T)
        res[res>=0.9]<-0.9
        res[res<=0.009]<-0.009
        wd <- c(12.5,rep(25,11),12.5)
    
    ###################
    ### 1 optical group: small phyto and diatoms are equal
    ###################    
        image2D(x=qy, y=ab, z=res, las=1, xlab="", ylab="a* (m2 mgChla)", zlim=c(0.009,0.9),mgp=c(3.1,0.8,0), col=mycol, colkey=F)
        contour(x=qy, y=ab, z=res, levels=c(0.14), col="darkcyan", add=TRUE, lty=1, lwd=3, labcex=2.0, font=2)
        contour(x=qy, y=ab, z=res, levels=c(0.19), col="firebrick", add=TRUE, lty=1, lwd=3, labcex=2.0, font=2)
        mtext(1,at=0.0002, line=2, outer=FALSE,text="QY (mmolC J-1)", cex=0.8)
        
        aph1 <- c(0.007,  0.015,  0.023)
        aph2 <- c(0.03,   0.023,  0.015)
        tabla_ap<-matrix(NA,ncol=2,nrow=8)
        tabla_ps<-matrix(NA,ncol=2,nrow=8)
        tabla_qy<-matrix(NA,ncol=2,nrow=8)
        
              ### Diatoms
              datos <- read.csv(paste(dir,"datosdiatoms_carbon.csv", sep=""))
              lambda <- datos[,2]
              alpha <- 0.19/86400
                  ## APS
                  dat <- datos[,4]
                  datos21 <- dat*(aph1[1]/(sum(wd*dat)/(700-400)))
                  media11 <- sum(wd*datos21)/(700-400)
                  QYm <- alpha / media11
                  points(x=QYm, y=media11, pch=21, cex=2.0,        bg="coral")    
                  tabla_ps[1:3,2]<-media11
                  tabla_qy[1:3,2]<-QYm
                  
                      # AP
                      dat <- datos[,3]
                      datos21 <- dat*(aph2[1]/(sum(wd*dat)/(700-400)))
                      media11 <- sum(wd*datos21)/(700-400)
                      points(x=QYm, y=media11, pch=22, cex=2.4,  col="coral", lwd=2)             
                      tabla_ap[3,2]<-media11
                      datos21 <- dat*(aph2[2]/(sum(wd*dat)/(700-400)))
                      media11 <- sum(wd*datos21)/(700-400)
                      points(x=QYm, y=media11, pch=22, cex=1.4,  col="coral", lwd=2)                      
                      tabla_ap[2,2]<-media11
                      datos21 <- dat*(aph2[3]/(sum(wd*dat)/(700-400)))
                      media11 <- sum(wd*datos21)/(700-400)
                      points(x=QYm, y=media11, pch=22, cex=1.4,  col="coral", lwd=2)              
                      tabla_ap[1,2]<-media11
                      
              dat <- datos[,4]
              datos21 <- dat*(aph1[2]/(sum(wd*dat)/(700-400)))
              media11 <- sum(wd*datos21)/(700-400)
              QYm <- alpha / media11
              points(x=QYm, y=media11, pch=21, cex=1.0,        bg="coral")              
              tabla_ps[4:6,2]<-media11
              tabla_qy[4:6,2]<-QYm
              
                      # AP
                      dat <- datos[,3]
                      datos21 <- dat*(aph2[1]/(sum(wd*dat)/(700-400)))
                      media11 <- sum(wd*datos21)/(700-400)
                      points(x=QYm, y=media11, pch=22, cex=1.4,  col="coral", lwd=2)             
                      tabla_ap[6,2]<-media11
                      datos21 <- dat*(aph2[2]/(sum(wd*dat)/(700-400)))
                      media11 <- sum(wd*datos21)/(700-400)
                      points(x=QYm, y=media11, pch=22, cex=1.4,  col="coral", lwd=2)              
                      tabla_ap[5,2]<-media11
                      datos21 <- dat*(aph2[3]/(sum(wd*dat)/(700-400)))
                      media11 <- sum(wd*datos21)/(700-400)
                      points(x=QYm, y=media11, pch=22, cex=1.4,  col="coral", lwd=2)                      
                      tabla_ap[4,2]<-media11
                      
              dat <- datos[,4]
              datos21 <- dat*(aph1[3]/(sum(wd*dat)/(700-400)))
              media11 <- sum(wd*datos21)/(700-400)
              QYm <- alpha / media11
              points(x=QYm, y=media11, pch=21, cex=1.0,        bg="coral")              
              tabla_ps[7:8,2]<-media11
              tabla_qy[7:8,2]<-QYm
              
                      # AP
                      dat <- datos[,3]
                      datos21 <- dat*(aph2[1]/(sum(wd*dat)/(700-400)))
                      media11 <- sum(wd*datos21)/(700-400)
                      points(x=QYm, y=media11, pch=22, cex=1.4,  col="coral", lwd=2)              
                      tabla_ap[8,2]<-media11
                      datos21 <- dat*(aph2[2]/(sum(wd*dat)/(700-400)))
                      media11 <- sum(wd*datos21)/(700-400)
                      points(x=QYm, y=media11, pch=22, cex=1.4,  col="coral", lwd=2)                         
                      tabla_ap[7,2]<-media11
              
              
       ### Small Phyto
       datos <- read.csv(paste(dir,"datosothers_carbon.csv", sep=""))
       lambda <- datos[,2]
       alpha <- 0.14/86400 
              
              ## APS
              dat <- datos[,4]
              datos21 <- dat*(aph1[1]/(sum(wd*dat)/(700-400)))
              media11 <- sum(wd*datos21)/(700-400)
              QYm <- alpha / media11
              points(x=QYm, y=media11, pch=21, cex=2.0,        bg="yellowgreen")    
              tabla_ps[1:3,1]<-media11
              tabla_qy[1:3,1]<-QYm
              
                      # AP
                      dat <- datos[,3]
                      datos21 <- dat*(aph2[1]/(sum(wd*dat)/(700-400)))
                      media11 <- sum(wd*datos21)/(700-400)
                      points(x=QYm, y=media11, pch=22, cex=2.4,  col="yellowgreen", lwd=2)             
                      tabla_ap[3,1]<-media11
                      datos21 <- dat*(aph2[2]/(sum(wd*dat)/(700-400)))
                      media11 <- sum(wd*datos21)/(700-400)
                      points(x=QYm, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)                      
                      tabla_ap[2,1]<-media11
                      datos21 <- dat*(aph2[3]/(sum(wd*dat)/(700-400)))
                      media11 <- sum(wd*datos21)/(700-400)
                      points(x=QYm, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)              
                      tabla_ap[1,1]<-media11
                      
              dat <- datos[,4]
              datos21 <- dat*(aph1[2]/(sum(wd*dat)/(700-400)))
              media11 <- sum(wd*datos21)/(700-400)
              QYm <- alpha / media11
              points(x=QYm, y=media11, pch=21, cex=1.0,        bg="yellowgreen")              
              tabla_ps[4:6,1]<-media11
              tabla_qy[4:6,1]<-QYm 
              
                    # AP
                    dat <- datos[,3]
                    datos21 <- dat*(aph2[1]/(sum(wd*dat)/(700-400)))
                    media11 <- sum(wd*datos21)/(700-400)
                    points(x=QYm, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)             
                    tabla_ap[6,1]<-media11
                    datos21 <- dat*(aph2[2]/(sum(wd*dat)/(700-400)))
                    media11 <- sum(wd*datos21)/(700-400)
                    points(x=QYm, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)              
                    tabla_ap[5,1]<-media11
                    datos21 <- dat*(aph2[3]/(sum(wd*dat)/(700-400)))
                    media11 <- sum(wd*datos21)/(700-400)
                    points(x=QYm, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)                      
                    tabla_ap[4,1]<-media11
                    
              dat <- datos[,4]
              datos21 <- dat*(aph1[3]/(sum(wd*dat)/(700-400)))
              media11 <- sum(wd*datos21)/(700-400)
              QYm <- alpha / media11
              points(x=QYm, y=media11, pch=21, cex=1.0,     bg="yellowgreen")              
              tabla_ps[7:8,1]<-media11
              tabla_qy[7:8,1]<-QYm
              
                    # AP
                    dat <- datos[,3]
                    datos21 <- dat*(aph2[1]/(sum(wd*dat)/(700-400)))
                    media11 <- sum(wd*datos21)/(700-400)
                    points(x=QYm, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)              
                    tabla_ap[8,1]<-media11
                    datos21 <- dat*(aph2[2]/(sum(wd*dat)/(700-400)))
                    media11 <- sum(wd*datos21)/(700-400)
                    points(x=QYm, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)                         
                    tabla_ap[7,1]<-media11
                    
            
         ## Segments 
         segments(x0=tabla_qy[1,1],x1=tabla_qy[1,2],y0=tabla_ps[1,1],y1=tabla_ps[1,2])         
         segments(x0=tabla_qy[4,1],x1=tabla_qy[4,2],y0=tabla_ps[4,1],y1=tabla_ps[4,2])         
         segments(x0=tabla_qy[7,1],x1=tabla_qy[7,2],y0=tabla_ps[7,1],y1=tabla_ps[7,2])         
         
         tabla_ap_1group<-tabla_ap
         tabla_ps_1group<-tabla_ps
         tabla_qy_1group<-tabla_qy  
         legend(x="bottomleft",legend=c("ap", "ap_ps"), col="black", bty = "n", pch=c(22,21), pt.cex=2)
         
         # Labels
         alpha <- 0.13/86400 
         numero<-c("a","b","c")
         QY1 <- alpha / aph1
         points(x=QY1-(QY1*0.05), y=aph1, pch=numero, cex=1.0, col="black", lwd=2)
         numero<-c("3","2","1")
         points(x=rep(3.5e-4,3), y=aph2, pch=numero, cex=1.0, col="black", lwd=2)
    ####################         
         
    ####################         
    #### 2 optical group: small phyto have higher aph and aps than diatoms
    ###################         
           image2D(x=qy, y=ab, z=res, las=1, xlab="", ylab="a* (m2 mgChla)",
                   zlim=c(0.009,0.9), mgp=c(3.1,0.8,0), col=mycol, colkey=F)
           contour(x=qy, y=ab, z=res, levels=c(0.14), col="darkcyan", add=TRUE, lty=1, lwd=3, labcex=2.0, font=2)
           contour(x=qy, y=ab, z=res, levels=c(0.19), col="firebrick", add=TRUE, lty=1, lwd=3, labcex=2.0, font=2)
           mtext(1,at=0.0002, line=2, outer=FALSE,text="QY (mmolC J-1)", cex=0.8)
                    
             aph1 <- c(0.007,   0.0143,   0.021)
             aph2 <- c(0.0143,  0.0106,   0.007,  0.018,  0.021, 0.025, 0.029)
             tabla_ap<-matrix(NA,ncol=2,nrow=9)
             tabla_ps<-matrix(NA,ncol=2,nrow=9)
             tabla_qy<-matrix(NA,ncol=2,nrow=9)
             
             ### Diatoms
             datos <- read.csv(paste(dir,"datosdiatoms_carbon.csv", sep=""))
             lambda <- datos[,2]
             alpha <- 0.19/86400
                      
                      ## APS
                      dat <- datos[,4]
                      # Initial
                      media1 <- sum(wd*dat)/(700-400)
                      QYi <- alpha / media1
                      points(x=QYi, y=media1, pch=21, cex=2.0,         bg="coral")                  
                      tabla_ps[4:6,2]<-media1
                      tabla_qy[4:6,2]<-QYi
                      
                            # AP
                            dat <- datos[,3]
                            datos21 <- dat*(aph2[1]/(sum(wd*dat)/(700-400)))
                            media11 <- sum(wd*datos21)/(700-400)
                            points(x=QYi, y=media11, pch=22, cex=2.4,  col="coral", lwd=2)             
                            tabla_ap[4,2]<-media11
                            datos21 <- dat*(aph2[4]/(sum(wd*dat)/(700-400)))
                            media11 <- sum(wd*datos21)/(700-400)
                            points(x=QYi, y=media11, pch=22, cex=1.4,  col="coral", lwd=2)                      
                            tabla_ap[5,2]<-media11
                            datos21 <- dat*(aph2[5]/(sum(wd*dat)/(700-400)))
                            media11 <- sum(wd*datos21)/(700-400)
                            points(x=QYi, y=media11, pch=22, cex=1.4,  col="coral", lwd=2)              
                            tabla_ap[6,2]<-media11
                            
                      dat <- datos[,4]
                      datos21 <- dat*(aph1[1]/(sum(wd*dat)/(700-400)))
                      media11 <- sum(wd*datos21)/(700-400)
                      QYm <- alpha / media11
                      points(x=QYm, y=media11, pch=21, cex=1.0,        bg="coral")              
                      tabla_ps[1:3,2]<-media11
                      tabla_qy[1:3,2]<-QYm  
                      
                            # AP
                            dat <- datos[,3]
                            datos21 <- dat*(aph2[1]/(sum(wd*dat)/(700-400)))
                            media11 <- sum(wd*datos21)/(700-400)
                            points(x=QYm, y=media11, pch=22, cex=1.4,  col="coral", lwd=2)             
                            tabla_ap[3,2]<-media11
                            datos21 <- dat*(aph2[2]/(sum(wd*dat)/(700-400)))
                            media11 <- sum(wd*datos21)/(700-400)
                            points(x=QYm, y=media11, pch=22, cex=1.4,  col="coral", lwd=2)              
                            tabla_ap[2,2]<-media11
                            datos21 <- dat*(aph2[3]/(sum(wd*dat)/(700-400)))
                            media11 <- sum(wd*datos21)/(700-400)
                            points(x=QYm, y=media11, pch=22, cex=1.4,  col="coral", lwd=2)                      
                            tabla_ap[1,2]<-media11
                            
                      dat <- datos[,4]
                      datos21 <- dat*(aph1[3]/(sum(wd*dat)/(700-400)))
                      media11 <- sum(wd*datos21)/(700-400)
                      QYm <- alpha / media11
                      points(x=QYm, y=media11, pch=21, cex=1.0,        bg="coral")              
                      tabla_ps[7:9,2]<-media11
                      tabla_qy[7:9,2]<-QYm  
                      
                            # AP
                            dat <- datos[,3]
                            datos21 <- dat*(aph2[5]/(sum(wd*dat)/(700-400)))
                            media11 <- sum(wd*datos21)/(700-400)
                            points(x=QYm, y=media11, pch=22, cex=1.4,  col="coral", lwd=2)              
                            tabla_ap[7,2]<-media11
                            datos21 <- dat*(aph2[6]/(sum(wd*dat)/(700-400)))
                            media11 <- sum(wd*datos21)/(700-400)
                            points(x=QYm, y=media11, pch=22, cex=1.4,  col="coral", lwd=2)                         
                            tabla_ap[8,2]<-media11
                            datos21 <- dat*(aph2[7]/(sum(wd*dat)/(700-400)))
                            media11 <- sum(wd*datos21)/(700-400)
                            points(x=QYm, y=media11, pch=22, cex=1.4,  col="coral", lwd=2)                        
                            tabla_ap[9,2]<-media11
                      
      ### Small Phyto
      datos <- read.csv(paste(dir,"datosothers_carbon.csv", sep=""))
      lambda <- datos[,2]
      alpha <- 0.14/86400 
      aph1 <- c(0.009,    0.0166,   0.025)
      aph2 <- c(0.02095,  0.0160,   0.012,  0.026,  0.030, 0.035, 0.040)
                           
              ## APS
              dat <- datos[,4]
              # Initial
              media1 <- sum(wd*dat)/(700-400)
              QYi <- alpha / media1
              points(x=QYi, y=media1, pch=21, cex=2.0,         bg="yellowgreen")                  
              tabla_ps[4:6,1]<-media1
              tabla_qy[4:6,1]<-QYi
              
                    # AP
                    dat <- datos[,3]
                    datos21 <- dat*(aph2[1]/(sum(wd*dat)/(700-400)))
                    media11 <- sum(wd*datos21)/(700-400)
                    points(x=QYi, y=media11, pch=22, cex=2.4,  col="yellowgreen", lwd=2)             
                    tabla_ap[4,1]<-media11
                    datos21 <- dat*(aph2[4]/(sum(wd*dat)/(700-400)))
                    media11 <- sum(wd*datos21)/(700-400)
                    points(x=QYi, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)                      
                    tabla_ap[5,1]<-media11
                    datos21 <- dat*(aph2[5]/(sum(wd*dat)/(700-400)))
                    media11 <- sum(wd*datos21)/(700-400)
                    points(x=QYi, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)              
                    tabla_ap[6,1]<-media11
                    
              dat <- datos[,4]
              datos21 <- dat*(aph1[1]/(sum(wd*dat)/(700-400)))
              media11 <- sum(wd*datos21)/(700-400)
              QYm <- alpha / media11
              points(x=QYm, y=media11, pch=21, cex=1.0,        bg="yellowgreen")              
              tabla_ps[1:3,1]<-media11
              tabla_qy[1:3,1]<-QYm  
              
                    # AP
                    dat <- datos[,3]
                    datos21 <- dat*(aph2[1]/(sum(wd*dat)/(700-400)))
                    media11 <- sum(wd*datos21)/(700-400)
                    points(x=QYm, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)             
                    tabla_ap[3,1]<-media11
                    datos21 <- dat*(aph2[2]/(sum(wd*dat)/(700-400)))
                    media11 <- sum(wd*datos21)/(700-400)
                    points(x=QYm, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)              
                    tabla_ap[2,1]<-media11
                    datos21 <- dat*(aph2[3]/(sum(wd*dat)/(700-400)))
                    media11 <- sum(wd*datos21)/(700-400)
                    points(x=QYm, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)                      
                    tabla_ap[1,1]<-media11
                    
              dat <- datos[,4]
              datos21 <- dat*(aph1[3]/(sum(wd*dat)/(700-400)))
              media11 <- sum(wd*datos21)/(700-400)
              QYm <- alpha / media11
              points(x=QYm, y=media11, pch=21, cex=1.0,        bg="yellowgreen")              
              tabla_ps[7:9,1]<-media11
              tabla_qy[7:9,1]<-QYm  
              
                      # AP
                      dat <- datos[,3]
                      datos21 <- dat*(aph2[5]/(sum(wd*dat)/(700-400)))
                      media11 <- sum(wd*datos21)/(700-400)
                      points(x=QYm, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)              
                      tabla_ap[7,1]<-media11
                      datos21 <- dat*(aph2[6]/(sum(wd*dat)/(700-400)))
                      media11 <- sum(wd*datos21)/(700-400)
                      points(x=QYm, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)                         
                      tabla_ap[8,1]<-media11
                      datos21 <- dat*(aph2[7]/(sum(wd*dat)/(700-400)))
                      media11 <- sum(wd*datos21)/(700-400)
                      points(x=QYm, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)                        
                      tabla_ap[9,1]<-media11      
                          
            
           ## Segments 
           segments(x0=tabla_qy[1,1],x1=tabla_qy[1,2],y0=tabla_ps[1,1],y1=tabla_ps[1,2])         
           segments(x0=tabla_qy[4,1],x1=tabla_qy[4,2],y0=tabla_ps[4,1],y1=tabla_ps[4,2])         
           segments(x0=tabla_qy[7,1],x1=tabla_qy[7,2],y0=tabla_ps[7,1],y1=tabla_ps[7,2])         
              
           tabla_ap_2group<-tabla_ap
           tabla_ps_2group<-tabla_ps
           tabla_qy_2group<-tabla_qy              

           legend(x="bottomleft",legend=c("ap", "ap_ps"), col="black", bty = "n", pch=c(22,21), pt.cex=2)
              
           # Labels
           alpha <- 0.13/86400 
           numero<-c("d","e","f")
           QY1 <- alpha / aph1
           points(x=QY1-(QY1*0.05), y=aph1, pch=numero, cex=1.0, col="black", lwd=2)
           #numero<-c("3","2","1")
           numero<-c("3","2","1")
           points(x=rep(3.5e-4,3), y=aph2[7:5], pch=numero, cex=1.0, col="black", lwd=2)
    ####################               
       
    #############       
    #### extremes: small phyto and diatoms have aph and aps in the extremes of the range
    #############           
           image2D(x=qy, y=ab, z=res, las=1, xlab="", ylab="a* (m2 mgChla)", zlim=c(0.009,0.9), mgp=c(3.1,0.8,0), col=mycol)
           contour(x=qy, y=ab, z=res, levels=c(0.14), col="darkcyan", add=TRUE, lty=1, lwd=3, labcex=2.0, font=2)
           contour(x=qy, y=ab, z=res, levels=c(0.19), col="firebrick", add=TRUE, lty=1, lwd=3, labcex=2.0, font=2)
           mtext(1,at=0.0002, line=2, outer=FALSE,text="QY (mmolC J-1)", cex=0.8)
           
           aph1 <- c(0.007,   0.0143,   0.021)
           aph2 <- c(0.0143,  0.0106,   0.007,  0.018,  0.021, 0.025, 0.029)
           tabla_ap<-matrix(NA,ncol=2,nrow=3)
           tabla_ps<-matrix(NA,ncol=2,nrow=3)
           tabla_qy<-matrix(NA,ncol=2,nrow=3)
           
             ### Diatoms
             datos <- read.csv(paste(dir,"datosdiatoms_carbon.csv", sep=""))
             lambda <- datos[,2]
             alpha <- 0.19/86400

                 ## APS
                 dat <- datos[,4]
                 # Initial
                 media1 <- sum(wd*dat)/(700-400)
                 QYi <- alpha / media1
                 points(x=QYi, y=media1, pch=21, cex=1.0,         bg="coral")                  
                 tabla_ps[1,2]<-media1
                 tabla_qy[1,2]<-QYi
                 
                         # AP
                         dat <- datos[,3]
                         datos21 <- dat*(aph2[1]/(sum(wd*dat)/(700-400)))
                         media11 <- sum(wd*datos21)/(700-400)
                         points(x=QYi, y=media11, pch=22, cex=1.4, col="coral", lwd=2)             
                         tabla_ap[1,2]<-media11
                         
                 dat <- datos[,4]
                 datos21 <- dat*(aph1[1]/(sum(wd*dat)/(700-400)))
                 media11 <- sum(wd*datos21)/(700-400)
                 QYm <- alpha / media11
                 points(x=QYm, y=media11, pch=21, cex=1.0,        bg="coral")              
                 tabla_ps[2:3,2]<-media11
                 tabla_qy[2:3,2]<-QYm  
                 
                         # AP
                         dat <- datos[,3]
                         datos21 <- dat*(aph2[3]/(sum(wd*dat)/(700-400)))
                         media11 <- sum(wd*datos21)/(700-400)
                         points(x=QYm, y=media11, pch=22, cex=1.4,  col="coral", lwd=2)             
                         tabla_ap[2:3,2]<-media11
                         
            
               ### Small Phyto
               datos <- read.csv(paste(dir,"datosothers_carbon.csv", sep=""))
               lambda <- datos[,2]
               alpha <- 0.14/86400 
               
               aph1 <- c(0.009,    0.0166,   0.025)
               aph2 <- c(0.02095,  0.0160,   0.012,  0.024,  0.030, 0.035, 0.040)
               
                 ## APS
                 dat <- datos[,4]
                 # Initial
                 media1 <- sum(wd*dat)/(700-400)
                 QYi <- alpha / media1
                 points(x=QYi, y=media1, pch=21, cex=1.0,         bg="yellowgreen")                  
                 tabla_ps[3,1]<-media1
                 tabla_qy[3,1]<-QYi
                 
                       # AP
                       dat <- datos[,3]
                       datos21 <- dat*(aph2[1]/(sum(wd*dat)/(700-400)))
                       media11 <- sum(wd*datos21)/(700-400)
                       points(x=QYi, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)             
                       tabla_ap[3,1]<-media11
                       
                 dat <- datos[,4]
                 datos21 <- dat*(aph1[3]/(sum(wd*dat)/(700-400)))
                 media11 <- sum(wd*datos21)/(700-400)
                 QYm <- alpha / media11
                 points(x=QYm, y=media11, pch=21, cex=1.0,        bg="yellowgreen")              
                 tabla_ps[1:2,1]<-media11
                 tabla_qy[1:2,1]<-QYm  
                 
                       # AP
                       dat <- datos[,3]
                       datos21 <- dat*(aph2[7]/(sum(wd*dat)/(700-400)))
                       media11 <- sum(wd*datos21)/(700-400)
                       points(x=QYm, y=media11, pch=22, cex=1.4,  col="yellowgreen", lwd=2)             
                       tabla_ap[1:2,1]<-media11
                       
       
           ## Segments 
           segments(x0=tabla_qy[1,1],x1=tabla_qy[1,2],y0=tabla_ps[1,1],y1=tabla_ps[1,2])         
           segments(x0=tabla_qy[2,1],x1=tabla_qy[2,2],y0=tabla_ps[2,1],y1=tabla_ps[2,2])         
           segments(x0=tabla_qy[3,1],x1=tabla_qy[3,2],y0=tabla_ps[3,1],y1=tabla_ps[3,2])         
           
           tabla_ap_3group<-tabla_ap
           tabla_ps_3group<-tabla_ps
           tabla_qy_3group<-tabla_qy              
           
           legend(x="bottomleft",legend=c("ap", "ap_ps"), col="black", bty = "n", pch=c(22,21), pt.cex=2)
           
           # Labels
           alpha <- 0.13/86400 
           numero<-c("g","h","i")
           QY1 <- alpha / c(0.023,0.018,0.014)
           QY1<-QY1-(QY1*0.05)
           QY1[2]<-2.0e-04
           points(x=QY1, y=c(0.023,0.018,0.014), pch=numero, cex=1.0, col="black", lwd=2)
           #numero<-c("3","2","1")
           #points(x=rep(3.5e-4,3), y=aph2, pch=numero, cex=1.0, col="black", lwd=2)
    #############           
  
###############################    
    
 
### Tabla de valores para sensitivity analisis    

    tabla_ap<-rbind(tabla_ap_1group,tabla_ap_2group,tabla_ap_3group)
    tabla_ps<-rbind(tabla_ps_1group,tabla_ps_2group,tabla_ps_3group)
    tabla_qy<-rbind(tabla_qy_1group,tabla_qy_2group,tabla_qy_3group)       

    optical_groups<-c(rep(1,8),rep(2,9),rep(2,3))
    ps_types <-c(rep(c("a","b","c"),each=3)[-9], rep(c("d","e","f"),each=3),  c("g","h","i"))
    ap_types <-c(rep(c("1","2","3"),times=3)[-7],rep(c("1","2","3"),times=3), c("1","1","1"))

    tabla_total<-data.frame(optical_groups,ps_types,ap_types,
                            round(tabla_ps,4),round(tabla_ap,4),round(tabla_qy,6))
    
    colnames(tabla_total)<-c("OG","ps_type","ap_type",
                             "ps_phy","ps_dia",
                             "ap_phy","ap_dia",
                             "qy_phy","qy_dia")
    
    write.csv(tabla_total,file=paste(dir,"tabla_total_sensitivity.csv", sep=""))
    save(tabla_total,file=paste(dir,"tabla_total_sensitivity.RData", sep=""))
########################
    
    
    
    
#############################################################    
## Create all phyto_optics.dat for sensitivity analysis    ##    
#############################################################    
    
    dir<-"/Users/ealvarez/Datos/Res_C20_radtrans/phyto_optics/"
    load(paste(dir,"tabla_total_sensitivity.RData", sep=""))
    
    png(file="/Users/ealvarez/Datos/Res_C20_radtrans/phyto_optics/sensitivity/Figure_sensitivity.png", width = 1200, height = 800, units = "px", pointsize = 15, bg = "white")  
    
      layout(matrix(c(1:20),ncol=4,byrow=T),widths=c(1,1,1,1)) 
      par(mar=c(1,4,0,1)) 
      par(oma=c(1,1,1,1)) 
      layout.show(n=20)   
      lamda_new <- c(400,425,450,475,500,525,550,575,600,625,650,675,700)
      
      matrix_smallph_ap<-matrix(NA,ncol=20,nrow=length(lamda_new))
      matrix_smallph_ps<-matrix(NA,ncol=20,nrow=length(lamda_new))
      matrix_diatoms_ap<-matrix(NA,ncol=20,nrow=length(lamda_new))
      matrix_diatoms_ps<-matrix(NA,ncol=20,nrow=length(lamda_new))
      
        
        ### Small Phyto
        datos <- read.csv(paste(dir,"datosothers_carbon.csv", sep=""))
        lambda <- datos[,2]
          # AP
          datAP_phy <- datos[,3]
          # AP_PS      
          datPS_phy <- datos[,4]
          datos13<-datos[,5]
          media3 <- sum(wd*datos13)/(700-400)
          datos14<-datos[,6]     
          
        ### Diatoms
        datos <- read.csv(paste(dir,"datosdiatoms_carbon.csv", sep=""))
          # AP
          datAP_dia <- datos[,3]
          # AP_PS      
          datPS_dia <- datos[,4]
          datos23<-datos[,5]
          media3 <- sum(wd*datos23)/(700-400)
          datos24<-datos[,6]       
 
    for (i in c(1:nrow(tabla_total))){  
      
      # Small phyto
      plot(lambda,datAP_phy, type="n", ylim=c(0,0.10), pch=19, lty=1, lwd=2, las=1,
          col="cyan2", ylab="", xlab=expression(lambda), mgp=c(3,1,0), yaxs="i") 
          # AP
          media1 <-  sum(wd*datAP_phy)/(700-400)
          datos11 <- datAP_phy*(tabla_total[i,6]/media1)
          # AP_PS
          media2 <-  sum(wd*datPS_phy)/(700-400)
          datos12 <- datPS_phy*(tabla_total[i,4]/media2)
          
          #datos11<-datos11/datos11[12]*datos12[12]
          points(lambda,datos11, type="l", pch=19, lty=1, col="darkcyan", lwd=2)
          media11 <- sum(wd*datos11)/(700-400)
          max11<-datos11[lambda==450]      
          points(lambda,datos12, type="l", pch=19, lty=1, col="cyan2", lwd=2)
          media12 <- sum(wd*datos12)/(700-400)
          max12<-datos12[lambda==450]
          
          matrix_smallph_ap[,i]<-datos11
          matrix_smallph_ps[,i]<-datos12
          
          legend(x="topright",legend=c(paste("<aPH> = ",round(media11,3),sep=""),paste("aPH(450) = ",round(max11,3),sep=""),
                                       paste("<aPS> = ",round(media12,3),sep=""),paste("aPS(450) = ",round(max12,3),sep="")),
                                       text.col=c("darkcyan","darkcyan","cyan2","cyan2"), bty = "n") 
          
      ## Diatoms
      media1 <- sum(wd*datAP_dia)/(700-400)
      datos21 <- datAP_dia*(tabla_total[i,7]/media1)
      # AP_PS
      media2 <- sum(wd*datPS_dia)/(700-400)
      datos22 <- datPS_dia*(tabla_total[i,5]/media2)
        #datos11<-datos11/datos11[12]*datos12[12]
        points(lambda,datos21, type="l", pch=19, lty=1, col="brown", lwd=2)
        media11 <- sum(wd*datos21)/(700-400)
        max11<-datos21[lambda==450]      
        points(lambda,datos22, type="l", pch=19, lty=1, col="orange", lwd=2)
        media12 <- sum(wd*datos22)/(700-400)
        max12<-datos22[lambda==450]
        
        matrix_diatoms_ap[,i]<-datos21
        matrix_diatoms_ps[,i]<-datos22      
        
        legend(x="topleft",legend=c(paste("<aPH> = ",round(media11,3),sep=""),paste("aPH(450) = ",round(max11,3),sep=""),
                                     paste("<aPS> = ",round(media12,3),sep=""),paste("aPS(450) = ",round(max12,3),sep="")),
                                     text.col=c("brown","brown","orange","orange"), bty = "n")       
        
      
      ## Save .dat
      datosothers  <- cbind(paste("",lambda, sep=" "),
                            format(round(datos11,4), nsmall = 4),
                            format(round(datos12,4), nsmall = 4),
                            format(round(datos13,4), nsmall = 4),
                            format(round(datos14,9), nsmall = 9))
      datosdiatoms <- cbind(paste("",lambda, sep=" "),
                            format(round(datos22,4), nsmall = 4),
                            format(round(datos22,4), nsmall = 4),
                            format(round(datos23,4), nsmall = 4),
                            format(round(datos24,9), nsmall = 9))    
      
      # SAVE phyto_optics.dat file
      pdir<-"/Users/ealvarez/Datos/Res_C20_radtrans/"
      sink(paste(pdir,"phyto_optics/sensitivity/optics_phyto_recom_carbon_",i,".dat",sep=""))
      # HEADER
      cat("Computed in 10_Optics2_Marshall_sensitivity_carbon.R\n")
      cat("# aPH compiled from literature (absorption_spectra.csv) (m2 mgChla-1)\n")    
      cat("# aPS reconstructed based on Bidigare/Babin/Hickman (m2 mgChla-1)\n")
      cat("# b converted to mol-specific from Dutkiewicz2015 (m2 mmolC-1)\n")
      cat("# non-spectral bb converted to mol-specific from Dutkiewicz2015 (m2 molC-1)\n")
      cat("Format I4,3F10.4,F20.14\n")
      
      cat("*** Others ***\n")    
      write.table(datosothers, file=paste(pdir,"phyto_optics/sensitivity/optics_phyto_recom_carbon_",i,".dat",sep=""),
                  sep = "    ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
      #sep = "\t"
      cat("*** Diatom ***\n",   file=paste(pdir,"phyto_optics/sensitivity/optics_phyto_recom_carbon_",i,".dat",sep=""), append = TRUE)
      write.table(datosdiatoms, file=paste(pdir,"phyto_optics/sensitivity/optics_phyto_recom_carbon_",i,".dat",sep=""),
                  sep = "    ", row.names = FALSE, col.names = FALSE, quote = FALSE, append=TRUE)
      sink()       
      
      
      }    
    dev.off()
    
    
#####################       
    