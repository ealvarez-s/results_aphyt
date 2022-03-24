source("00_Summary_EDITME.R")
source("functions_misc.R")

###############################
## BIO-OPTICAL RELATIONSHIPS
###############################

  o_dir<-paste(global_path,"Res_model/",sep="")        
  p_dir<-paste(global_path,"Res_model/interpolated/",sep="")
  path_figures<-paste(global_path,"Figures/",sep="")
  # List of simulations
  experimentos<-read.csv(paste(o_dir,"run_log_marshall_PPC_PS_SA_G100.csv",sep=""), sep=",")
  s_dir <- paste(global_path,"Dat_observations/satellite/",sep="")
  is_dir <- paste(global_path,"Dat_observations/optics/grids_optics",sep="")

          modelnames<-experimentos$name[c(21,9)]
          ruta<-experimentos$path[c(21,9)]
          nombre <- c(expression(paste({a^{"*"}}[PH],"(", lambda, ") constant",sep="")),
                      expression(paste({a^{"*"}}[PH],"(", lambda, ") variable",sep="")))
          nombre <- c(expression(paste("EXP-1: ",{a^{"*"}}[PH],"(", lambda, ") constant",sep="")),
                      expression(paste("EXP-2: ",{a^{"*"}}[PH],"(", lambda, ") variable",sep="")))
          


############# 
### FIGURE 9:  Chla vs Aph(450), with obs from Valente+Bracher, and Bricaud (1995, 2004) relships
############

png(file=paste(path_figures,"Figure9_chla_vs_Aph.png", sep=""), width = 1320, height = 450, units = "px",pointsize = 19, bg = "white")
            layout(matrix(c(1,2,3),ncol=3,byrow=T))
            par(family="")
            par(mar=c(4,4.5,2,2))  
            par(oma=c(1,1,1,0))                 
            layout.show(n=3)
            lettersize=1.8
            textinner=1.4
            cexlab=1.3
            cexaxis=1.2
            cextitle=1.8
            
# OBSERVATIONS
##############     
            modelo_lambda <- c(400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0,600.0, 625.0, 650.0, 675.0, 700.0)
            lambda_extremos<-sort(unique(c(modelo_lambda[-length(modelo_lambda)]+(diff(modelo_lambda)/2),400,700)))
            s_dir <- paste(global_path,"Dat_observations/optics/total_tables/",sep="")
            
        #### Valente  -(awi&nomad) 
            load(paste(s_dir,"Valente2_ChlaRssIOP_modbands12.RData", sep=""))       
            #names(tabla_total)
            equis<-log10(tabla_total$Chla_HPLC.mg.m..3.)
            ies<-log10(tabla_total$aph_450)  
            flag1<-tabla_total[,which(colnames(tabla_total)=="QFChla")] 
            equis[flag1==1]<-NA              
            
            plot(equis,ies,pch=19,col="grey30", las=1,xlim=c(-2.5,2), ylim=c(-3.5,0.2),xaxt="n", yaxt="n", xaxs="i",yaxs="i",
                 xlab=expression(paste("TChla (mg ",m^-3,")",sep='')),ylab=expression(paste(a[PH]," (450nm) [",m^-1,"]",sep='')),
                 cex.lab=cexlab, cex.main=cextitle, main=expression(italic("In situ"))) 
                      axis(1,at=seq(-2,2,by=1),labels=10^seq(-2,2,by=1), las=1, cex.axis=cexaxis)
                      axis(2,at=seq(-3,0,by=1),labels=10^seq(-3,0,by=1), las=1, cex.axis=cexaxis)
                      
            mtext(3,at=-2.5, adj=0, line=0.5, text="a)", cex=lettersize, font=1)
                      tabla<-cbind(equis,ies)
                      tabla[is.nan(tabla)]<-NA
                      tabla[tabla=="Inf"]<-NA
                      tabla[tabla=="-Inf"]<-NA
                      tabla<-tabla[complete.cases(tabla),]
               # My FIT
               newx <- seq(2e-2, 20, length=50)
               logcloro <- tabla[,1]
               logX     <- tabla[,2]
               fit5 <- lmodel2(logX ~ logcloro, data=data.frame(tabla)) 
               newc <- (10^fit5$regression.results[1,2])*(newx^fit5$regression.results[1,3])
               points(log10(newx), log10(newc), type="l", col="black", lwd=2)
            
          #### Bracher         
            load(paste(s_dir,"Astrid6_IopHplc2nm_modbands12.RData", sep=""))       
            equis<-log10(tabla_total$MTChla)
            ies<-log10(tabla_total$aph_450)
            flag3<-tabla_total[,which(colnames(tabla_total)=="MYFLAG3")]
            flag5<-tabla_total[,which(colnames(tabla_total)=="MYFLAG5")]
            equis[flag3==3]<-NA 

            points(equis,ies, pch=19,col="grey70", las=1, xlim=c(-3,1), ylim=c(-3,0)) 
                tabla<-cbind(equis,ies)
                tabla[is.nan(tabla)]<-NA
                tabla[tabla=="Inf"]<-NA
                tabla[tabla=="-Inf"]<-NA
                tabla<-tabla[complete.cases(tabla),]
            
      points(log10(newx), log10(newc), type="l", col="black", lwd=2) # this line just overplots Valente fit again, so it is visible
              
              # My FIT
              logcloro <- tabla[,1]
              logX     <- tabla[,2]
              fit2 <- lmodel2(logX ~ logcloro, data=data.frame(tabla)) 
              newc <- (10^fit2$regression.results[1,2])*(newx^fit2$regression.results[1,3])
              points(log10(newx), log10(newc), type="l", col="black", lwd=2, lty=2)
            
            
      ## Bricaud 1995
      newc <- 0.0383*(newx^0.651)
      #points(log10(newx), log10(newc), type="l", col="grey50", lty=2, lwd=2)
      a=round(0.0383,4)
      b=round(0.651,3)
      legend(x=-1.6, y=-2.9,cex=1.1,#lty=c(2), lwd=2, col="grey50",
             legend=substitute(paste("Bricaud et al. 1995: ",a ,"+" , TChla^b,sep=""),list(a=a, b=b)), bty="n")  

      ## Bricaud 2004
      newc <- 0.0654*(newx^0.728)
      #points(log10(newx), log10(newc), type="l", col="grey50", lty=1, lwd=2)
      a=round(0.0654,4)
      b=round(0.728,3)
      legend(x=-1.6, y=-3.2,cex=1.1,#lty=c(1), lwd=2, col="grey50",
             legend=substitute(paste("Bricaud et al. 2004: ",a ,"+" , TChla^b,sep=""),list(a=a, b=b)), bty="n")  
      
    # Valente fit
     a=round(10^fit5$regression.results[1,2],3)
     b=round(fit5$regression.results[1,3] ,3)
     legend(x=-2.5,y=0.25, lty=c(1), lwd=2, col="black",cex=1.15,
            legend=substitute(paste("Valente: ",a ,"+" , TChla^b,sep=""),list(a=a, b=b)), bty="n") 

     # Bracher fit
     a=round(10^fit2$regression.results[1,2],3)
     b=round(fit2$regression.results[1,3] ,3) 
     legend(x=-2.5, y=-0.05, lty=c(2), lwd=2, col="black",cex=1.15,
            legend=substitute(paste("Bracher: ",a ,"+" , TChla^b,sep=""),list(a=a, b=b)), bty="n")  

##############            
            
            
# MODEL
############       
            modelnames<-experimentos$name[c(21,9)]
            ruta<-experimentos$path[c(21,9)]
                mycol1<-usecol(pal_bordeaux,n=4)
                mycol2<-usecol(pal_petrol,n=4)
                colores_puntos<-c(mycol1[3],mycol2[3])
                colores_bins<-c(mycol1[1],mycol2[1])            
             
                    for (w in 1:length(modelnames)){  
                        modelname<-modelnames[w]
                        # CHLA
                        load(paste(p_dir,"2cloro_lin/",modelname,".RData",sep=""))
                        cloro <- res
                        cloro[cloro<1e-4] <- 1e-4
                        # Aph
                        load(paste(p_dir,"2iop_aph/",modelname,".RData",sep=""))
                        res[cloro<=1e-4] <- NA
                        X <- res[,,which(dimnames(res)$l=="450"),]

                 plot(log10(cloro), log10(X), type="n", las=1, pch=19, cex=0.08, ylim=c(-3.5,0.2), xlim=c(-2.5, 2), xaxt="n", yaxt="n", xaxs="i",yaxs="i", main=nombre[w],
                      xlab=expression(paste("TChla (mg ",m^-3,")",sep='')),ylab=expression(paste(a[PH]," (450nm) [",m^-1,"]",sep='')), cex.lab=cexlab, cex.main=cextitle)
                   if (w==1){mtext(3,at=-2.5,adj=0,line=0.5, text="b)", cex=lettersize, font=1)}       
                   if (w==2){mtext(3,at=-2.5,adj=0,line=0.5, text="c)", cex=lettersize, font=1)}  
                   
                 points(log10(cloro), log10(X), las=1, pch=19, cex=0.08, ylim=c(-3.5,0.2), xlim=c(-2.5, 2), xaxt="n", yaxt="n", col=colores_puntos[w])
                        axis(1,at=seq(-2,2,by=1),labels=10^seq(-2,2,by=1), las=1, cex.axis=cexaxis)
                        axis(2,at=seq(-3,0,by=1),labels=10^seq(-3,0,by=1), las=1, cex.axis=cexaxis)
                        
                 # Bin aph(450) in TChla bins
                 clases_chla     <- seq(-2,  3, length=101)[seq(1,101,by=2)]
                 clases_chla_mid <- seq(-2,  3, length=101)[seq(2,100,by=2)]
                 cloro  <-log10(cloro)
                 abs  <-log10(X)
                 particulas<-cut(cloro, breaks=clases_chla, include.lowest = TRUE,right = TRUE)
                 grupos<-split(abs,particulas)
                 nombres<-names(grupos)
                 niveles<-levels(particulas)
                 na<-sapply(grupos,mean,na.rm=T)
                 points(clases_chla_mid,na,pch=19, col=colores_bins[w],cex=1.2)           
  
                 # My FIT
                 logcloro <- clases_chla_mid
                 logX <- na
                 newx <- seq(1e-2, 16, length=50)
                 fit4 <- lmodel2(logX ~ logcloro, data=data.frame(cbind(logX,logcloro))) 
                 newc <- (10^fit4$regression.results[1,2])*(newx^fit4$regression.results[1,3])
              points(log10(newx), log10(newc), type="l", col="black", lwd=2) 
                 a=round(10^fit4$regression.results[1,2],3)
                 b=round(fit4$regression.results[1,3],3)                
              text(x=-2,y=-0.5, pos=4,cex=textinner, labels=substitute(paste(" ",a ,"+" , TChla^b,sep=""),list(a=a, b=b)))

              # Show TChla bins in the plot
              axis(1,at=clases_chla[1:33],labels=F,tcl=0.8)
         }  # end loop modelnames
############            
            
dev.off()
#######################



                
                
                
            
            
############# 
### FIGURE 10   TChla vs A*ph(450), with obs from Bracher
############
               
     png(file=paste(path_figures,"Figure10_chl_specific_aph.png", sep=""), width = 1220, height = 850, units = "px",pointsize = 18.5, bg = "white")
            layout(matrix(c(1:6),ncol=3,byrow=T))
            par(family="")
            par(mar=c(4,4.7,2,2))  
            par(oma=c(10,0.5,1,0))                 
            layout.show(n=6)
            lettersize=1.8
            textinner=1.4
            cexlab=1.3
            cexaxis=1.2            
            cextitle=1.7
            
# OBSERVATIONS
############## 
            modelo_lambda <- c(400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0,600.0, 625.0, 650.0, 675.0, 700.0)
            lambda_extremos<-sort(unique(c(modelo_lambda[-length(modelo_lambda)]+(diff(modelo_lambda)/2),400,700)))
            s_dir <- paste(global_path,"Dat_observations/optics/total_tables/",sep="")
            
            ### Valente
            load(paste(s_dir,"Valente2_ChlaRssIOP_modbands12.RData", sep=""))       
            equis<-log10(tabla_total$Chla_HPLC.mg.m..3.)
            ies<-tabla_total$aph_450/tabla_total$Chla_HPLC.mg.m..3.
            ies[ies==0]<-NA
            
            plot(equis,ies,pch=19,col="grey30", las=1,xlim=c(-2.5,2), ylim=c(0,0.3),xaxt="n", yaxt="n", xaxs="i",yaxs="i",
              xlab=expression(paste("TChla (mg ",m^-3,")",sep='')),
              ylab=expression(paste({a^{"*"}}[PH]," (450nm) [",m^2," mg ",Chl^-1,"]",sep='')),
              cex.lab=cexlab, cex.main=cextitle, main=expression(italic("In situ")),)
              axis(1,at=seq(-3,3,by=1),labels=10^seq(-3,3,by=1), las=1, cex.axis=cexaxis)
              axis(2,at=seq(0,0.3,by=0.05), las=1, cex.axis=cexaxis)
              mtext(3,at=-2.5,adj=0, line=0.5, text="a)", cex=lettersize, font=1)
            
              tabla<-cbind(equis,ies)
              tabla[is.nan(tabla)]<-NA
              tabla[tabla=="Inf"]<-NA
              tabla[tabla=="-Inf"]<-NA
              tabla<-tabla[complete.cases(tabla),]
              
            # My FIT
            newx <- seq(-2, 2, length=50)
            logcloro <- tabla[,1]
            logX     <- tabla[,2]
            fit4 <- lmodel2(logX ~ logcloro, data=data.frame(tabla))
            newc <- fit4$regression.results[1,2]+(newx*fit4$regression.results[1,3])

            ### Bracher
            load(paste(s_dir,"Astrid6_IopHplc2nm_modbands12.RData", sep=""))       
            equis<-log10(tabla_total$MTChla)
            ies<-tabla_total$aph_450/tabla_total$MTChla
            equis[equis<=-2]<-NA
            ies[ies==0]<-NA
            points(equis,ies, pch=19,col="grey70", las=1) 
              tabla<-cbind(equis,ies)
              tabla[is.nan(tabla)]<-NA
              tabla[tabla=="Inf"]<-NA
              tabla[tabla=="-Inf"]<-NA
              tabla<-tabla[complete.cases(tabla),]
              
            # overplot Valente fit again, so it is visible
            points(newx, newc, type="l", col="black", lwd=2)
            
            # My FIT
            logcloro <- tabla[,1]
            logX     <- tabla[,2]
            fit3 <- lmodel2(logX ~ logcloro, data=data.frame(tabla))
            newc <- fit3$regression.results[1,2]+(newx*fit3$regression.results[1,3])
            points(newx, newc, type="l", col="black", lwd=2, lty=2)
            
            # Valente fit
            a=round(fit4$regression.results[1,2],3)
            b=round(fit4$regression.results[1,3],3)
            legend(x=-2.3,y=0.30, lty=c(1), lwd=2, col="black",cex=1.15,
                   legend=substitute(paste("Valente ",b ," " ,log[10],"(TChla)", " + ", a,sep=""),list(a=a, b=b)), bty="n") 
            # Bracher fit            
            a=round(fit3$regression.results[1,2],3)
            b=round(fit3$regression.results[1,3],3)
            legend(x=-2.3,y=0.27, lty=c(2), lwd=2, col="black",cex=1.15,
                   legend=substitute(paste("Bracher ",b ," " ,log[10],"(TChla)", " + ", a,sep=""),list(a=a, b=b)), bty="n")              
            
##############            
            
            
# MODEL
############       
            mycol1<-usecol(pal_bordeaux,n=4)
            mycol2<-usecol(pal_petrol,n=4)
            colores_puntos<-c(mycol1[3],mycol2[3])
            colores_bins<-c(mycol1[1],mycol2[1])
            
            for (w in 1:length(modelnames)){  
              modelname<-modelnames[w]
              
              # CHLA
              load(paste(p_dir,"2cloro_lin/",modelname,".RData",sep=""))
              cloro <- res
              cloro[cloro<1e-4] <- 1e-4
              # Aph
              load(paste(p_dir,"2iop_aph/",modelname,".RData",sep=""))
              res[cloro<=1e-4] <- NA
              X <- res[,,which(dimnames(res)$l=="450"),]
              equis<-log10(cloro)
              ies<-X/cloro

              plot(equis, ies, type="n", las=1, pch=19, cex=0.08, ylim=c(0,0.3), xlim=c(-2.5, 2),
                   xaxt="n", yaxt="n", xaxs="i",yaxs="i", main=nombre[w],  cex.main=cextitle,
                   xlab=expression(paste("TChla (mg ",m^-3,")",sep='')),
                   ylab=expression(paste({a^{"*"}}[PH]," (450nm) [",m^2," mg ",Chl^-1,"]",sep='')), cex.lab=cexlab)
              if (w==1){mtext(3,at=-2.5,adj=0, line=0.5, text="b)", cex=lettersize, font=1)}       
              if (w==2){mtext(3,at=-2.5,adj=0, line=0.5, text="c)", cex=lettersize, font=1)}  
              
              points(equis, ies, las=1, pch=19, cex=0.08, ylim=c(-3.5,0.2), xlim=c(-2.5,2),xaxt="n", yaxt="n", main=nombre[w], col=colores_puntos[w])
              axis(1,at=seq(-3,3,by=1),labels=10^seq(-3,3,by=1), las=1, cex.axis=cexaxis)
              axis(2,at=seq(0,0.3,by=0.05), las=1, cex.axis=cexaxis)
              
              # Bin a*ph(450) in TChla classes
              clases_chla     <- seq(-2, 3, length=101)[seq(1,101,by=2)]
              clases_chla_mid <- seq(-2, 3, length=101)[seq(2,100,by=2)]              
              cloro <- equis
              abs <- ies
              particulas<-cut(cloro, breaks=clases_chla, include.lowest = TRUE,right = TRUE)
              grupos<-split(abs,particulas)
              nombres<-names(grupos)
              niveles<-levels(particulas)
              na<-sapply(grupos,mean,na.rm=T)
              points(clases_chla_mid,na,pch=19,col=colores_bins[w],cex=1.1)           
              
              # My FIT
              log10(1e-2)
              log10(16)
              newx <- seq(-2, 1.2, length=50)
              logcloro <- clases_chla_mid
              logX <- na
              fit4 <- lmodel2(logX ~ logcloro, data=data.frame(cbind(logX,logcloro)))
              newc <- fit4$regression.results[1,2]+(newx*fit4$regression.results[1,3])
              points(newx, newc, type="l", col="black", lwd=2)
              # Bracher fit
              a=round(fit4$regression.results[1,2],3)
              b=round(fit4$regression.results[1,3],3)
              legend(x="topright",cex=1.15,
                     legend=substitute(paste("",b ," " ,log[10],"(TChla)", " + ", a,sep=""),list(a=a, b=b)), bty="n")              
        
              # Show TChla bins in the plot
              axis(1,at=clases_chla[1:33],labels=F,tcl=0.8)                
        
      }  # end loop modelnames
############  
            
            limsurface=0      
            
# OBSERVATIONS in depth
#################       
      o_dir <-paste(global_path,"Res_model", sep="")
      archivo<-"/grid.nc"
      filename <- paste(o_dir,archivo, sep="")
      filen <- open.nc(filename)
      filedepth <- read.nc(filen)       
      profun<-abs(filedepth$Z)
      prof_model <- profun 
      profunU<-abs(filedepth$Zu)        
      profunL<-abs(filedepth$Zl)    
      clases_pro  <-c(-1,profunU)
            
            modelo_lambda <- c(400.0, 425.0, 450.0, 475.0, 500.0, 525.0, 550.0, 575.0,600.0, 625.0, 650.0, 675.0, 700.0)
            lambda_extremos<-sort(unique(c(modelo_lambda[-length(modelo_lambda)]+(diff(modelo_lambda)/2),400,700)))
            
            ### Bracher
            s_dir <- paste(global_path,"Dat_observations/optics/total_tables/",sep="")
            load(paste(s_dir,"Astrid6_IopHplc2nm_modbands12.RData", sep=""))       
            #names(tabla_total)
            equis  <-tabla_total$latitude
            ies    <-tabla_total$depth
            zetas  <-tabla_total$aph_450/tabla_total$MTChla
            zetas2  <- tabla_total$MTChla
            zetas3  <-tabla_total$MPPC/tabla_total$MTChla
            cruise <-tabla_total$lista
                zetas[zetas<0]<-NA
                zetas2[zetas2<0]<-NA
                zetas3[zetas3<0]<-NA
            niveles<-cut(ies,breaks=clases_pro)
            profundidad<-match(niveles,levels(niveles))
            experiment<-cbind(zetas, niveles, zetas2, zetas3)
            experiment[experiment=="Inf"]<-NA
            experiment[experiment==0]<-NA
            experiment<-experiment[complete.cases(experiment),]
            
                # a*PH
                cuantos <- aggregate(zetas~niveles,experiment, length, drop=F, simplify=FALSE)
                medias <- aggregate(zetas~niveles,experiment, mean, na.rm=TRUE, drop=F, simplify=FALSE)
                medias[unlist(cuantos[,2])<5,2]<-NA
                #cajas <- aggregate(zetas~niveles,experiment, boxplot, plot=F, drop=F, simplify=FALSE) 
                Q3 <- aggregate(zetas~niveles,experiment, quantile, probs=0.75, na.rm=TRUE, drop=F, simplify=FALSE)
                Q3[unlist(cuantos[,2])<10,2]<-NA
                Q2 <- aggregate(zetas~niveles,experiment, quantile, probs=0.25, na.rm=TRUE, drop=F, simplify=FALSE)
                Q2[unlist(cuantos[,2])<10,2]<-NA            
                Q4 <- aggregate(zetas~niveles,experiment, quantile, probs=0.95, na.rm=TRUE, drop=F, simplify=FALSE)
                Q4[unlist(cuantos[,2])<10,2]<-NA
                Q1 <- aggregate(zetas~niveles,experiment, quantile, probs=0.05, na.rm=TRUE, drop=F, simplify=FALSE)
                Q1[unlist(cuantos[,2])<10,2]<-NA 
            
                # TChla   
                media_cloro <- aggregate(zetas2~niveles,experiment, mean, na.rm=TRUE, drop=F, simplify=FALSE)
                media_cloro[unlist(cuantos[,2])<5,2]<-NA            
                Q3_cloro <- aggregate(zetas2~niveles,experiment, quantile, probs=0.75, na.rm=TRUE, drop=F, simplify=FALSE)
                Q3_cloro[unlist(cuantos[,2])<5,2]<-NA
                Q2_cloro <- aggregate(zetas2~niveles,experiment, quantile, probs=0.25, na.rm=TRUE, drop=F, simplify=FALSE)
                Q2_cloro[unlist(cuantos[,2])<5,2]<-NA            
                sd_cloro <- aggregate(zetas2~niveles,experiment, sd, na.rm=TRUE, drop=F, simplify=FALSE)
                Q3_cloro<-media_cloro+sd_cloro
                Q2_cloro<-media_cloro-sd_cloro
                
                # PPC:TChla
                media_ppc <- aggregate(zetas3~niveles,experiment, mean, na.rm=TRUE, drop=F, simplify=FALSE)
                media_ppc[unlist(cuantos[,2])<5,2]<-NA 
                Q3_ppc <- aggregate(zetas3~niveles,experiment, quantile, probs=0.75, na.rm=TRUE, drop=F, simplify=FALSE)
                Q3_ppc[unlist(cuantos[,2])<5,2]<-NA
                Q2_ppc <- aggregate(zetas3~niveles,experiment, quantile, probs=0.25, na.rm=TRUE, drop=F, simplify=FALSE)
                Q2_ppc[unlist(cuantos[,2])<5,2]<-NA
                sd_ppc <- aggregate(zetas3~niveles,experiment, sd, na.rm=TRUE, drop=F, simplify=FALSE)
                Q3_ppc<-media_ppc+sd_ppc
                Q2_ppc<-media_ppc-sd_ppc                
                
       #### plot TChla 
            equis<-media_cloro[,2]
            ies<-prof_model[media_cloro[,1]]
            plot(equis[!is.na(equis)],ies[!is.na(equis)], col="darkseagreen4", type="l",
                 xlim=c(0,1.2), ylim=c(200,limsurface),
                 xaxt="n",yaxt="n", xlab="",ylab="",xaxs="i",yaxs="i", lwd=2,
                 bty="l")
            points(equis[!is.na(equis)],ies[!is.na(equis)], pch=19, col="darkseagreen4", type="p")
            axis(1,line=3.5,at=seq(0,1.2,by=0.2), col="darkseagreen4", col.axis="darkseagreen4", cex.axis=cexaxis)
            mtext(1,line=5.5, at=0.6,text=expression(paste("TChla (mg ",m^-3,")",sep='')), cex=0.8, col="darkseagreen4")            
            # Range
            mycol<-usecol(pal_seegruen,n=4,alpha=0.2)
            ies<-c(prof_model[Q3_cloro[,1]],
                   prof_model[Q3_cloro[,1]][length(prof_model[Q3_cloro[,1]]):1])
            equis<-c(unlist(Q3_cloro[,2]),unlist(Q2_cloro[,2])[length(unlist(Q2_cloro[,2])):1])
            polygon(equis[!is.na(ies)], ies[!is.na(ies)], col=mycol[4], border=NA)
            
       #### add PPC:TChla
            par(new=T)
            equis<-media_ppc[,2]
            ies<-prof_model[media_ppc[,1]]
            plot(equis[!is.na(equis)],ies[!is.na(equis)], pch=21, col="indianred4", type="l",
                 xlim=c(0,0.5), ylim=c(200,limsurface),
                 xaxt="n",yaxt="n", xlab="", ylab="", xaxs="i",yaxs="i", lwd=2, bty="n")
            points(equis[!is.na(equis)],ies[!is.na(equis)], pch=19, col="indianred4", type="p")
            axis(1,line=7,at=seq(0,0.5,by=0.1), col="indianred4", col.axis="indianred4", cex.axis=cexaxis)
            mtext(1,line=9, at=0.25,text="PPC:TChla (g:g)", cex=0.8,col="indianred4")              
            # Range
            mycol<-usecol(pal_bordeaux,n=4,alpha=0.2)
            ies<-c(prof_model[Q3_ppc[,1]],
                   prof_model[Q3_ppc[,1]][length(prof_model[Q3_ppc[,1]]):1])
            equis<-c(unlist(Q3_ppc[,2]),unlist(Q2_ppc[,2])[length(unlist(Q2_ppc[,2])):1])
            polygon(equis[!is.na(ies)], ies[!is.na(ies)], col=mycol[4], border=NA)
            
        #### add a*ph
            par(new=T)      
            plot(1, 1, type="n",xlim=c(0,0.30), ylim=c(200,limsurface), las=1, xlab="",ylab="depth (m)",
               xaxs="i",yaxs="i",
               cex.lab=cexlab, cex.axis=cexaxis, bty="n")
                mtext(3,at=0,adj=0, line=0.5, text="d)", cex=lettersize, font=1)
                mtext(1,line=2.2, at=0.150, text=expression(paste({a^{"*"}}[PH]," (450nm) [",m^2," mg ",Chl^-1,"]",sep='')), cex=0.8)
                
                for (i in c(1:10)){
                segments(x0=unlist(Q1[i,2]),x1=unlist(Q4[i,2]),y0=prof_model[i], lwd=2, lend="square")
                segments(x0=unlist(Q2[i,2]),x1=unlist(Q3[i,2]),y0=prof_model[i], lwd=10,col="grey60", lend="square")
                points(medias[i,2],prof_model[i], pch=22, col="white", bg="white")
                #boxplot(zetas[profundidad==i], add=T, at=prof_model[i], horizontal=T) 
                }
     # n           
     text(x=unlist(Q4[,2])[1:8],y=prof_model[1:8], pos=4,labels=paste("n = ",unlist(cuantos[,2])[1:8],sep=""))
            
#################              
  
                      
## MODEL in depth
#################      
                  
            for (w in 1:length(modelnames)){  
              modelname<-modelnames[w]
              
              # CHLA
              load(paste(p_dir,"1cloro_global/",modelname,".RData",sep=""))
              cloro <- res
              cloro[cloro<1e-4] <- 1e-4
              cloro[cloro==0]<-NA
              
              # CHLA:C
              load(paste(p_dir,"1ratio_global/",modelname,".RData",sep=""))
              ratio <- res
              ratio[cloro[,,1:10]<1e-4] <- 1e-4
              ratio[ratio==0]<-NA
              
              # PPC:CHLA
              load(paste(p_dir,"1ppc_global/",modelname,".RData",sep=""))
              ppc <- res
              ppc[cloro[,,1:10]<1e-4] <- 1e-4
              range(ppc,na.rm=T)
              ppc[ppc==0]<-NA
              
              # Aph
              load(paste(p_dir,"2iop_abphy_wb/",modelname,"_450.RData",sep=""))
              res<-apply(res, MARGIN=c("x","y","z"), FUN=mean, na.rm=T)
              res[cloro<=1e-4] <- NA
              X <- res/cloro

        ### Empty formated plot      
              plot(1, 1, type="n",xlim=c(0,0.30), ylim=c(200,limsurface), las=1, xlab="",ylab="depth (m)",xaxs="i",yaxs="i", cex.lab=cexlab, cex.axis=cexaxis, bty="l")
              clorofila<-rep(NA,length=30)
                  Q3_clorofila<-rep(NA,length=30)
                  Q2_clorofila<-rep(NA,length=30)
              ratioC<-rep(NA,length=30)
                  Q3_ratioC<-rep(NA,length=30)
                  Q2_ratioC<-rep(NA,length=30)
              ratioP<-rep(NA,length=30)
                  Q3_ratioP<-rep(NA,length=30)
                  Q2_ratioP<-rep(NA,length=30)

              if (w==1){mtext(3,at=0,adj=0, line=0.5, text="e)", cex=lettersize, font=1)}       
              if (w==2){mtext(3,at=0,adj=0, line=0.5, text="f)", cex=lettersize, font=1)}  
              
              
              for (i in c(1:10)){
                Q4<-quantile(X[,,i],probs=0.95,na.rm=T)
                Q3<-quantile(X[,,i],probs=0.75,na.rm=T)
                Q2<-quantile(X[,,i],probs=0.25,na.rm=T)
                Q1<-quantile(X[,,i],probs=0.05,na.rm=T)
                media<-mean(X[,,i],na.rm=T)
                segments(x0=Q1,x1=Q4,y0=prof_model[i], lwd=2, lend="square")
                segments(x0=Q2,x1=Q3,y0=prof_model[i], lwd=10, col="grey70", lend="square")
                points(media,prof_model[i], pch=22, col="white", bg="white")
                #boxplot(zetas[profundidad==i], add=T, at=prof_model[i], horizontal=T) 
                clorofila[i]<-mean(cloro[,,i],na.rm=T)
                       Q3_clorofila[i]<-mean(cloro[,,i],na.rm=T)+sd(cloro[,,i],na.rm=T)#quantile(cloro[,,i],probs=0.75,na.rm=T)
                       Q2_clorofila[i]<-mean(cloro[,,i],na.rm=T)-sd(cloro[,,i],na.rm=T)#quantile(cloro[,,i],probs=0.25,na.rm=T)
                ratioC[i]<-mean(ratio[,,i],na.rm=T)
                       Q3_ratioC[i]<-mean(ratio[,,i],na.rm=T)+sd(ratio[,,i],na.rm=T)#quantile(ratio[,,i],probs=0.75,na.rm=T)
                       Q2_ratioC[i]<-mean(ratio[,,i],na.rm=T)-sd(ratio[,,i],na.rm=T)#quantile(ratio[,,i],probs=0.25,na.rm=T)
                ratioP[i]<-mean(ppc[,,i],na.rm=T)
                       Q3_ratioP[i]<-mean(ppc[,,i],na.rm=T)+sd(ppc[,,i],na.rm=T)#quantile(ppc[,,i],probs=0.75,na.rm=T)
                       Q2_ratioP[i]<-mean(ppc[,,i],na.rm=T)-sd(ppc[,,i],na.rm=T)#quantile(ppc[,,i],probs=0.25,na.rm=T)                
              }
              mtext(1,line=2.2, at=0.150,text=expression(paste({a^{"*"}}[PH]," (450nm) [",m^2," mg ",Chl^-1,"]",sep='')), cex=0.8)
              
              
        #### Add TChla
              par(new=T)
              plot(clorofila,prof_model, pch=21,col="darkseagreen4", type="l",
                   xlim=c(0,1.2), ylim=c(200,limsurface),
                   xaxt="n",yaxt="n", xlab="",ylab="", xaxs="i",yaxs="i", lwd=2, bty="n")
                   points(clorofila,prof_model, pch=19,col="darkseagreen4", type="p")
              axis(1,line=3.5,at=seq(0,1.2,by=0.2), col="darkseagreen4", col.axis="darkseagreen4", cex.axis=cexaxis)
              mtext(1,line=5.5, at=0.6,text=expression(paste("TChla (mg ",m^-3,")",sep='')), cex=0.8, col="darkseagreen4")
              # Range
              mycol<-usecol(pal_seegruen,n=4,alpha=0.2)
              ies<-c(prof_model, prof_model[length(prof_model):1])
              equis<-c(Q2_clorofila,Q3_clorofila[length(Q3_clorofila):1])
              #cbind(equis, ies)
              polygon(equis[!is.na(equis)], ies[!is.na(equis)], col=mycol[4], border=NA)              
              
         #### Add Chla:C
              par(new=T)
              plot(ratioC,prof_model, pch=21, col="dodgerblue4", type="l",
                   xlim=c(0,0.08), ylim=c(200,limsurface),
                   xaxt="n",yaxt="n", xlab="",ylab="", xaxs="i",yaxs="i",lwd=2, bty="n") 
                   points(ratioC,prof_model, pch=19,col="dodgerblue4", type="p")
              axis(1,line=7,at=seq(0,0.08,by=0.02), col="dodgerblue4", col.axis="dodgerblue4", cex.axis=cexaxis)
              mtext(1,line=9, at=0.04,text="TChla:C (g:g)", cex=0.8, col="dodgerblue4")              
              # Range
              mycol<-usecol(pal_karpfenblau,n=4,alpha=0.2)
              ies<-c(prof_model, prof_model[length(prof_model):1])
              equis<-c(Q2_ratioC,Q3_ratioC[length(Q3_ratioC):1])
              polygon(equis[!is.na(equis)], ies[!is.na(equis)], col=mycol[4], border=NA)            
              

         #### Add PPC:TChla
              par(new=T)
              plot(ratioP,prof_model, pch=21, col="indianred4", type="l",
                   xlim=c(0,0.5), ylim=c(200,limsurface),
                   xaxt="n",yaxt="n", xlab="",ylab="", xaxs="i",yaxs="i", lwd=2, bty="n")
                   points(ratioP,prof_model, pch=19, col="indianred4", type="p")
              axis(1,line=10.5,at=seq(0,0.5,by=0.1), col="indianred4", col.axis="indianred4", cex.axis=cexaxis)
              mtext(1,line=12.5, at=0.25,text="PPC:TChla (g:g)", cex=0.8, col="indianred4")               
              # Range
              mycol<-usecol(pal_bordeaux,n=4,alpha=0.2)
              ies<-c(prof_model, prof_model[length(prof_model):1])
              equis<-c(Q2_ratioP,Q3_ratioP[length(Q3_ratioP):1])
              polygon(equis[!is.na(equis)], ies[!is.na(equis)], col=mycol[4], border=NA)            
              
        ### Replot aph on top
              par(new=T)
              plot(1, 1, type="n",xlim=c(0,0.30), ylim=c(200,limsurface),
                   las=1, xlab="",ylab="",xaxs="i",yaxs="i",xaxt="n",yaxt="n",
                   cex.lab=cexlab, cex.axis=cexaxis, bty="n")
              for (i in c(1:10)){
                Q4<-quantile(X[,,i],probs=0.95,na.rm=T)
                Q3<-quantile(X[,,i],probs=0.75,na.rm=T)
                Q2<-quantile(X[,,i],probs=0.25,na.rm=T)
                Q1<-quantile(X[,,i],probs=0.05,na.rm=T)
                media<-mean(X[,,i],na.rm=T)
                segments(x0=Q1,x1=Q4,y0=prof_model[i], lwd=2, lend="square")
                segments(x0=Q2,x1=Q3,y0=prof_model[i], lwd=10, col="grey60", lend="square")
                points(media,prof_model[i], pch=22, col="white", bg="white")
             }              
 
            } # end loop modelnames
#################             

            dev.off()
            
#######################            
            