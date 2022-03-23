if (Sys.info()['sysname']=="Windows") {OS<-"C:/"} else {OS<-paste("/Users/",Sys.info()['user'],"/", sep="")}
library(chron)

#############################################
## Make aPH variable based on PPC:TChla ratio
#############################################

# Steps:
# Start with best constant curve: optics_phyto_recom_carbon.dat == optics_phyto_recom_carbon_12.dat
# Find the best shape descriptor (using the NOMADv2 dataset)
# Re-construct aPH based on shape descriptor and put it into REcoM2 code (recom_aphyt.F90)


###########
## NOMAD v2 : In situ bio-optical data + HPLC ##  Werdell & Bailey 2005
###########
# Phytoplankton absorption spectra and associated HPLC pigment measurements
# were obtained from the NASA NOMAD-SeaBASS data archive [Werdell and Bailey, 2005]. 
    p_dir <- paste(OS,"Datos/Dat_SGlobal_depth/Valente-BO/NASA_NOMAD/",sep="")
    nomadv2<-read.csv(paste(p_dir,"NOMAD_v2.csv", sep=""))

    # total particulate = phyto + detritus
    #total <- (nomadv2$ap465+nomadv2$ag465)

    # Components
    sp_ap <- nomadv2[!is.na(nomadv2$TPPC),72:91]
    sp_ad <- nomadv2[!is.na(nomadv2$TPPC),92:111]
    sp_ag <- nomadv2[!is.na(nomadv2$TPPC),112:131] 
    sp_aph <- sp_ap - sp_ad
    TChla <- nomadv2$Tchla[!is.na(nomadv2$TPPC)]
    lambda <- as.numeric(substring(colnames(sp_ap),3,5))

### Normalize with aph670
          normalize <- function(x){
                       #media <- weighted.mean(x,w=lambda)
                       media<-x[lambda==670]
                       res <- x/media
                       return(res)}
          sp_aphN <-apply(sp_aph, MARGIN=1, FUN=normalize) 
          #dim(sp_aphN)
          sp_aphN <- t(sp_aphN)
                    
###  Compute fdiatoms: DP and fdChla        
       #names(nomadv2)         
       Fuco <- nomadv2$fuco
       NDP <- (nomadv2$fuco*1.41 +
               nomadv2$perid*1.41 +
               nomadv2$hex.fuco*1.27 +
               nomadv2$allo*0.6 +
               nomadv2$but.fuco*0.35 +
               nomadv2$chl_b*1.01 +
               nomadv2$zea*0.86)

      NfdChla <- (1.41*Fuco)/NDP  
      fdChla <- NfdChla[!is.na(nomadv2$TPPC)] 
                              

### Characterize spectral shape of aph(l)
    library(abind)
        # Chla-specific ap
        specific <- function(x){
            media <- x[21]
            res <- x[1:20]/media
            return(res)}
        
        resu <- abind(sp_aph,TChla,along=2)
        sp_aphC <-apply(resu, MARGIN=1, FUN=specific) 
        sp_aphC <- t(sp_aphC)    
        
#############################
    
        
        
        
################     
### FIGURE 3 ###        
################
        
    path_figures<-paste(OS,"Documentos/5_Trabajos/20_Radtrans_aph/reviews_coauthors2/para_enviar_JAMES/Figuras/",sep="") 
    png(file=paste(path_figures,"Figure3_PPC_aph_2groups.png",sep=""),width = 1210, height = 730, units = "px", pointsize = 22, bg = "white")

################            
        layout(matrix(c(1:4),ncol=2, byrow=F),widths=c(0.6,1))
        par(mar=c(4,5,2,2))
        layout.show(n=4)
        lettersize=2
        lettersizeS=1.2
        
        equis2 <- nomadv2$TPPC[!is.na(nomadv2$TPPC)]/nomadv2$Tchla[!is.na(nomadv2$TPPC)]
        #sum(complete.cases(cbind(equis2,sp_aphC,fdChla)))
        
#### A Small phyto
##################             
        equis4  <- equis2[fdChla<=0.5]
        sp_aphS <- sp_aphC[fdChla<=0.5,] 
        #sum(complete.cases(cbind(equis4,sp_aphS)))

        # Eisner 2003
        ies <- (sp_aphS[,6]-sp_aphS[,9])/((lambda[6]-lambda[9]))
        # Hirata 2008  # aPH443-aPH510/(443-510)
        #ies <- (sp_aphS[,3]-sp_aphS[,7])/((lambda[3]-lambda[7]))
        # RATIOS
        #ies <- (sp_aphS[,6]/sp_aphS[,9])   # aPH489:aPH530
        #ies <- (sp_aphS[,3]/sp_aphS[,7])   # aPH443:aPH510
        #ies <- (sp_aphS[,4]/sp_aphS[,6])   # aPH455:aPH489
        
        plot(equis4, ies, las=1, main="",cex.main=0.8, xlab="",
             ylab="slope (489-532) nm",mgp=c(3.5,1,0), yaxt="n",
             xlim=c(0,0.6), xaxs="i", pch=21, bg="grey30", cex=0.8)   # ylim=c(-0.0008,0.0), 
        axis(2, at=c(-2e-4, -4e-4,-6e-4,-8e-4, -1e-3,-1.2e-3,-1.4e-3),
             labels=c("-2e-4","-4e-4","-6e-4","-8e-4", "-1e-3","-1.2e-3","-1.4e-3"), las=1)
        regre1 <- lm(ies~equis4)
        summ <- summary(regre1)
        abline(reg=regre1, col="black", lwd=2)
        a=round(summ$r.squared,2)
        text(x=0.6, y=-1e-4, labels=paste("y = ",round(regre1$coefficients[2],5),"x",round(regre1$coefficients[1],5),sep=""),pos=2)
        text(x=0.6, y=-2.5e-4, labels=expression(paste(R^2, " = 0.76",sep="")), pos=2)        
        text(x=0.6, y=-4e-4, labels=paste("n = ",sum(complete.cases(cbind(equis4,ies))),sep=""), pos=2)
        
        mtext(3,at=0,adj=0, line=0.4,text="a) fDiatoms<0.5", cex=lettersizeS, font=1)
        mtext(1,at=0.3,line=2.4, text="PPC:TChla", cex=0.8, font=1)
#############        
        
        
#### B Diatoms
##############        
        equis3  <- equis2[fdChla>=0.5]
        sp_aphD <- sp_aphC[fdChla>=0.5,] 
        #sum(complete.cases(cbind(equis3,sp_aphD)))

        # Eisner 2003
        ies <- (sp_aphD[,6]-sp_aphD[,9])/((lambda[6]-lambda[9]))
        # Hirata 2008  # aPH443-aPH510/(443-510)
        #ies <- (sp_aphD[,3]-sp_aphD[,7])/((lambda[3]-lambda[7]))
        # RATIOS
        #ies <- (sp_aphD[,6]/sp_aphD[,9])   # aPH489:aPH530
        #ies <- (sp_aphD[,3]/sp_aphD[,7])   # aPH443:aPH510
        #ies <- (sp_aphD[,4]/sp_aphD[,6])   # aPH455:aPH489        
        
        plot(equis3, ies, las=1, main="",cex.main=0.8, xlab="", ylab="slope (489-532) nm",mgp=c(3.5,1,0),
            xlim=c(0,0.4), xaxs="i", pch=21, bg="grey30", cex=0.8, yaxt="n")  # xlim=c(0,0.5), ylim=c(-0.0008,0.0), 
            axis(2, at=c(-0.5e-4,-1e-4,-1.5e-4,-2e-4,-2.5e-4,-3e-4,-3.5e-4,-4e-4),
                 labels=c("5e-5","-1e-4","-1.5e-4","-2e-4","-2.5e-4","-3e-4","-3.5e-4","-4e-4"), las=1)
            regre2 <- lm(ies~equis3)
            summ <- summary(regre2)
            abline(reg=regre2, col="black", lwd=2)        
            #legend(x="bottomleft", legend=c(paste("y = ",round(regre2$coefficients[1],5),"x",round(regre2$coefficients[2],5),sep=""),
            a=round(summ$r.squared,2)
            substitute(paste("R2 = ",a,sep=""), list(a))      
            
            text(x=0.4, y=-0.00005, labels=paste("y = ",round(regre2$coefficients[2],5),"x",round(regre2$coefficients[1],5),sep=""),pos=2)
            text(x=0.4, y=-0.00008, labels=expression(paste(R^2, " = 0.39",sep="")), pos=2) 
            text(x=0.4, y=-0.00011, labels=paste("n = ",sum(complete.cases(cbind(equis3,ies))),sep=""), pos=2)
            
          mtext(3,at=0,adj=0,line=0.4,text="b) fDiatoms>0.5", cex=lettersizeS, font=1)
          mtext(1,at=0.2,line=2.4, text="PPC:TChla", cex=0.8, font=1)
#################        
     save(regre1,regre2,file=paste(OS,"Datos/Res_C20_radtrans/phyto_optics/regre1.RData",sep=""))
        
        
    
#### Re-construct aph(l)
    dir <-paste(OS,"Datos/Res_C20_radtrans/phyto_optics/",sep="")   
    pdir<-paste(OS,"Datos/Res_C20_radtrans/",sep="")
    paleta<-wes_palette(6, name = "Zissou1", type = "continuous")
    colores<-c(paleta[3],paleta[5])    
    colores<-c("grey60","black")  
    
    i=12  #  optics_phyto_recom_carbon_12.dat", initial constant
        
        
### C Small Phyto
#################        
        alpha <- 0.14/86400
        filename<-paste(pdir,"phyto_optics/sensitivity/optics_phyto_recom_carbon_",i,".dat",sep="")
        datos <- read.table(filename, skip=7, nrows=13)
        lambda <- datos[,1]
          # APH
          dat3 <- datos[,2]
          alpha <- 0.14/86400 
          aph1 <- c(0.008, 0.0143, 0.025)
          QY1 <- alpha / aph1
          plot(lambda, dat3, type="l", ylim=c(0,0.15), pch=19, lty=1, lwd=3, las=1,
               col=colores[2], ylab=expression(paste({a^{"*"}}," (",m^2," mg ",Chl^-1,")",sep='')), xlab="", mgp=c(3,1,0), yaxs="i") 
               #media2 <- sum(25*datos[,34])/(700-400)
          # AP_PS
          dat <- datos[,3]
          points(lambda,dat, type="l", ylim=c(0,0.08), pch=19, lty=1, lwd=3, las=1,
               col=colores[1], ylab="", xlab="", mgp=c(3,1,0), yaxs="i") 
               #media2 <- sum(25*datos[,34])/(700-400)
      
          load(paste(dir,"regre1.RData",sep=""))
          new_PPCs <- c(0.2,0.4,0.6,0.8,1.0)
          for (k in c(1:5)){
                    new_PPC <- new_PPCs[k]
                    slope1 <- new_PPC*coefficients(regre1)[2]+coefficients(regre1)[1]
                    #slope2 <- coefficients(regre2)[1]+new_PPC*coefficients(regre2)[2]
                    dat2 <- dat
                    dat2[lambda==500] <- dat2[lambda==525]+(-slope1)*(525-500)
                    dat2[lambda==475] <- dat2[lambda==525]+(-slope1)*(525-475)
                    dat2[lambda==450] <- dat2[lambda==450]+ (dat2[lambda==475]-dat[lambda==475])
                    dat2[lambda==425] <- dat2[lambda==425]+ (dat2[lambda==475]-dat[lambda==475])
                    dat2[lambda==400] <- dat2[lambda==400]+ (dat2[lambda==475]-dat[lambda==475])
                points(lambda,dat2, type="l", ylim=c(0,0.08), pch=19, lty=2, lwd=2, las=1,col=colores[2], mgp=c(3,1,0), yaxs="i") 
                text(x=450,y=dat2[lambda==450]-0.005,pos=3,labels=paste("PPC:TChla = ",new_PPC, sep=""), cex=0.7)
                    } # end loop k
          mtext(3,at=392,adj=0, line=0.4, text="c) Small phytoplankton", cex=lettersizeS, font=1)
          mtext(1,at=550,line=2.4, text=expression(lambda), cex=0.8, font=1)

          legend(x="topright",bty="n",lty=c(1,1,2),lwd=2, col=c(colores[1],colores[2],colores[2]),
          legend=c(expression(paste({a^{"*"}}[PS],"(", lambda, ")",sep="")),
                      expression(paste({a^{"*"}}[PH],"(", lambda, ") constant",sep="")),
                      expression(paste({a^{"*"}}[PH],"(", lambda, ") variable",sep=""))) )         
          
#################
                        
### D Diatoms
##################      
        alpha <- 0.19/86400
        filename<-paste(pdir,"phyto_optics/sensitivity/optics_phyto_recom_carbon_",i,".dat",sep="")
        datos <- read.table(filename, skip=21)
        lambda <- datos[,1]
        # APH        
        dat3 <- datos[,2]
        alpha <- 0.19/86400 
        aph2 <- c(0.007, 0.0115, 0.020)
               QY2 <- alpha / aph2    
               plot(lambda,dat3, las=1, type="l", ylim=c(0,0.15), pch=19, lty=1, lwd=3,
                    col=colores[2], ylab=expression(paste({a^{"*"}}," (",m^2," mg ",Chl^-1,")",sep='')), xlab="", mgp=c(3,1,0), yaxs="i")
               #media1 <- sum(25*datos[,30])/(700-400)
         # AP_PS
         dat <- datos[,3]
               points(lambda,dat, las=1, type="l", ylim=c(0,0.14), pch=19, lty=1, lwd=3,
                    col=colores[1], ylab="", xlab="", mgp=c(3,1,0), yaxs="i")
               #media1 <- sum(25*datos[,30])/(700-400)           
          
          load(paste(dir,"regre1.RData",sep=""))     
          new_PPCs <- c(0.2,0.4,0.6,0.8,1.0)
          for (k in c(1,3,5)){
          new_PPC <- new_PPCs[k]
                      slope1 <- new_PPC*coefficients(regre2)[2]+coefficients(regre2)[1]
                      dat2 <- dat
                      dat2[lambda==500] <- dat2[lambda==525]+(-slope1)*(525-500)
                      dat2[lambda==475] <- dat2[lambda==525]+(-slope1)*(525-475)
                      dat2[lambda==450] <- dat2[lambda==450] + (dat2[lambda==475]-dat[lambda==475])
                      dat2[lambda==425] <- dat2[lambda==425]+ (dat2[lambda==475]-dat[lambda==475])
                      dat2[lambda==400] <- dat2[lambda==400]+ (dat2[lambda==475]-dat[lambda==475])
                      
                  points(lambda,dat2, type="l", ylim=c(0,0.08), pch=19, lty=2, lwd=2, las=1,col=colores[2], mgp=c(3,1,0), yaxs="i") 
                  text(x=450,y=dat2[lambda==450]-0.005,pos=3,labels=paste("PPC:TChla = ",new_PPC, sep=""), cex=0.7)
                      } # end loop k   
          mtext(3,at=392,adj=0, line=0.4,text="d) Diatoms", cex=lettersizeS, font=1)
          mtext(1,at=550,line=2.4, text=expression(lambda), cex=0.8, font=1)
 
          legend(x="topright",bty="n",lty=c(1,1,2),lwd=2, col=c(colores[1],colores[2],colores[2]),
                 legend=c(expression(paste({a^{"*"}}[PS],"(", lambda, ")",sep="")),
                          expression(paste({a^{"*"}}[PH],"(", lambda, ") constant",sep="")),
                          expression(paste({a^{"*"}}[PH],"(", lambda, ") variable",sep=""))) )         
          
##################          
          
                 
    dev.off()
          