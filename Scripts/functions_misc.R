
flip<-function(vector)
        {
            vector<-vector[length(vector):1];
            return(vector);
        }

fliplr<-function(matrix)
        {
            matrix<-matrix[nrow(matrix):1,];
            return(matrix);
        }
flipud<-function(matrix)
        {
            matrix<-matrix[,ncol(matrix):1];
            return(matrix);
        }


jet<-function(n)
      {
       m<-n
       n <- max(c(round(m/4),1));
       x = (1:n)/n;
       y = ((n/2):n)/n;
       e = rep(1,times=length(x));
       r = c(0*y, 0*e, x, e, flip(y));
       g = c(0*y, x, e, flip(x), 0*y);
       b = c(y, e, flip(x), 0*e, 0*y);
       J = cbind(r,g,b);
      while (nrow(J) > m)
      {
       J <- J[-c(1),];
      if(nrow(J) > m) J <- J[-c(nrow(J)),];
      }
      newpalette<-rgb(J[,1],J[,2],J[,3])
      return(newpalette);
      }


interp3D<-function(longitud,latitud,profundidad,valor,largo,ancho,fondo){
  
  library(lattice)
  library(akima)
  
  b<-as.vector(interaction(longitud,latitud, sep=":",drop=TRUE))
  pon<-which(!is.na(valor))
  longitud<-longitud[pon]
  latitud<-latitud[pon]
  profundidad<-profundidad[pon]
  valor<-valor[pon]
  b<-b[pon]
  xout<-seq(from=max(profundidad), to=min(profundidad), length=fondo)
  times<-matrix(NA,nrow=length(unique(b)),ncol=length(xout))
  latitude<-matrix(NA,nrow=length(unique(b)),ncol=length(xout))
  depth<-matrix(NA,nrow=length(unique(b)),ncol=length(xout))
  value<-matrix(NA,nrow=length(unique(b)),ncol=length(xout))
  
  d<-unique(b)
  for ( i in 1:length(d)) {
    colv<-rep(NA,length(xout))
    colz<-xout
    y<-valor[b==d[i]]
    x<-profundidad[b==d[i]]
    if(length(y)>1) {
      a<-approx(x,y,xout=xout)
      colz<-as.vector(a$x)
      colv<-as.vector(a$y)
    }
    colx<-rep(unique(longitud[b==d[i]]),times=length(xout))
    coly<-rep(unique(latitud[b==d[i]]),times=length(xout))
    
    times[i,]<-rbind(colx)
    latitude[i,]<-rbind(coly)
    depth[i,]<-rbind(colz)
    value[i,]<-rbind(colv)
  }
  
  longitud2<-rbind(as.vector(times))
  latitud2<-rbind(as.vector(latitude))
  prof2<-rbind(as.vector(depth))
  sinechos2<-rbind(as.vector(value))
  
  quita<-which(!is.na(sinechos2))
  longitud2<-longitud2[quita]
  latitud2<-latitud2[quita]
  prof2<-prof2[quita]
  sinechos2<-sinechos2[quita]
  p<-as.factor(prof2)
  k<-unique(p)
  
  xo<-seq(from=min(longitud), to=max(longitud), length=largo)
  yo<-seq(from=min(latitud), to=max(latitud), length=ancho)
  RESULTADOS<-array(NA,dim=c(length(xo),length(yo),length(k)))
  dimnames(RESULTADOS)<-list(xo,yo,k)
  
  for ( i in 1:length(k)) {
    z<-sinechos2[p==k[i]]
    y<-latitud2[p==k[i]]
    x<-longitud2[p==k[i]]
    j<-matrix(NA,nrow=length(xo), ncol=length(yo))
    if(length(z[!is.na(z)])>1) {
      g<-interp(x[!is.na(z)],y[!is.na(z)],z[!is.na(z)],xo,yo, duplicate="mean")
      j<-g$z
    }
    RESULTADOS[,,i]<-j  }
  return(RESULTADOS)
  }


interp3Drecom<-function(x, lon_x0=346, lon_x1=0, lat_y0=40, lat_y1=50, depth=120){     
          arr<-x
          #arr<-m[,,,1]
          #mu<-abind(arr[91:180,,],arr[1:90,,], along=1)   #  180 126 30
          #X<-c(filerc$X[91:180],filerc$X[1:90])         
          mu<-arr
          X<-filerc$X
    
          # Bay of Biscay #
          mbob<-mu[(X>=lon_x0 & X<=lon_x1),(filerc$Y>=lat_y0 & filerc$Y<=lat_y1),profun>=-(depth)]  # 6 10 6
          valor<-c(mbob, recursive=TRUE)
          longi<-    rep(dimnames(mbob)[[1]], times=dim(mbob)[2]*dim(mbob)[3])
          lati <-rep(rep(dimnames(mbob)[[2]], each =dim(mbob)[1]),times=dim(mbob)[3]) 
          profi<-    rep(dimnames(mbob)[[3]], each =dim(mbob)[1]*dim(mbob)[2])
  
          # resolution
                 boxes_depth<-profun[profun>=-(depth)]
          res_depth<-((max(boxes_depth, na.rm=TRUE))-(min(boxes_depth, na.rm=TRUE)))+1
                 boxes_long<-X[(X>=lon_x0 & X<=lon_x1)]
          res_long <-((max(boxes_long, na.rm=TRUE)+1)-(min(boxes_long, na.rm=TRUE)-1))*10 
                 boxes_lat<-filerc$Y[(filerc$Y>=lat_y0 & filerc$Y<=lat_y1)]
          res_lat  <-((max(boxes_lat, na.rm=TRUE)+1)-(min(boxes_lat, na.rm=TRUE)-1))*10 
  
          # Interpolate in 3D  
          valor[is.nan(valor)]<-NA
          longi <-as.numeric(longi[!is.na(valor)])
          lati  <-as.numeric(lati[!is.na(valor)])
          profi <-as.numeric(profi[!is.na(valor)])
          valor <-valor[!is.na(valor)]        
          library(akima)
          resul<-interp3D(longitud=longi,latitud=lati,profundidad=profi,valor=valor, largo=res_long, ancho=res_lat, fondo=res_depth)
          return(list(resul)) 
          }


interp2Drecom<-function(x, lon_x0=346, lon_x1=0, lat_y0=40, lat_y1=50){     
  mu<-x                            #  180 126 30
  X<-as.numeric(dimnames(mu)$x)
  Y<-as.numeric(dimnames(mu)$y)
  mu[mu==0]<-NA

  mbob <-mu[(X>=lon_x0 & X<=lon_x1),(Y>=lat_y0 & Y<=lat_y1)]   # dim(mbob) 
  valor<-c(mbob, recursive=TRUE)
  longi<-    rep(rownames(mbob), times=dim(mbob)[2])
  lati <-    rep(colnames(mbob), each =dim(mbob)[1]) 

  # resolution
  boxes_long<-X[(X>=lon_x0 & X<=lon_x1)]
  res_long <-((max(boxes_long, na.rm=TRUE)+1)-(min(boxes_long, na.rm=TRUE)-1))*10 
  boxes_lat<-Y[(Y>=lat_y0 & Y<=lat_y1)]
  res_lat  <-((max(boxes_lat, na.rm=TRUE)+1)-(min(boxes_lat, na.rm=TRUE)-1))*10 
  
  # Interpolate in 2D  
  valor[is.nan(valor)]<-NA
  longi <-as.numeric(longi[!is.na(valor)])
  lati  <-as.numeric(lati[!is.na(valor)])
  valor <-valor[!is.na(valor)]        
  library(akima)
  resul<-interp(x=longi,y=lati,z=valor,
          xo=seq(min(longi,na.rm=TRUE), max(longi,na.rm=TRUE), length = res_long),
          yo=seq(min(lati,na.rm=TRUE), max(lati,na.rm=TRUE), length = res_lat))
  res<-resul[[3]]   # dim(res)
  rownames(res)<-resul[[1]]
  colnames(res)<-resul[[2]]
  return(list(res)) 
}


interpSURFACErecom<-function(x, lon_x0=346, lon_x1=0, lat_y0=40, lat_y1=50){     
  arr<-x
  mu<-abind(arr[91:180,,],arr[1:90,,], along=1)   #  180 126 30
  X<-c(filerc$X[91:180],filerc$X[1:90])         
  mu[mu==0]<-NA
  #filled.contour(mbob,zlim=c(0.0,500),nlevels=20,color.palette=topo.colors, xaxt="n", yaxt="n",
  #                            plot.axes={axis(2, at=seq(40,50,by=2), labels=FALSE, las=1)
  #                             axis(1, at=seq(-20,0,by=2), labels=FALSE, las=1)})   
  # Bay of Biscay #
  mbob <-mu[(X>=lon_x0 | X<=lon_x1),(filerc$Y>=lat_y0 & filerc$Y<=lat_y1),profun>=-(10)]  # 7 5 7
  valor<-c(mbob, recursive=TRUE)
  longi<-    rep(rownames(mbob), times=dim(mbob)[2])
  lati <-    rep(colnames(mbob), each =dim(mbob)[1]) 
  
  # resolution
  boxes_long<-X[(X>=lon_x0 | X<=lon_x1)]-360
  res_long <-((max(boxes_long, na.rm=TRUE)+1)-(min(boxes_long, na.rm=TRUE)-1))*10 
  boxes_lat<-filerc$Y[(filerc$Y>=lat_y0 & filerc$Y<=lat_y1)]
  res_lat  <-((max(boxes_lat, na.rm=TRUE)+1)-(min(boxes_lat, na.rm=TRUE)-1))*10 
  
  # Interpolate in 2D  
  valor[is.nan(valor)]<-NA
  longi <-as.numeric(longi[!is.na(valor)])
  lati  <-as.numeric(lati[!is.na(valor)])
  valor <-valor[!is.na(valor)]        
  library(akima)
  resul<-interp(x=longi,y=lati,z=valor,
                xo=seq(min(longi,na.rm=TRUE), max(longi,na.rm=TRUE), length = res_long),
                yo=seq(min(lati,na.rm=TRUE), max(lati,na.rm=TRUE), length = res_lat))
  
  res<-resul[[3]]
  rownames(res)<-resul[[1]]
  colnames(res)<-resul[[2]]
  return(list(res)) 
}

interpolateSurface <- function(m,depth=10){
  md<-array(data = NA, dim=c(180,90,dim(m)[4]-1), dimnames = NULL)   
  for (k in 2:dim(m)[4]) {
    z<-m[,,profun>=-(depth),k]
    z[which(z==0)] <- NA       
    ZETA<-z 
    X<-filerc$X
    
    # Interpolate in 2D: para que el mapa no se deforme en el sur 
    nuevo<-seq(-79,79,by=2)
    viejas<-as.numeric(rownames(ZETA))
    ZETA_NEW<-matrix(NA,ncol=length(nuevo), nrow=nrow(ZETA))
    rownames(ZETA_NEW)<-viejas
    colnames(ZETA_NEW)<-nuevo
    indices<-as.numeric(colnames(ZETA))
    
    for (i in 1:(length(nuevo)-1)){
      cuales <- which(indices>=nuevo[i] & indices<nuevo[i+1])
      mat<-ZETA[,cuales]
      mat <- matrix(mat,nrow=nrow(ZETA_NEW))
      ZETA_NEW[,i]<-apply(mat, MARGIN=1, mean, na.rm=TRUE)
    }  # end loop i
    
    range(ZETA_NEW, na.rm=TRUE)
    # Anadir bloque NAS artico     
    cajasur<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
    colnames(cajasur)<-seq(-89,-81,by=2)
    cajanorte<-matrix(NA,ncol=5,nrow=nrow(ZETA_NEW))
    colnames(cajanorte)<-seq(81,89,by=2)
    res<-cbind(cajasur, ZETA_NEW, cajanorte)   
    rownames(res)<-rownames(ZETA_NEW)          
    md[,,k-1] <- res  }  # end loop k
  
  dimnames(md) <- list(x=X, y=colnames(res), z=fechas_recom[-1])
  range(md, na.rm=TRUE)
  return(md)}     


taylor.diagram.mine <- function (ref, model, add = FALSE, col = "red", pch = 19, pos.cor = TRUE, 
          xlab = "", ylab = "", main = "Taylor Diagram", show.gamma = TRUE, 
          ngamma = 3, gamma.col = 8, gamma.cex=1, sd.arcs = 0, ref.sd = FALSE, sd.method = "sample", 
          grad.corr.lines = c(0.2, 0.4, 0.6, 0.8, 0.9, 0.99), pcex = 1, cex.axis = 1, 
          normalize = FALSE, mar = c(5, 4, 6, 6), maxsd=2, maxray=2, minsd=0.0, ...) 
{
  grad.corr.full <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.99)
  R <- cor(ref, model, use = "pairwise")
  if (is.list(ref)) 
    ref <- unlist(ref)
  if (is.list(model)) 
    ref <- unlist(model)
  SD <- function(x, subn) {
    meanx <- mean(x, na.rm = TRUE)
    devx <- x - meanx
    ssd <- sqrt(sum(devx * devx, na.rm = TRUE)/(length(x[!is.na(x)]) - subn))
    return(ssd)
  }
  subn <- sd.method != "sample"
  sd.r <- SD(ref, subn)
  sd.f <- SD(model, subn)
  if (normalize) {
    sd.f <- sd.f/sd.r
    sd.r <- 1
  }
  #maxsd <- 1.5 * max(sd.f, sd.r)
  oldpar <- par("mar", "xpd", "xaxs", "yaxs")
  if (!add) {
    if (pos.cor) {
      if (nchar(ylab) == 0) 
        ylab = "sd model / sd ref"
      par(mar = mar)
      plot(0, xlim = c(minsd, maxsd), ylim = c(minsd, maxsd), xaxs = "i", 
           yaxs = "i", axes = FALSE, main = main, xlab = xlab, 
           ylab = "", type = "n", cex = cex.axis, cex.lab = cex.axis)
      if (grad.corr.lines[1]) {
        tramo <-(maxsd-minsd) 
        for (gcl in grad.corr.lines){
          lines(c(minsd, minsd+(tramo*gcl)), c(minsd, minsd+(tramo*sqrt(1 - gcl^2))), lty = 3)
      }}
      segments(c(minsd, minsd), c(minsd, minsd), c(minsd, maxsd), c(maxsd,minsd))
      axis.ticks <- seq(minsd, maxsd, by=0.2)
      axis.ticks <- axis.ticks[axis.ticks <= maxsd]
      axis(1, at = axis.ticks, cex.axis = cex.axis, las=1)
      axis(2, at = axis.ticks, cex.axis = cex.axis, las=1)
      if (sd.arcs[1]) {
        if (length(sd.arcs) == 1) 
          sd.arcs <- axis.ticks
        for (sdarc in sd.arcs) {
          xcurve <- cos(seq(0, pi/2, by = 0.03)) * sdarc
          ycurve <- sin(seq(0, pi/2, by = 0.03)) * sdarc
          lines(xcurve, ycurve, lty = 3)
        }
      }
      if (show.gamma[1]) {
        if (length(show.gamma) > 1) 
          gamma <- show.gamma
        else gamma <- pretty(c(0, maxsd), n = ngamma)[-1]
        if (gamma[length(gamma)] > maxsd) 
          gamma <- gamma[-length(gamma)]
        labelpos <- seq(45, 70, length.out = length(gamma))
        for (gindex in 1:length(gamma)) {
          xcurve <- cos(seq(0, pi, by = 0.03)) * gamma[gindex] + 
            sd.r
          endcurve <- which(xcurve < 0)
          endcurve <- ifelse(length(endcurve), min(endcurve) - 1, 105)
          ycurve <- sin(seq(0, pi, by = 0.03)) * gamma[gindex]
          maxcurve <- xcurve * xcurve + ycurve * ycurve
          startcurve <- which(maxcurve > maxsd * maxsd)
          startcurve <- ifelse(length(startcurve), max(startcurve) + 
                                 1, 0)
          lines(xcurve[startcurve:endcurve], ycurve[startcurve:endcurve], 
                col = gamma.col, lty=2)
          if (xcurve[labelpos[gindex]] > 0) 
            boxed.labels(xcurve[labelpos[gindex]], ycurve[labelpos[gindex]], 
                         gamma[gindex], border = FALSE, col = gamma.col, cex=gamma.cex)
        }
      }
      tramo <-(maxsd-minsd)
      xcurve <- minsd+(cos(seq(0, pi/2, by = 0.01)) * tramo)
      ycurve <- minsd+(sin(seq(0, pi/2, by = 0.01)) * tramo)
      lines(xcurve, ycurve)
      bigtickangles <- acos(seq(0.1, 0.9, by = 0.1))
      medtickangles <- acos(seq(0.05, 0.95, by = 0.1))
      smltickangles <- acos(seq(0.91, 0.99, by = 0.01))
      segments(minsd+(cos(bigtickangles) * tramo),
               minsd+(sin(bigtickangles) * tramo),
               minsd+(cos(bigtickangles) * 0.97 * tramo),
               minsd+(sin(bigtickangles) * 0.97 * tramo))
     
      bigtickangles <- bigtickangles[c(seq(2,8,by=2),9)]
      
      par(xpd = TRUE)
      if (ref.sd) {
        tramo2 <- sd.r-minsd
        xcurve <- minsd+(cos(seq(0, pi/2, by = 0.01)) * tramo2)
        ycurve <- minsd+(sin(seq(0, pi/2, by = 0.01)) * tramo2)
        lines(xcurve, ycurve)
      }
      points(sd.r, minsd, cex = pcex)
      text(minsd+(cos(c(bigtickangles, acos(c(0.99)))) * 1.05 * tramo),
           minsd+(sin(c(bigtickangles, acos(c(0.99)))) * 1.05 * tramo),
           c(seq(0.2, 0.8, by = 0.2), 0.9, 0.99), cex = cex.axis)
      
      #text(minsd+(tramo * 0.9), minsd+(tramo * 0.8), "Correlation", srt = 305 , cex = cex.axis)
      segments(minsd+(cos(medtickangles) * tramo),
               minsd+(sin(medtickangles) * tramo),
               minsd+(cos(medtickangles) * 0.98 * tramo),
               minsd+(sin(medtickangles) * 0.98 * tramo))
      
      segments(minsd+(cos(smltickangles) * tramo),
               minsd+(sin(smltickangles) * tramo),
               minsd+(cos(smltickangles) * 0.99 * tramo),
               minsd+(sin(smltickangles) * 0.99 * tramo))
    }
    else {
      x <- ref
      y <- model
      R <- cor(x, y, use = "pairwise.complete.obs")
      E <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
      xprime <- x - mean(x, na.rm = TRUE)
      yprime <- y - mean(y, na.rm = TRUE)
      sumofsquares <- (xprime - yprime)^2
      Eprime <- sqrt(sum(sumofsquares)/length(complete.cases(x)))
      E2 <- E^2 + Eprime^2
      if (add == FALSE) {
        #maxray <- 1.5 * max(sd.f, sd.r)
        par(mar = mar)
        lim <- maxsd+(maxsd*0.4)
        plot(c(-maxray, maxray), c(0, maxray), type = "n", 
             asp = 1.00, bty = "n", xaxt = "n", yaxt = "n", 
             xlab = xlab, ylab = ylab, main = main,
             cex = cex.axis, cex.lab = cex.axis, las=1,
             xlim=c(-lim,lim), xaxs="i")
        discrete <- seq(180, 0, by = -1)
        listepoints <- NULL
        for (i in discrete) {
          listepoints <- cbind(listepoints, maxray * 
                                 cos(i * pi/180), maxray * sin(i * pi/180))
        }
        listepoints <- matrix(listepoints, 2, length(listepoints)/2)
        listepoints <- t(listepoints)
        lines(listepoints[, 1], listepoints[, 2])
        lines(c(-maxray, maxray), c(0, 0))
        lines(c(0, 0), c(0, maxray))
        for (i in grad.corr.lines) {
          lines(c(0, maxray * i), c(0, maxray * sqrt(1 - i^2)), lty = 3)
          lines(c(0, -maxray * i), c(0, maxray * sqrt(1 - i^2)), lty = 3)
        }
        for (i in grad.corr.full) {
          text(1.05 * maxray * i, 1.05 * maxray * sqrt(1 - i^2), i, cex = cex.axis)
          text(-1.05 * maxray * i, 1.05 * maxray * sqrt(1 - i^2), -i, cex = cex.axis)
        }
        
        #seq.sd <- seq.int(0, 2 * maxray, by = (maxray/10))[-1]
        #for (i in seq.sd) {
        #  xcircle <- sd.r + (cos(discrete * pi/180) *i)
        #  ycircle <- sin(discrete * pi/180) * i
        #  for (j in 1:length(xcircle)) {
        #    if ((xcircle[j]^2 + ycircle[j]^2) < (maxray^2)) {
        #      points(xcircle[j], ycircle[j], col = "darkgreen", pch = ".")
        #      if (j == 10) 
        #        text(xcircle[j], ycircle[j], signif(i, 2), cex = 0.5, col = "darkgreen")
        #    }
        #  }
        # }
        seq.sd <- seq.int(0, maxray, length.out = 5)
        for (i in seq.sd) {
          xcircle <- (cos(discrete * pi/180) * i)
          ycircle <- sin(discrete * pi/180) * i
          if (i) 
            lines(xcircle, ycircle, lty = 3, col = "black")
          #text(min(xcircle), -0.03 * maxray, signif(i, 2), cex = cex.axis, col = "black")
          #text(max(xcircle), -0.03 * maxray, signif(i, 2), cex = cex.axis, col = "black")
        }
        #text(0, -0.08 * maxray, "Standard Deviation", cex = cex.axis, col = "black")
        #text(0, -0.12 * maxray, "Centered RMS Difference", 
        #cex = 0.7, col = "darkgreen")
        #points(sd.r, 0, pch = 22, bg = "darkgreen", cex = 1.1)
        #text(0, 1.1 * maxray, "Correlation Coefficient", cex = cex.axis, col = "black")
        #text(maxray * 0.8, maxray * 0.8, "Correlation", srt = 315 , cex = cex.axis)
        #text(maxray * 0.9, maxray * 0.8, "Correlation", srt = 305 , cex = cex.axis)
        }
      S <- (2 * (1 + R))/(sd.f + (1/sd.f))^2
      if (ref.sd) {
        xcurve <- cos(seq(0, pi/2, by = 0.01)) * sd.r
        ycurve <- sin(seq(0, pi/2, by = 0.01)) * sd.r
        lines(xcurve, ycurve)
        points(sd.r, 0, cex = pcex)
      }
    }
  }
  points(sd.f * R, sd.f * sin(acos(R)), pch = pch, col = col,cex = pcex)
  invisible(oldpar)
}

find_eudepth <- function(x, prof, dif=0.01){
    l_sur <- x[1]
    light <- x[-1]
    light[light<=1e-08] <- 0
    light[light==0] <- NA
    depth <- prof[!is.na(light)]
    light <- light[!is.na(light)]
    if (length(light)>=2){
        #plot(light, depth, las=1, type="b", xlim=c(0,40), ylim=c(-300,0))
        #points(l_sur,0, las=1, pch=19) 
        valor <- (light/l_sur)-0.01
        index <- which(abs(valor)==min(abs(valor)))
        #points(light[index],  depth[index], pch=19, col="red")
        minimos <- sort(abs(valor),na.last=NA)[1:4]
        cuales <- match(minimos,abs(valor))
        cuales <- cuales[!is.na(cuales)]
        #points(light[cuales], depth[cuales], pch=21, col="green", lwd=2)
        valorcitos <- (light/l_sur)[cuales]
        # plot((light/l_sur),depth) 
        # points(valorcitos,depth[cuales], pch=19)
        if (sum(!is.na(valorcitos))>=2){
            res <- approx(x=valorcitos, y=depth[cuales], xout=0.01, rule=2)
            res <- res$y
        } else if (length(index)==1){res <- depth[index]}else{res<-NA}
    } else{res<-NA}
    # abline(h=res)
    return(res)}

find_opticaldepth<-function(x){
  light <- x
  light[light<=5.5e-5] <- 0
  light[depth>100] <- 0
  light[light==0] <- NA
  depth <- depth[!is.na(light)]
  light <- light[!is.na(light)]
  if (length(light)>=2){
    regre<-lm(log(light)~depth)
    Eo<-exp(regre$coefficients[1])
    kd<-(-regre$coefficients[2])
    #equis<-seq(0,100,by=1)
    #ies<-Eo*exp(-kd*equis)
    #plot(depth, light, las=1,type="b", col="red")
    #points(equis, ies, col="green")
    res<-1/kd
  } else{res<-NA}
  # abline(h=res)
  return(res)} 

find_kd<-function(x){
  light <- x
  light[light<=5.5e-5] <- 0
  light[depth>100] <- 0
  light[light==0] <- NA
  depth <- depth[!is.na(light)]
  light <- light[!is.na(light)]
  if (length(light)>=2){
    regre<-lm(log(light)~depth)
    Eo<-exp(regre$coefficients[1])
    res<-(-regre$coefficients[2])
  } else{res<-NA}
  return(res)}