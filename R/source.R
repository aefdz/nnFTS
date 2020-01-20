#' @export
#Algorithm

banddpeeling<-function(data, focal, dist, plot){
  P<-data$x[data$y[,focal]!=matrixStats::rowMaxs(data$y) & data$y[,focal]!=matrixStats::rowMins(data$y)]
  length_P<-length(P)
  if(length_P>0){
    Subsample<-vector(mode='numeric', length = 0)
    I <- colnames(data$y)[!colnames(data$y) %in% focal]
    iter<-c(1)
    iterDepth<-c(0)
    while(length(I)>1){
      x <- as.character(P)
      #Primera etapa cogemos dos
      candidates<-I
      if(dist=='supremum'){supDist<-colMaxs(abs(matrix(data$y[,I]-data$y[,focal], ncol=length(I), nrow=length(data$x)))); names(supDist)<-c(I)}
      if(dist=='l2'){supDist <- colSums((data$y[,I]-data$y[,focal])^2); names(supDist)<-c(I)}
      d <- sort(supDist, decreasing = FALSE)

      SubsampleIter<-names(d[1])
      candidates=candidates[!candidates %in% names(d[1])]

      while(length(x)!=0 & 0<length(candidates)){
        #supDist<-colMaxs(abs(matrix(data$y[,candidates]-data$y[,focal], ncol=length(candidates), nrow=length(data$x)))); names(supDist)<-c(candidates)
        #supDist <- colSums((as.matrix(data$y[,candidates]-data$y[,focal])^2)); names(supDist) <- c(candidates)
        dAUX<-sort(d[candidates], decreasing = FALSE)

        Ji<-which(abs(rowSums(sign(as.matrix(data$y[,c(SubsampleIter, names(dAUX[1]))]-data$y[,focal]))))<length(c(SubsampleIter, names(d[1]))))

        # Envuelven algo las dos primeras?
        if(length(x)<=length(x[!x %in% names(Ji)])){ #No envuelven
          candidates=candidates[!candidates %in% names(dAUX[1])]
        }else{ # si Envuelven
          x<-x[!x %in% names(Ji)] #remaining points to cover
          SubsampleIter<-union(SubsampleIter,names(dAUX[1]))
          candidates<-candidates[!candidates %in% SubsampleIter]
          I=I[!I %in% SubsampleIter]
        }
      }
      iter<-c(iter, iter[length(iter)]+1)
      aux<-fMBD(data$y[,c(focal,Subsample,SubsampleIter)]);
      iterDepthAux<-depthPercentile(names(aux$x),focal)
      iterDepth<-c(iterDepth,iterDepthAux)
      if(max(iterDepth[1:(length(iter)-1)])<=iterDepth[iter[length(iter)]]){
        Subsample<-c(Subsample, SubsampleIter)
        if(plot==TRUE){
          #print(c('X for covering', length(x)))
          #print(c('Remaining Sample',length(I)))
          #print(c('In this envelope',length(SubsampleIter)))
          #print(iterDepthAux)
          fdata<-reshape::melt(data$y, id='x')
          aux<-matrix(c(data$y[,SubsampleIter]), nrow=length(data$x))
          low<-  vector(mode="numeric", length=length(data$x))
          high<-  vector(mode="numeric", length=length(data$x))
          for (j in c(1:length(data$x))){
            low[j]<-min(aux[j,])
            high[j]<-max(aux[j,])
          }
          dataBand<-data.frame(X1=c(data$x, rev(data$x)), X2='NULL', value=c(high, rev(low)))
          pl1 <- ggplot2::ggplot(data=fdata, ggplot2::aes(x=X1, y=value,colour=as.factor(X2), group=as.factor(X2))) +
            ggplot2::geom_line(color='grey50')+
            ggplot2::theme(legend.position="none")+
            ggplot2::geom_line(data=fdata[fdata$X2 %in% Subsample,], color='black', cex=1.25)+
            ggplot2::geom_line(data=fdata[fdata[,2]==focal,],color='red', cex=1.25)+
            ggplot2::geom_polygon(data=dataBand, ggplot2::aes(x=X1,y=value), color='grey25')+
            ggplot2::geom_line(data=fdata[fdata$X2 %in% Subsample,], color='black', cex=1)+
            ggplot2::geom_line(data=fdata[fdata$X2 %in% SubsampleIter,], color='grey', cex=1)+
            ggplot2::geom_line(data=fdata[fdata$X2==focal,], color='red', cex=1.2)
          print(pl1)

          readline(prompt = "Press <Enter> to continue...")
        }

      }
      if(length(candidates)==0){break}
    }
  }else{Subsample=c('No Subsample')}
  if(dist=='supremum'){supDist<-colMaxs(abs(matrix(data$y-data$y[,focal], ncol=length(I), nrow=length(data$x)))); names(supDist)<-colnames(data$y)}
  if(dist=='l2'){supDist <- colSums((data$y-data$y[,focal])^2); names(supDist)<-colnames(data$y)}

  if(Subsample[1]!='No Subsample'){
    auxDist <- sort(supDist, decreasing = FALSE, index.return=TRUE)
    auxDepth<-fMBD(data$y[,c(focal,Subsample)]);
    return(list(subsample=names(auxDepth$x)))
  }else{
    return(list(subsample='No Subsample'))
  }
}

#Forecasting
extension <- function(data, focal, tcut, dist){
  datayshort<-data$y[1:tcut,]; dataxshort<-data$x[1:tcut]
  rownames(datayshort) <- dataxshort
  dataShort<-rainbow::fds(dataxshort, datayshort, xname='x', yname='y')
  sco <-  PCAproj(t(dataShort$y), center = median)$scores; rownames(sco) <-  colnames(dataShort$y)

  SubsampleAll <- banddpeeling(dataShort, focal, dist, plot = FALSE)
  Subsample<-SubsampleAll$subsample
  if (length(Subsample) - 1 < 2) { #if the focal is a missing case
    dimJ <- 'No Subsample'
    width <- 'No Subsample'
    depth <- c('No Subsample','No Subsample','No Subsample')
    stats <- c(dimJ, width, depth);
    output<-list(Jordered='No Subsample', bag=cbind('No Subsample','No Subsample'), point='No Subsample', stat=stats)
    return(output)

    Jordered <- orderByMBD(dataShort, Subsample) #Ordered with focal
  }

  depth <- matrix(data = 0, nrow = 1, ncol = 4); colnames(depth) <- c('TukeyObserved', 'MBDObserved','TukeyPredicted', 'MBDPredicted')

  dimJ <- length(Subsample) - 1

  #Now we have the curves alphaJ that build the bag
  #Order them by tukey and by MBD
  orderMBDObserved<-orderByMBD(dataShort, Subsample)
  depth[1,1] <- depthPercentile(names(orderMBDObserved$x), focal)

  orderMBDPredicted<-orderByMBD(extract(data, direction = "x", xorder = data$x[(tcut+1):length(data$x)]), Subsample)
  depth[1,2] <- depthPercentile(names(orderMBDPredicted$x), focal)

  Jordered <- orderMBDObserved
  #Now compute the statistics
  bagId <- names(Jordered$x)[!names(Jordered$x) %in% focal]
  bag <- band(data,bagId); low<-bag[,1]; high<-bag[,2]
  bagTotal <- band(data,colnames(data$y)); lowTotal<-bagTotal[,1]; highTotal<-bagTotal[,2]

  coverageObserved<-coverageStat(data$y[,focal][1:tcut], high[1:tcut], low[1:tcut])/possibleCoverage(rainbow::fds(data$x[1:tcut],data$y[1:tcut,]),focal)
  coveragePredicted<-coverageStat(data$y[,focal][(tcut+1):length(data$x)], high[(tcut+1):length(data$x)],low[(tcut+1):length(data$x)])/possibleCoverage(rainbow::fds(data$x[(tcut+1):length(data$x)],data$y[(tcut+1):length(data$x),]),focal)

  widthObserved <- widthStat(high[1:tcut], low[1:tcut])/widthStat(highTotal[1:tcut], lowTotal[1:tcut])
  widthPredicted <- widthStat(high[(tcut+1):length(data$x)], low[(tcut+1):length(data$x)])/widthStat(highTotal[(tcut+1):length(data$x)], lowTotal[(tcut+1):length(data$x)])

  stats <- cbind(dimJ, widthObserved, widthPredicted, depth, coverageObserved, coveragePredicted); rownames(stats) <- focal
  output<-list(Jordered=names(Jordered$x)[!names(Jordered$x) %in% focal], stat=stats)
  return(output)
}

#Quit modified band depth computation
fMBD<-function(data){
  p=dim(data)[1]
  n=dim(data)[2]
  rmat=apply(data,1,rank)
  down=rmat-1
  up=n-rmat
  depth<-(rowSums(up*down)/p+n-1)/combinat(n,2)
  fMBD<-sort(depth, decreasing=TRUE,index.return=TRUE)
  return(fMBD)
}
combinat<-function(n,p){
  if (n<p){combinat=0}
  else {combinat=exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))}
}

#Order by depth measure
orderByMBD <- function(data, alphaJ){
  y <- t(data$y);
  mbd <- fMBD(t(y[alphaJ,])); names(mbd$x) <- alphaJ[mbd$ix]
  return(mbd)
}
depthPercentile<-function(depthOrder,focal){
  depthP<-1-(which(depthOrder==focal)-1)/(length(depthOrder)-1)
  return(depthP)
}

#Compute a band given the data and the indices of the functions that will build the band
band <- function(data, bagId){
  aux <- matrix(c(data$y[,bagId]), nrow = length(data$x))
  low <-  vector(mode = "numeric", length = length(data$x))
  up <-  vector(mode = "numeric", length = length(data$x))
  for (j in c(1:length(data$x))) {
    low[j] <- min(aux[j,], na.rm = TRUE)
    up[j] <- max(aux[j,], na.rm = TRUE)
  }
  band<-cbind(low, up); colnames(band)<-c('Lower', 'Upper')
  return(band=cbind(low, up))
}
coverageStat <- function(focaly, high, low){
  coverage <- sum(1*(low <= focaly & focaly <= high))/length(focaly)
  return(coverage)
}
possibleCoverage <- function(data, focal){
  aux <- matrix(c(data$y[,colnames(data$y)[!colnames(data$y) %in% focal]]), nrow = length(data$x))
  low <-  vector(mode = "numeric", length = length(data$x))
  high <-  vector(mode = "numeric", length = length(data$x))
  focaly<-data$y[,focal]
  for (j in c(1:length(data$x))) {
    low[j] <- min(aux[j,], na.rm = TRUE)
    high[j] <- max(aux[j,], na.rm = TRUE)
  }
  band<-cbind(low, high); colnames(band)<-c('Lower', 'Upper')
  coverage <- sum(1*(low <= focaly & focaly <= high))/length(focaly)
  return(coverage)
}
widthStat <- function(high, low){
  width <- sum(abs(high - low))/length(high)
  return(width)
}

plotBand<-function(data, cut, Jordered, kcurves, focal){
  dataShort<- rainbow::fds(data$x[1:cut], data$y[1:cut,], xname=data$xname, yname = data$yname) # , cex.axis= 0.5, cex.lab=0.5)

  #Rule, take de min of normJordered and k
  J<-min(length(Jordered[!Jordered %in% focal]),kcurves)
  bagId <- Jordered[!Jordered %in% focal][1:J];

  fdata<-reshape::melt(data$y, id='x')
  fdatashort<-reshape::melt(dataShort$y, id='x')

  #Now we compute the bag
  bag <- band(data,bagId); low<-bag[,1]; high<-bag[,2]
  dataBand<-data.frame(X1=c(data$x, rev(data$x)), X2='NULL', value=c(high, rev(low)))
  pl1 <- ggplot2::ggplot(data=fdata, ggplot2::aes(x=X1, y=value,colour=as.factor(X2), group=as.factor(X2))) +
    ggplot2::geom_line(color='grey50')+theme(legend.position="none")+
    ggplot2::geom_polygon(data=dataBand, ggplot2::aes(x=X1,y=value), color='grey25')+
    ggplot2::geom_line(data=fdata[fdata[,2]==focal,],color='red', cex=1.25,linetype = 2)+
    ggplot2::geom_line(data=fdatashort[fdatashort$X2==focal,], color='red', cex=1.2)+
    ggplot2::geom_vline(xintercept = fdata$X1[cut],cex=1)

  print(pl1)
}

