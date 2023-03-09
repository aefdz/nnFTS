#' Envelope Algorithm
#'
#' This function implement the algorithm to optain neighbourhoods as in the paper.
#'
#' @param data a fts or fds object see R package rainbow or ftsa.
#' @param focal the name of the curve to envelope. In the article, the most recent curve.
#' @param distance a character string naming the distance. Available options are "supremum" and "l2".
#' @param plot if plotting the selected curves in the different iterations of the algorithm.
#'
#' @return a list containing the selected curves.
#'
#' @import ggplot2
#' @import matrixStats
#' @importFrom reshape2 melt
#' @importFrom rainbow fts
#'
#' @examples
#' data(sinedata)
#' focal <- '1'
#' resultsBand <- envelope(sinedata, focal, distance = "l2", plot = FALSE)
#'
#' @export
envelope <- function(data, focal, distance, plot, max_iter){
  P<-data$x[data$y[,focal]!=matrixStats::rowMaxs(data$y) & data$y[,focal]!=matrixStats::rowMins(data$y)]
  length_P<-length(P)
  if(length_P>0){
    Subsample<-vector(mode='numeric', length = 0)
    I <- colnames(data$y)[!colnames(data$y) %in% focal]
    iter<-c(1)
    iterDepth<-c(0)
    while(length(I)>1){
      x <- as.character(P)

      candidates<-I
      if(distance=='supremum'){
        dist<-matrixStats::colMaxs(abs(matrix(data$y[,I]-data$y[,focal], ncol=length(I), nrow=length(data$x))))
        names(dist)<-c(I)}
      if(distance=='l2'){
        dist <- colSums((data$y[,I]-data$y[,focal])^2)
        names(dist)<-c(I)}

      d <- sort(dist, decreasing = FALSE)

      SubsampleIter <- names(d[1])
      candidates <- candidates[!candidates %in% names(d[1])]

      while(length(x)!=0 & 0<length(candidates)){
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
          fdata<-reshape2::melt(data$y, id='x')
          colnames(fdata) <- c("X1", "X2", "value")

          aux<-matrix(c(data$y[,SubsampleIter]), nrow=length(data$x))
          low<-  vector(mode="numeric", length=length(data$x))
          high<-  vector(mode="numeric", length=length(data$x))
          for (j in c(1:length(data$x))){
            low[j]<-min(aux[j,])
            high[j]<-max(aux[j,])
          }
          dataBand<-data.frame(X1=c(data$x, rev(data$x)), X2='NULL', value=c(high, rev(low)))
          pl1 <- ggplot2::ggplot(data=fdata, ggplot2::aes(x=.data$X1, y=.data$value,colour=as.factor(.data$X2), group=as.factor(.data$X2))) +
            ggplot2::geom_line(color='grey50')+
            ggplot2::theme(legend.position="none")+
            ggplot2::geom_line(data=fdata[fdata$X2 %in% Subsample,], color='black', cex=1.25)+
            ggplot2::geom_line(data=fdata[fdata[,2]==focal,],color='red', cex=1.25)+
            ggplot2::geom_polygon(data=dataBand, ggplot2::aes(x=.data$X1,y=.data$value), color='grey25')+
            ggplot2::geom_line(data=fdata[fdata$X2 %in% Subsample,], color='black', cex=1)+
            ggplot2::geom_line(data=fdata[fdata$X2 %in% SubsampleIter,], color='grey', cex=1)+
            ggplot2::geom_line(data=fdata[fdata$X2==focal,], color='red', cex=1.2)
          print(pl1)

          readline(prompt = "Press <Enter> to continue...")
        }

      }
      if(length(candidates)==0){break}
      if(!missing(max_iter)){if(max_iter == max(iter)-1){break}}
    }
  }else{Subsample=c('No Subsample')}
  if(distance=='supremum'){dist <- matrixStats::colMaxs(abs(matrix(data$y-data$y[,focal], ncol=length(I), nrow=length(data$x)))); names(dist)<-colnames(data$y)}
  if(distance=='l2'){dist <- colSums((data$y-data$y[,focal])^2); names(dist)<-colnames(data$y)}

  if(Subsample[1]!='No Subsample'){
    auxDist <- sort(dist, decreasing = FALSE, index.return=TRUE)
    auxDepth <- fMBD(data$y[,c(focal,Subsample)]);
    return(list(Jordered=names(auxDepth$x)))
  }else{
    return(list(Jordered='No Subsample'))
  }
}

#Quit modified band depth computation
fMBD <- function(data){
  p=dim(data)[1]
  n=dim(data)[2]
  rmat=apply(data,1,rank)
  down=rmat-1
  up=n-rmat
  depth<-(rowSums(up*down)/p+n-1)/combinat(n,2)
  fMBD<-sort(depth, decreasing=TRUE,index.return=TRUE)
  return(fMBD)
}

combinat <- function(n,p){
  if (n<p){combinat=0}
  else {combinat=exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))}
}

#Order by depth measure
orderByMBD <- function(data, alphaJ){
  y <- t(data$y);
  mbd <- fMBD(t(y[alphaJ,])); names(mbd$x) <- alphaJ[mbd$ix]
  return(mbd)
}

depthPercentile <- function(depthOrder,focal){
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

plotBand <- function(data, cut, Jordered, kcurves, focal){
  dataShort<- rainbow::fts(data$x[1:cut], data$y[1:cut,], xname=data$xname, yname = data$yname) # , cex.axis= 0.5, cex.lab=0.5)

  #Rule, take de min of normJordered and k
  J<-min(length(Jordered[!Jordered %in% focal]),kcurves)
  bagId <- Jordered[!Jordered %in% focal][1:J];

  fdata<-reshape2::melt(data$y, id='x')
  colnames(fdata) <- c("X1", "X2", "value")
  fdatashort<-reshape2::melt(dataShort$y, id='x')
  colnames(fdatashort) <- c("X1", "X2", "value")

  #Now we compute the bag
  bag <- band(data,bagId); low<-bag[,1]; high<-bag[,2]
  dataBand<-data.frame(X1=c(data$x, rev(data$x)), X2='NULL', value=c(high, rev(low)))
  pl1 <- ggplot2::ggplot(data=fdata, ggplot2::aes(x=.data$X1, y=.data$value,colour=as.factor(.data$X2), group=as.factor(.data$X2))) +
    ggplot2::geom_line(color='grey50')+
    ggplot2::theme(legend.position="none")+
    ggplot2::geom_polygon(data=dataBand, ggplot2::aes(x=.data$X1,y=.data$value), color='grey25')+
    ggplot2::geom_line(data=fdata[fdata[,2]==focal,],color='red', cex=1.25,linetype = 2)+
    ggplot2::geom_line(data=fdatashort[fdatashort$X2==focal,], color='red', cex=1.2)+
    ggplot2::geom_vline(xintercept = fdata$X1[cut],cex=1)

  print(pl1)
}
