#' Plot Bands
#'
#' This function plot the band given by the k deepest or closest functions.
#'
#' @param data a fts or fds object see R package rainbow or ftsa.
#' @param focal the name of the curve to envelope. In the article, the most recent curve.
#' @param Joredered set of curves ordered by distance or depth.
#' @param kcurves number of curves to build the band.
#'
#' @return a list containing the selected curves.
#'
#' @import ggplot2
#' @import matrixStats
#' @importFrom reshape2 melt
#' @importFrom rainbow fts
#'
#' @export
plotBand <- function(data, focal, Jordered, kcurves){
  cut <- sum(!is.na(data$y[,focal]))

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
