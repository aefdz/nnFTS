#' Envelope Functional Time Series Forecasting
#'
#' @param data matrix PXN being N the total number of functions and
#' P the total number of observed values for each function. The columns are
#' ordered by time from 0 to T.
#' @param focal the name of the curve to envelope. In the article, the most recent curve.
#' @param h forecasting horizon.
#' @param distance vector of distances $D(i,N)$, i.e. distance between
#' the last observed curve and the others.
#' @param typePoint type of point estimate: "w" for weighted mean and "expw" for exponentially weighted mean.
#' @param theta parameter for the "expw" point type. The default value is 1.
#'
#' @return a matrix containing the point forecast.
#' @examples
#'
#' # One-period-ahead
#' data(electricityDemand)
#' focal <- "saturday/29/12/2018"
#' data <- rainbow::fts(electricityDemand$x, electricityDemand$y[,1:1825])
#' point <- envelope.forecast(data, focal, h = 1, distance = "l2", typePoint = "expw", theta = 1)
#'
#' # DU half day
#' focal <- "monday/31/12/2018"
#' data <- electricityDemand
#' data$y[72:144, focal] <- NA
#' point <- envelope.forecast(data, focal, h = 1, distance = "l2", typePoint = "expw", theta = 1)
#'
#' @import matrixStats
#' @importFrom rainbow fts
#'
#' @export
envelope.forecast <- function(data, focal, h, distance, typePoint, theta = 1, max_iter){
  if(anyNA(data$y[,focal])){ #DU
    tcut <- sum(!is.na(data$y[,focal]))
    results <- envelope.DU(data, focal, tcut, distance)
    selected <- results$Jordered[!results$Jordered %in% focal]

    if(distance=='supremum'){
      dist<-matrixStats::colMaxs(abs(matrix(data$y[1:tcut,selected]-data$y[1:tcut,focal], ncol=length(selected), nrow=length(data$x))))
      names(dist)<-c(selected)}
    if(distance=='l2'){
      dist <- colSums((data$y[1:tcut,selected]-data$y[1:tcut,focal])^2)
      names(dist)<-c(selected)}

    envelopeAhead <- selected

    w <- wDistance(dist, theta = theta, typePoint = typePoint)

    pointPrediction <- t(as.matrix(w[selected]))%*%t(data$y[,envelopeAhead])

  }else{
    results <- envelope(data, focal, distance = distance, plot = FALSE, ...)
    selected <- results$Jordered[!results$Jordered %in% focal]

    if(distance=='supremum'){
      dist<-matrixStats::colMaxs(abs(matrix(data$y[,selected]-data$y[,focal], ncol=length(selected), nrow=length(data$x))))
      names(dist)<-c(selected)}
    if(distance=='l2'){
      dist <- colSums((data$y[,selected]-data$y[,focal])^2)
      names(dist)<-c(selected)}

    envelopeAhead <- neighbours_H_Ahead(data, selected, h = h)

    w <- wDistance(distancesToFocal = dist, theta = theta, typePoint = typePoint)

    pointPrediction <- t(as.matrix(w[selected]))%*%t(data$y[,envelopeAhead])
  }

  return(pointPrediction)
}

envelope.DU <- function(data, focal, tcut, dist){
  datayshort<-data$y[1:tcut,]; dataxshort<-data$x[1:tcut]
  rownames(datayshort) <- dataxshort
  dataShort<-rainbow::fts(dataxshort, datayshort, xname='x', yname='y')

  SubsampleAll <- envelope(dataShort, focal, dist, plot = FALSE)
  Subsample<-SubsampleAll$Jordered
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
  orderMBDObserved <- orderByMBD(dataShort, Subsample)
  Jordered <- orderMBDObserved

  output<-list(Jordered = names(Jordered$x)[!names(Jordered$x) %in% focal])
  return(output)
}

wDistance<-function(distancesToFocal, theta = 1, typePoint){
  if(typePoint == 'w'){
    w <- (1/distancesToFocal)/sum(1/distancesToFocal)
    names(w) <- names(distancesToFocal)
  }
  if(typePoint == 'expw'){
    minDist <- min(distancesToFocal)
    expDistances <- exp(-theta*distancesToFocal/minDist)
    w<- expDistances/sum(expDistances)
    names(w)<-names(distancesToFocal)
  }
  return(w)
}

neighbours_H_Ahead <- function(data, neighbours, h = 1){
  neighboursAhead <- vector(mode='numeric', length = 0)
  idAll <- colnames(data$y)
  for(j in c(1:length(neighbours))){
    oneAhead <- which(idAll == neighbours[j]) + h
    if(oneAhead <= length(idAll)){
      neighboursAhead <- c(neighboursAhead, c(idAll)[oneAhead])
    }
  }
  return(neighboursAhead)
}
