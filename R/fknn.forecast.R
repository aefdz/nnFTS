#' Functional knn forecasting
#'
#' @param data matrix PXN being N the total number of functions and
#' P the total number of observed values for each function. The columns are
#' ordered by time from 1 to T.
#' @param focal the name of the curve to envelope. In the article, the most recent curve.
#' @param k number of nearest neighbourds.
#' @param h forecasting horizon. If Dynamic Updating the forecasting horizon, h, will be such that
#' the partially observed part it completed.
#' @param distance vector of distances $D(i,N)$, i.e. distance between
#' the last observed curve and the others.
#' @param typePoint type of point estimate: "w" for weighted mean and "exp" for exponentially weighted mean.
#' @param theta parameter for the "exp" point type. The default value is 1.
#'
#' @return a matrix containing the point forecast.
#' @examples
#'
#'# One-period-ahead
#' data(electricityDemand)
#' focal <- "saturday/29/12/2018"
#' data <- rainbow::fts(electricityDemand$x, electricityDemand$y[,1:1825])
#' point <- fknn.forecast(data, focal, k = 5, h = 1, distance = "l2", typePoint = "expw", theta = 1)
#'
#' # DU half day
#' focal <- "monday/31/12/2018"
#' data <- electricityDemand
#' data$y[72:144, focal] <- NA
#' point <- fknn.forecast(data, focal, k = 5, h = 1, distance = "l2", typePoint = "expw", theta = 1)
#'
#' @import matrixStats
#' @importFrom rainbow fts
#'
#' @export

fknn.forecast <- function(data, focal, k, h, distance, typePoint, theta = 1){
  if(anyNA(data$y[,focal])){ #DU
    tcut <- sum(!is.na(data$y[,focal]))

    selected <- fknn.distance(data, focal, distance)[1:k]
    envelopeAhead <- names(selected)

    w <- wDistance(selected, theta = theta, typePoint = typePoint)

    pointPrediction <- t(as.matrix(w))%*%t(data$y[,envelopeAhead])

  }else{
    selected <- fknn.distance(data, focal, distance)[1:k]
    envelopeAhead <- neighbours_H_Ahead(data, names(selected), h = h)

    w <- wDistance(selected, theta = theta, typePoint = typePoint)

    pointPrediction <- t(as.matrix(w))%*%t(data$y[,envelopeAhead])
  }
  return(pointPrediction)
}

fknn.distance <- function(data, focal, distance){
  if(anyNA(data$y[,focal])){ #DU
    tcut <- sum(!is.na(data$y[,focal]))

    if(distance=='supremum'){
      dist <- matrixStats::colMaxs(abs(matrix(data$y[1:tcut,]-data$y[1:tcut,focal],
                                              ncol = ncol(data$y),
                                              nrow = length(1:tcut))))
      names(dist)<- colnames(data$y)
    }
    if(distance=='l2'){
      dist <- colSums((data$y[1:tcut,]-data$y[1:tcut,focal])^2)
      names(dist)<-colnames(data$y)
    }
  }else{
    if(distance=='supremum'){
      dist <- matrixStats::colMaxs(abs(matrix(data$y-data$y[,focal],
                                              ncol = ncol(data$y),
                                              nrow = length(data$x))))
      names(dist)<- colnames(data$y)
    }
    if(distance=='l2'){
      dist <- colSums((data$y-data$y[,focal])^2)
      names(dist)<-colnames(data$y)
    }
  }

  dist <- dist[-which(names(dist) == focal)]
  return(sort(dist))
}
