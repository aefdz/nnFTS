#' extension
#'
#' This function provides a prediction band.
#' @param data An object of class fds. data$y is a matrix of dimension \code{c(p,n)},
#' being \code{c(p)} the number of observed points and \code{c(n)} the number of functions.
#' data$x is the grid where the functions are recorded with length equal to \code{c(p)}.
#' @param focal Column name of the function to envelope.
#' @param tcut Order of the last observed point of the partially observed curve.
#' It is a number between 1 and \code{c(p)}.
#' @param dist Distance between functions ("l2" or "supremum"). Defaults is L2 distance.
#' @keywords prediction
#' @export
#' @examples
#' extension()


extension <- function(data, focal, tcut, dist) {
  datayshort <- data$y[1:tcut, ]
  dataxshort <- data$x[1:tcut]
  rownames(datayshort) <- dataxshort
  dataShort <- fds(dataxshort, datayshort, xname = 'x', yname = 'y')
  sco <-
    PCAproj(t(dataShort$y), center = median)$scores
  rownames(sco) <-  colnames(dataShort$y)

  SubsampleAll <- envelope(dataShort, focal, dist, plot = FALSE)
  Subsample <- SubsampleAll$subsample
  if (length(Subsample) - 1 < 2) {
    #if the focal is a missing case
    dimJ <- 'No Subsample'
    width <- 'No Subsample'
    depth <- c('No Subsample', 'No Subsample', 'No Subsample')
    stats <- c(dimJ, width, depth)

    output <-
      list(
        Jordered = 'No Subsample',
        bag = cbind('No Subsample', 'No Subsample'),
        point = 'No Subsample',
        stat = stats
      )
    return(output)

    Jordered <- orderByMBD(dataShort, Subsample) #Ordered with focal
  }

  depth <-
    matrix(data = 0,
           nrow = 1,
           ncol = 4)
  colnames(depth) <-
    c('TukeyObserved',
      'MBDObserved',
      'TukeyPredicted',
      'MBDPredicted')

  dimJ <- length(Subsample) - 1

  #Now we have the curves alphaJ that build the bag
  #Order them by tukey and by MBD
  orderMBDObserved <- orderByMBD(dataShort, Subsample)
  depth[1, 1] <- depthPercentile(names(orderMBDObserved$x), focal)

  orderMBDPredicted <-
    orderByMBD(extract(data, direction = "x", xorder = data$x[(tcut + 1):length(data$x)]), Subsample)
  depth[1, 2] <- depthPercentile(names(orderMBDPredicted$x), focal)

  Jordered <- orderMBDObserved
  #Now compute the statistics
  bagId <- names(Jordered$x)[!names(Jordered$x) %in% focal]
  bag <- band(data, bagId)
  low <- bag[, 1]
  high <- bag[, 2]
  bagTotal <-
    band(data, colnames(data$y))
  lowTotal <- bagTotal[, 1]
  highTotal <- bagTotal[, 2]

  coverageObserved <-
    coverageStat(data$y[, focal][1:tcut], high[1:tcut], low[1:tcut]) / possibleCoverage(fds(data$x[1:tcut], data$y[1:tcut, ]), focal)
  coveragePredicted <-
    coverageStat(data$y[, focal][(tcut + 1):length(data$x)], high[(tcut + 1):length(data$x)], low[(tcut +
                                                                                                     1):length(data$x)]) / possibleCoverage(fds(data$x[(tcut + 1):length(data$x)], data$y[(tcut +
                                                                                                                                                                                             1):length(data$x), ]), focal)

  widthObserved <-
    widthStat(high[1:tcut], low[1:tcut]) / widthStat(highTotal[1:tcut], lowTotal[1:tcut])
  widthPredicted <-
    widthStat(high[(tcut + 1):length(data$x)], low[(tcut + 1):length(data$x)]) /
    widthStat(highTotal[(tcut + 1):length(data$x)], lowTotal[(tcut + 1):length(data$x)])

  stats <-
    cbind(dimJ,
          widthObserved,
          widthPredicted,
          depth,
          coverageObserved,
          coveragePredicted)
  rownames(stats) <- focal
  output <-
    list(Jordered = names(Jordered$x)[!names(Jordered$x) %in% focal], stat =
           stats)
  return(output)
}
