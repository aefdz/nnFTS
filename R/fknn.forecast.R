#' Functional knn forecasting
#'
#' @param data matrix PXN being N the total number of functions and
#' P the total number of observed values for each function. The columns are
#' ordered by time from 0 to T.
#' @param h forecasting horizon.
#' @param distance vector of distances $D(i,N)$, i.e. distance between
#' the last observed curve and the others.
#' @param point type of point estimate: "w" for weighted mean and "exp" for exponentially weighted mean.
#' @param theta parameter for the "exp" point type. The default value is 1.
#'
#' @return
#' @export
#'
#' @examples
fknn.forecast <- function(data, h, distance, point, theta = NULL){

}
