#' envelope
#'
#' This function allows to get the envelopes.
#' @param data Do you love cats? Defaults to TRUE.
#' @param focal Column name of the function to envelope.
#' @param dist Distance between functions ("l2" or "supremum"). Defaults is L2 distance.
#' @param plot If plotting or not the iterations. Defaults to TRUE.
#' @keywords envelope
#' @export
#' @examples
#' envelope()

envelope <- function(data,
                         focal,
                         dist = "l2",
                         plot = TRUE) {
  P <-
    data$x[data$y[, focal] != rowMaxs(data$y) &
             data$y[, focal] != rowMins(data$y)]

  length_P <- length(P)

  if (length_P > 0) {
    Subsample <- vector(mode = 'numeric', length = 0)

    I <- colnames(data$y)[!colnames(data$y) %in% focal]

    iter <- c(1)

    iterDepth <- c(0)

    while (length(I) > 1) {
      x <- as.character(P)

      # First step

      candidates <- I

      if (dist == 'supremum') {
        supDist <-
          colMaxs(abs(matrix(
            data$y[, I] - data$y[, focal],
            ncol = length(I),
            nrow = length(data$x)
          )))
        names(supDist) <- c(I)
      }

      if (dist == 'l2') {
        supDist <-
          colSums((data$y[, I] - data$y[, focal]) ^ 2)
        names(supDist) <- c(I)
      }

      d <- sort(supDist, decreasing = FALSE)

      SubsampleIter <- names(d[1])

      candidates = candidates[!candidates %in% names(d[1])]

      while (length(x) != 0 & 0 < length(candidates)) {
        dAUX <- sort(d[candidates], decreasing = FALSE)

        Ji <-
          which(abs(rowSums(sign(
            as.matrix(data$y[, c(SubsampleIter, names(dAUX[1]))] - data$y[, focal])
          ))) < length(c(SubsampleIter, names(d[1]))))

        # Do the first two envelope the focal?
        if (length(x) <= length(x[!x %in% names(Ji)])) {
          #No
          candidates = candidates[!candidates %in% names(dAUX[1])]
        } else{
          # si Envuelven
          x <- x[!x %in% names(Ji)] #remaining points to cover
          SubsampleIter <- union(SubsampleIter, names(dAUX[1]))
          candidates <- candidates[!candidates %in% SubsampleIter]
          I = I[!I %in% SubsampleIter]
        }
      }

      iter <- c(iter, iter[length(iter)] + 1)
      aux <- fMBD(data$y[, c(focal, Subsample, SubsampleIter)])

      iterDepthAux <- depthPercentile(names(aux$x), focal)
      iterDepth <- c(iterDepth, iterDepthAux)
      if (max(iterDepth[1:(length(iter) - 1)]) <= iterDepth[iter[length(iter)]]) {
        Subsample <- c(Subsample, SubsampleIter)
        if (plot == TRUE) {
          #print(c('X for covering', length(x)))
          #print(c('Remaining Sample',length(I)))
          #print(c('In this envelope',length(SubsampleIter)))
          #print(iterDepthAux)
          fdata <- melt(data$y, id = 'x')
          aux <-
            matrix(c(data$y[, SubsampleIter]), nrow = length(data$x))
          low <-   vector(mode = "numeric", length = length(data$x))
          high <-   vector(mode = "numeric", length = length(data$x))
          for (j in c(1:length(data$x))) {
            low[j] <- min(aux[j, ])
            high[j] <- max(aux[j, ])
          }
          dataBand <-
            data.frame(
              X1 = c(data$x, rev(data$x)),
              X2 = 'NULL',
              value = c(high, rev(low))
            )
          pl1  <-
            ggplot(data = fdata,
                   ggplot2::aes(
                     x = X1,
                     y = value,
                     colour = as.factor(X2),
                     group = as.factor(X2)
                   )) +
            geom_line(color = 'grey50') + theme(legend.position = "none") +
            geom_line(data = fdata[fdata$X2 %in% Subsample, ],
                      color = 'black',
                      cex = 1.25) +
            geom_line(data = fdata[fdata[, 2] == focal, ],
                      color = 'red',
                      cex = 1.25) +
            geom_polygon(data = dataBand,
                         ggplot2::aes(x = X1, y = value),
                         color = 'grey25') +
            geom_line(data = fdata[fdata$X2 %in% Subsample, ],
                      color = 'black',
                      cex = 1) +
            geom_line(data = fdata[fdata$X2 %in% SubsampleIter, ],
                      color = 'grey',
                      cex = 1) +
            geom_line(data = fdata[fdata$X2 == focal, ],
                      color = 'red',
                      cex = 1.2)
          print(pl1)

          readline(prompt = "Press <Enter> to continue...")
        }

      }
      if (length(candidates) == 0) {
        break
      }
    }
  } else{
    Subsample = c('No Subsample')
  }
  if (dist == 'supremum') {
    supDist <-
      colMaxs(abs(matrix(
        data$y - data$y[, focal],
        ncol = length(I),
        nrow = length(data$x)
      )))
    names(supDist) <- colnames(data$y)
  }
  if (dist == 'l2') {
    supDist  <-
      colSums((data$y - data$y[, focal]) ^ 2)
    names(supDist) <- colnames(data$y)
  }

  if (Subsample[1] != 'No Subsample') {
    auxDist  <-  sort(supDist, decreasing = FALSE, index.return = TRUE)
    auxDepth <- fMBD(data$y[, c(focal, Subsample)])

    return(list(subsample = names(auxDepth$x)))
  } else{
    return(list(subsample = 'No Subsample'))
  }
}
