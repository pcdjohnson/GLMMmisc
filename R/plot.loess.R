#' Scatter plot with LOESS line
#'
#' This function makes an x-y scatter plot and adds a LOESS line
#'
#' @param x Either a data frame with two columns, which will be taken as x and y, or x coordinates
#' @param y If x is not a data frame, y coordinates
#' @export
#' @examples
#' plot.loess(cars)

plot.loess <-
  function(x, y = NULL, loess.col = "red", ...)
  {
    if(is.data.frame(x) && ncol(x) == 2) {
      plot.data <- x
      xy.labels <- names(plot.data)
      names(plot.data) <- c("x", "y")
    } else {
      plot.data <- data.frame(x = x, y = y)
      xy.labels <- names(plot.data)
    }
    plot.data$loess.line <- predict(loess(y ~ x, data = plot.data))
    plot.data <- plot.data[order(plot.data$x), ]
    plot(plot.data[,c("x", "y")], xlab = xy.labels[1], ylab = xy.labels[2], ...)
    abline(h = 0)
    points(plot.data[,c("x", "loess.line")], type = "l", col = loess.col)
  }
