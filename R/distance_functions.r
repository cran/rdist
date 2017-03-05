available_metrics <- c("euclidean", "minkowski", "manhattan", 
                       "chebyshev", "maximum", "canberra", 
                       "angular", "correlation", "absolute_correlation", 
                       "hamming", "jaccard")

#' @rdname rdist
#' @export
rdist <- function(X, 
                  metric = "euclidean", 
                  p = 2){
  # make sure input is well-defined
  metric <- match.arg(metric, available_metrics)
  X <- as.matrix(X)
  # use metric
  ans <- switch(metric, 
                "euclidean" = euclidean_rdist(X), 
                "minkowski" = minkowski_rdist(X, p = p), 
                "manhattan" = manhattan_rdist(X), 
                "chebyshev" = maximum_rdist(X), 
                "maximum" = maximum_rdist(X), 
                "canberra" = canberra_rdist(X), 
                "angular" = angular_rdist(X), 
                "correlation" = correlation_rdist(X), 
                "absolute_correlation" = absolute_correlation_rdist(X), 
                "hamming" = hamming_rdist(X), 
                "jaccard" = jaccard_rdist(X))
  # change attributes
  attributes(ans) <- NULL
  attr(ans, "Size") <- nrow(X)
  attr(ans, "call") <- match.call()
  attr(ans, "method") <- metric
  class(ans) <- "dist"
  return(ans)
}

#' @rdname rdist
#' @export
pdist <- function(X, 
                  metric = "euclidean", 
                  p = 2){
  # make sure input is well-defined
  metric <- match.arg(metric, available_metrics)
  X <- as.matrix(X)
  # use metric
  switch(metric, 
         "euclidean" = euclidean_pdist(X), 
         "minkowski" = minkowski_pdist(X, p = p), 
         "manhattan" = manhattan_pdist(X), 
         "chebyshev" = maximum_pdist(X), 
         "maximum" = maximum_pdist(X), 
         "canberra" = canberra_pdist(X), 
         "angular" = angular_pdist(X), 
         "correlation" = correlation_pdist(X), 
         "absolute_correlation" = absolute_correlation_pdist(X), 
         "hamming" = hamming_pdist(X),
         "jaccard" = jaccard_pdist(X))
}

#' @rdname rdist
#' @export
cdist <- function(X, Y, 
                  metric = "euclidean",
                  p = 2){
  # make sure input is well-defined
  metric <- match.arg(metric, available_metrics)
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  stopifnot(ncol(X) == ncol(Y))
  # use metric
  switch(metric, 
         "euclidean" = euclidean_cdist(X, Y), 
         "minkowski" = minkowski_cdist(X, Y, p = p), 
         "manhattan" = manhattan_cdist(X, Y),
         "chebyshev" = maximum_cdist(X, Y),
         "maximum" = maximum_cdist(X, Y), 
         "canberra" = canberra_cdist(X, Y), 
         "angular" = angular_cdist(X, Y), 
         "correlation" = correlation_cdist(X, Y), 
         "absolute_correlation" = absolute_correlation_cdist(X, Y), 
         "hamming" = hamming_cdist(X, Y),
         "jaccard" = jaccard_cdist(X, Y))
}