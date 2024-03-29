#' rdist: an R package for distances
#'
#' \code{rdist} provide a common framework to calculate distances. There are three main functions: 
#' \itemize{
#' \item \code{rdist} computes the pairwise distances between observations in one matrix and returns a \code{dist} object,
#' \item \code{pdist} computes the pairwise distances between observations in one matrix and returns a \code{matrix}, and
#' \item \code{cdist} computes the distances between observations in two matrices and returns a \code{matrix}. 
#' }
#' In particular the \code{cdist} function is often missing in other distance functions. All 
#' calculations involving \code{NA} values will consistently return \code{NA}. 
#' 
#' @details Available distance measures are (written for two vectors v and w):
#' \itemize{
#' \item \code{"euclidean"}: \eqn{\sqrt{\sum_i(v_i - w_i)^2}}{sqrt(sum_i((v_i - w_i)^2))}
#' \item \code{"minkowski"}: \eqn{(\sum_i|v_i - w_i|^p)^{1/p}}{(sum_i(|v_i - w_i|^p))^{1/p}}
#' \item \code{"manhattan"}: \eqn{\sum_i(|v_i-w_i|)}{sum_i(|v_i-w_i|)}
#' \item \code{"maximum"} or \code{"chebyshev"}: \eqn{\max_i(|v_i-w_i|)}{max_i(|v_i-w_i|)}
#' \item \code{"canberra"}: \eqn{\sum_i(\frac{|v_i-w_i|}{|v_i|+|w_i|})}{sum_i(|v_i-w_i|/(|v_i|+|w_i|))}
#' \item \code{"angular"}: \eqn{\cos^{-1}(cor(v, w))}{arccos(cor(v, w))}
#' \item \code{"correlation"}: \eqn{\sqrt{\frac{1-cor(v, w)}{2}}}{sqrt((1-cor(v, w))/2)}
#' \item \code{"absolute_correlation"}: \eqn{\sqrt{1-|cor(v, w)|^2}}{sqrt((1-|cor(v, w)|^2))}
#' \item \code{"hamming"}: \eqn{(\sum_i v_i \neq w_i) / \sum_i 1}{sum_i(v_i != w_i)/sum_i(1)}
#' \item \code{"jaccard"}: \eqn{(\sum_i v_i \neq w_i) / \sum_i 1_{v_i \neq 0 \cup w_i \neq 0}}{sum_i(v_i != w_i)/sum_i(v_i != 0 or w_i != 0)}
#' \item Any function that defines a distance between two vectors. 
#' }
#' @param X,Y A matrix
#' @param metric The distance metric to use
#' @param p The power of the Minkowski distance
#' @name rdist
#' @docType package
#' @useDynLib rdist, .registration = TRUE 
#' @importFrom Rcpp sourceCpp
#' @importFrom stats as.dist cor
NULL

