# Laplacian following Ng et al formulation
.LaplacianNg <- function(mat)
{
  D <- rowSums(mat)
  sqriD <- diag(1/sqrt(D))
  return(sqriD %*% mat %*% sqriD)
}

# Gives the suggested sigma to spectral cluster given a distance matrix
# neighbour parameter:
# NULL -> select avg distance according to Luxburg criterium (avg distance to
# the log(number of samples)th neighbour);
# Integer > 1: avg distance to the neighbour-th closest sample
# Real number 0<neighbour<=1: avg distance to the (neighbour-th*number of samples)
# closest sample
.suggestedSigma <- function(distances, neighbour=NULL)
{
  if(is.null(neighbour))
    n <- ceiling(log(nrow(distances)))
  else #if(neighbour > 1)
    n <- neighbour
  # else
  #   n <- ceiling(neighbour*nrow(distances))

  dist.ord <- apply(distances, 2, sort)
  # Compute the mean removing NAs and infinite values just in case
  # If it is 0 then return 1 (or we will get errors)
  result <- mean(dist.ord[n,is.finite(dist.ord[n,])], na.rm=TRUE)
  if(result == 0.0)
    result <- 1.0
  return(result)
}

# Given a distance matrix, compute the gaussian similarity of its values
.distanceGaussianSimilarity <- function(distances, sigma)
{
  myfactor <- -1/(2*sigma^2)
  result   <- exp(distances^2*myfactor)

  # 0's are dangerous afterwards, they should be replaced by sth safer
  result[result == 0] <- 1e-16

  return(result)
}



#' Multiview spectral clustering on a list of matrices or \code{dist} objects.
#'
#' Computes the multiview spectral clustering of data on a list of matrices
#' or dist objects (or a mix of both), supposed to be different views of the same data.
#'
#' @param x A list of feature matrices or \code{dist} objects (or a mix of both).
#' @param k Number of desired clusters.
#' @param sigmas Either \code{NULL}, a single real value or a vector of real values.
#'               They correspond
#'               to the \code{sigma} parameter in the Gaussian radial basis function.
#'               If it is NULL then the default sigma computation is used (average
#'               distance to the log(n)-th neighbour, with n = number of samples), unless
#'               \code{neighbours} has a value different from \code{NULL}.
#'               If it is a single number then the same sigma is applied to all input
#'               views. If it is a vector each value in it is applied to the corresponding
#'               input view.
#' @param neighbours Either \code{NULL}, a single integer value or a vector of integers.
#'                   They correspond
#'                   to the expected number of neighbours per point, used to estimate the
#'                   \code{sigma} values of the Gaussian radial basis function.
#'                   If it is NULL then the default sigma computation is used (average
#'                   distance to the log(n)-th neighbour, with n = number of samples).
#'                   If it is a single value then the same number of neighbours is used
#'                   on all input views, else each value in the vector is applied to the
#'                   corresponding input view.
#'                   Does not have effect if \code{sigma} is different from \code{NULL}.
#' @param clustering Tells \code{mvsc} if it has to perform the clustering on the projection
#'                   (if TRUE) or to skip the clustering step of the algorithm.
#'
#' @return A list with four elements: \code{clustering} is a vector of integers with the
#'   clustering assignment of each sample (not included if \code{clustering = FALSE}),
#'   \code{evalues} is a matrix with the eigenvalues
#'   of the common principal components (CPC) step, \code{evectors} is a matrix with the
#'   eigenvectors of the CPC step, and \code{sigmas} is a vector with the sigmas used on
#'   the Gaussian radial basis function of each input view.
#'
#' @note All input views must have the same number of samples (rows).
#'
#' @examples
#' m1 <- iris[, 1:2]
#' m2 <- iris[, 3:4]
#' mvsc(list(m1, m2), k = 3)
#'
#' @export
mvsc <- function(x, k, sigmas = NULL, neighbours = NULL, clustering = TRUE)
{
  # Prepare sigmas and neighbours parameters: if sigmas/neighbours have only
  # 1 element, repeat them nviews times (length(NULL) == 0)
  nviews <- length(x)
  if(length(sigmas)     == 1) sigmas     <- rep(sigmas,     nviews)
  if(length(neighbours) == 1) neighbours <- rep(neighbours, nviews)

  # Placeholder to store the actual sigmas used
  mysigmas <- rep(0, nviews)

  # Compute the joint diagonal matrix of the similarity matrices
  # First we have to create a p x p x n array with the laplacian matrices
  # p = number of samples, n = number of views

  # as.matrix is required as the first view can be a "dist"
  numPoints <- nrow(as.matrix(x[[1]]))
  lapMatrix <- array(dim=c(numPoints, numPoints, nviews))

  for(i in 1:length(x))
  {
      # Compute distance if necessary, then the kernel
      if(class(x[[i]]) != "dist")
        view.dist <- as.matrix(stats::dist(x[[i]]))
      else
        view.dist <- as.matrix(x[[i]])

      if(!is.null(sigmas)) {
        mysigmas[i] <- sigmas[i]
      } else if(!is.null(neighbours)) {
        mysigmas[i] <- .suggestedSigma(view.dist, neighbours[i])
      } else {
        mysigmas[i] <- .suggestedSigma(view.dist, NULL)
      }
      view.grbf        <- .distanceGaussianSimilarity(view.dist, mysigmas[i])
      #diag(view.grbf)  <- 0
      view.lsym        <- .LaplacianNg(view.grbf)

      lapMatrix[, , i] <- view.lsym
  }

  # Now we have to compute CPC on the laplacians to get the eigen values and vectors
  cpc.result <- .cpc(lapMatrix, k)

  # Normalize evectors by row (so each has unit length)
  cpc.result$evectors <- t(apply(cpc.result$evectors, 1, function(r) r/sqrt(sum(r^2))))

  # Run KMeans on the eigenvectors (only first K columns are computed) and return everything
  kmeans.clust <- stats::kmeans(cpc.result$evectors[,1:k], k, iter.max = 10 + round(log(numPoints)), nstart=10)$cluster


  # Projection is the projection of the samples received on the CPC space
  return(list(clustering = kmeans.clust, evalues = cpc.result$evalues,
              evectors   = cpc.result$evectors, sigmas = mysigmas))

}
