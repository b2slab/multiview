###################################################################################

# Preprocess of a distance matrix for mvMDS (square + double center).
#
# @param values The distance matrix (not a "dist") to be preprocessed
# @return The preprocessed matrix (same dimensions)
.preprocess.mvmds <- function(values)
{
  m <- values^2

  row.means <- apply(m, 1, mean)
  col.means <- apply(m, 2, mean)
  all.means <- mean(m)


  # TODO: optimize this code
  for(i in 1:nrow(m))
    for(j in 1:nrow(m))
      m[i,j] <- m[i,j] + all.means - row.means[i] - col.means[j]

  return(m)
}

#' Multiview MDS on a list of matrices or \code{dist} objects.
#'
#' Multiview multidimensional scaling (\code{mvmds}) receives two or more
#' feature matrices or \code{dist} objects
#' (or any combination of both) and produces a low-dimensional representation
#' of the samples according to the combined information in all the input data.
#'
#' @param x A list of data matrices or dist objects. Both types can be mixed.
#'   In the case of plain data matrices, euclidean distance will be used to
#'   generate a \code{dist} object for that data view.
#' @param k Number of desired dimensions of the low-dimensional projection.
#' @return A \code{n x k} matrix with the k-dimensional projection, where \code{n}
#'    is the number of samples in the dataset.
#'
#' @note All input views must have the same number of samples (rows).
#'
#' @examples
#' x1 <- iris[, 1:2]
#' x2 <- iris[, 3:4]
#' mvmds(list(x1, x2), k = 2)
#'
#' @export
mvmds <- function(x, k = 2)
{
  # Set up a n x n x nviews matrix for the preprocessed results
  # Compute dist of feature matrices.
  # Convert "dist" to matrix.
  # Preprocess the matrices

  num.obs <- nrow(as.matrix(x[[1]]))
  my.mat  <- array(dim = c(num.obs, num.obs, length(x)))
  for(i in 1:length(x))
  {
    if(class(x[[i]]) != "dist")
      my.view <- as.matrix(stats::dist(x[[i]]))
    else
      my.view <- as.matrix(x[[i]])

    my.view2 <- .preprocess.mvmds(my.view)

    my.mat[, , i] <- - my.view2 / 2
  }

  # Apply CPC to the 3D matrix of preprocessed data views
  common <- .cpc(my.mat, k)

  # Return first k columns
  return(common$evectors[,1:k])
}
