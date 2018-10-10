# MV-tSNE with log-linear opinion pooling, using GAs to optimize the weights of each view (opinion)
# The function is split in two parts, one that computes P and other that performs tsne on a given P
# plus an "interface" function that does everything



#' Multiview tSNE using an expert opinion pooling on the input probability matrices
#'
#' Given a list of of input views and other parameters, \code{mvtsne} computes a neighbouring probability matrix
#' for each input view, then finds the optimal set of weights to combine these matrices using a log-linear pool,
#' and applies the pooled probability matrix as input to the standard tSNE procedure, where the probability matrix
#' of the output space is adjusted to the pooled probability matrix using Kullback-Liebler divergence.
#'
#' @param X A list of R matrices or "dist" objects, where each matrix/dist is one of the views of the dataset.
#' @param k The desired dimension of the resulting embedding.
#' @param initial_dims Number of dimensions to use in the reduction method.
#' @param perplexity This perplexity parameter is roughly equivalent to the optimal number of neighbours.
#' @param max_iter Maximum number of iterations to perform.
#' @param min_cost The minimum cost value (error) to stop iterations.
#' @param epoch_callback A callback function to be called after each epoch (which is a number of iterations controlled
#'                       parameter \code{epoch}, see next).
#' @param whiten A boolean value indicating if the data matrices should be whitened.
#' @param epoch The number of iterations between update messages.
#'
#' @return A list with two elements: \code{embedding} with the k-dimensional embedding of the input samples, and
#'         \code{weights} with the weights associated to each input data view.
#'
#' @note All input views must have the same number of samples (rows).
#'
#' @examples
#' m1 <- iris[, 1:2]
#' m2 <- iris[, 3:4]
#' mvtsne(list(m1, m2), k = 2)
#'
#' @export
mvtsne <- function(X, k = 2, initial_dims = 30, perplexity = 30, max_iter = 1000, min_cost = 0, epoch_callback = NULL, whiten = TRUE, epoch = 100)
{
    pooling   <- find.pooling(X, initial_dims=initial_dims, perplexity=perplexity, whiten=whiten)

    embedding <- tsne.p(pooling$P, k = k, initial_dims, max_iter, min_cost, epoch_callback, epoch)

    return(list(embedding = embedding, weights = pooling$weights))
}



# Compute optimal pooling of P's for MV-tSNE
# X-> list of matrices/dist
# method -> best, log, linear (best computes both and selects the best according to Abbas 2009)
# returns the weights and the pooled probability matrix
find.pooling <- function(X, initial_dims=30, perplexity=30, whiten=TRUE, method='best')
{
  nviews <- length(X)
  eps = 2^(-52) # typical machine precision

  # Preprocess input matrices
  for(i in 1:nviews)
  {
    if ("dist" %in% class(X[[i]])) {
      n = attr(X[[i]], "Size")
    }
    else {
      X[[i]] = as.matrix(X[[i]])
      X[[i]] = X[[i]] - min(X[[i]])
      X[[i]] = X[[i]]/max(X[[i]])
      initial_dims = min(initial_dims, ncol(X[[i]]))
      if (whiten)
        X[[i]] <- .whiten(as.matrix(X[[i]]), n.comp = initial_dims)
      n = nrow(X[[i]])
    }
  }


  # Compute the probability distribution of each input view
  # Now P is an array of n x n x nviews (instead of n x n)
  P = array(0, dim=c(n, n, nviews))
  for(i in 1:nviews)
  {
    P[, , i] = .x2p(X[[i]], perplexity, 1e-05)$P
    P[, , i] = 0.5 * (P[, , i] + t(P[, , i]))
    P[, , i][P[, , i] < eps] <- eps
    P[, , i] = P[, , i]/sum(P[, , i])
  }


  ###################################################################################
  # Find optimum weight combination

  # Computes log.linear pooled opinion, returns pooled matrix (norm) + reg constant (1/sum)
  log.linear.pooling <- function(P, weights)
  {
    P.exp.w <- sweep(P, 3, weights, FUN='^')
    pooled  <- apply(P.exp.w, c(1,2), prod)
    reg.const <- 1/sum(pooled)
    pooled <- pooled * reg.const

    result <- {}
    result$pooled <- pooled
    result$reg.const <- reg.const
    return(result)
  }

  # General objective function for log-linear pooling (Abbas 2009 (9))
  objective.log.linear <- function(weights)
  {
    # Compute log-linear pooled prob with given weights
    pooling <- log.linear.pooling(P, weights)

    # Compute log-linear payoff (Abbas (9)) (here higher is worse)
    kls <- apply(P, 3, function(q.i) entropy::KL.plugin(pooling$pooled, q.i))
    payoff <- sum(kls*weights)

    # Introduce constraint sum(weights)=1 through a penalty
    penalty <- abs(1 - sum(weights))
    goal    <- payoff + penalty
    return(-goal)
  }

  # Gradient of the objective
  gradient.log.linear <- function(weights)
  {
    pooling <- log.linear.pooling(P, weights)

    log.pooling <- log(pooling$pooled)

    res <- sapply(1:nviews, function(v) sum(weights[v] * log.pooling/P[, , v]))

    return(res)
  }

  # Run the optimizer to find the best set of weights. Start point at (1/v, ...)
  # Using bounds to limit the weights to 0..1
  opt <- stats::optim(runif(nviews), objective.log.linear, gradient.log.linear, method='L-BFGS-B',
                      lower = rep(0, nviews), upper = rep(1, nviews))

  pooled.log  <- log.linear.pooling(P, opt$par)
  reg.const   <- pooled.log$reg.const
  pooled.log  <- pooled.log$pooled

  ###################################################################################


  # Compute KL between the different P's
  # Yes I compute every combination as it is not symmetrical
  #message('KL on Ps')
  kl.ps <- matrix(0, nrow=nviews, ncol=nviews)
  for(i in 1:nviews)
    for(j in 1:nviews)
      kl.ps[i, j] <- sum(P[,,i] * log((P[,,i]+eps)/(P[,,j]+eps)))

  # Compute KL between each input P and the consensus one
  kl.central <- sapply(1:nviews, function(i) sum(P[,,i] * log((P[,,i]+eps)/(pooled.log+eps))))

  # Return info
  result            <- {}
  result$P          <- pooled.log
  result$weights    <- opt$par
  result$kl.interps <- kl.ps
  result$kl.central <- kl.central
  return(result)
}


tsne.p <- function(P, k=2, initial_dims=30, max_iter = 1000, min_cost=0, epoch_callback=NULL, epoch=100 )
{
  n <- nrow(P)

	momentum = .5
	final_momentum = .8
	mom_switch_iter = 250

	epsilon = 500
	min_gain = .01
	initial_P_gain = 4

	eps = 2^(-52) # typical machine precision

	ydata = matrix(rnorm(k * n),n)

	# Prepare for standard tSNE
	P <- P * initial_P_gain

	grads =  matrix(0,nrow(ydata),ncol(ydata))
	incs =  matrix(0,nrow(ydata),ncol(ydata))
	gains = matrix(1,nrow(ydata),ncol(ydata))

	for (iter in 1:max_iter){
		if (iter %% epoch == 0) { # epoch
			cost =  sum(apply(P * log((P+eps)/(Q+eps)),1,sum))
			#message("Epoch: Iteration #",iter," error is: ",cost)

			if (cost < min_cost) break
			if (!is.null(epoch_callback)) epoch_callback(ydata)

		}

		sum_ydata = apply(ydata^2, 1, sum)
		num =  1/(1 + sum_ydata +    sweep(-2 * ydata %*% t(ydata),2, -t(sum_ydata)))
		diag(num)=0
		Q = num / sum(num)
		if (any(is.nan(num))) message ('NaN in grad. descent')
		Q[Q < eps] = eps
		stiffnesses = 4 * (P-Q) * num
		for (i in 1:n){
			grads[i,] = apply(sweep(-ydata, 2, -ydata[i,]) * stiffnesses[,i],2,sum)
		}

		gains = ((gains + .2) * abs(sign(grads) != sign(incs)) +
				 gains * .8 * abs(sign(grads) == sign(incs)))

		gains[gains < min_gain] = min_gain

		incs = momentum * incs - epsilon * (gains * grads)
		ydata = ydata + incs
		ydata = sweep(ydata,2,apply(ydata,2,mean))
		if (iter == mom_switch_iter) momentum = final_momentum

		if (iter == 100) P = P/4

	}
	ydata
}

