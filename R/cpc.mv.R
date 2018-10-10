#require(RSpectra)

# Common principal components of a set of matrices.
#
# \code{cpc} uses a variation of Trendafilov (2010) method to compute
# the k first common principal components of a set of matrices in an
# efficient way.
#
# @param x A set of \code{n} matrices of dimension \code{pxp} given as a \code{p x p x n} matrix.
# @param k Number of components to extract (0 means all \code{p} components).
# @return A list with two elements: \code{evalues} (the eigenvalues) and \evalues{evectors} (the common eigenvectors).
.cpc <- function(x, k = 0)
{
  # Adapt parameter names to those used by Trendafilov on his code
  n_g <- rep(dim(x)[1], dim(x)[3])

  #p=size(x,1); mcas=size(x,3);
  p <- dim(x)[1]
  # If k is 0 then retrieve all the components
  if(k == 0) k <- p
  mcas <- dim(x)[3]

  #iter=15; n = n_g./sum(n_g);
  iter <- 15
  n <- n_g/sum(n_g)

  #D=zeros(p,mcas); CPC=zeros(p); Qw=eye(p);
  D <- array(0, dim=c(k, mcas))
  CPC <- array(0, dim=c(p, k))
  Qw <- diag(1, p)

  #s=zeros(p);
  #for m=1:mcas;
  #  s = s + double(n(m))*x(:,:,m);
  #end
  s <- array(0, dim=c(p,p))
  for(m in 1:mcas) {
    s <- s + n[m]*x[, , m]
  }

  #[q0,d0]=eig(s);
  #if d0(1,1)<d0(p,p)
  #  q0=q0(:,p:-1:1);
  #end
  #res <- eigen(s)  ######
  res <- RSpectra::eigs_sym(s, k)
  q0 <- res$vectors
  #print(dim(q0))
  #print(res$values)
  # d0 <- diag(res$values, p)
  # if(d0[1, 1] < d0[p, p]) {
  #   q0 <- q0[, ncol(q0):1]
  # }

  # loop 'for ncomp=1:p'
  # EDIT: extract only the first k components
  for(ncomp in 1:k) {
    #q = q0(:,ncomp);
    #d = zeros(1,mcas);
    #for m=1:mcas
    #    d(m) = q'*x(:,:,m)*q;
    #end
    #print(ncomp)
    q <- q0[, ncomp]
    d <- array(0, dim=c(1, mcas))
    for(m in 1:mcas) {
      d[, m] <- t(q) %*% x[, , m] %*% q
    }

    # loop 'for i=1:iter'
    for(i in 1:iter) {
      #s=zeros(p);
      #for m=1:mcas
      #  s = s + double(n_g(m))*x(:,:,m)/d(m);
      #end

      s <- array(0, dim=c(p, p))
      for(m in 1:mcas) {
        s <- s + n_g[m] * x[, , m] / d[, m]
      }

      #w = s*q;
      #if ncomp~=1;
      #  w = Qw*w;
      #end;

      w <- s %*% q
      if( ncomp != 1) {
        w <- Qw %*% w
      }

      #q = w/((w'*w)^.5);
      #for m=1:mcas
      #  d(m) = q'*x(:,:,m)*q;
      #end

      q <- w / as.numeric(sqrt((t(w) %*% w)))
      for(m in 1:mcas) {
        d[, m]  <- t(q) %*% x[, , m] %*% q
      }

    }
    # end of loop 'for i=1:iter'

    #D(ncomp,:) = d; CPC(:,ncomp) = q;
    #Qw=Qw-q*q';

    D[ncomp, ] <- d
    CPC[, ncomp] <- q
    Qw <- Qw - q %*% t(q)
  }
  # end of loop 'for ncomp=1:k'

  return(list(evalues=D, evectors=CPC));
}
