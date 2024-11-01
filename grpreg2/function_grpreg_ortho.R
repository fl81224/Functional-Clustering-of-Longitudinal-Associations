orthogonalize <- function(X, group) {
  n <- nrow(X)
  J <- max(group)
  T <- vector("list", J)
  XX <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  XX[, which(group == 0)] <- X[, which(group == 0)]
  for (j in seq_along(integer(J))) {
    ind <- which(group == j)
    if (length(ind) == 0)
      next
    SVD <- svd(X[, ind, drop = FALSE], nu = 0)
    r <- which(SVD$d > 1e-10)
    T[[j]] <-
      sweep(SVD$v[, r, drop = FALSE], 2, sqrt(n) / SVD$d[r], "*")
    XX[, ind[r]] <- X[, ind] %*% T[[j]]
  }
  nz <- !apply(XX == 0, 2, all)
  XX <- XX[, nz, drop = FALSE]
  attr(XX, "T") <- T
  attr(XX, "group") <- group[nz]
  XX
}
unorthogonalize <- function(b, XX, group, intercept = TRUE) {
  ind <- !sapply(attr(XX, "T"), is.null)
  T <- bdiag(attr(XX, "T")[ind])
  if (intercept) {
    ind0 <- c(1, 1 + which(group == 0))
    val <-
      Matrix::as.matrix(rbind(b[ind0, , drop = FALSE], T %*% b[-ind0, , drop =
                                                                 FALSE]))
  } else if (sum(group == 0)) {
    ind0 <- which(group == 0)
    val <-
      Matrix::as.matrix(rbind(b[ind0, , drop = FALSE], T %*% b[-ind0, , drop =
                                                                 FALSE]))
  } else {
    val <- as.matrix(T %*% b)
  }
}
##############################################################################################
setupLambda <-
  function(X,
           y,
           group,
           family,
           penalty,
           alpha,
           lambda.min,
           log.lambda,
           nlambda,
           group.multiplier) {
    # Fit to unpenalized covariates
    n <- length(y)
    K <- table(group)
    K1 <- if (min(group) == 0)
      cumsum(K)
    else
      c(0, cumsum(K))
    storage.mode(K1) <- "integer"
    if (K1[1] != 0) {
      fit <- glm(y ~ X[, group == 0], family = family)
    } else {
      fit <- glm(y ~ 1, family = family)
    }
    
    ## Determine lambda.max
    if (family == "gaussian") {
      r <- fit$residuals
    } else {
      w <- fit$weights
      if (max(w) < 1e-4)
        stop("Unpenalized portion of model is already saturated; exiting...",
             call. = FALSE)
      r <- residuals(fit, "working") * w
    }
    if (strtrim(penalty, 2) == "gr") {
      zmax <- .Call("maxgrad", X, r, K1, as.double(group.multiplier)) / n
    } else {
      zmax <- .Call("maxprod", X, r, K1, as.double(group.multiplier)) / n
    }
    lambda.max <- zmax / alpha
    
    if (log.lambda) {
      # lambda sequence on log-scale
      if (lambda.min == 0) {
        lambda <-
          c(exp(seq(
            log(lambda.max), log(.001 * lambda.max), length = nlambda - 1
          )), 0)
      } else {
        lambda <-
          exp(seq(log(lambda.max), log(lambda.min * lambda.max), length = nlambda))
      }
    } else {
      # lambda sequence on linear-scale
      if (lambda.min == 0) {
        lambda <-
          c(seq(lambda.max, 0.001 * lambda.max, length = nlambda - 1), 0)
      } else {
        lambda <- seq(lambda.max, lambda.min * lambda.max, length = nlambda)
      }
    }
    lambda
  }


unstandardize <- function(b, XG) {
  beta <- matrix(0, nrow = 1 + length(XG$scale), ncol = ncol(b))
  beta[1 + XG$nz, ] <- b[-1, ] / XG$scale[XG$nz]
  beta[1, ] <- b[1, ] - crossprod(XG$center, beta[-1, , drop = FALSE])
  beta
}
