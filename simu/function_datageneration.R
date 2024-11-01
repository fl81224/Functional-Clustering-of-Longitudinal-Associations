beta_generation = function(beta_reso, print = F) {
  c = 3
  Ts = seq(0, 1, , beta_reso)
  betatrue = list()
  for (i in 1:3)
    betatrue[[i]] = matrix(0, beta_reso, 3)
  
  
  
  betatrue[[1]][, 2] = sin(0.5 * pi * Ts + 1.5 * pi) + (Ts)
  betatrue[[1]][, 1] = sin(0.5 * pi * Ts + pi) - Ts
  betatrue[[1]][, 3] = (Ts - 0.5) ^ 2
  
  betatrue[[2]][, 2] = sin(2 * pi * Ts) - (Ts)
  betatrue[[2]][, 1] = (cos(2 * pi * Ts) - 1) ^ 2
  betatrue[[2]][, 3] = -1 * sin(2 * pi * Ts) + (-1)
  
  
  betatrue[[3]][, 1] = -1 * sin(0.5 * pi * Ts + 1.5 * pi) - (Ts - 0.5) -
    1
  betatrue[[3]][, 2] = -1 * (cos(2 * pi * Ts) - 1) ^ 2
  betatrue[[3]][, 3] = -1 * (-1 * sin(0.5 * pi * Ts + 1.5 * pi)) + 1
  
  betatrue[[1]] = apply(betatrue[[1]], 2,
                        function(x)
                          (x / sqrt(sum(x ^ 2)) / beta_reso))
  betatrue[[2]] = apply(betatrue[[2]], 2,
                        function(x)
                          (x / sqrt(sum(x ^ 2)) / beta_reso))
  betatrue[[3]] = apply(betatrue[[3]], 2,
                        function(x)
                          (x / sqrt(sum(x ^ 2)) / beta_reso))
  
  if (print == T) {
    par(mfrow = c(3, 3))
    ylim = range(as.vector(sapply(betatrue, unlist)))
    for (i in 1:3) {
      apply(betatrue[[i]], 2,
            function(x)
              plot(
                Ts,
                x,
                type = 'l',
                main = paste0("var", i),
                ylim = ylim
              ))
    }
  }
  
  return(betatrue)
}
#########################################################################################
#########################################################################################
HQ_generation = function(time_reso,
                         nbeta_basis,
                         N,
                         X_i,
                         lambda_RP = 0.5,
                         norder = 4) {
  Ts = seq(0, 1, , time_reso)
  
  p = ncol(X_i[[1]])
  beta_basis = create.bspline.basis(c(Ts[1], Ts[length(Ts)]), nbeta_basis, norder =
                                      norder)
  Phi0 = (eval.basis(Ts, beta_basis))
  
  RP_betamat = bsplinepen(
    beta_basis,
    Lfdobj = 2,
    rng = beta_basis$rangeval,
    returnMatrix = FALSE
  )
  L2_betamat = bsplinepen(
    beta_basis,
    Lfdobj = 0,
    rng = beta_basis$rangeval,
    returnMatrix = FALSE
  )
  L = chol(lambda_RP * RP_betamat + (1 - lambda_RP) * L2_betamat)
  Linv = solve(L)
  
  H = list()
  Q = list()
  for (i in 1:N) {
    H[[i]] = matrix(0, time_reso, nbeta_basis * p)
    Q[[i]] = matrix(0, time_reso, nbeta_basis * p)
    for (j in 1:p) {
      for (t in 1:time_reso) {
        H[[i]][t, (1 + nbeta_basis * (j - 1)):(j * nbeta_basis)] = X_i[[i]][t, j] *
          Phi0[t, ]
      }
      Q[[i]][, (1 + nbeta_basis * (j - 1)):(j * nbeta_basis)] =
        H[[i]][, (1 + nbeta_basis * (j - 1)):(j * nbeta_basis)] %*% Linv
    }
  }
  
  res = list(
    H = H,
    Q = Q,
    beta_basis = beta_basis,
    L = L,
    Linv = Linv
  )
  
  return(res)
}
#########################################################################################
X_generation = function(Rseed,
                        N = 200,
                        rho = 0,
                        p,
                        time_reso) {
  set.seed(Rseed)
  M = 4
  Ts = seq(1, 10, , time_reso)
  maxvar = 500
  
  psi.basis <-
    fda::create.fourier.basis(rangeval = c(1, 10),
                              nbasis = M + 1,
                              period = 10)
  psi.true <- eval.basis(Ts, psi.basis)[, 2:(M + 1)]
  
  Sigma = matrix(0, M * maxvar, M * maxvar)
  for (i in 1:M) {
    temp = matrix(0, nrow = maxvar, ncol = maxvar)
    for (j1 in 1:maxvar) {
      for (j2 in 1:maxvar) {
        temp[j1, j2] = rho ^ (abs(j1 - j2)) / (i ^ 2)
      }
    }
    Sigma[(1 + maxvar * (i - 1)):(maxvar * i), (1 + maxvar * (i - 1)):(maxvar *
                                                                         i)] = temp
  }
  X_i = list()
  
  xivec = mvrnorm(N, rep(0, maxvar * M), Sigma)
  for (i in 1:N) {
    coef = NULL
    for (m in 1:M) {
      coef = cbind(coef, xivec[i, (maxvar * (m - 1) + 1):(maxvar * m)])
    }
    X_i[[i]] = psi.true %*% t(coef[1:p, ])
    
  }
  return(X_i)
}
#########################################################################################



data_generation = function(Rseed,
                           noise = c(0.05, 0.05, 0.05),
                           N,
                           numbers,
                           p,
                           time_reso = 10,
                           nbeta_basis = 30,
                           lambda_RP,
                           depTX,
                           depFX = 0,
                           vs = NULL) {
  set.seed(seed = Rseed)
  Ts = seq(1, 10, , time_reso)
  
  
  K = 3
  Ztrue = c(c(rep(1, numbers[1]), rep(2, numbers[2]), rep(3, numbers[3])))
  
  X_i = X_generation(
    Rseed = Rseed,
    N = N,
    rho = depTX,
    p = p,
    time_reso = time_reso
  )
  
  
  betatrue = beta_generation(beta_reso = time_reso, print = T)
  Y = list()
  for (i in 1:N) {
    Y[[i]] = X_i[[i]][, 1] * betatrue[[1]][, Ztrue[i]] +
      X_i[[i]][, 2] * betatrue[[1]][, Ztrue[i]] +
      X_i[[i]][, 3] * betatrue[[2]][, Ztrue[i]] +
      X_i[[i]][, 4] * betatrue[[2]][, Ztrue[i]] +
      X_i[[i]][, 5] * betatrue[[3]][, Ztrue[i]] +
      X_i[[i]][, 6] * betatrue[[3]][, Ztrue[i]]
    
    Y[[i]] = as.vector(Y[[i]])
  }
  
  Y_nonoise = Y
  
  Ts = seq(1, 10, , time_reso)
  error_basis = create.fourier.basis(c(0, 1), nbasis = 4)
  Phierror = (eval.basis(seq(0, 1, , 10), error_basis))
  Phierror = Phierror[, -1]
  
  epsilon_var = sum(unlist(lapply(Y, function(x)
    sum(x ^ 2)))) / N / time_reso * noise
  
  epsilon_cor = matrix(0, 10, 10)
  for (i in 1:10) {
    for (j in 1:10) {
      epsilon_cor[i, j] = exp(-1 * abs(i - j))
    }
  }
  
  epsilon = mvrnorm(N, rep(0, 10), epsilon_cor)
  
  for (i in 1:N) {
    Y[[i]] = Y[[i]]    +
      epsilon[i, ] * sqrt(epsilon_var[Ztrue[i]] / 2) +
      rnorm(10, 0, sqrt(epsilon_var[Ztrue[i]] / 2))
  }
  Ydist = list()
  Ydist[[1]] = Ydist[[2]] = Ydist[[3]] = matrix(0, time_reso, N)
  for (i in 1:N) {
    for (k in 1:3) {
      Ydist[[k]][, i] = X_i[[i]][, 1] * betatrue[[1]][, k] +
        X_i[[i]][, 2] * betatrue[[2]][, k] +
        X_i[[i]][, 3] * betatrue[[3]][, k]
    }
  }
  
  AVGdist = matrix(0, K, K)
  for (k1 in 1:3)
    for(k2 in 1:3) {
      AVGdist[k1, k2] = sqrt((mean((Ydist[[k1]] - Ydist[[k2]]) ^ 2)) /
                               (epsilon_var[k1] ^
                                  2 + epsilon_var[k2] ^ 2))
    }
  AVGdist
  if (!is.null(vs)) {
    X_i = lapply(X_i, function(x)
      matrix(x[, vs], , length(vs)))
    p = length(vs)
  }
  
  
  res = list()
  ind = 1
  for (lRP in lambda_RP) {
    temp = HQ_generation(
      time_reso = time_reso,
      nbeta_basis = nbeta_basis,
      X_i = X_i,
      lambda_RP = lRP,
      N = N
    )
    H = temp$H
    Q = temp$Q
    beta_basis = temp$beta_basis
    L = temp$L
    Linv = temp$Linv
    
    res[[ind]] = list(
      N = N,
      K = K,
      numbers = numbers,
      Ts = Ts,
      Y = Y,
      Y_nonoise = Y_nonoise,
      X_i = X_i,
      Ztrue = Ztrue,
      p = p,
      betatrue = betatrue,
      H = H,
      Q = Q,
      beta_basis = beta_basis,
      L = L,
      Linv = Linv,
      nbeta_basis = nbeta_basis,
      lRP = lRP
    )
    ind = ind + 1
    
  }
  return(res)
}
