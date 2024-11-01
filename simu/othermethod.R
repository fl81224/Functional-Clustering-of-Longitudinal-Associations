flexmix_simu = function(Rseed, noise, N, numbers, p, depTX, depFX = 0) {
  inputlist = data_generation(
    Rseed = Rseed,
    noise = noise,
    N = N,
    numbers = numbers,
    p = p,
    time_reso = 10,
    nbeta_basis = 4,
    lambda_RP = c(0.5),
    depTX = depTX
  )
  
  inputlist = inputlist[[1]]
  
  Y = inputlist$Y
  X = inputlist$X_i
  Ts = inputlist$Ts
  p = inputlist$p
  N = inputlist$N
  
  Y_vec = unlist(Y)
  X_mat = do.call(rbind, X)
  
  fit0 = cv.glmnet(x = X_mat, y = Y_vec)
  ind = which.min(fit0$cvm)
  
  vs = which(fit0$glmnet.fit$beta[, ind] != 0)
  
  X_mat = matrix(X_mat[, vs], , length(vs))
  
  data = cbind(Y_vec, X_mat)
  data = data.frame(data)
  
  
  names(data)[1] = "Y"
  for (i in 2:(length(vs) + 1))
    names(data)[i] = paste0("x", i - 1)
  data$ind = as.factor(rep(1:N, each = length(Ts)))
  weights = matrix(0, N * 10, 3)
  weights[1:600, 1] = weights[601:1200, 2] = weights[1200:1800, 3] = 1
  
  
  form = paste0("Y~")
  for (i in 1:(length(vs) - 1))
    form = paste0(form, "x", i, "+")
  form = paste0(form, "x", length(vs), "|ind")
  form = as.formula(form)
  
  if (length(vs) == 1) {
    form = "Y~x1|ind"
    form = as.formula(form)
  }
  
  
  fit = flexmix(formula = form,
                data = data,
                k = 3)#,cluster=weights)
  Z = apply(fit@posterior$scaled, 1, which.max)
  Z = Z[1 + 10 * (0:(N - 1))]
  
  coef = matrix(0, 3, p)
  for (k in 1:3) {
    comp = as.numeric(names(sort(table(Z[(1 + (k - 1) * 60):(k * 60)]))[length(table(Z[(1 +
                                                                                          (k - 1) * 60):(k * 60)]))]))
    coef[k, vs] = fit@components[[comp]][[1]]@parameters$coef[-1]
  }
  
  return(list(vs = vs, Z = as.vector(Z), coef = coef))
}
