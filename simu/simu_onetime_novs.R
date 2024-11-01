simu_onetime_novs = function(task,
                             Rseed,
                             noise,
                             N,
                             numbers,
                             p,
                             samediff = "same",
                             inflation,
                             save = TRUE,
                             depTX,
                             depFX = 0,
                             max_outer_iter) {
  inputlist = data_generation(
    Rseed = Rseed,
    noise = noise,
    N = N,
    numbers = numbers,
    p = p,
    time_reso = 10,
    nbeta_basis = 10,
    lambda_RP = c(0.999),
    depTX = depTX
  )
  weights_init = matrix(0, N, 3)
  weights_init[1:(N / 3), 1] = weights_init[(N / 3 + 1):(N / 3 * 2), 2] =
    weights_init[(N / 3 * 2 + 1):N, 3] = 1
  
  
  res = EM_novs(
    weights_init = weights_init,
    inputlist = inputlist[[1]],
    Rseed = Rseed,
    K = 3,
    task = Rseed,
    numbers = numbers,
    plot_M_step = FALSE,
    max_outer_iter = 30,
    inflation = 2,
    save = FALSE
  )
  return(res)
}
