.pad_vec <- function(v, p) {
  if (length(v) >= p) v[seq_len(p)] else c(v, rep(0, p - length(v)))
}

Cas2_deriv <- function(nc = c(100, 100, 100),
                       len_t = 60,
                       delta = c(1, 1, 1),
                       p_z = 10L,
                       coef_mu_base = list(
                         # Classe 1
                         c(  1.0, -0.8, -1.2,  0.3, rep(0, 6)),
                         # Classe 2
                         c( -1.0,  0.8,  1.2, -0.3, rep(0, 6)),
                         # Classe 3
                         c(  0.0,  0.0,  0.0,  1.2, rep(0, 6))
                       ),
                       lambda_list_10 = list(
                         0.5^((1:10) - 1),
                         0.5^((1:10) - 1),
                         0.5^((1:10) - 1)
                       ),
                       m_rho_base = list(
                         c(-0.8,  0.0, -0.6,  0.0, rep(0, 6)),  # classe 1
                         c( 0.8,  0.0,  0.6,  0.4, rep(0, 6)),  # classe 2
                         c( 0.0,  0.0,  0.0, -0.4, rep(0, 6))   # classe 3
                       ),
                       muY_base = NULL,
                       sigma2 = 0.2,
                       tau2   = 0.2) {

  p_z <- as.integer(p_z)
  stopifnot(p_z >= 2L, is.finite(p_z))

  if (is.null(muY_base)) {
    muY_base <- list(
      .pad_vec(rep(0.8, 10), p_z),
      .pad_vec(rep(-0.8, 10), p_z),
      .pad_vec(rep(0.0, 10), p_z)
    )
  }

  stopifnot(length(nc) == 3,
            length(lambda_list_10) == 3,
            length(coef_mu_base) == 3,
            length(m_rho_base) == 3,
            length(muY_base) == 3,
            length(delta) == 3)
  
  D <- 10
  stopifnot(all(vapply(lambda_list_10, length, integer(1)) == D))
  stopifnot(all(vapply(coef_mu_base,   length, integer(1)) == D))
  stopifnot(all(vapply(m_rho_base,     length, integer(1)) == D))
  stopifnot(all(vapply(muY_base,     length, integer(1)) == p_z))

  t_grid <- seq(0, 1, length.out = len_t)
  
  Psi <- matrix(NA, nrow = len_t, ncol = D)
  for (j in 1:5) {
    Psi[, 2*j - 1] <- sin(2 * j * pi * t_grid)
    Psi[, 2*j]     <- cos(2 * j * pi * t_grid)
  }
  
  idx_curve   <- 1:2
  idx_deriv   <- 3:4
  idx_neutral <- 5:10
  
  scale_coef <- function(a, delta1, delta2) {
    a1 <- a
    a1[idx_curve]   <- delta1 * a[idx_curve]
    a1[idx_deriv]   <- delta2 * a[idx_deriv]
    a1[idx_neutral] <- a[idx_neutral]
    a1
  }
  
  coef_mu <- lapply(coef_mu_base, function(a) scale_coef(a, delta[1], delta[2]))
  m_rho   <- lapply(m_rho_base,   function(a) scale_coef(a, delta[1], delta[2]))
  
  mu_fun_list <- lapply(coef_mu, function(a) as.vector(Psi %*% a))
  
  sim_one_m_rho <- function(nc_c, lambda_vec, sigma2, mu_rho, mu_fun) {
    Z   <- matrix(rnorm(nc_c * D), nc_c, D)
    Rho <- sweep(Z, 2, sqrt(lambda_vec), "*")
    Rho <- sweep(Rho, 2, mu_rho, "+")
    X_rand <- Rho %*% t(Psi)
    X_true <- sweep(X_rand, 2, mu_fun, "+")
    X_obs  <- X_true + matrix(rnorm(nc_c * length(mu_fun), sd = sqrt(sigma2)),
                              nc_c, length(mu_fun))
    list(X_true = X_true, X_obs = X_obs, Rho = Rho)
  }
  
  sim_list <- lapply(1:3, function(c)
    sim_one_m_rho(nc[c],
                  lambda_list_10[[c]],
                  sigma2,
                  m_rho[[c]],
                  mu_fun_list[[c]])
  )
  
  X_true <- do.call(rbind, lapply(sim_list, `[[`, "X_true"))
  X_obs  <- do.call(rbind, lapply(sim_list, `[[`, "X_obs"))
  Rho    <- do.call(rbind, lapply(sim_list, `[[`, "Rho"))
  class  <- rep(1:3, times = nc)
  
  p <- p_z
  R <- matrix(0.2, nrow = p, ncol = p); diag(R) <- 1
  Theta <- eigen(R, symmetric = TRUE)$vectors[, 1:D, drop = FALSE]
  
  muY_list <- lapply(muY_base, function(m) delta[3] * m)
  
  muY_mat <- do.call(rbind, lapply(1:3, function(k)
    matrix(muY_list[[k]], nrow = nc[k], ncol = p, byrow = TRUE)
  ))
  
  Y_signal <- Rho %*% t(Theta) + muY_mat
  Y <- Y_signal + matrix(rnorm(sum(nc) * p, sd = sqrt(tau2)), sum(nc), p)
  
  return(list(
    t_grid      = t_grid,
    class       = class,
    Psi         = Psi,
    coef_mu     = coef_mu,
    mu_fun_list = mu_fun_list,
    Rho         = Rho,
    X_true      = X_true,
    X_obs       = X_obs,
    Theta       = Theta,
    Y           = Y,
    sigma2      = sigma2,
    tau2        = tau2,
    delta       = delta,
    idx_curve   = idx_curve,
    idx_deriv   = idx_deriv,
    p_z         = p_z
  ))
}



