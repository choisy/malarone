# This version has a no patient effect.

fit_oocyst_model_no_patient <- function(
    counts, sizes,
    control = "H00",
    n_quad = 25,
    d50_init = 25,
    b_init = 3
) {
  
  # 1. Prepare data -----------------------------------
  
  stopifnot(
    all(c("time_point", "count") %in% names(counts)),
    all(c("patient_id", "mosquito_id",
          "time_point", "area") %in% names(sizes))
  )
  
  time_points <- c(
    control,
    sort(setdiff(unique(counts$time_point), control))
  )
  
  G <- length(time_points)
  n_beta <- G - 1
  
  stopifnot(
    G >= 2,
    control %in% counts$time_point,
    setequal(unique(sizes$time_point), time_points),
    !anyNA(counts),
    !anyNA(sizes),
    all(is.finite(counts$count)),
    all(counts$count >= 0),
    all(counts$count == floor(counts$count)),
    all(is.finite(sizes$area)),
    all(sizes$area > 0)
  )
  
  counts <- counts |>
    mutate(g = match(time_point, time_points))
  
  sizes <- sizes |>
    mutate(
      g = match(time_point, time_points),
      log_area = log(area)
    )
  
  # Group by patient_id as well to avoid merging
  # mosquitoes with identical IDs across patients
  size_data <- sizes |>
    group_by(patient_id, g, mosquito_id) |>
    summarise(
      n = n(),
      sum_x = sum(log_area),
      sum_x2 = sum(log_area^2),
      .groups = "drop"
    )
  
  count_data <- counts |>
    count(g, count, name = "frequency")
  
  # 2. Gaussian quadrature ----------------------------
  
  gh <- gauss.quad.prob(n_quad, dist = "normal")
  
  z <- gh$nodes
  w <- gh$weights
  log_w <- log(w)
  
  logsumexp <- function(x) {
    m <- max(x)
    if (!is.finite(m)) return(-Inf)
    m + log(sum(exp(x - m)))
  }
  
  # 3. Parameters -------------------------------------
  
  # Parameter order:
  # lambda, k, a,
  # beta_2 ... beta_G,
  # sigma, tau, d50, b
  
  unpack <- function(par) {
    
    beta <- c(0, par[4:(3 + n_beta)])
    
    j <- 4 + n_beta
    
    list(
      lambda = par[1],
      k = par[2],
      a = par[3],
      beta = beta,
      sigma = par[j],
      tau = par[j + 1],
      d50 = par[j + 2],
      b = par[j + 3]
    )
  }
  
  # 4. Negative log composite likelihood -------------
  
  nll <- function(par) {
    
    p <- unpack(par)
    
    if (!all(is.finite(unlist(p)))) return(1e100)
    
    if (p$lambda <= 0 || p$k <= 0 ||
        p$sigma <= 0 || p$tau <= 0 ||
        p$d50 <= 0 || p$b <= 0) {
      return(1e100)
    }
    
    total <- 0
    
    # Mean detection probability for each
    # mosquito-specific mean log-area
    get_q <- function(mu) {
      
      log_x <- outer(
        mu,
        p$sigma * z,
        "+"
      )
      
      detection <- plogis(
        p$b * (log_x - log(p$d50))
      )
      
      pmax(
        as.vector(detection %*% w),
        .Machine$double.xmin
      )
    }
    
    # A. Marginal visible oocyst counts ---------------
    
    for (g in seq_len(G)) {
      
      mu <- p$a + p$beta[g] + p$tau * z
      
      q <- get_q(mu)
      
      mean_visible <- p$lambda * q
      
      dat <- count_data |>
        filter(.data$g == !!g)
      
      for (r in seq_len(nrow(dat))) {
        
        log_prob <- logsumexp(
          log_w +
            dnbinom(
              dat$count[r],
              size = p$k,
              mu = mean_visible,
              log = TRUE
            )
        )
        
        total <- total +
          dat$frequency[r] * log_prob
      }
    }
    
    # B. Measured oocyst areas ------------------------
    
    log_det <- plogis(
      p$b * (sizes$log_area - log(p$d50)),
      log.p = TRUE
    )
    
    total <- total + sum(log_det)
    
    for (i in seq_len(nrow(size_data))) {
      
      dat <- size_data[i, ]
      
      g <- dat$g
      
      mu <- p$a + p$beta[g] + p$tau * z
      
      q <- get_q(mu)
      
      ss <- dat$sum_x2 -
        2 * mu * dat$sum_x +
        dat$n * mu^2
      
      log_density <-
        -dat$n * log(p$sigma) -
        dat$n * log(2 * pi) / 2 -
        ss / (2 * p$sigma^2) -
        dat$n * log(q)
      
      # Probability mosquito has >= 1 visible oocyst
      mean_visible <- p$lambda * q
      
      log_prob_zero <- -p$k *
        log1p(mean_visible / p$k)
      
      log_selection <- log(
        pmax(
          -expm1(log_prob_zero),
          .Machine$double.xmin
        )
      )
      
      log_numerator <- logsumexp(
        log_w + log_selection + log_density
      )
      
      log_denominator <- logsumexp(
        log_w + log_selection
      )
      
      total <- total +
        log_numerator - log_denominator
    }
    
    if (!is.finite(total)) return(1e100)
    
    -total
  }
  
  # 5. Initial values ---------------------------------
  
  x <- sizes$log_area
  
  sigma0 <- sqrt(log(1 + 0.25^2))
  
  tau0 <- sqrt(max(
    var(x) - sigma0^2,
    0.01^2
  ))
  
  # Control time-point counts
  y0 <- counts$count[
    counts$time_point == control
  ]
  
  # Initial detection probability = 0.8
  lambda0 <- max(mean(y0) / 0.8, 0.1)
  
  # Negative binomial dispersion
  m <- mean(y0)
  v <- var(y0)
  
  k0 <- if (is.finite(v) && v > m && m > 0) {
    m^2 / (v - m)
  } else {
    100
  }
  
  k0 <- pmax(0.01, pmin(k0, 1e3))
  
  # Initial time-point effects
  beta0 <- sapply(time_points[-1], function(g) {
    
    xg <- sizes$log_area[sizes$time_point == g]
    xc <- sizes$log_area[sizes$time_point == control]
    
    if (length(xg) > 0 && length(xc) > 0) {
      mean(xg) - mean(xc)
    } else {
      0
    }
  })
  
  beta0 <- pmax(
    log(0.01),
    pmin(beta0, log(2))
  )
  
  initial <- c(
    lambda0,
    k0,
    mean(x),
    beta0,
    sigma0,
    tau0,
    d50_init,
    b_init
  )
  
  # 6. Bounds -----------------------------------------
  
  lower <- c(
    0.01,                  # lambda
    0.01,                  # k
    -Inf,                  # a
    rep(log(0.01), G - 1), # beta
    0.01,                  # sigma
    0.01,                  # tau
    min(sizes$area) / 100, # d50
    0.1                    # b
  )
  
  upper <- c(
    1e5,                   # lambda
    1e3,                   # k
    Inf,                   # a
    rep(log(2), G - 1),    # beta
    3,                     # sigma
    3,                     # tau
    max(sizes$area) * 10,  # d50
    100                    # b
  )
  
  if (any(!is.finite(initial)) ||
      any(initial < lower | initial > upper)) {
    stop("Initial values must lie within bounds.")
  }
  
  # 7. Optimization -----------------------------------
  
  fit <- optim(
    par = initial,
    fn = nll,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(
      maxit = 2000,
      factr = 1e7
    )
  )
  
  p <- unpack(fit$par)
  
  # 8. Results ----------------------------------------
  
  effects <- data.frame(
    time_point = time_points,
    beta = p$beta,
    size_ratio = exp(p$beta),
    size_reduction = 1 - exp(p$beta),
    reduction_percent = 100 * (1 - exp(p$beta))
  )
  
  list(
    effects = effects,
    d50 = p$d50,
    b = p$b,
    sigma = p$sigma,
    tau = p$tau,
    k = p$k,
    lambda = p$lambda,
    logLik = -fit$value,
    convergence = fit$convergence,
    optim = fit
  )
}


#######################################################################################

# Takes 7"  
system.time(
  fit_no_patient_effect <- fit_oocyst_model_no_patient(counts = oocyst_count2,
                                                       sizes = oocyst_size2)
)

saveRDScache(fit_no_patient_effect, "fit_no_patient_effect.rds")

#######################################################################################

fit <- fit_no_patient_effect

#library(ggplot2)

