# This version has a patient × time-point random-effect interaction while time_point is
# fixed effect.
# I reduced the default quadrature nodes to 15 and 9, respectively. Once you obtain a
# stable fit, I recommend repeating it with n_quad = 25 and n_quad_patient = 15 to
# check numerical accuracy.

fit_oocyst_model_patient_x_timepoint_random <- function(
    counts, sizes,
    control = "H00",
    n_quad = 15,
    n_quad_patient = 9,
    d50_init = 25,
    b_init = 3
) {
  
  # 1. Prepare data -----------------------------------
  
  stopifnot(
    all(c("patient_id", "time_point", "count") %in%
          names(counts)),
    all(c("patient_id", "mosquito_id",
          "time_point", "area") %in% names(sizes))
  )
  
  time_points <- c(
    control,
    sort(setdiff(unique(counts$time_point), control))
  )
  
  patients <- sort(unique(counts$patient_id))
  
  G <- length(time_points)
  P <- length(patients)
  
  stopifnot(
    G >= 2,
    control %in% counts$time_point,
    setequal(unique(sizes$time_point), time_points),
    all(sizes$patient_id %in% patients),
    !anyNA(counts),
    !anyNA(sizes),
    all(is.finite(counts$count)),
    all(counts$count >= 0),
    all(counts$count == floor(counts$count)),
    all(is.finite(sizes$area)),
    all(sizes$area > 0),
    n_quad >= 3,
    n_quad_patient >= 3
  )
  
  counts <- counts |>
    mutate(
      g = match(time_point, time_points),
      p = match(patient_id, patients)
    )
  
  sizes <- sizes |>
    mutate(
      g = match(time_point, time_points),
      p = match(patient_id, patients),
      log_area = log(area)
    )
  
  size_data <- sizes |>
    group_by(p, g, mosquito_id) |>
    summarise(
      n = n(),
      sum_x = sum(log_area),
      sum_x2 = sum(log_area^2),
      .groups = "drop"
    )
  
  count_data <- counts |>
    count(p, g, count, name = "frequency")
  
  # Precompute data for each patient/time-point group
  groups <- counts |>
    distinct(p, g) |>
    arrange(p, g)
  
  group_data <- lapply(seq_len(nrow(groups)), function(h) {
    
    pp <- groups$p[h]
    gg <- groups$g[h]
    
    cd <- count_data |>
      filter(p == pp, g == gg)
    
    sd <- size_data |>
      filter(p == pp, g == gg)
    
    observed_log_areas <- sizes$log_area[
      sizes$p == pp & sizes$g == gg
    ]
    
    list(
      p = pp,
      g = gg,
      y = cd$count,
      frequency = cd$frequency,
      n = sd$n,
      sum_x = sd$sum_x,
      sum_x2 = sd$sum_x2,
      observed_log_areas = observed_log_areas
    )
  })
  
  # 2. Gaussian quadrature ----------------------------
  
  gh_m <- gauss.quad.prob(
    n_quad, dist = "normal"
  )
  
  gh_v <- gauss.quad.prob(
    n_quad_patient, dist = "normal"
  )
  
  z <- gh_m$nodes
  w <- gh_m$weights
  log_w <- log(w)
  
  zv <- gh_v$nodes
  log_wv <- log(gh_v$weights)
  
  logsumexp <- function(x) {
    m <- max(x)
    if (!is.finite(m)) return(-Inf)
    m + log(sum(exp(x - m)))
  }
  
  # 3. Parameters -------------------------------------
  
  # lambda_1 ... lambda_P
  # k, a
  # beta_2 ... beta_G
  # sigma, tau, omega, d50, b
  
  unpack <- function(par) {
    
    j <- 1
    
    lambda <- par[j:(j + P - 1)]
    j <- j + P
    
    k <- par[j]
    j <- j + 1
    
    a <- par[j]
    j <- j + 1
    
    beta <- c(0, par[j:(j + G - 2)])
    j <- j + G - 1
    
    list(
      lambda = lambda,
      k = k,
      a = a,
      beta = beta,
      sigma = par[j],
      tau = par[j + 1],
      omega = par[j + 2],
      d50 = par[j + 3],
      b = par[j + 4]
    )
  }
  
  # 4. Negative log composite likelihood -------------
  
  nll <- function(par) {
    
    p <- unpack(par)
    
    if (!all(is.finite(unlist(p)))) return(1e100)
    
    if (any(p$lambda <= 0) ||
        p$k <= 0 ||
        p$sigma <= 0 ||
        p$tau <= 0 ||
        p$omega < 0 ||
        p$d50 <= 0 ||
        p$b <= 0) {
      return(1e100)
    }
    
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
    
    total <- 0
    
    # Each group shares one patient x time-point
    # random effect, integrated over below.
    
    for (dat in group_data) {
      
      patient <- dat$p
      g <- dat$g
      
      # Detection contribution for measured areas
      log_det <- sum(plogis(
        p$b * (dat$observed_log_areas -
                 log(p$d50)),
        log.p = TRUE
      ))
      
      log_lik_v <- numeric(length(zv))
      
      for (h in seq_along(zv)) {
        
        # Patient x time-point random effect
        v <- p$omega * zv[h]
        
        # Mosquito-specific mean log-area
        mu <- p$a + p$beta[g] + v +
          p$tau * z
        
        q <- get_q(mu)
        
        mean_visible <- p$lambda[patient] * q
        
        # A. Marginal count likelihood
        log_count <- 0
        
        for (r in seq_along(dat$y)) {
          
          log_prob <- logsumexp(
            log_w +
              dnbinom(
                dat$y[r],
                size = p$k,
                mu = mean_visible,
                log = TRUE
              )
          )
          
          log_count <- log_count +
            dat$frequency[r] * log_prob
        }
        
        # B. Conditional size likelihood
        log_size <- log_det
        
        if (length(dat$n) > 0) {
          
          # Probability of a positive mosquito
          log_prob_zero <- -p$k *
            log1p(mean_visible / p$k)
          
          log_selection <- log(pmax(
            -expm1(log_prob_zero),
            .Machine$double.xmin
          ))
          
          log_denominator <- logsumexp(
            log_w + log_selection
          )
          
          # Matrix: mosquito quadrature nodes x
          # mosquitoes with measured areas
          
          ss <- outer(mu^2, dat$n, "*") -
            2 * outer(mu, dat$sum_x, "*") +
            matrix(
              dat$sum_x2,
              nrow = length(mu),
              ncol = length(dat$n),
              byrow = TRUE
            )
          
          log_density <-
            -ss / (2 * p$sigma^2) -
            matrix(
              dat$n * (
                log(p$sigma) +
                  0.5 * log(2 * pi)
              ),
              nrow = length(mu),
              ncol = length(dat$n),
              byrow = TRUE
            ) -
            outer(log(q), dat$n, "*")
          
          log_integrand <- sweep(
            log_density,
            1,
            log_w + log_selection,
            "+"
          )
          
          log_numerators <- apply(
            log_integrand,
            2,
            logsumexp
          )
          
          log_size <- log_size +
            sum(log_numerators - log_denominator)
        }
        
        log_lik_v[h] <- log_count + log_size
      }
      
      # Integrate over the patient x time-point
      # random effect
      total <- total + logsumexp(
        log_wv + log_lik_v
      )
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
  
  lambda0 <- sapply(patients, function(patient) {
    
    y <- counts$count[
      counts$patient_id == patient &
        counts$time_point == control
    ]
    
    if (length(y) == 0) {
      y <- counts$count[
        counts$patient_id == patient
      ]
    }
    
    max(mean(y) / 0.8, 0.1)
  })
  
  lambda0 <- pmin(lambda0, 1e5)
  
  m <- mean(counts$count)
  v <- var(counts$count)
  
  k0 <- if (is.finite(v) && v > m && m > 0) {
    m^2 / (v - m)
  } else {
    100
  }
  
  k0 <- pmax(0.01, pmin(k0, 1e3))
  
  # Initial fixed time-point effects
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
    0.1,          # omega
    d50_init,
    b_init
  )
  
  # 6. Bounds -----------------------------------------
  
  lower <- c(
    rep(0.01, P),
    0.01,                  # k
    -Inf,                  # a
    rep(log(0.01), G - 1), # beta
    0.01,                  # sigma
    0.01,                  # tau
    0,                     # omega
    min(sizes$area) / 100, # d50
    0.1                    # b
  )
  
  upper <- c(
    rep(1e5, P),
    1e3,                   # k
    Inf,                   # a
    rep(log(2), G - 1),    # beta
    3,                     # sigma
    3,                     # tau
    3,                     # omega
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
    omega = p$omega,
    k = p$k,
    lambda = setNames(p$lambda, patients),
    logLik = -fit$value,
    convergence = fit$convergence,
    optim = fit
  )
}

#######################################################################################

# Takes 2h:4':15"

system.time(
  fit_patient_x_timepoint_random <- fit_oocyst_model_patient_x_timepoint_random(
    counts = oocyst_count2,
    sizes = oocyst_size2)
)

saveRDScache(fit_patient_x_timepoint_random, "fit_patient_x_timepoint_random.rds")
fit_patient_x_timepoint_random <- readRDScache("fit_patient_x_timepoint_random.rds")
