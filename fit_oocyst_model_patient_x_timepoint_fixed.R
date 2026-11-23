# This version has a patient × time-point fixed-effect interaction.

fit_oocyst_model <- function(
    counts, sizes,
    control = "H00",
    n_quad = 25,
    d50_init = 25,
    b_init = 3
) {
  
  # 1. Prepare data -----------------------------------
  
  stopifnot(
    all(c("patient_id", "time_point", "count") %in% names(counts)),
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
    P >= 2,
    control %in% counts$time_point,
    setequal(unique(sizes$time_point), time_points),
    all(sizes$patient_id %in% patients),
    !anyNA(counts),
    !anyNA(sizes),
    all(is.finite(counts$count)),
    all(counts$count >= 0),
    all(counts$count == floor(counts$count)),
    all(is.finite(sizes$area)),
    all(sizes$area > 0)
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
  
  # Check complete patient x time-point design
  combinations <- expand.grid(
    p = seq_len(P),
    g = seq_len(G)
  )
  
  observed <- counts |>
    distinct(p, g)
  
  if (nrow(anti_join(combinations, observed,
                     by = c("p", "g"))) > 0) {
    stop("Every patient must have counts at every time point.")
  }
  
  # Mosquito-level size summaries
  size_data <- sizes |>
    group_by(p, g, mosquito_id) |>
    summarise(
      n = n(),
      sum_x = sum(log_area),
      sum_x2 = sum(log_area^2),
      .groups = "drop"
    )
  
  # Count frequencies
  count_data <- counts |>
    count(p, g, count, name = "frequency")
  
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
  # lambda_1 ... lambda_P
  # k
  # alpha_(p,g), P*G values
  # sigma, tau, d50, b
  
  unpack <- function(par) {
    
    j <- 1
    
    lambda <- par[j:(j + P - 1)]
    j <- j + P
    
    k <- par[j]
    j <- j + 1
    
    alpha <- matrix(
      par[j:(j + P * G - 1)],
      nrow = P,
      ncol = G
    )
    
    j <- j + P * G
    
    list(
      lambda = lambda,
      k = k,
      alpha = alpha,
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
    
    if (any(p$lambda <= 0) ||
        p$k <= 0 ||
        p$sigma <= 0 ||
        p$tau <= 0 ||
        p$d50 <= 0 ||
        p$b <= 0) {
      return(1e100)
    }
    
    get_q <- function(mu) {
      
      log_x <- outer(mu, p$sigma * z, "+")
      
      detection <- plogis(
        p$b * (log_x - log(p$d50))
      )
      
      pmax(
        as.vector(detection %*% w),
        .Machine$double.xmin
      )
    }
    
    total <- 0
    
    # A. Marginal visible counts ----------------------
    
    groups <- unique(count_data[c("p", "g")])
    
    for (h in seq_len(nrow(groups))) {
      
      patient <- groups$p[h]
      g <- groups$g[h]
      
      mu <- p$alpha[patient, g] + p$tau * z
      
      q <- get_q(mu)
      
      mean_visible <- p$lambda[patient] * q
      
      dat <- count_data[
        count_data$p == patient &
          count_data$g == g, ,
        drop = FALSE
      ]
      
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
      
      patient <- dat$p
      g <- dat$g
      
      mu <- p$alpha[patient, g] + p$tau * z
      
      q <- get_q(mu)
      
      ss <- dat$sum_x2 -
        2 * mu * dat$sum_x +
        dat$n * mu^2
      
      log_density <-
        -dat$n * log(p$sigma) -
        dat$n * log(2 * pi) / 2 -
        ss / (2 * p$sigma^2) -
        dat$n * log(q)
      
      mean_visible <- p$lambda[patient] * q
      
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
  
  lambda0 <- sapply(patients, function(patient) {
    
    y <- counts$count[
      counts$patient_id == patient &
        counts$time_point == control
    ]
    
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
  
  # Initial alpha for each patient x time point
  alpha0 <- matrix(
    mean(x),
    nrow = P,
    ncol = G
  )
  
  for (patient in seq_len(P)) {
    for (g in seq_len(G)) {
      
      xx <- sizes$log_area[
        sizes$p == patient & sizes$g == g
      ]
      
      if (length(xx) > 0) {
        alpha0[patient, g] <- mean(xx)
      }
    }
  }
  
  initial <- c(
    lambda0,
    k0,
    as.vector(alpha0),
    sigma0,
    tau0,
    d50_init,
    b_init
  )
  
  # 6. Bounds -----------------------------------------
  
  lower <- c(
    rep(0.01, P),
    0.01,
    rep(-Inf, P * G),
    0.01,
    0.01,
    min(sizes$area) / 100,
    0.1
  )
  
  upper <- c(
    rep(1e5, P),
    1e3,
    rep(Inf, P * G),
    3,
    3,
    max(sizes$area) * 10,
    100
  )
  
  if (any(!is.finite(initial)) ||
      any(initial < lower | initial > upper)) {
    stop("Initial values must lie within optimization bounds.")
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
  
  patient_effects <- expand.grid(
    patient_id = patients,
    time_point = time_points
  )
  
  patient_effects$alpha <- as.vector(p$alpha)
  
  patient_effects <- patient_effects |>
    group_by(patient_id) |>
    mutate(
      beta = alpha - alpha[time_point == control],
      size_ratio = exp(beta),
      size_reduction = 1 - size_ratio,
      reduction_percent = 100 * size_reduction
    ) |>
    ungroup()
  
  patient_effects$lambda <- p$lambda[
    match(patient_effects$patient_id, patients)
  ]
  
  # Average time-point effects across patients
  effects <- patient_effects |>
    group_by(time_point) |>
    summarise(
      mean_size_ratio = mean(size_ratio),
      mean_reduction_percent = 100 *
        (1 - mean(size_ratio)),
      .groups = "drop"
    )
  
  list(
    effects = effects,
    patient_effects = patient_effects,
    d50 = p$d50,
    b = p$b,
    sigma = p$sigma,
    tau = p$tau,
    k = p$k,
    lambda = setNames(p$lambda, patients),
    logLik = -fit$value,
    convergence = fit$convergence,
    optim = fit
  )
}

################################

# This takes 3h:36':40"
system.time(
  fit_patient_x_timepoint_fixed <- fit_oocyst_model(counts = oocyst_count2,
                                                    sizes = oocyst_size2)
)

saveRDScache(fit_patient_x_timepoint_fixed, "fit_patient_x_timepoint_fixed.rds")

fit_patient_x_timepoint_fixed <- readRDScache("fit_patient_x_timepoint_fixed.rds")
