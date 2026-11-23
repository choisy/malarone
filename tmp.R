set.seed(123)

user <- Sys.getenv("USER")

if (user %in% c("MarcChoisy", "marc")) {  
  path2data <- paste0(Sys.getenv("HOME"), "/Library/CloudStorage/",
                      "OneDrive-OxfordUniversityClinicalResearchUnit/",
                      "GitHub/choisy/malarone/")
}

path2cache <- paste0(path2data, "cache/")
if (! dir.exists(path2cache)) dir.create(path2cache)

make_path <- function(x) paste0(path2cache, x)
file_exists <- function(x) file.exists(make_path(x))
readRDScache <- function(x) readr::read_rds(make_path(x))
saveRDScache <- function(object, file) readr::write_rds(object, make_path(file))
remove_cache <- function(x) {
  file.remove(grep(x, dir(path2cache, full.names = TRUE), value = TRUE))
}



lwd_val <- 2

YlOrRd11 <- hcl.colors(11, "YlOrRd", rev = TRUE)

n_workers <- 12

library(statmod)
library(stringr)
library(furrr)
library(tidyr)
library(purrr)
library(dplyr)
future::plan(future::multisession, workers = n_workers)

# library(statmod)
# library(ggplot2)

tune_reader <- function(f) function(file, ...) f(paste0(path2data, file), ...)
read_excel2 <- tune_reader(readxl::read_excel)
read_table2 <- tune_reader(readr::read_table)
read_rds2 <- tune_reader(readr::read_rds)

polygon2 <- function(x, y1, y2, ...) {
  polygon(c(x, rev(x)), c(y1, rev(y2)), border = NA, ...)
}

reading_exp_data <- function(file) {
  hours <- c(0, 5, 48, 72)
  days <- c(5, 7, 14, 28)
  hash <- c(setNames(hours, paste0("H", hours)),
            setNames(24 * days, paste0("D", days)))
  
  file |> 
    read_rds2() |> 
    filter(treatment == "AL+AP") |> 
    mutate(t = hash[stringr::str_remove(time_points, "0")],
           value = Inhibition / 100) |> 
    select(t, value)
}

stratified_bootstrap <- function(df, N) {
  future_map(1:ceiling(N), ~ df |> 
                               group_by(t) |> 
                               group_split() |> 
                               map(slice_sample, prop = 1, replace = TRUE) |> 
                               bind_rows(), .options = furrr_options(seed = TRUE))
}

plot2 <- function(...) plot(..., xlab = "days post infection",
                                 ylab = "proportion of infectious mosquitoes")

plot3 <- function(...) plot2(..., type = "l", lwd = lwd_val)

predict_glm <- function(object, newdata, ci = .95) {
  q_val <- qt(1 - (1 - ci) / 2, object$df.residual)
  f <- object$family$linkinv
  newdata <- as_tibble(newdata)
  newdata |> 
    bind_cols(head(predict(object, newdata, se.fit = TRUE), -1)) |> 
    mutate(lwr = f(fit - q_val * se.fit),
           upr = f(fit + q_val * se.fit),
           fit = f(fit)) |> 
    select(-se.fit)
}

arrows2 <- function(...) arrows(..., length = .1, angle = 90, code = 3, lwd = lwd_val)

predict2 <- function(...) predict(..., type = "response")

integrate2 <- function(...) integrate(...)$value

f4dpi <- function(t, A, k, d) A * (1 - exp(-k * t)) * exp(-d * t)

initial_values <- function(x, y) {
  df <- tibble(x, y) |> 
    arrange(x)

  df2 <- df |> 
    group_by(x) |> 
    summarise(y = mean(y))
  
  A <- df2 |> 
    pull(y) |> 
    max()
  
  xs <- df2 |> 
    filter(y == A) |> 
    pull(x)

  k <- df |> 
    filter(x < first(xs)) |> 
    lm(y ~ -1 + x, data = _) |> 
    coef() |> 
    last()
  
  d <- df |>
    filter(x > last(xs)) |> 
    filter(y > 0) |> 
    mutate(across(y, log)) |> 
    lm(y ~ x, data = _) |> 
    coef() |> 
    last()
  
  c(A = A, k = k / A, d = -d)
}

calibrate4dpi <- function(x, y) {
  optim(initial_values(x, y), function(p) sum((y - f4dpi(x, p[1], p[2], p[3]))^2))$par
}

alpha_data <- reading_exp_data("figure2a.rds")

alpha_fit4dpi <- with(alpha_data, calibrate4dpi(t, value))

make_effect <- function(A, k, d, r) {
  A <- A * exp(4 * 24 * r)
  function(t, tau) pmin(1, A * (1 - exp(-k * t)) * exp(-d * t - r * tau))
}

alpha1 <- with(as.list(alpha_fit4dpi), make_effect(A, k, d, d))

alpha2 <- with(as.list(alpha_fit4dpi), make_effect(A, k, d, d / 3))

n_tau <- 10
n_t <- 100
t <- seq(0, 350, le = n_t)
tau <- seq(0, 350, le = n_tau)
tau_val <- 4 * 24

plot3d <- function(f) {
  opar <- par(plt = c(.05, .95, .05, .95))
  pmat <- expand_grid(t, tau) |> 
    mutate(value = f(t, tau)) |> 
    pull(value) |>
    matrix(n_tau) |> 
    persp(tau, t, z = _, theta = 45, border = 4, zlim = 0:1,
          xlab = "\u03C4 days post infection",
          ylab = "t days post drug administration",
          zlab = "TBA \u03B1")
  
  trans3d(rep(tau_val, n_t), t, f(t, tau_val), pmat) |> 
    with(lines(x, y, col = "red", lwd = 3))
  par(opar)
}

plot3d(alpha2)
plot3d(alpha1)

t_max <- 700
t_val <- 0:t_max

plot2d <- function(f) {
  plot(t_val, f(t_val, tau = 0), type = "n", ylim = 0:1,
       xlab = "time post drug administration (hours)",
       ylab = bquote("reduction" ~ alpha ~ "of mosquitoes with oocysts"))
  walk(seq(0, 700, 24), ~ lines(t_val, f(t_val, tau = .x), col = 4))
  lines(t_val, f(t_val, tau = 4 * 24), col = 2, lwd = lwd_val)
}

plot2d(alpha2)
plot2d(alpha1)

boot_calibration_4dpi <- function(data, N = 1e4, eps = 1.1) {
  model <- data |> 
    stratified_bootstrap(eps * N) |> 
    future_map(~ with(.x, safely(calibrate4dpi)(t, value)))

  model[map_lgl(model, ~ is.null(.x$error))] |>
    map("result") |> 
    sample(N)
}



if (file_exists("alpha_model.rds")) {
  alpha_model_4dpi <- readRDScache("alpha_model.rds")
} else {
  alpha_model_4dpi <- boot_calibration_4dpi(alpha_data)
  saveRDScache(alpha_model_4dpi, "alpha_model.rds")
}

t_val2 <- seq(0, 1.05 * t_max)

trajectories <- function(model, xs = t_val2) {
  model |> 
    set_names() |> 
    map_dfc(~ with(as.list(.x), f4dpi(xs, A, k, d)))
}

alpha_trajectories <- trajectories(alpha_model_4dpi)

plot_data_model <- function(data, simulations, ylabel, col_data = 4, col_model = 2) {
  with(data, plot(t, value, type = "n",
                  xlim = c(0, t_max), ylim = c(min(data$value), 1),
                  xlab = "time post drug administration (hours)",
                  ylab = ylabel))
  abline(h = 0)
  
  simulations |> 
    t() |> 
    as.data.frame() |> 
    map_dfr(quantile, c(.025, .5, .975)) |> 
    with({
      polygon2(t_val2, `2.5%`, `97.5%`, col = adjustcolor(col_model, .1))
      lines(t_val2, `50%`, lwd = lwd_val, col = col_model)
    })
  
  with(data, points(jitter(t), value, col = col_data, lwd = lwd_val, pch = 3))
}

plot_data_model(alpha_data, alpha_trajectories,
                bquote("reduction" ~ alpha ~ "of mosquitoes with oocysts"))

lines(t_val2, with(as.list(alpha_fit4dpi), f4dpi(t_val2, A, k, d)))

beta_data <- reading_exp_data("figure2c.rds")

beta_fit4dpi <- with(beta_data, calibrate4dpi(t, value))

beta1 <- with(as.list(beta_fit4dpi), make_effect(A, k, d, d))

beta2 <- with(as.list(beta_fit4dpi), make_effect(A, k, d, d / 3))



if (file_exists("beta_model.rds")) {
  beta_model_4dpi <- readRDScache("beta_model.rds")
} else {
  beta_model_4dpi <- boot_calibration_4dpi(beta_data)
  saveRDScache(beta_model_4dpi, "beta_model.rds")
}

beta_trajectories <- trajectories(beta_model_4dpi)

plot_data_model(beta_data, beta_trajectories,
                bquote("reduction" ~ beta ~ "of the number of oocysts"))

lines(t_val2, with(as.list(beta_fit4dpi), f4dpi(t_val2, A, k, d)))

bompard2020 <- "doi_10_5061_dryad_9p8cz8wfh__v20210123/act2datasetRpcr.txt" |>
  read_table2()

bompard2020 |> 
  filter(0 < Total, Total < 50) |> 
  with(hist(Total, 0:50, xlab = "oocyst load L", ylab = "frequency",
            main = NA, col = 4))

bompard2020 |> 
  group_by(Total) |> 
  tally() |> 
  mutate(cum_n = cumsum(n) / sum(n)) |> 
  with(plot(Total, cum_n, type = "h", col = 4, xlim = c(0, 100), ylim = c(.85, 1),
            xlab = "oocyst load L", ylab = "cumulative proportion of samples"))

weights <- bompard2020 |> 
  filter(Total > 0) |> 
  group_by(Total) |> 
  tally() |> 
  mutate(cum_n = cumsum(n) / sum(n)) |> 
  with(approx(Total, cum_n, 1:max(Total))) |> 
  setNames(c("load", "weight")) |> 
  as_tibble() |>
  mutate(across(weight, ~ c(.x[1], diff(.x))))

aleshnick2020 <- "ppat.1008181.s006.xlsx" |> 
  read_excel2() |> 
  mutate(load = `Salivary Gland Load` / 1000) |> 
  rename(infection = `Malaria Infection`) |> 
  select(load, infection)

g <- function(L, B, s) B * (1 - exp(-s * L))

th <- 25

B <- aleshnick2020 |> 
  filter(load > th) |> 
  pull(infection) |> 
  mean()

s <- aleshnick2020 |> 
  filter(load < th) |> 
  arrange(load) |> 
  map(~ zoo::rollmean(.x, 10)) |> 
  as_tibble() |> 
  lm(infection ~ -1 + load, data = _) |> 
  coef()

mLL <- function(x) {
  -sum(with(aleshnick2020, dbinom(infection, 1, g(load, exp(x[1]), exp(x[2])), TRUE)))
}

params <- exp(optim(c(B = log(B), s = log(s)), mLL)$par)
g_calibrated <- function(x) g(x, params[1], params[2])

max_load <- 250
load_vals <- 1:max_load

plot(load_vals, g_calibrated(load_vals), type = "h", lwd = lwd_val, col = 4,
     xlim = c(0, max_load), ylim = c(0, .5),
     xlab = "oocyst load L", ylab = "probability of transmission")

zeta <- function(lambda, w = weights) {
  lambda |> 
    map_dbl(~ with(w, sum(weight * g_calibrated(.x * load) / g_calibrated(load))))
}

lambda_vals <- seq(0, 1, le = 512)
zeta_vals <- zeta(lambda_vals)

opar <- par(pty = "s", xpd = TRUE, plt = c(.15, .95, .15, .95))
plot(lambda_vals, zeta_vals, type = "l", lwd = lwd_val, col = 4,
     xaxs = "i", yaxs = "i",
     xlab = bquote("reducing factor of oocyst load"),
     ylab = "reducing factor of transmission")
segments(0, 0, 1, 1)
box(bty = "o")
par(opar)

delta_data <- reading_exp_data("percentage_inhibition_size2.rds")

delta_fit4dpi <- with(delta_data, calibrate4dpi(t, value))

delta1 <- with(as.list(delta_fit4dpi), make_effect(A, k, d, d))

delta2 <- with(as.list(delta_fit4dpi), make_effect(A, k, d, d / 3))



if (file_exists("delta_model.rds")) {
  delta_model_4dpi <- readRDScache("delta_model.rds")
} else {
  delta_model_4dpi <- boot_calibration_4dpi(delta_data)
  saveRDScache(delta_model_4dpi, "delta_model.rds")
}

delta_trajectories <- trajectories(delta_model_4dpi)

plot_data_model(delta_data, delta_trajectories,
                bquote("reduction" ~ delta ~ "of the size of oocysts"))

shaw_fig1 <- read_excel2("journal.ppat.1009131.s013.xlsx", "Fig1BC")
shaw_fig3 <- read_excel2("journal.ppat.1009131.s013.xlsx", "Fig3 & FigS2")

(oo_size <- shaw_fig1 |> 
  select(Group, `Blood Feeds`, `Mean oocyst size`) |> 
  mutate(across(`Blood Feeds`, as.factor)) |> 
  na.exclude())

colors <- oo_size |> 
  select(-`Mean oocyst size`) |> 
  unique() |> 
  mutate(color = c(yellow = rgb(.944, .745, .510),
                   blue   = rgb(.636, .657, .812),
                   red    = rgb(.929, .669, .651),
                   purple = rgb(.790, .679, .822)))

oo_size <- left_join(oo_size, colors, c("Group", "Blood Feeds"))

lty_val <- setNames(rep(1:2, each = 2), colors$color)

with(oo_size, table(Group, `Blood Feeds`))

oosise_lab <- expression(mean ~ oocyst ~ size ~ at ~ 7 ~ dpi ~ (µm^2))

oo_size |> 
  mutate(xlab = paste(Group, `Blood Feeds`)) |> 
  group_split(Group, `Blood Feeds`) |>
  bind_rows(.id = "x") |> 
  mutate(across(x, as.integer)) |> 
  with({
    plot(jitter(x), `Mean oocyst size`, col = color, pch = 19, ylim = c(0, 2e3),
         xaxt = "n", xlab = NA, ylab = oosise_lab)
    axis(1, 1:4, unique(xlab), tick = FALSE, line = -1, lty = 0)
  })

(mean_oosize <- oo_size |> 
  mutate(bf = as.factor(`Blood Feeds`)) |> 
  summarise(mod = list(lm(`Mean oocyst size` ~ 1)), .by = c(Group, bf, color)) |> 
  mutate(oosize_M = map_dbl(mod, coef),
         ci = map(mod,
                  ~ setNames(as.vector(confint(.x)), c("oosize_L", "oosize_U")))) |> 
  unnest_wider(ci) |> 
  select(-mod))

(EIP <- shaw_fig3 |> 
  mutate(infected = `Infected?` == "Yes") |> 
  select(Group, `Blood Feeds`, `Days pIBM`, infected) |> 
  mutate(across(`Blood Feeds`, as.factor)) |> 
  na.exclude() |> 
  left_join(colors, c("Group", "Blood Feeds")))

second_meal <- 3

models <- EIP |> 
  mutate(across(`Blood Feeds`, as.factor),
         time = `Days pIBM` - second_meal) |>
  summarise(mod = list(glm(infected ~ time + time:`Blood Feeds`, binomial)),
            .by = Group)

add_data <- function(x) {
  x |> 
    summarise(prop = list(prop.test(sum(infected), length(infected))),
              .by = c(`Days pIBM`, color)) |> 
    mutate(p  = map_dbl(prop, "estimate"),
           ci = map(prop, ~ setNames(.x$conf.int, c("lwr", "upr")))) |> 
    unnest_wider(ci) |> 
    with({
      col <- unique(color)
      points(`Days pIBM`, p, pch = 19, col = col)
      arrows2(`Days pIBM`, lwr, `Days pIBM`, upr, col = col)
    })
}

find_EIP <- function(mod, bf, val = "fit") {
  optimise(
    function(x) abs(predict_glm(mod, list(time = x, `Blood Feeds` = bf))[[val]] - .5),
    c(5, 20)
  )$min
}

plot1combination <- function(gp, bf) {
  col <- colors |>
    filter(Group == gp, `Blood Feeds` == bf) |>
    pull(color)
  
  mod <- models |> 
    filter(Group == gp) |> 
    pull(mod) |> 
    first()
  
  dpis <- seq(0, 15, le = 512)
  plot3(dpis + second_meal,
        predict2(mod, tibble(time = dpis, `Blood Feeds` = bf)), col = col)
    
  EIP |> 
    filter(Group == gp, `Blood Feeds` == bf) |> 
    add_data()
  
  abline(h = .5)
  abline(v = second_meal + c(0, find_EIP(mod, bf)))
}

opar <- par(plt = c(.15, .95, .21, .89))
with(colors, walk2(Group, `Blood Feeds`, ~ plot1combination(.x, .y)))
par(opar)

plot2(NA, xlim = c(3, 18), ylim = 0:1)

dpis <- seq(0, 15, le = 512)
colors |> 
  left_join(models, "Group") |> 
  mutate(predictions = map2(mod, `Blood Feeds`,
                            ~ predict2(.x, tibble(time          = dpis,
                                                  `Blood Feeds` = .y)))) |> 
  with(walk2(predictions, color,
             ~ lines(dpis + second_meal, .x,
                     lwd = lwd_val, col = .y, lty = lty_val[.y])))

abline(h = .5)
abline(v = c(3, 7))

rate_val <- models |> 
  pull(mod) |> 
  map(~ cumsum(coef(.x)[-1])) |> 
  unlist()

placebo <- c(7:9, 13, 17, 19, 20)

oocyst_size <- "oocyst_size2.txt" |>
  read_table2() |> 
  mutate(across(Area, ~ .x * 1e6)) |> 
  rename(patient_id = patient,
         time_point = time_points,
         mosquito_id = Mosq_no,
         area = Area) |> 
  mutate(treatment = case_when(patient_id %in% placebo ~ "AL+Placebo",
                               .default = "AL+AP")) |> 
  select(treatment, patient_id, mosquito_id, time_point, area, Total_oocysts)

oocyst_size_crtl <- oocyst_size |> 
  filter(time_point %in% c("Unfed", "Control")) |> 
  pull(area)

hist(oocyst_size_crtl, 20, main = NA, col = 4,
     xlab = "oocyst size", ylab = "frequency")

mean_oosize_crtl <- mean(oocyst_size_crtl)
abline(v = mean_oosize_crtl, col = 2, lwd = lwd_val)

mos <- mean_oosize$oosize_M
mod <- lm(rate_val ~ mos)
compute_rate <- function(mean_oo_size) predict(mod, list(mos = mean_oo_size))

plot(mean_oosize$oosize_M, rate_val, col = colors$color, pch = 19, cex = 1.5,
     xlim = c(0, 1100), ylim = c(0, 1.1),     
     xlab = oosise_lab, ylab = "oocyst maturity rate")

abline(mod, col = 4)
abline(v = mean_oosize_crtl)
abline(h = compute_rate(mean_oosize_crtl))

m_oo_size_val <- mean(oocyst_size_crtl)

size2rate_reduction <- function(reduction_in_size, ref = m_oo_size_val) {
  compute_rate(reduction_in_size * ref) / compute_rate(ref)
}

#size2rate_reduction(c(.8, .7))
#size2rate_reduction(0)

reduction1 <- function(alpha, beta, delta, n, mu, a, p, m, infinity) {
  l <- function(tau, m) {
    map_dbl(tau,
            ~ integrate2(function(x) (1 - alpha(x, .x)) * (1 - beta(x, .x)), 0, m) / m)
  }
  
  qf <- function(tau, m) {
    map_dbl(tau,
            ~ integrate2(function(x) size2rate_reduction(1 - delta(x, .x)), 0, m) / m)
  }
  
  lambda <- function(tau, a, p, m) {
    map_dbl(tau, ~ exp(-a * p * integrate2(function(x) 1 - l(x, m), 0, .x)))
  }
  
  kappaf <- function(tau, a, p, m) {
    map_dbl(tau, ~ exp(-a * p * integrate2(function(x) 1 - qf(x, m), 0, .x)))
  }
  
  theta <- function(mu, a, p, m, infinity) {
    integrate2(function(x) dexp(x, mu) * zeta(lambda(x, a, p, m)), 0, infinity)
  }
  
  EIPf <- function(n, a, p, m) {
    optimise(function(x) abs(x * kappaf(x, a, p, m) - n), c(0, 24 * 100))$min
  }
  
  theta(mu, a, p, m, infinity) * exp((n - EIPf(n, a, p, m)) * mu)  
}

mass_treatment <- function(n = 10, mu = 1 / 14, a = 1 / 3, T_min = 3, T_max = 8,
                           p_vals = 0:10 / 10, infinity = 250,
                           alpha = alpha2, beta = beta2, delta = delta2) {
  n <- n * 24
  mu <- mu / 24
  a <- a / 24
  infinity <- infinity * 24
  
  expand_grid(p_vals, m_vals = (T_min * 7):(T_max * 7) * 24) |> 
    mutate(reduction =
             future_map2_dbl(p_vals, m_vals,
                             ~ reduction1(alpha, beta, delta, n, mu, a, p = .x, m = .y,
                                          infinity)))
}

add_legend <- function(where, title = title) {
  legend(where, legend = paste(1:10 * 10, "%"), col = tail(YlOrRd11, -1), lty = 1,
         lwd = lwd_val, title = title, adj = c(1, .5), x.intersp = 4, y.intersp = .85,
         xjust = .5, ncol = 2, inset = .03)
}

plot_mass_treatment <- function(data) {
  with(data, plot(m_vals, reduction, type = "n", xaxt = "n", xaxs = "i", yaxs = "i",
                  xlab = "intake frequency (/week)",
                  ylab = bquote(R[0] ~ "reduction factor"),
                  xlim = rev(range(m_vals)), ylim = 0:1))
  
  hours_per_week <- 7 * 24
  T_range <- range(data$m_vals) / hours_per_week
  labels <- T_range[1]:T_range[2]
  ats <- labels * hours_per_week
  axis(1, ats, labels)
  abline(v = ats, col = "grey")
  abline(h = 0:10 / 10, col = "grey")
  
  data |> 
    group_split(p_vals) |> 
    walk2(YlOrRd11, ~ with(.x, lines(m_vals, reduction, col = .y, lwd = lwd_val)))
  
  add_legend("bottomleft", title = "proportion of population treated")
  box(bty = "o")
}


####################################################################################### OK


if (file_exists("mtreat1.rds")) {
  mtreat1 <- readRDScache("mtreat1.rds")
} else {
  mtreat1 <- mass_treatment()
  saveRDScache(mtreat1, "mtreat1.rds")
}

plot_mass_treatment(mtreat1)

susceptible <- function(R0, c, a, mu) {
  tmp <- a * c / mu
  (1 + tmp) / (R0 + tmp)
}

proba_infection <- function(R0, r, m, a, c, mu) {
  1 - exp(-r * m * (R0 - 1) / (1 + a * c / mu))
}

incidence <- function(R0, c, a, mu, r, m)  {
  tmp1 <- a * c / mu
  tmp2 <- tmp1 + 1
  tmp2 * (1 - exp(-r * m * (R0 - 1) / tmp2)) / (R0 + tmp1)
}

clinical_cases <- function(a = 1 / 3, mu = 1 / 14, r = 1 / 200, m = 30, c = .3) {
  R0s <- seq(1, 100, le = 512)
  prop_susc <- map_dbl(R0s, ~ susceptible(.x, c, a, mu))
  p_inf <- map_dbl(R0s, ~ proba_infection(.x, r, m, a, c, mu))
  plot(R0s, prop_susc, type = "l", lwd = lwd_val, col = 4,
       xlab = bquote(R[0]), ylab = "proportion or probability")
  lines(R0s, p_inf, lwd = lwd_val, col = 3)
  lines(R0s, p_inf * prop_susc, lwd = lwd_val, col = 2)
  legend("right", legend = c("proportion of susceptibles",
                             "probability of getting severely infected",
                             "proportion of clinical cases"),
         lty = 1, lwd = lwd_val, col = 4:2, bty = "n")
}

clinical_cases(m = 30)

R0s <- seq(1, 100, le = 512)
m <- 30
inc_vals <- map_dbl(R0s, ~ incidence(.x, c = .3, a = 1 / 3, mu = 1 / 14,
                                     r = 1 / 200, m = m))
plot(R0s, inc_vals, type = "l", lwd = lwd_val, col = 2,
     xlab = bquote(R[0]), ylab = paste("new clinical cases over the past", m, "days"))

dynamics <- function(R0i = 40, n = 10, mu = 1 / 14, a = 1 / 3, prop = 1, m = 50,
                     cv = .3, r = 1 / 200, infinity = 250, nb = 10,
                     alpha = alpha2, beta = beta2, delta = delta2) {
  R0 <- p <- rep(0, nb + 1)
  R0 <- R0_tmp <- R0i
  for (i in 2:(nb + 1)) {
    p[i - 1] <- p_tmp <- incidence(R0_tmp, cv, a, mu, r, m)
    R0[i] <- R0_tmp <- reduction1(alpha, beta, delta,
                                  n, mu, a, prop * p_tmp, m, infinity) * R0i
  }
  list(R0 = R0, p = head(p, -1))
}

plot_dynamics <- function(x) {
  plott <- function(...) {
    plot(..., type = "o", xlab = "iteration", col = 4, lwd = lwd_val)
  }
  with(x, {
    plott(R0, ylab = bquote(R[0]))
    plott(p,  ylab = "incidence")
  })
}



if (file_exists("out40.rds")) {
  out40 <- readRDScache("out40.rds")
} else {
  out40 <- dynamics()
  saveRDScache(out40, "out40.rds")
}

if (file_exists("out30.rds")) {
  out30 <- readRDScache("out30.rds")
} else {
  out30 <- dynamics()
  saveRDScache(out30, "out30.rds")
}

plot_dynamics(out40)
plot_dynamics(out30)

finalR0nl <- function(R0i, n, mu, a, prop, m, cv, r, infinity = 250,
                      alpha = alpha2, beta = beta2, delta = delta2) {
  nleqslv::nleqslv(c(R0i, incidence(R0i, cv, a, mu, r, m)),
                   function(vars) {
                      R0 <- vars[1]
                      p  <- vars[2]
                      c(R0 - reduction1(alpha, beta, delta,
                                        n, mu, a, prop * p, m, infinity) * R0i,
                        p  - incidence(R0, cv, a, mu, r, m)) })$x[1]
}

finalR0rs <- function(R0i, n, mu, a, prop, m, cv, r, infinity = 250,
                      alpha = alpha2, beta = beta2, delta = delta2) {
  rootSolve::multiroot(function(vars) {
                         R0 <- vars[1]
                         p  <- vars[2]
                         c(R0 - reduction1(alpha, beta, delta,
                                           n, mu, a, prop * p, m, infinity) * R0i,
                           p  - incidence(R0, cv, a, mu, r, m)) },
                       c(R0i, incidence(R0i, cv, a, mu, r, m)))$root[1]
}



m_vals <- c(4:7 * 10, seq(100, 250, 50))
if (file_exists("R0red.rds")) {
  R0red <- readRDScache("R0red.rds")
} else {
  R0red <- expand_grid(m  = m_vals,
                       R0 = 5:15 * 10) |> 
    mutate(reduction = future_map2_dbl(R0, m,
                                       ~ finalR0nl(R0i = .x, n = 10, mu = 1 / 14,
                                                   a = 1 / 3, prop = 1, m = .y,
                                                   cv = .3, r = 1 / 200) / .x))
  saveRDScache(R0red, "R0red.rds")
}

with(R0red, plot(R0, reduction, col = 4))

R0_new <- seq(15, 150, le = 512)

R0red_Hill <- R0red |> 
  summarise(fit = list(nls(reduction ~ R0^h / (K^h + R0^h),
                           start = c(K = 5, h = .7))), .by = m) |> 
  mutate(reduction = map(fit, ~ predict(.x, newdata = data.frame(R0 = R0_new)))) |> 
  select(-fit) |> 
  left_join(expand_grid(m = m_vals, R0 = list(R0_new)), "m") |> 
  unnest(c(reduction, R0))

plot(NA, xlim = c(0, 150), ylim = c(.6, 1),
     xlab = bquote(initial ~ value ~ of ~ R[0]), ylab = bquote(R[0] ~ reduction))
R0red_Hill |> 
  group_by(m) |> 
  group_walk(~ with(.x, lines(R0, reduction, col = 4)))

R0red_Hill |> 
  filter(R0 == 15) |> 
  with(plot(m, reduction, type = "o", col = 4, lwd = lwd_val,
            xlab = "time period of incidence definition (days)",
            ylab = bquote(R[0] ~ reduction ~ at ~ R[0] == 15)))



m_vals2 <- 64:75
if (file_exists("reduc_val.rds")) {
  reduc_val <- readRDScache("reduc_val.rds")
} else {
  reduc_val <- future_map_dbl(m_vals2,
                              ~ finalR0nl(R0i = 40, n = 10, mu = 1 / 14, a = 1 / 3,
                                          prop = 1, m = .x, cv = .3, r = 1 / 200) / 40)
  saveRDScache(reduc_val, "reduc_val.rds")
}

plot(m_vals2, reduc_val, type = "o", col = 4, lwd = lwd_val,
     xlab = "time period of incidence definition (days)",
     ylab = bquote(R[0] ~ reduction ~ at ~ R[0] == 40))

search_m <- function(from = 15, R0i = 40, n = 10, mu = 1 / 14, a = 1 / 3, prop = 1,
                     cv = .3, r = 1 / 200, infinity = 250,
                     alpha = alpha2, beta = beta2, delta = delta2) {
  f <- function(x) {
    finalR0nl(R0i, n, mu, a, prop, m = x, cv, r, infinity, alpha, beta, delta)
  }
  
  out <- from |>
    seq(by = n_workers, le = n_workers) |> 
    tibble(x = _) |> 
    mutate(y = future_map_dbl(x, f))
  
  idx <- out |>
    mutate(row_id = row_number()) |> 
    filter(y == min(y)) |> 
    pull(row_id)
  
  from <- out |>
    slice(idx - 1) |> 
    pull(x)

  out <- from |>
    seq(by = 2, le = n_workers) |> 
    tibble(x = _) |> 
    mutate(y = future_map_dbl(x, f)) |> 
    filter(y == min(y))
  
  x <- out$x
  
  c(x - 1, x + 1) |> 
    tibble(x = _) |> 
    mutate(y = future_map_dbl(x, f)) |> 
    bind_rows(out) |> 
    filter(y == min(y)) |> 
    pull(x)
}



if (file_exists("m_min.rds")) {
  readRDScache("m_min.rds")
} else {
  (m_min <- search_m(15))
  saveRDScache(m_min, "m_min.rds")
}



if (file_exists("R0red2.rds")) {
  R0red2 <- readRDScache("R0red2.rds")
} else {
  R0red2 <- expand_grid(prop = 1:10 / 10,
                        R0   = 5:15 * 10) |> 
    mutate(reduction = future_map2_dbl(R0, prop,
                                       ~ finalR0nl(R0i = .x, n = 10, mu = 1 / 14,
                                                   a = 1 / 3, prop = .y, m = 69,
                                                   cv = .3, r = 1 / 200) / .x))
  saveRDScache(R0red2, "R0red2.rds")
}

extrapolate <- function(x, R0_new = seq(15, 150, le = 512), K0 = 5, h0 = .7) {
# Trying to fit a Hill model on reduction as a function of R0 for each value of prop:
  first_try <- x |> 
    summarise(fit = list(safely(nls)(reduction ~ R0^h / (K^h + R0^h),
                                     data = pick(everything()),
                                     start = c(K = K0, h = h0))), .by = prop) |> 
    filter_out(map_lgl(map(fit, "result"), is.null)) |> 
    mutate(across(fit, ~ map(.x, "result")),
           coeff = map(fit, coef)) |>
    unnest_wider(coeff)

# Modelling the values of the parameters of the Hill model as a function of prop:
  modK <- first_try |> 
    select(prop, K) |> 
    bind_rows(tibble(prop = 0, K = 0)) |> 
    lm(K ~ prop + I(prop^2), data = _)
  
  modh <- lm(h ~ prop + I(prop^2), first_try)

# Using the output of the second modelling as starting values of the first modelling:
  x |> 
    nest(.by = prop) |>
    mutate(K0  = predict(modK, list(prop = prop)),
           h0  = predict(modh, list(prop = prop)),
           fit = pmap(list(data, K0, h0),
                      ~ nls(reduction ~ R0^h / (K^h + R0^h),
                            data = ..1, start = c(K = ..2, h = ..3))),
           reduction = map(fit,
                           ~ tibble(R0 = R0_new,
                                    prediction =
                                      predict(.x, newdata = list(R0 = R0_new))))) |> 
    select(-c(K0, h0, fit))
}

out <- extrapolate(R0red2)

plot(NA, xlim = c(0, 150), ylim = 0:1,
     xlab = bquote(initial ~ value ~ of ~ R[0]), ylab = bquote(R[0] ~ reduction))

out |> 
  select(prop, reduction) |> 
  unnest(reduction) |> 
  group_by(prop) |> 
  group_walk(~ with(.x, lines(R0, prediction, col = 4)))

out |> 
  select(prop, data) |> 
  unnest(data) |> 
  group_by(prop) |> 
  group_walk(~ with(.x, points(R0, reduction, col = 4)))

treatment_infected <- function(R0 = 5:15 * 10, prop = 1:10 / 10, from = 15, n = 10,
                               mu = 1 / 14, a = 1 / 3, c = .3, r = 1 / 200,
                               R0_new = seq(15, 150, le = 512), K0 = 5, h0 = .7,
                               infinity = 250,
                               alpha = alpha2, beta = beta2, delta = delta2) {
  m_opt <- search_m(from, first(R0), n, mu, a, prop = 1, c, r, infinity,
                    alpha, beta, delta)
  expand_grid(prop, R0) |> 
    mutate(reduction = future_map2(R0, prop,
                                   ~ safely(finalR0nl)(R0i = .x, n, mu, a, prop = .y,
                                                       m = m_opt, c, r, infinity,
                                                       alpha, beta, delta))) |> 
    filter(map_lgl(reduction, ~ is.null(.x$error))) |> 
    mutate(across(reduction, ~ map_dbl(.x, "result") / R0)) |> 
    extrapolate(R0_new, K0, h0)
}

plot_treatment_infected <- function(data) {
  plot(NA, xlim = c(15, 150), ylim = 0:1, xaxs = "i", yaxs = "i",
     xlab = bquote(initial ~ value ~ of ~ R[0]), ylab = bquote(R[0] ~ reduction))
  
  abline(v = 0:15 * 10, col = "grey")
  abline(h = 0:10 / 10, col = "grey")
  
  data |> 
    select(prop, reduction) |> 
    unnest(reduction) |> 
    group_split(prop) |> 
    walk2(tail(YlOrRd11, -1),
          ~ with(.x, lines(R0, prediction, col = .y, lwd = lwd_val)))
  
  add_legend("bottomright", title = "proportion of cases treated")
  
  box(bty = "o")
}



if (file_exists("treati1.rds")) {
  treati1 <- readRDScache("treati1.rds")
} else {
  treati1 <- treatment_infected()
  saveRDScache(treati1, "treati1.rds")
}

plot_treatment_infected(treati1)

tmp1 <- unique(oocyst_size$time_point)
tmp2 <- c(tmp1[str_starts(tmp1, "H")], tmp1[str_starts(tmp1, "D")])

color_dpa <- c("black", "grey", rev(tail(YlOrRd11, length(tmp2)))) |> 
  setNames(c(setdiff(tmp1, tmp2), tmp2))

oocyst_size <- oocyst_size |> 
  mutate(color = color_dpa[time_point])

densities <- oocyst_size |> 
  mutate(across(area, log10)) |> 
  select(time_point, area, color) |> 
  nest(data = area) |> 
  mutate(density = map(data, ~ density(.x$area)))

densities_plot <- densities |>
  mutate(across(density, ~ map(.x, function(x) {x$x <- 10^x$x; x})))

x_max <- densities_plot |> 
  pull(density) |> 
  map("x") |> 
  unlist() |> 
  max()

x_min <- densities_plot |> 
  pull(density) |> 
  map("x") |> 
  unlist() |> 
  min()

y_max <- densities_plot |> 
  pull(density) |> 
  map("y") |> 
  unlist() |> 
  max()

plot(NA, xlim = c(x_min, x_max), ylim = c(0, y_max), log = "x",
     xlab = "oocsyt size at 7 dpi", ylab = "density")
with(densities_plot, walk2(density, color, ~ lines(.x, col = .y, lwd = lwd_val)))

legend("topleft", legend = names(color_dpa), col = color_dpa, lty = 1, lwd = lwd_val,
       bty = "n")

xs <- seq(1, log10(5000), le = 512)

time_points <- names(color_dpa) 

densities |> 
  arrange(match(time_point, time_points)) |> 
  mutate(den_val = map(density, ~ with(.x, tibble(approx(x, y, xs)$y)))) |> 
  pull(den_val) |> 
  bind_cols() |> 
  as.matrix() |> 
  t() |> 
  image(1:10, xs, z = _, xaxt = "n", yaxt = "n", xlab = NA, ylab = "oocyst size at 7 dpi")

labels <- 10^(1:3) |> 
  (\(x) sort(c(x, 5 * x)))()

axis(1, 1:10, time_points)
axis(2, log10(labels), labels)
box(bty = "o")

H00 <- oocyst_size |> 
  filter(time_point == "H00") |> 
  pull(area)

dst <- density(H00, from = 0)

xs <- seq(0, max(H00), le = 512)
fitdstr <- as.list(fitdistrplus::fitdist(H00, "lnorm")$estimate)
ys <- with(fitdstr, dlnorm(xs, meanlog, sdlog))

plot(dst, lwd = lwd_val, col = 4, main = NA, ylim = c(0, max(dst$y, ys)),
     xlab = "mean oocyst sizes at 7 dpi and H00", ylab = "density")
lines(xs, ys, lwd = lwd_val, col = 2)

legend("topright", legend = c("kernel density", "log-normal fit"), lty = 1,
       lwd = lwd_val, col = c(4, 2), bty = "n")

abline(v = 50)

ys2 <- with(fitdstr, dlnorm(dst$x, meanlog, sdlog))
plot(dst$y, ys2, xlab = "kernel density", ylab = "log-normal fit", col = 4)
abline(0, 1)

densities <- oocyst_size |> 
  filter(time_point == "H00") |> 
  select(patient_id, mosquito_id, area) |> 
  nest(data = area) |> 
  filter(map_int(data, nrow) > 1) |> 
  mutate(density = map(data, ~ density(.x$area, from = 0))) |> 
  pull(density)

plot(NA, xlim = c(0, max(unlist(map(densities, "x")))), ylim = c(0, .01),
     xlab = "oocyst size at 7 dpi and H00", ylab = "density")
walk(densities, ~ lines(.x, col = adjustcolor(1, .2)))

oocyst_size |> 
  filter(time_point == "H00") |> 
  pull(area) |> 
  density() |> 
  lines(col = 2, lwd = lwd_val)

abline(v = 50, lwd = lwd_val, col = 4)

oo_size_mu_sd <- oocyst_size |> 
  filter(time_point == "H00") |> 
  select(patient_id, mosquito_id, area) |> 
  summarise(mean = mean(area),
            sd   = sd(area),
            n    = n(), .by = c(patient_id, mosquito_id))

with(oo_size_mu_sd, plot(mean, sd, xlim = c(0, max(mean)), cex = sqrt(n), col = 4,
     xlab = "mean oocyst size", ylab = "standard deviation of oocyst size"))
mod <- lm(sd ~ -1 + mean, oo_size_mu_sd, weights = n)
abline(mod)

(gamma_shape <- 1 / coef(mod)^2)
(log_sd <- sqrt(log(1 + coef(mod)^2)))

oo_size_mu_sd <- oocyst_size |> 
  select(time_point, patient_id, mosquito_id, area) |> 
  summarise(mean = mean(area),
            sd   = sd(area),
            n    = n(), .by = c(time_point, patient_id, mosquito_id))

with(oo_size_mu_sd, plot(mean, sd, xlim = c(0, max(mean)), cex = sqrt(n),
                         col = adjustcolor(4, .15),
                         xlab = "mean oocyst size",
                         ylab = "standard deviation of oocyst size"))
mod <- lm(sd ~ -1 + mean, oo_size_mu_sd, weights = n)
abline(mod)

(gamma_shape <- 1 / coef(mod)^2)
(log_sd <- sqrt(log(1 + coef(mod)^2)))

oocyst_count <- "ASAAP4b_oocyst2.txt" |>
  read_table2() |> 
  rename(time_point = time_points,
         count      = total_oocysts) |> 
  mutate(across(time_point, ~ case_when(.x == "H0" ~ "H00",
                                        .x == "H5" ~ "H05",
                                        .x == "D5" ~ "D05",
                                        .x == "D7" ~ "D07", .default = .x)),
         treatment = case_when(patient_id %in% placebo ~ "AL+Placebo",
                               .default = "AL+AP")) |> 
  select(treatment, patient_id, time_point, day_of_dissection, count)

oocyst_size2 <- oocyst_size |> 
  filter_out(time_point %in% c("Unfed", "Control")) |> 
  filter(treatment == "AL+AP") |> 
  select(patient_id, mosquito_id, time_point, area) |> 
  na.exclude()

oocyst_count2 <- oocyst_count |> 
  filter_out(time_point %in% c("Unfed", "Control")) |> 
  filter(day_of_dissection == "7dpi", treatment == "AL+AP") |> 
  select(patient_id, time_point, count) |> 
  na.exclude()

####################################################################################### NOT OK 

fit_patient_x_timepoint_fixed <- readRDScache("fit_patient_x_timepoint_fixed.rds")

hours <- c(0, 5, 48, 72)
days <- c(5, 7, 14, 28)
hash <- c(setNames(hours, sprintf("H%02d", hours)),
          setNames(24 * days, sprintf("D%02d", days)))


#######################################################################################

delta_data2 <- fit_patient_x_timepoint_fixed |> 
  magrittr::use_series(patient_effects) |> 
  mutate(t = hash[as.character(time_point)]) |> 
  select(t, size_reduction) |> 
  filter(t > 0) |> 
  rename(value = size_reduction) |> 
  arrange(t)

delta_fit4dpi2 <- with(delta_data2, calibrate4dpi(t, value))

delta3 <- with(as.list(delta_fit4dpi2), make_effect(A, k, d, d / 3))

# mtreat1 <- mass_treatment()

#######################################################################################
bootstrap_functions <- function(x, f = function(x) x / 3) {
  map(x, ~ with(as.list(.x), make_effect(A, k, d, f(d))))
}

alpha2s <- bootstrap_functions(alpha_model_4dpi)
beta2s  <- bootstrap_functions(beta_model_4dpi)
delta2s <- bootstrap_functions(delta_model_4dpi)



N_target <- 1e3
N <- n_workers * ceiling(1.05 * N_target / n_workers)
if (file_exists("conv.rds")) {
  conv <- readRDScache("conv.rds")
} else {
  conv <- future_pmap(list(sample(alpha2s, N), sample(beta2s, N), sample(delta2s, N)),
                      ~ safely(reduction1)(..1, ..2, ..3, 10 * 24, 1 / (14 * 24),
                                           1 / (3 * 24), .1, 50 * 24, 50 * 24)) |> 
    keep(~ is.null(.x$error)) |> 
    map_dbl("result")
  saveRDScache(conv, "conv.rds")
}

cumconv <- cummean(conv)
plot(cumconv, ylim = c(.9, .95), type = "n")
abline(v = 0:15 * 100, col = "grey")
abline(h = 91:95 / 100, col = "grey")
lines(cumconv, col = 4, lwd = lwd_val)

reduction2 <- function(n, mu, a, p, m, infinity, alpha, beta, delta) {
  list(alpha, beta, delta) |> 
    pmap(~ safely(reduction1)(..1, ..2, ..3, n, mu, a, p, m, infinity)) |> 
    keep(~ is.null(.x$error)) |> 
    map_dbl("result") |> 
    mean()
}

#reduction2(10 * 24, 1 / (14 * 24), 1 / (3 * 24), .1, 50 * 24, 50 * 24,
#           alpha2s[1:2], beta2s[1:2], delta2s[1:2])

mass_treatment2 <- function(n = 10, mu = 1 / 14, a = 1 / 3, T_min = 3, T_max = 8,
                            p_vals = 0:10 / 10, infinity = 250, alpha, beta, gamma) {
  n <- n * 24
  mu <- mu / 24
  a <- a / 24
  infinity <- infinity * 24
  
  expand_grid(p_vals, m_vals = (T_min * 7):(T_max * 7) * 24) |> 
    mutate(reduction =
             future_map2_dbl(p_vals, m_vals,
                             ~ reduction2(n, mu, a, p = .x, m = .y, infinity,
                                          alpha, beta, gamma)))
}



N <- 50
if (file_exists("mtreat2_50.rds")) {
  mtreat2_50 <- readRDScache("mtreat2_50.rds")
} else {
  mtreat2_50 <- mass_treatment2(alpha = sample(alpha2s, N),
                                beta  = sample(beta2s,  N),
                                gamma = sample(delta2s, N))
  saveRDScache(mtreat2_50, "mtreat2_50.rds")
}

plot_mass_treatment(mtreat2_50)

finalR0nl2 <- function(R0i, n, mu, a, prop, m, cv, r, infinity = 250,
                       alpha, beta, delta) {
  nleqslv::nleqslv(c(R0i, incidence(R0i, cv, a, mu, r, m)),
                   function(vars) {
                      R0 <- vars[1]
                      p  <- vars[2]
                      c(R0 - reduction2(n, mu, a, prop * p, m, infinity,
                                        alpha, beta, delta) * R0i,
                        p  - incidence(R0, cv, a, mu, r, m)) })$x[1]
}

#m = 63:74
#R0 = 40
#n = 10
#mu = 1 / 14
#a = 1 / 3
#prop = 1
#c = .3
#r = 1 / 200

#finalR0nl2(R0i, n, mu, a, prop, m, cv, r, infinity = 250, alpha, beta, delta)
#finalR0nl2(R0, n, mu, a, prop, m = .x, c, r,                 alpha, beta, delta)

search_m2 <- function(m = 63:74, R0 = 40, n = 10, mu = 1 / 14, a = 1 / 3, prop = 1,
                      c = .3, r = 1 / 200, infinity = 250, alpha, beta, delta) {
  m |>
    tibble() |> 
    mutate(reduction =
             future_map_dbl(m, ~ finalR0nl2(R0, n, mu, a, prop, m = .x, c, r, infinity,
                                            alpha, beta, delta) / R0)) |> 
    filter(reduction == min(reduction)) |> 
    pull(m)
}

treatment_infected2 <- function(R0 = 5:15 * 10, prop = 1:10 / 10, m = 63:74, n = 10,
                                mu = 1 / 14, a = 1 / 3, c = .3, r = 1 / 200,
                                R0_new = seq(15, 150, le = 512), K0 = 5, h0 = .7,
                                infinity = 250, alpha, beta, delta) {
  m_opt <- search_m2(m, first(R0), n, mu, a, prop = 1, c, r, infinity,
                     alpha, beta, delta)
  expand_grid(prop, R0) |> 
    mutate(reduction = future_map2(R0, prop,
                                   ~ safely(finalR0nl2)(R0i = .x, n, mu, a, prop = .y,
                                                        m = m_opt, c, r, infinity,
                                                        alpha, beta, delta))) |> 
    filter(map_lgl(reduction, ~ is.null(.x$error))) |> 
    mutate(across(reduction, ~ map_dbl(.x, "result") / R0)) |> 
    extrapolate(R0_new, K0, h0)
}



#R0 = 5:15 * 10
#prop = 1:10 / 10
#m = 63:74
#n = 10
#mu = 1 / 14
#a = 1 / 3
#c = .3
#r = 1 / 200
#R0_new = seq(15, 150, le = 512)
#K0 = 5
#h0 = .7

#N <- 1
#alpha = sample(alpha2s, N)
#beta  = sample(beta2s,  N)
#delta = sample(delta2s, N)

# N <- 2
# if (file_exists("treatinf2_50.rds")) {
#   treatinf2_50 <- readRDS("treatinf2_50.rds")
# } else {
#   system.time(
#   treatinf2_50 <- treatment_infected2(alpha = sample(alpha2s, N),
#                                       beta  = sample(beta2s,  N),
#                                       delta = sample(delta2s, N))
#   )
#   saveRDScache(treatinf2_50, "treatinf2_50.rds")
# }

#reduction1(alpha, beta, delta, n, mu, a, p, m, infinity)

reduction_dpi <- function(dpi, n, mu, a, p, m, infinity) {
  reduction(with(as.list(alpha_fit4dpi), make_effect(A, k, d, d / dpi)),
            with(as.list(beta_fit4dpi),  make_effect(A, k, d, d / dpi)),
            with(as.list(delta_fit4dpi), make_effect(A, k, d, d / dpi)),
            n, mu, a, p, m, infinity)
}

#reduction_dpi(2, 10 * 24, 1 / (14 * 24), 1 / (3 * 24), .1, 50 * 24, 50 * 24)

mass_treatment_dpi <- function(dpi = c(1 / 3, .5, 1, 2, 3),
                               n = 10, mu = 1 / 14, a = 1 / 3,
                               T_min = 3, T_max = 8, p = 1, infinity = 250) {
  n <- n * 24
  mu <- mu / 24
  a <- a / 24
  infinity <- infinity * 24
  
  expand_grid(dpi, m_vals = (T_min * 7):(T_max * 7) * 24) |> 
    mutate(reduction =
             future_map2_dbl(dpi, m_vals,
                             ~ reduction_dpi(dpi = .x, n, mu, a, p, m = .y, infinity)))
}

# system.time(
# mtreat_dpi <- mass_treatment_dpi()
# )

# with(mtreat_dpi, plot(m_vals, reduction))
# with(filter(mtreat_dpi, dpi < .4), lines(m_vals, reduction))
# with(filter(mtreat_dpi, dpi == .5), lines(m_vals, reduction))
# with(filter(mtreat_dpi, dpi == 1), lines(m_vals, reduction))
# with(filter(mtreat_dpi, dpi == 2), lines(m_vals, reduction))
# with(filter(mtreat_dpi, dpi == 3), lines(m_vals, reduction))

future::plan(future::sequential)
