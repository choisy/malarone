### infections prevented per treated

From a given $R_0$, a given prevalence $\bar{X}$ and a given proportion $\tau$ of
infected people that are treated, we can compute the number of human infections 
averted per person treated as so:
  
  $$
  \Delta(\bar{X}, R_0, \tau |\vec\theta) = \frac{(1 - \bar{X})(1 - \rho(\bar{X}, \tau|\vec\theta))R_0}{((1-(1-\rho(\bar{X}, \tau|\vec\theta))\bar{X})R_0 - 1)\tau}
$$
  
  where $\vec\theta = [g,n,L,r,D]^T$. Let's now visualize $\Delta$ as a function of
$R_0$, $\bar{X}$ and $\tau$ for given values of $g$, $n$, $L$, $r$ and $D$. First the
fixed parameters values:

```{r}
# mosquito death rate (inverse of mosquito life expectancy):
g <- 1 / 14  # same unit as n (here in days)
# EPI:
n <- 7  # same unit as g (here in days)
# mean sporozoite load in absence of malarone:
L <- 6e4
# recovery rate in humans (inverse of the mean duration of infection in humans):
r <- 1 / (14 * 24)  # in hours
# duration over which to perform the integration of the effect of malarone:
D <- 28 * 24  # in hours (note: this parameter is not a biological parameter)
```

Let's now compute the $\Delta$ values for a 3D figure:
  
  ```{r message = FALSE}
tau <- c(.001, seq(.1, 1, .1))
X3D <- c(.01, .05, seq(.1, .9, .1))
R0 <- c(5, 10, 15)

lambda_val <- lambda(D, L, g, n, nb = 1e6)  # that's the bit that takes the most time

rho2 <- function(D, X, tau, r) 1 - tau * X * r * D * (1 - lambda_val)

Delta2 <- function(X, R0, D, tau, r) {
  rho_val <- rho2(D, X, tau, r)
  ((1 - X) * (1 - rho_val) * R0) / (((1 - (1 - rho_val) * X) * R0 - 1) * tau)
}

Delta3D <- map(R0, ~ map_dfc(X3D, Delta2, .x, D, tau, r))
```

The 3D figure:
  
  ```{r}
persp2 <- function(x, y, z, ...) {
  persp(x, y, as.matrix(z), theta = 65, col = NA,
        zlim = c(0, max(unlist(Delta3D))), ...)
}

persp3 <- function(...) {
  persp2(..., axes = FALSE, box = FALSE)
  par(new = TRUE)
}

epsilon <- .05
opar <- par(plt = c(epsilon, 1 - epsilon, epsilon, 1 - epsilon))
persp3(tau, X3D, Delta3D[[1]], border = 2)
persp3(tau, X3D, Delta3D[[2]], border = 3)
persp2(tau, X3D, Delta3D[[3]], border = 4,
       xlab = "treatment rate",
       ylab = "prevalence",
       zlab = "inf. prevented per treated")
par(opar)
```

Computing the $\Delta$ values for a 2D figure (takes 17"):

```{r eval = FALSE}
X2D <- seq2(.1, .9)
Delta2D <- map(R0,
               ~ map_dfc(X2D, Delta2, .x, D, tau, r) |>
                 t() |>
                 as.data.frame())
```

```{r include = FALSE}
X2D <- seq2(.1, .9)

if (file_exists("Delta2D.rds")) {
  Delta2D <- readRDS3("Delta2D.rds")
} else {
  Delta2D <- map(R0,
                 ~ map_dfc(X2D, Delta2, .x, D, tau, r) |>
                   t() |>
                   as.data.frame())
  saveRDS2(Delta2D, "Delta2D.rds")
}
```

The 2D figure:

```{r fig.width = 8, fig.height = 3, margin1 = FALSE, margin2 = TRUE}
plot2D <- function(x, ...) {
  plot(rep(X2D, length(x)), unlist(x), type = "n",
       xlim = 0:1, ylim = c(0, max(unlist(Delta2D))),
       xlab = "prevalence", ylab = "infections prevented per treated", ...)
  walk(x, lines, x = X2D, ...)
}

mcx <- 3 / 2
opar <- par(mfrow = c(1, 3), mex = mcx, cex.axis = mcx, cex.lab = mcx, cex.main = mcx)
for (i in 1:3) {
  plot2D(Delta2D[[i]], col = i + 1)
  title(bquote(R[0] ~ "=" ~ .(R0[i])), line = 0)
}
par(opar)
```
