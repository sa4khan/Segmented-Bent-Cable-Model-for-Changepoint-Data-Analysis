#######################################################
#######################################################
# STAN Codes
#######################################################
#######################################################

# STAN Model for the Segmented Bent-Cable Regression: 5-Phase

write("
functions {
vector f(vector t, real b0, real b1, real b2, real b3, real b4, 
    real gam1, real tau1, real tau2, real gam2, int n) {
vector[n] fun;
vector[n] ph;
for (i in 1:n){
if(t[i]>gam1 && t[i]<=tau1)
ph[i] = b2*(t[i]-gam1);
else if(t[i]>tau1 && t[i]<=tau2)
ph[i] = b2*(square(t[i]-tau1)+(tau1-gam1));
else if(t[i]>tau2 && t[i]<=gam2)
ph[i] = b3*((t[i]-tau2)+b2*(square(tau2-tau1)+(tau1-gam1))/b3);
else if(t[i]>gam2)
ph[i] = b4*((t[i]-gam2)+b3*(gam2-tau2)/b4
    + b2*(square(tau2-tau1)+(tau1-gam1))/b4);
else
ph[i] = 0.0;
}
fun = b0+b1*t+ph;
return fun;
}
}
data {
int<lower=1> n;
real b0_mu;
real<lower=0> b0_sd;
real b1_mu;
real<lower=0> b1_sd;
real b2_mu;
real<lower=0> b2_sd;
real b3_mu;
real<lower=0> b3_sd;
real b4_mu;
real<lower=0> b4_sd;
real<lower=0> a_gam1;
real<lower=0> b_gam1;
real<lower=0> a_tau1;
real<lower=0> b_tau1;
real<lower=0> a_tau2;
real<lower=0> b_tau2;
real<lower=0> a_gam2;
real<lower=0> b_gam2;
real<lower=0> a_sigma;
real<lower=0> b_sigma;
real lower_gam1;
vector[n] y;
vector[n] t;
}
parameters {
real b0;
real b1;
real b2;
real b3;
real b4;
real<lower=lower_gam1> gam1;
real<lower=gam1> tau1;
real<lower=tau1> tau2;
real<lower=tau2,upper=1.0*n> gam2;
real<lower=0> sigma;
}
transformed parameters {
real alpha1;
real alpha2;
real alpha3;
real alpha4;
real incoming_diff;
real outgoing_diff;
real<lower=0> inv_sigma2;
alpha1 = b1;
alpha2 = b1 + b2;
alpha3 = b1 + b3;
alpha4 = b1 + b4;
incoming_diff = alpha1 - alpha2;
outgoing_diff = alpha3 - alpha4;
inv_sigma2 = inv_square(sigma);
}
model {
target += normal_lpdf(y | f(t,b0,b1,b2,b3,b4,gam1,
    tau1,tau2,gam2,n), sigma);
target += gamma_lpdf(gam1 | a_gam1, b_gam1);
target += gamma_lpdf(tau1 | a_tau1, b_tau1);
target += gamma_lpdf(tau2 | a_tau2, b_tau2);
target += gamma_lpdf(gam2 | a_gam2, b_gam2);
target += normal_lpdf(b0 | b0_mu, b0_sd);
target += normal_lpdf(b1 | b1_mu, b1_sd);
target += normal_lpdf(b2 | b2_mu, b2_sd);
target += normal_lpdf(b3 | b3_mu, b3_sd);
target += normal_lpdf(b4 | b4_mu, b4_sd);
target += gamma_lpdf(inv_sigma2 | a_sigma, b_sigma);
}

",
"stan_bcable_5phase.stan")

#######################################################

# STAN Model for the Bent-Cable Regression: 3-Phase

write("
functions {
vector f(vector t, real b0, real b1, real b2, 
    real tau1, real tau2, int n) {
vector[n] fun;
vector[n] ph;
for (i in 1:n){
if(t[i]>tau1 && t[i]<=tau2)
ph[i] = square(t[i]-tau1)/(2*(tau2-tau1));
else if(t[i]>tau2)
ph[i] = t[i]-(tau1+tau2)/2;
else
ph[i] = 0.0;
}
fun = b0+b1*t+b2*ph;
return fun;
}
}
data {
int<lower=1> n;
real b0_mu;
real<lower=0> b0_sd;
real b1_mu;
real<lower=0> b1_sd;
real b2_mu;
real<lower=0> b2_sd;
real<lower=0> a_tau1;
real<lower=0> b_tau1;
real<lower=0> a_tau2;
real<lower=0> b_tau2;
real<lower=0> a_sigma;
real<lower=0> b_sigma;
real lower_tau1;
vector[n] y;
vector[n] t;
}
parameters {
real b0;
real b1;
real b2;
real<lower=lower_tau1> tau1;
real<lower=tau1,upper=1.0*n> tau2;
real<lower=0> sigma;
}
transformed parameters {
real alpha1;
real alpha2;
real<lower=0> inv_sigma2;
alpha1 = b1;
alpha2 = b1 + b2;
inv_sigma2 = inv_square(sigma);
}
model {
target += normal_lpdf(y | f(t,b0,b1,b2,tau1,tau2,n), sigma);
target += gamma_lpdf(tau1 | a_tau1, b_tau1);
target += gamma_lpdf(tau2 | a_tau2, b_tau2);
target += normal_lpdf(b0 | b0_mu, b0_sd);
target += normal_lpdf(b1 | b1_mu, b1_sd);
target += normal_lpdf(b2 | b2_mu, b2_sd);
target += gamma_lpdf(inv_sigma2 | a_sigma, b_sigma);
}

",
"stan_bcable_3phase.stan")


#######################################################
#######################################################
# R Functions
#######################################################
#######################################################

# f for the segmented bent-cable model
fcable <- function(t, b0, b1, b2, b3, b4, gam1, tau1, tau2, gam2) {
  ph2 <- ifelse(t > gam1 & t <= tau1, b2 * (t - gam1), 0)
  ph3 <- ifelse(t > tau1 & t <= tau2, b2 * ((t - tau1) ^ 2 + 
    (tau1 - gam1)), 0)
  ph4 <- ifelse(t > tau2 & t <= gam2, b3 * ((t - tau2) + 
    b2 * ((tau2 - tau1) ^ 2 + (tau1 - gam1)) / b3), 0)
  ph5 <- ifelse(t > gam2, b4 * ((t - gam2) + b3 * (gam2 - tau2) / b4 + 
    b2 * ((tau2 - tau1) ^ 2 + (tau1 - gam1)) / b4), 0)
  f <- b0 + b1 * t + (ph2 + ph3 + ph4 + ph5)
  f
}

#######################################################

# f for the bent-cable model
bfcable <- function(t, b0, b1, b2, tau1, tau2) {
  ph2 <- ifelse(t > tau1 & t <= tau2, (t - tau1) ^ 2 / 
    (2 * (tau2 - tau1)), 0)
  ph3 <- ifelse(t > tau2, t - (tau1 + tau2) / 2, 0)
  f <- b0 + b1 * t + b2 * (ph2 + ph3)
  f
}

#######################################################

# Likelihood function
fcable.llik <- function(par, y, t, phase) {
  b0 <- par[1]
  b1 <- par[2]
  b2 <- par[3]
  if (phase == 5) {
    b3 <- par[4]
    b4 <- par[5]
    gam1 <- par[6]
    tau1 <- par[7]
    tau2 <- par[8]
    gam2 <- par[9]
    sigma <- par[10]
    f <- fcable(t, b0, b1, b2, b3, b4, gam1, tau1, tau2, gam2)
  }
  if (phase == 3) {
    tau1 <- par[4]
    tau2 <- par[5]
    sigma <- par[6]
    f <- bfcable(t, b0, b1, b2, tau1, tau2)
  }
  llik <- dnorm(y, mean = f, sd = sigma, log = TRUE)
  return(llik)
}

#######################################################

# Calculate WAIC
waic <- function(fit, y, t, phase) {
  if (phase == 5) {
    mcmc.sample <- as.matrix(fit)[, c("b0", "b1", "b2", "b3", "b4", 
        "gam1", "tau1", "tau2", "gam2", "sigma")]
  }
  if (phase == 3) {
    mcmc.sample <-
      as.matrix(fit)[, c("b0", "b1", "b2", "tau1", "tau2", "sigma")]
  }
  llik <- apply(mcmc.sample, 1, fcable.llik, y = y, 
    t = t, phase = phase)
  pwaic <- sum(apply(llik, 1, var, na.rm = TRUE), na.rm = TRUE)
  lppd <-  sum(log(rowSums(exp(llik), na.rm = TRUE) / 
    dim(llik)[2]), na.rm = TRUE)
  w <- -2 * (lppd - pwaic)
  return(w)
}

#######################################################

# Calculate DIC

dic <- function(fit, y, t, phase) {
  if (phase == 5) {
    mcmc.sample <- as.matrix(fit)[, c("b0", "b1", "b2", "b3", "b4",
         "gam1", "tau1", "tau2", "gam2", "sigma")]
    p.mean <- rstan::summary(fit)$summary[c("b0", "b1", "b2", "b3", "b4",
         "gam1", "tau1", "tau2", "gam2", "sigma"), 1]
  }
  if (phase == 3) {
    mcmc.sample <- as.matrix(fit)[, c("b0", "b1", "b2", "tau1", "tau2", "sigma")]
    p.mean <- rstan::summary(fit)$summary[c("b0", "b1", "b2", "tau1", "tau2", "sigma"), 1]
  }
  llik0 <- fcable.llik(par = p.mean, y = y, t = t, phase = phase)
  llik <- apply(mcmc.sample, 1, fcable.llik, y = y, t = t, phase = phase)
  llik.bar <- mean(colSums(llik, na.rm = TRUE))
  llik.hat <- sum(llik0, na.rm = TRUE)
  p.dic <- 2 * (llik.hat - llik.bar)
  d <- -2 * (llik.hat - p.dic)
  dic.results <- cbind(llik.hat, p.dic, d)
  colnames(dic.results) <- c("lpd(posterior mean)", "p", "DIC")
  return(dic.results)
}


#######################################################
#######################################################
# Example: Segmented Bent-Cable Fit to the Housing data
#######################################################
#######################################################

library(rstan)
library(Quandl)                    # If Quandl does not work, you can download data
housing <- Quandl("FRED/COMPUTSA") # from https://fred.stlouisfed.org/series/COMPUTSA.
dat1 <- housing[190:1, ]           # We consider time-series data (in 100,000’s of units) from
y <- dat1[, 2] / 100               # March of 2006 to December of 2021. 
y

#  [1] 22.45 20.71 18.97 20.50 19.34 18.77 20.11 19.18 18.93 18.88 18.22 16.40
# [13] 16.23 15.39 15.36 14.81 15.34 15.12 13.56 14.05 13.90 13.28 13.31 12.74
# [25] 11.95 10.22 11.42 11.42 10.87 10.17 11.60 10.55 10.76 10.21  7.77  8.19
# [37]  8.39  8.46  8.18  7.97  7.97  7.92  7.21  7.46  8.44  7.50  6.89  6.70
# [49]  6.35  7.37  7.02  8.94  5.72  5.92  6.32  6.05  5.52  5.65  5.20  6.15
# [61]  5.91  5.49  5.49  5.79  6.34  6.14  6.03  5.66  5.85  6.10  5.45  5.61
# [73]  5.82  6.71  6.26  6.25  6.76  6.84  6.52  7.30  6.68  6.77  7.14  7.29
# [85]  8.16  7.12  7.08  7.36  7.84  7.65  7.76  8.01  8.39  7.78  8.38  8.76
# [97]  8.93  8.34  8.91  8.00  8.39  9.08  9.75  9.20  8.63  9.55  9.54  8.63
#[109]  7.81 10.01 10.12  9.62  9.88  9.74 10.40  9.85 10.00 10.23 10.56 10.37
#[121] 10.03  9.45 10.00 11.09 10.72 10.54 10.35 10.70 12.44 11.01 10.73 11.52
#[133] 11.62 10.66 11.70 12.38 11.74 11.00 11.08 11.91 11.80 12.03 12.04 12.90
#[145] 11.89 12.13 12.63 12.16 11.69 12.46 11.74 11.16 11.42 10.52 12.46 13.37
#[157] 13.16 13.12 12.35 11.70 12.42 12.72 11.36 12.79 12.52 13.08 12.88 12.94
#[169] 12.67 11.91 11.78 12.43 13.40 12.16 14.26 13.56 12.44 13.86 13.28 13.47
#[181] 14.97 14.17 13.50 13.12 13.80 12.91 12.35 12.55 14.18 12.95

n <- length(y)
t <- 1:n

# 5-phase analysis
data.stan <- list(n = n, t = t, y = y, b0_mu = max(y), b0_sd = 1, b1_mu = 0,
    b1_sd = 1, b2_mu = 0, b2_sd = 1, b3_mu = 0, b3_sd = 1, b4_mu = 0, b4_sd = 1,
    a_gam1 = 1, b_gam1 = 0.01, a_tau1 = 1, b_tau1 = 0.01, a_tau2 = 1,
    b_tau2 = 0.005, a_gam2 = 1, b_gam2 = 0.001, a_sigma = 0.01,
    b_sigma = 0.01, lower_gam1 = 1)
chains <- 2
inits <- rep(list(list(b0 = max(y), b1 = -0.8, b2 = -0.2, b3 = 0.7,
      b4 = 0.5, gam1 = 10, tau1 = 40, tau2 = 90, gam2 = 130, 
      sigma = 1)), chains)
par.stan <- c("b0", "b1", "b2", "b3", "b4", "alpha1", "alpha2",
    "alpha3", "alpha4", "incoming_diff", "outgoing_diff",
    "gam1", "tau1", "tau2", "gam2", "sigma")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
stanc("stan_bcable_5phase.stan")
bcable_model <- stan_model("stan_bcable_5phase.stan")
fit <- rstan::sampling(object = bcable_model, data = data.stan,
  init = inits, pars = par.stan, warmup = 2000, iter = 5000,
  chains = chains, thin = 1, 
  control = list(adapt_delta = 0.95, max_treedepth = 15))
fit

#Inference for Stan model: anon_model.
#2 chains, each with iter=5000; warmup=2000; thin=1; 
#post-warmup draws per chain=3000, total post-warmup draws=6000.
#
#                 mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
#b0              21.91    0.01 0.28   21.36   21.72   21.91   22.10   22.47  1613    1
#b1              -0.41    0.00 0.03   -0.47   -0.43   -0.41   -0.39   -0.37  1258    1
#b2               0.00    0.00 0.00    0.00    0.00    0.00    0.00    0.00  4909    1
#b3               0.50    0.00 0.03    0.45    0.48    0.49    0.51    0.56  1261    1
#b4               0.45    0.00 0.03    0.40    0.43    0.45    0.47    0.51  1254    1
#alpha1          -0.41    0.00 0.03   -0.47   -0.43   -0.41   -0.39   -0.37  1258    1
#alpha2          -0.41    0.00 0.03   -0.47   -0.42   -0.41   -0.39   -0.36  1259    1
#alpha3           0.08    0.00 0.01    0.07    0.08    0.08    0.09    0.09  4938    1
#alpha4           0.04    0.00 0.01    0.02    0.03    0.04    0.04    0.05  5182    1
#incoming_diff    0.00    0.00 0.00    0.00    0.00    0.00    0.00    0.00  4909    1
#outgoing_diff    0.05    0.00 0.01    0.03    0.04    0.05    0.05    0.06  5899    1
#gam1             7.51    0.09 4.50    1.26    3.76    6.85   10.65   17.35  2613    1
#tau1            14.06    0.11 4.16    5.03   11.39   14.28   16.90   21.68  1443    1
#tau2            75.48    0.15 7.52   62.58   69.14   75.44   81.87   88.58  2481    1
#gam2           134.81    0.09 5.78  122.79  131.23  134.98  138.60  145.47  3859    1
#sigma            0.64    0.00 0.03    0.58    0.62    0.64    0.66    0.71  3252    1
#lp__          -205.43    0.06 2.36 -210.88 -206.77 -205.09 -203.70 -201.88  1408    1
#
#Samples were drawn using NUTS(diag_e) at Wed May 29 10:24:36 2024.
#For each parameter, n_eff is a crude measure of effective sample size,
#and Rhat is the potential scale reduction factor on split chains (at 
#convergence, Rhat=1).

# Find WAIC
waic(fit = fit, y = y, t = t, phase = 5)
#[1] 374.7714

# Find DIC
dic(fit=fit,y=y,t=t,phase=5)
#     lpd(posterior mean)         p      DIC
#[1,]           -186.8349 -8.419979 356.8298

# Pairs Plots
pairs(fit,pars = c("b0", "alpha1", "alpha2", "alpha3",
    "alpha4", "gam1", "tau1", "tau2", "gam2", "sigma"))

# Trace plots
stan_trace(fit, pars = c("b0", "alpha1", "alpha2", "alpha3",
    "alpha4", "gam1", "tau1", "tau2", "gam2", "sigma"))

# Convert numbers into dates
sum.fit <- rstan::summary(fit, pars = c("b0", "b1", "b2",
    "b3", "b4", "gam1", "tau1", "tau2",
    "gam2", "sigma"))$summary
start.date <- "2006-03-01"
end.date <- "2021-12-01"
data.frame(gam1 = c(as.Date(ceiling(sum.fit[6, 1] / 12 * 365), origin = start.date),
    as.Date(ceiling(sum.fit[6, 4] / 12 * 365), origin = start.date),
    as.Date(ceiling(sum.fit[6, 8] / 12 * 365), origin = start.date)),
  tau1 = c(as.Date(ceiling(sum.fit[7, 1] / 12 * 365), origin = start.date),
    as.Date(ceiling(sum.fit[7, 4] / 12 * 365), origin = start.date),
    as.Date(ceiling(sum.fit[7, 8] / 12 * 365), origin = start.date)),
  tau2 = c(as.Date(ceiling(sum.fit[8, 1] / 12 * 365), origin = start.date),
    as.Date(ceiling(sum.fit[8, 4] / 12 * 365), origin = start.date),
    as.Date(ceiling(sum.fit[8, 8] / 12 * 365), origin = start.date)),
  gam2 = c(as.Date(ceiling(sum.fit[9, 1] / 12 * 365), origin = start.date),
    as.Date(ceiling(sum.fit[9, 4] / 12 * 365), origin = start.date),
    as.Date(ceiling(sum.fit[9, 8] / 12 * 365), origin = start.date)))

        gam1       tau1       tau2       gam2
1 2006-10-16 2007-05-03 2012-06-13 2017-05-23
2 2006-04-09 2006-08-02 2011-05-18 2016-05-22
3 2007-08-11 2007-12-21 2013-07-17 2018-04-12

# Plot fitted curve
b0 <- sum.fit[1, 1]
b1 <- sum.fit[2, 1]
b2 <- sum.fit[3, 1]
b3 <- sum.fit[4, 1]
b4 <- sum.fit[5, 1]
gam1 <- sum.fit[6, 1]
tau1 <- sum.fit[7, 1]
tau2 <- sum.fit[8, 1]
gam2 <- sum.fit[9, 1]
start.y <- as.numeric(substr(start.date, 1, 4))
start.m <- as.numeric(substr(start.date, 6, 7))
end.y <- as.numeric(substr(end.date, 1, 4))
end.m <- as.numeric(substr(end.date, 6, 7))
f <- fcable(t, b0, b1, b2, b3, b4, gam1, tau1, tau2, gam2)
plot(ts(y, frequency = 12, start = c(start.y, start.m),
    end = c(end.y, end.m)), type = "p",
    ylab = "", xlab = "", col = "grey")
lines(ts(f, frequency = 12, start = c(start.y, start.m),
    end = c(end.y, end.m)), type = "l",
    ylab = "", xlab = "")
abline(v = (start.y + start.m / 12 + gam1 / 12), lty = 2)
abline(v = (start.y + start.m / 12 + tau1 / 12), lty = 2)
abline(v = (start.y + start.m / 12 + tau2 / 12), lty = 2)
abline(v = (start.y + start.m / 12 + gam2 / 12), lty = 2)

#######################################################
#######################################################
# Example: Bent-Cable Fit to the Housing data (3-Phase)
#######################################################
#######################################################

# 3-phase bent-cable
data.stan <- list(n = n, t = t, y = y, b0_mu = max(y), b0_sd = 1,
    b1_mu = 0, b1_sd = 1, b2_mu = 0, b2_sd = 1, a_tau1 = 1,
    b_tau1 = 0.005, a_tau2 = 1, b_tau2 = 0.005, a_sigma = 0.01,
    b_sigma = 0.01, lower_tau1 = 1)
chains <- 2
inits <- rep(list(list(b0 = max(y), b1 = -0.8, b2 = 1, tau1 = 40,
    tau2 = 90, sigma = 1)), chains)
par.stan <- c("b0", "b1", "b2", "alpha1", "alpha2", "tau1", "tau2", "sigma")
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
stanc("stan_bcable_3phase.stan")
bcable_model3 <- stan_model("stan_bcable_3phase.stan")
fit <- rstan::sampling(object = bcable_model3, data = data.stan,
  init = inits, pars = par.stan, warmup = 2000, iter = 5000, chains = chains, 
  thin = 1, control = list(adapt_delta = 0.95, max_treedepth = 15))
fit

#Inference for Stan model: anon_model.
#2 chains, each with iter=5000; warmup=2000; thin=1; 
#post-warmup draws per chain=3000, total post-warmup draws=6000.
#
#          mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
#b0       21.99    0.01 0.35   21.34   21.75   21.98   22.22   22.69  1226    1
#b1       -0.42    0.00 0.03   -0.50   -0.44   -0.41   -0.39   -0.36   891    1
#b2        0.48    0.00 0.04    0.43    0.46    0.48    0.50    0.57   885    1
#alpha1   -0.42    0.00 0.03   -0.50   -0.44   -0.41   -0.39   -0.36   891    1
#alpha2    0.07    0.00 0.00    0.06    0.06    0.07    0.07    0.07  5525    1
#tau1     15.49    0.18 5.41    3.44   11.95   16.05   19.43   24.85   930    1
#tau2     67.83    0.04 1.72   64.35   66.67   67.89   68.99   71.07  2389    1
#sigma     0.74    0.00 0.04    0.67    0.71    0.74    0.76    0.82  1990    1
#lp__   -223.56    0.05 1.90 -228.09 -224.62 -223.17 -222.14 -220.92  1326    1
#
#Samples were drawn using NUTS(diag_e) at Wed May 29 10:44:47 2024.
#For each parameter, n_eff is a crude measure of effective sample size,
#and Rhat is the potential scale reduction factor on split chains (at 
#convergence, Rhat=1).

# WAIC
waic(fit = fit, y = y, t = t, phase = 3)
# [1] 426.5846

# DIC
dic(fit=fit,y=y,t=t,phase=3)
#     lpd(posterior mean)        p      DIC
#[1,]           -208.1755 3.631975 423.6149


# Pairs plots
pairs(fit, pars = c("b0", "alpha1", "alpha2", "tau1", "tau2", "sigma"))

# Trace plots
stan_trace(fit, pars = c("b0", "alpha1", "alpha2", "tau1", "tau2", "sigma"))

# Convert numbers into dates
sum.fit <- rstan::summary(fit, pars = c("b0", "b1", "b2", "tau1", "tau2", "sigma"))$summary
start.date <- "2006-03-01"
end.date <- "2021-12-01"
data.frame(tau1 = c(as.Date(ceiling(sum.fit[4, 1] / 12 * 365), origin = start.date),
    as.Date(ceiling(sum.fit[4, 4] / 12 * 365), origin = start.date),
    as.Date(ceiling(sum.fit[4, 8] / 12 * 365), origin = start.date)),
  tau2 = c(as.Date(ceiling(sum.fit[5, 1] / 12 * 365), origin = start.date),
    as.Date(ceiling(sum.fit[5, 4] / 12 * 365), origin = start.date),
    as.Date(ceiling(sum.fit[5, 8] / 12 * 365), origin = start.date)))

        tau1       tau2
1 2007-06-16 2011-10-25
2 2006-06-14 2011-07-11
3 2008-03-26 2012-01-31

# Plot fitted curve
b0 <- sum.fit[1, 1]
b1 <- sum.fit[2, 1]
b2 <- sum.fit[3, 1]
tau1 <- sum.fit[4, 1]
tau2 <- sum.fit[5, 1]
start.y <- as.numeric(substr(start.date, 1, 4))
start.m <- as.numeric(substr(start.date, 6, 7))
end.y <- as.numeric(substr(end.date, 1, 4))
end.m <- as.numeric(substr(end.date, 6, 7))
f <- bfcable(t, b0, b1, b2, tau1, tau2)
plot(ts(y, frequency = 12, start = c(start.y, start.m),
    end = c(end.y, end.m)), type = "p",
    ylab = "", xlab = "", col = "grey")
lines(ts(f, frequency = 12, start = c(start.y, start.m),
    end = c(end.y, end.m)), type = "l",
    ylab = "", xlab = "")
abline(v = (start.y + start.m / 12 + tau1 / 12), lty = 2)
abline(v = (start.y + start.m / 12 + tau2 / 12), lty = 2)









