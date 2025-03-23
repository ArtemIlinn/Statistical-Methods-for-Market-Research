####################################################################
###### SMMR by Ilya S. Slabolitskiy and Ekaterina A. Simonova ######
############################# Class 2 ##############################
####################################################################



# install.packages('mfx')
library('mfx') # margibal effects
# disable scientific notation
options(scipen = 999)


# problem 1 from class 2
# data simulation
n_obs <- 1000
beta_0 <- 2
beta_1 <- -3
epsilon <- rlogis(n = n_obs)
x <- rbinom(n = n_obs, size = 1, prob = 0.75)
y <- as.numeric(beta_0 + beta_1 * x + epsilon > 0)
df <- data.frame('y' = y, 'x' = x)

# logit model estimation
logit_1 <- glm(formula = y ~ x, data = df, family = binomial(link = 'logit'))
logit_1$coefficients

# numerical answer vs analytical answer
n_00 <- sum((y == 0) & (x == 0))
n_01 <- sum((y == 0) & (x == 1))
n_10 <- sum((y == 1) & (x == 0))
n_11 <- sum((y == 1) & (x == 1))
log(n_10 / n_00)
log(n_11 * n_00 / (n_01 * n_10))

# clear environment
rm(list = ls())


# marginar effects and interpretation
# data simulation
n_obs <- 1000
beta_0 <- 2
beta_1 <- -3
beta_2 <- 3
epsilon <- rlogis(n = n_obs)
x_1 <- rexp(n = n_obs, rate = 1)
x_2 <- rbinom(n = n_obs, size = 1, prob = 0.75)
y <- as.numeric(beta_0 + beta_1 * x_1 + beta_2 * x_2 + epsilon > 0)
df <- data.frame('y' = y, 'x_1' = x_1, 'x_2' = x_2)

# logit model estimation
logit_2 <- glm(y ~ x_1 + x_2, data = df, family = binomial(link = 'logit'))
logit_2$coefficients
# marginal effects
logitmfx(formula = y ~ x_1 + x_2, data = df, atmean = FALSE)


# custom log-likelihood function and numerical optimization
LogL <- function(beta)
{
  beta0 <- beta[1]
  beta1 <- beta[2]
  beta2 <- beta[3]
  
  cond0 <- y == 0
  cond1 <- y == 1
  
  y_star <- beta0 + beta1 * x_1 + beta2 * x_2
  log_lik <- sum(log(1 - plogis(q = y_star))[cond0]) +
             sum(log(plogis(q = y_star))[cond1])
  
  return(log_lik)
}

logL_optim <- optim(par = rep(0, 3),
                    fn = LogL,
                    method = 'Nelder-Mead',
                    hessian = TRUE,
                    lower = rep(-Inf, 3),
                    upper = rep(Inf, 3),
                    control = list(maxit = 10000,
                                   abstol = 1e-10,
                                   fnscale = -1))
logL_optim$par
logL_optim$value
logL_optim$hessian
-solve(logL_optim$hessian)