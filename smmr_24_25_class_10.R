####################################################################
###### SMMR by Ilya S. Slabolitskiy and Ekaterina A. Simonova ######
############################ Class 10 ##############################
####################################################################



# install.packages('GA')
library('GA') # genetic algorithm


# task 1 (from the lecture). MSM for a (simple) sample
n_obs <- 2000

# Markov chain simulation
S <- 1:2
P <- matrix(data = c(0.95, 0.05,
                     0.15, 0.85),
            nrow = 2,
            ncol = 2,
            byrow = T)
M <- sample(x = S, size = 1)
for (i in 2:n_obs){
  M[i] <- sample(x = S, size = 1, prob = P[M[i - 1], ])
}

# sample simulation
mu_1 <- 2
mu_2 <- -2
sigma_1 <- 1
sigma_2 <- 2
X <- c()
for (i in 1:n_obs){
  if (M[i] == 1){
    X[i] <- rnorm(n = 1, mean = mu_1, sd = sigma_1)
  } else{
    X[i] <- rnorm(n = 1, mean = mu_2, sd = sigma_2)
  }
}

# log-likelihood function
logL_SMM <- function(par){
  mu1 <- par[1]
  mu2 <- par[2]
  sigma1 <- par[3]
  sigma2 <- par[4]
  p11 <- par[5]
  p22 <- par[6]
  
  P_t <- matrix(data = c(p11, 1 - p11,
                         1 - p22, p22),
                nrow = 2,
                ncol = 2,
                byrow = T)
  
  xi <- c(0.5, 0.5)
  
  xi_all <- c()
  L <- c()
  for(i in 1:n_obs) {
    f1 <- dnorm(x = X[i], mean = mu1, sd = sigma1)
    f2 <- dnorm(x = X[i], mean = mu2, sd = sigma2)
    f <- c(f1, f2)
    xi_next <- t(P_t) %*% xi
    L_t <- t(f) %*% xi_next
    L <- c(L, L_t)
    xi <- f * xi_next / as.vector(L_t)
    xi_all <- rbind(xi_all, t(xi))
  }
  log_lik <- sum(log(L))
  list(xi_all = xi_all, log_lik = log_lik)
}

# constraints for constraint optimization (A * par >= b)
A <- matrix(c(0, 0, 1, 0, 0, 0,
              0, 0, 0, 1, 0, 0,
              0, 0, 0, 0, 1, 0,
              0, 0, 0, 0, 0, 1,
              0, 0, 0, 0, -1, 0,
              0, 0, 0, 0, 0, -1),
            nrow = 6,
            ncol = 6,
            byrow = T)
b <- c(0.01, 0.01, 0.01, 0.01, -0.99, -0.99)

logL_SMM_log_lik <- function(par){
  -logL_SMM(par)$log_lik
}

# optimization and results
logL_SMM_optim <- constrOptim(theta = c(0, 0, 1, 1, 0.5, 0.5),
                              f = logL_SMM_log_lik, grad = NULL,
                              ui = A, ci = b)
logL_SMM_optim$par
logL_SMM(par = logL_SMM_optim$par)$xi_all

# transition matrix estimate
P_est <- matrix(data = c(logL_SMM_optim$par[5], 1 - logL_SMM_optim$par[5],
                         1 - logL_SMM_optim$par[6], logL_SMM_optim$par[6]),
                nrow = 2,
                ncol = 2,
                byrow = T)
# stationary distribution of Markov chain
pi_stationary <- c(0.5, 0.5)
for (j in 1:100){
  pi_stationary <- pi_stationary %*% P_est
}


# task 2. MSM for an AR(1) process
T_obs <- 1000

# Markov chain simulation
S <- 1:2
P <- matrix(data = c(0.95, 0.05,
                     0.15, 0.85),
            nrow = 2,
            ncol = 2,
            byrow = T)
M <- sample(x = S, size = 1)
for (t in 2:T_obs){
  M[t] <- sample(x = S, size = 1, prob = P[M[t - 1], ])
}

# AR(1) process simulation
alpha_1 <- 0.1
alpha_2 <- 0.9
beta_1 <- 0.8
beta_2 <- 0.3
sigma_1 <- 2
sigma_2 <- 5
X <- 0
for (t in 2:T_obs){
  if (M[t] == 1){
    X[t] <- alpha_1 + beta_1 * X[t - 1] + rnorm(n = 1, mean = 0, sd = sigma_1)
  } else{
    X[t] <- alpha_2 + beta_2 * X[t - 1] + rnorm(n = 1, mean = 0, sd = sigma_2)
  }
}
plot(x = X[1:100], type = 'l')

# log-likelihood function
LogL_SMMAR <- function(par){
  alpha1 <- par[1]
  alpha2 <- par[2]
  beta1 <- par[3]
  beta2 <- par[4]
  sigma1 <- par[5]
  sigma2 <- par[6]
  p11 <- par[7]
  p22 <- par[8]
  
  P_t <- matrix(data = c(p11, 1 - p11,
                         1 - p22, p22),
                nrow = 2,
                ncol = 2,
                byrow = T)
  
  xi <- c(0.5, 0.5)
  
  xi_all <- c()
  L <- c()
  for(t in 1:(T_obs - 1)) {
    f1 <- dnorm(x = X[t + 1], mean = alpha1 + beta1 * X[t], sd = sigma1)
    f2 <- dnorm(x = X[t + 1], mean = alpha2 + beta2 * X[t], sd = sigma2)
    f <- c(f1, f2)
    xi_next <- t(P_t) %*% xi
    L_t <- t(f) %*% xi_next
    L <- c(L, L_t)
    xi <- (f * xi_next) / as.vector(L_t)
    xi_all <- rbind(xi_all, t(xi))
  }
  log_lik <- sum(log(L))
  return(log_lik)
}

# optimization and results
LogL_SMMAR_optim <- ga(type = 'real-valued',
                       fitness = LogL_SMMAR,
                       lower = c(-10, -10, -0.99, -0.99, 0.01, 0.01, 0.01, 0.01),
                       upper = c(10, 10, 0.99, 0.99, 10, 10, 0.99, 0.99),
                       popSize = 150, maxiter = 200, run = 100,
                       optim = TRUE,
                       optimArgs = list(method = 'Nelder-Mead', 
                                        poptim = 0.05,
                                        pressel = 0.5,
                                        control = list(maxit = 10000,
                                                       abstol = 1e-10,
                                                       fnscale = -1)))
summary(LogL_SMMAR_optim)
