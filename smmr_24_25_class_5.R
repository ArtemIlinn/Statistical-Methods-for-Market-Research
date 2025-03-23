####################################################################
###### SMMR by Ilya S. Slabolitskiy and Ekaterina A. Simonova ######
############################# Class 5 ##############################
####################################################################



# install.packages('dplyr')
library('dplyr')
# install.packages('rio')
library('rio')


# problem 2 from class 5
# data simulation
n_obs <- 1000
beta_0 <- 2
beta_1 <- -3
x <- rexp(n = n_obs, rate = 1)
epsilon <- rlogis(n = n_obs)
y <- as.numeric(beta_0 + beta_1 * x + epsilon > 0)
df <- data.frame('y' = y, 'x' = x)

# custom log-likelihood function and numerical optimization
LogL <- function(beta, data)
{
  beta0 <- beta[1]
  beta1 <- beta[2]
  
  cond0 <- data$y == 0
  cond1 <- data$y == 1
  
  y_star <- beta0 + beta1 * data$x
  log_lik <- sum(log(1 - plogis(q = y_star))[cond0]) +
             sum(log(plogis(q = y_star))[cond1])
  
  return(log_lik)
}

logit <- optim(par = rep(0, 2),
               fn = LogL,
               method = 'Nelder-Mead',
               hessian = TRUE,
               lower = rep(-Inf, 2),
               upper = rep(Inf, 2),
               control = list(maxit = 10000,
                              abstol = 1e-10,
                              fnscale = -1),
               data = df)
hat_beta_0 <- logit$par[1]
hat_beta_1 <- logit$par[2]
cov_mat <- -solve(logit$hessian)
hat_var_beta_0 <- cov_mat[1, 1]
hat_var_beta_1 <- cov_mat[2, 2]
hat_cov_beta <- cov_mat[1, 2]

# asymptotic variance of hat_beta_0 + hat_beta_1 * mean(x)
as_var <- hat_var_beta_0 + hat_var_beta_1 * mean(x) ^ 2 + 2 * mean(x) * hat_cov_beta

# CI via delta method (dm)
CI_dm_lower_bound <- plogis(hat_beta_0 + hat_beta_1 * mean(x)) -
                     qnorm(p = 0.975) * 
                     sqrt(as_var) * dlogis(hat_beta_0 + hat_beta_1 * mean(x))
CI_dm_upper_bound <- plogis(hat_beta_0 + hat_beta_1 * mean(x)) +
                     qnorm(p = 0.975) * 
                     sqrt(as_var) * dlogis(hat_beta_0 + hat_beta_1 * mean(x))

# CI via monotonic transform (mt)
CI_mt_lower_bound <- plogis(hat_beta_0 + hat_beta_1 * mean(x) - qnorm(p = 0.975) * 
                            sqrt(as_var))
CI_mt_upper_bound <- plogis(hat_beta_0 + hat_beta_1 * mean(x) + qnorm(p = 0.975) * 
                            sqrt(as_var))

# CI via bootstrap (BI)
prob_boot <- c()
for (b in 1:10000){
  df_boot <- sample_n(tbl = df, size = nrow(df), replace = T)
  logit_boot <- optim(par = rep(0, 2),
                      fn = LogL,
                      method = 'Nelder-Mead',
                      hessian = TRUE,
                      lower = rep(-Inf, 2),
                      upper = rep(Inf, 2),
                      control = list(maxit = 10000,
                                     abstol = 1e-10,
                                     fnscale = -1),
                      data = df_boot)
  prob_boot[b] <- plogis(logit_boot$par[1] + logit_boot$par[2] * mean(df_boot$x))
}

BI_lower_bound <- quantile(x = prob_boot, probs = 0.025)
BI_upper_bound <- quantile(x = prob_boot, probs = 0.975)

# comparison
CI <- data.frame('CI_dm' = c(CI_dm_lower_bound, CI_dm_upper_bound),
                 'CI_mt' = c(CI_mt_lower_bound, CI_mt_upper_bound),
                 'BI' = c(BI_lower_bound, BI_upper_bound))
CI

# clear environment
rm(list = ls())

# problem 3 from class 5
df_vk_wise <- import('vk_wise.csv')
hist(df_vk_wise$likes[df_vk_wise$I == TRUE], breaks = 40)
hist(df_vk_wise$likes[df_vk_wise$I == FALSE], breaks = 40)
median(df_vk_wise$likes[df_vk_wise$I == T])
median(df_vk_wise$likes[df_vk_wise$I == F])

# bootstrap interval
diff_median_likes_boot <- c()
for (b in 1:10000){
  df_vk_wise_boot <- sample_n(tbl = df_vk_wise, size = nrow(df_vk_wise),
                              replace = T)
  median_T_boot <- median(df_vk_wise_boot$likes[df_vk_wise_boot$I == T])
  median_F_boot <- median(df_vk_wise_boot$likes[df_vk_wise_boot$I == F])
  diff_median_likes_boot[b] <- median_T_boot - median_F_boot
}
hist(diff_median_likes_boot)
quantile(x = diff_median_likes_boot, probs = 0.025)
quantile(x = diff_median_likes_boot, probs = 0.975)