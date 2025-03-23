####################################################################
###### SMMR by Ilya S. Slabolitskiy and Ekaterina A. Simonova ######
############################ Class 13 ##############################
####################################################################



# install.packages('EnvStats')
library('EnvStats') # Pareto distribution
# install.packages('VGAM')
library('VGAM') # extreme value distributions



# problem 1 (based on problem 1, class 12)
n_obs <- 1000
n_sim <- 1000
theta <- 3
pareto_sam_max <- rep(NA, n_sim)
for (m in 1:n_sim)
{
  pareto_sam_aux <- rpareto(n = n_obs,
                            location = theta,
                            shape = 1)
  pareto_sam_max[m] <- max(pareto_sam_aux)
}
c_n <- theta * n_obs
alpha <- 1
frechet_sam <- rfrechet(n = n_obs, location = 0, scale = 1, shape = alpha)
ks.test(pareto_sam_max / c_n, frechet_sam)

hist(x = pareto_sam_max / c_n, breaks = 100, probability = T)
hist(x = frechet_sam, add = TRUE, col = rgb(0, 1, 0, 0.5), breaks = 100, probability = T)

# problem 2 (based on problem 5, class 12)
n_obs <- 100000
n_sim <- 1000
norm_sam_max <- rep(NA, n_sim)
for (m in 1:n_sim)
{
  norm_sam_aux <- rnorm(n = n_obs,
                        mean = 0,
                        sd = 1)
  norm_sam_max[m] <- max(norm_sam_aux)
}
c_n <- (2 * log(n_obs)) ^ (-1 / 2)
d_n <- sqrt(2 * log(n_obs)) - (log(4 * pi) + log(log(n_obs))) / (2 * sqrt(2 * log(n_obs)))
gumbel_sam <- rgumbel(n = n_obs,
                      location = 0,
                      scale = 1)
ks.test((norm_sam_max - d_n) / c_n, gumbel_sam)

hist(x = (norm_sam_max - d_n) / c_n, breaks = 100, probability = T)
hist(x = gumbel_sam, add = TRUE, col = rgb(0, 1, 0, 0.5), breaks = 100, probability = T)
hist(x = norm_sam_max, add = TRUE, col = rgb(0, 0, 1, 0.5), breaks = 10, probability = T)