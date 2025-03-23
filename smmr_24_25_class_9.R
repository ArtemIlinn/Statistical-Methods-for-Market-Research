####################################################################
###### SMMR by Ilya S. Slabolitskiy and Ekaterina A. Simonova ######
############################# Class 9 ##############################
####################################################################



# Markov chain simulation
S <- 1:3 # state space
P <- matrix(data = c(1/2, 1/4, 1/4,
                     1/3, 0, 2/3,
                     1/2, 1/2, 0),
            nrow = 3, ncol = 3, byrow = T) # transition matrix
X <- sample(x = S, size = 1)
for (t in 2:1000){
  X[t] <- sample(x = S, size = 1, prob = P[X[t - 1], ])
}
plot(x = X[1:50], type = 'l', xlab = '', ylab = 'time')

# test for homogeneity
P1_est <- matrix(0, nrow = 3, ncol = 3)
for (t in 2:500){
  if ((X[t] == 1) & (X[t - 1] == 1)){
    P1_est[1, 1] <- P1_est[1, 1] + 1
  } else if ((X[t] == 2) & (X[t - 1] == 1)){
    P1_est[1, 2] <- P1_est[1, 2] + 1
  } else if ((X[t] == 3) & (X[t - 1] == 1)){
    P1_est[1, 3] <- P1_est[1, 3] + 1
  } else if ((X[t] == 1) & (X[t - 1] == 2)){
    P1_est[2, 1] <- P1_est[2, 1] + 1
  } else if ((X[t] == 2) & (X[t - 1] == 2)){
    P1_est[2, 2] <- P1_est[2, 2] + 1
  } else if ((X[t] == 3) & (X[t - 1] == 2)){
    P1_est[2, 3] <- P1_est[2, 3] + 1
  } else if ((X[t] == 1) & (X[t - 1] == 3)){
    P1_est[3, 1] <- P1_est[3, 1] + 1
  } else if ((X[t] == 2) & (X[t - 1] == 3)){
    P1_est[3, 2] <- P1_est[3, 2] + 1
  } else{
    P1_est[3, 3] <- P1_est[3, 3] + 1
  }
}
P1_est

P2_est <- matrix(0, nrow = 3, ncol = 3)
for (t in 501:1000){
  if ((X[t] == 1) & (X[t - 1] == 1)){
    P2_est[1, 1] <- P2_est[1, 1] + 1
  } else if ((X[t] == 2) & (X[t - 1] == 1)){
    P2_est[1, 2] <- P2_est[1, 2] + 1
  } else if ((X[t] == 3) & (X[t - 1] == 1)){
    P2_est[1, 3] <- P2_est[1, 3] + 1
  } else if ((X[t] == 1) & (X[t - 1] == 2)){
    P2_est[2, 1] <- P2_est[2, 1] + 1
  } else if ((X[t] == 2) & (X[t - 1] == 2)){
    P2_est[2, 2] <- P2_est[2, 2] + 1
  } else if ((X[t] == 3) & (X[t - 1] == 2)){
    P2_est[2, 3] <- P2_est[2, 3] + 1
  } else if ((X[t] == 1) & (X[t - 1] == 3)){
    P2_est[3, 1] <- P2_est[3, 1] + 1
  } else if ((X[t] == 2) & (X[t - 1] == 3)){
    P2_est[3, 2] <- P2_est[3, 2] + 1
  } else{
    P2_est[3, 3] <- P2_est[3, 3] + 1
  }
}
P2_est