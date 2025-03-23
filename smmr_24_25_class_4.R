####################################################################
###### SMMR by Ilya S. Slabolitskiy and Ekaterina A. Simonova ######
############################# Class 4 ##############################
####################################################################



y_T_0 <- c(10, 10, 9, 11)
y_T_1 <- c(8, 8, 9, 7)
y_C_0 <- c(11, 12, 13, 13, 11)
y_C_1 <- c(7, 7, 7, 8, 9)

DiD <- mean(y_T_1) - mean(y_T_0) - (mean(y_C_1) - mean(y_C_0))

y <- c(y_T_0, y_T_1, y_C_0, y_C_1)
x <- c(rep(1, length(c(y_T_0, y_T_1))), rep(0, length(c(y_C_0, y_C_1))))
z <- c(rep(0, length(y_T_0)), rep(1, length(y_T_1)),
       rep(0, length(y_C_0)), rep(1, length(y_C_1)))

df <- data.frame('y' = y, 'x' = x, 'z' = z)

model <- lm(formula = y ~ x + z + x * z, data = df)
summary(model)

# BONUS. Estimate this model using fixed effect estimation technique