## Simulating rodenticide analysis

# Reproducible
set.seed(123)

M <- 350 # Number of liver samples (spatial locations)
J <- 3   # Number of years sampled over

# Generate values for covariates that are scaled to range (-1, 1)
wui <- runif(n=M, -1, 1)
baa <- runif(n=M, -1, 1)
bmi <- c(rep(1, 117), rep(2, 111), rep(3, 122))
bmi <- ordered(bmi, c(2,1,3), levels=c(2,1,3))

# Simulating ecological process and its outcome
mean.lambda <- 2
beta0 <- log(mean.lambda)
beta1 <- -2
beta2 <- 1
beta3 <- 1
beta4 <- 2

# Apply the linear model and obtain the logarithm of the expected outcome
