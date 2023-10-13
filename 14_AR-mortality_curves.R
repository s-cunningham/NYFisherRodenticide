
## Attempting to learn something about mortality from AR-age curves
## 2023-09-14

#### Elie's inital experimenting ####

library(magrittr)

# Exposure occurs as a negative binomial process with rate $\gamma$, so $P(C_{i,t+1} = k) = NegBinom(1 - \gamma)$, where $C_{i,t}$ refers to the number of compounds individual $i$ is exposed to by time $t$.  Just to give an idea, with a probability of annual contact $\gamma = 0.2$
par(mfrow=c(1,1))
par(bty = "l", mar = c(3,3,1,1), bty = "l", mgp= c(1,.2,0), 
    tck = 0.01, cex.axis = .8)
gamma <- 0.5
plot(0:6, 
     dnbinom(0:6, prob = 1-gamma, size = 1), type = "h", lwd = 5, ylab = "probability",
     xlab = "number of exposures per year")

# Mortality depends on exposure.  INTERESTINGLY, it also has to depend on an interaction between exposure and age(!)  I have a relatively simple two parameter model, but the AGE bit should probably be included as another parameter.
par(bty = "l", mar = c(3,3,1,1), bty = "l", mgp= c(1,.2,0), tck = 0.01, cex.axis = .8)
age.max <- 10; c.max <- 5
mu <- 0.2; beta <- 0.2

expit <- function(x) exp(x)/(1+exp(x))
exposure.matrix <- matrix(0:(c.max-1), nrow = age.max, ncol = c.max, byrow = TRUE)
age.matrix <- matrix(0:(age.max-1), nrow = age.max, ncol = c.max, byrow = FALSE)
mortality.matrix <- expit( log(mu) + beta * exposure.matrix * age.matrix)

matplot(mortality.matrix, type = "l", pch = 1, col = 1:5, lty = 1, lwd = 2, xlab = "Age", ylab = "probability of mortality", ylim = c(0,1))
legend("topleft", col = 1:c.max, legend = 1:c.max-1, lwd = 2, title = "RCs")


# In this model, animals with no exposure have a constant annual mortality of 10% through their adult life (actually not a terrible assumption). 
# But the more rodenticide they've been exposed to, the higher the mortality rate at higher ages. 

## Set up simulation
# basic paramaters
T <- seq(0, 10)
c.max <- 5
age.max <- 9
births <- 100

# rate of exposure
# gamma <- .75
gamma <- .83

# parameters of mortality 
beta <- 0.2
mu <- .32

expit <- function(x) exp(x)/(1+exp(x))
exposure.matrix <- matrix(0:(c.max-1), nrow = age.max, ncol = c.max, byrow = TRUE)
age.matrix <- matrix(0:(age.max-1), nrow = age.max, ncol = c.max, byrow = FALSE)
mortality.matrix <- expit( log(mu) + beta * exposure.matrix * age.matrix)

# construct exposure transition probability matrix
gamma.vector <- dnbinom(1:c.max, 1, prob = 1 - gamma)
gamma.vector <- gamma.vector/(sum(gamma.vector))
gamma.matrix <- matrix(0, nrow = c.max, ncol = c.max)
for(i in 1:c.max) gamma.matrix[i,i:c.max] <- gamma.vector[1:(c.max - i  + 1)]
gamma.matrix <- gamma.matrix %>% apply(1, function(x) x/sum(x)) %>% t

#
M.total <- array(0, dim= c(length(T), age.max, c.max))
M.total <- provideDimnames(M.total, sep = "", base = list("t", "a", "c"))
M.total[1,1,1] <- births

i <- 2

for(i in 2:length(T)){
  M.old <- M.total[i-1,,]
  
  # First expose them (with transition matrix)    
  M.exposed <- apply(M.old, 1, '%*%', gamma.matrix) %>% t
  
  # next, age + die
  M.aged <- M.exposed
  M.aged[2:age.max,] <- M.exposed[1:(age.max-1),] * (1-mortality.matrix[-1,])
  
  # put into array
  M.total[i,,] <- M.aged
  
  # reset births
  M.total[i,1,] <- c(births, rep(0, c.max-1))
}

M.final <- M.total[length(T),,]
CompoundsByAge <- apply(M.final, 1, 
                        function(x) sum( (x/sum(x)) * (1:c.max-1)))

par(mfrow = c(2,2), mar = c(3,3,1,1), bty = "l", mgp= c(1,.25,0), 
    tck = 0.01, cex.axis = .8)
plot(T, apply(M.total, 1, sum), type = "l", lty = 1, main = "total population", 
     ylab = "population")

plot(0:(age.max-1), CompoundsByAge, 
     type="l", main="mean number of compounds", xlab = "age")
points(0:(age.max-1), CompoundsByAge, pch=20, cex=1.8)
abline(v=1,col="gray60")

matplot(0:(age.max-1), M.final, type = "l", lty = 1, lwd = 2, 
        col = gplots::rich.colors(c.max), xlab = "age", 
        main = "compounds by age group")

legend("topright", legend = 0:(c.max - 1), 
       col = gplots::rich.colors(c.max), lty = 1, lwd = 2,  ncol = 3, 
       bty = "n")

plot(0:(age.max-1), apply(M.final, 1, sum), type = "h", main = "age distribution", lwd = 2, xlab = "age", ylab = "")


