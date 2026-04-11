library(cmna)
library(actuar)
library(univariateML)

# Step 1: Preparation of Lifetime-Data

t <- c(0.35, 0.59, 0.96, 0.99, 1.69,
       1.97, 2.07, 2.58, 2.71, 2.90,
       3.67, 3.99, 5.35, 13.77, 25.50)

# Step 2: Preparation of Functions

H <- function(beta){
  y <- (length(t)*sum((t^(-beta))*(log(t))))/(sum(t^(-beta)))
}

G <- function(beta){
  y <- (length(t)/beta)-sum(log(t))+H(beta)
}

# Step 3: Computation of beta and alpha

beta_fin  <- bisection(G,0.000001,5.000000,tol=0.000001)

alpha_fin <- ((1/length(t))*sum(t^(-beta_fin)))^(1/beta_fin)

test      <- mlinvweibull(t)
test2     <- mlweibull(t)

shape_par <- 1.0275
scale_par <- 1/0.6822

shape_wei <- 0.889
scale_wei <- 4.292

# Step 4: Preparation of Plotting

x2     <- seq(0, 30, length = 301)

eq1    <- function(x){
  (beta_fin/(alpha_fin^(beta_fin)))*(x^(-beta_fin-1))*exp(-(alpha_fin*x)^(-beta_fin))
}

eq2    <- function(x){
  (shape_par/x)*((scale_par/x)^(shape_par))*exp(-(scale_par/x)^(shape_par))
}

eq_wei <- function(x){
  (shape_wei/scale_wei)*(x/scale_wei)^(shape_wei-1)*exp(-(x/scale_wei)^(shape_wei))
}

y2     <- eq1(x2)
y2[1]  <- 0

y3     <- eq2(x2)
y3[1]  <- 0

y4     <- eq_wei(x2)
y4[1]  <- 0

hist(t,prob = TRUE, col = "white",
     breaks = seq(0,30,l=31), ylim = c(0,0.4),
     main = "Example 3: Histogram and Estimated PDF",
     xlab = "Time to Breakdown",
     ylab = "Relative Frequencies")
lines(x2,y2,col="darkorange3",lty="solid",lwd="3")
lines(x2,y3,col="darkorange3",lty="solid",lwd="3")
lines(x2,y4,col="green4",lty="solid",lwd="3")

# Step 5: R^2 Measure for 'Goodness-of-fit'-Measure for Inverse Weibull

x_hist        <- c(0.5, 1.5, 2.5, 3.5, 5.5, 13.5, 25.5)
y_hist        <- c(4/15, 2/15, 4/15, 2/15, 1/15, 1/15, 1/15)
y_hist_mean   <- 1/7
zaehler_iwei  <- sum((y_hist - eq2(x_hist))^2)
nenner_iwei   <- sum((y_hist - y_hist_mean)^2)
gof_meas_iwei <- zaehler_iwei/nenner_iwei

# Step 6: R^2 Measure for 'Goodness-of-fit'-Measure for Weibull

zaehler_wei   <- sum((y_hist - eq_wei(x_hist))^2)
nenner_wei    <- sum((y_hist - y_hist_mean)^2)
gof_meas_wei  <- zaehler_wei/nenner_wei