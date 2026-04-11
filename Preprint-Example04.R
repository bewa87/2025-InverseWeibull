library(cmna)
library(actuar)
library(univariateML)

# Step 1: Preparation of Wind Speed Data

# Please fill out your complete folder for the investigated data

A                    <- read.table('.../produkt_ff_stunde_01639.txt', header=TRUE, sep=";")
B                    <- A[,4]
B                    <- B[which(B>=-0.1)]
C                    <- B
C                    <- C[C>=0.05]
#D                    <- C %*% C

# Step 2: Preparation of Functions

H <- function(beta){
  y <- (length(C)*sum((C^(-beta))*(log(C))))/(sum(C^(-beta)))
}

G <- function(beta){
  y <- (length(C)/beta)-sum(log(C))+H(beta)
}

# Step 3: Computation of beta and alpha

beta_fin  <- bisection(G,0.5,2.5,tol=0.00000001)

alpha_fin <- ((1/length(C))*sum(C^(-beta_fin)))^(1/beta_fin)

eq <- function(x){
  beta_fin*((alpha_fin)^(-beta_fin))*x^(-(beta_fin+1))*exp(-(alpha_fin*x)^(-beta_fin))
}

test_wei <- mlweibull(C)

x2    <- seq(0, 35, length = 351)
y2    <- eq(x2)
y2[1] <- 0

# Step 4: Check

AA <- mlinvweibull(C)

eq3 <- function(x){
  1.05*((0.659)^(-1.05))*(x)^(-(1.05+1))*exp(-((0.659)*x)^(-1.05))
}

y4     <- eq3(x2)
y4[1]  <- 0

shape_wei <- 1.51
scale_wei <- 2.93

eq_wei <- function(x){
  (shape_wei/scale_wei)*(x/scale_wei)^(shape_wei-1)*exp(-(x/scale_wei)^(shape_wei))
}

y5     <- eq_wei(x2)
y5[1]  <- 0

hist(C,prob = TRUE, col = "white",
     breaks=seq(min(C),max(C),l=35),
     main="Example 4: Histogram and Estimated PDF",
     xlab="Wind Speed", ylab="Relative Frequencies", 
     xlim = c(0,20), ylim = c(0,0.45))
lines(x2,y2,col="darkorange3",lty="solid",lwd="3")
lines(x2,y5,col="green4",lty="solid",lwd="3")