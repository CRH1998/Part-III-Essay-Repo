############################################################
#                                                          #
#   The following script includes simulation functions     #
#   for the test of heterogenous effects simulation study  #
#                                                          #
############################################################


############################################################
library(MASS)
library(lava)

############################################################


############################################################
#                 Sánchez-Becerra simulations              #
############################################################


#' Title
#'
#' @param n sample size
#' @param p number of covariates in X0 and X1
#' @param rho correlation between X0j and X1j
#' @param mu mean of X0 and X1
#' @param Dp Treatment assignment probability
#' @param tau Average effect size
#' @param c Hyperparameter
#' @param sigma_d Hyperparameter
#' @param tilde_sigma_d Hyperparameter
#' @param lambda Exponential decay of coefficients
#' @param V_mu Baseline variance (variance of X0)
#' @param V_tau VCATE
#'
#' @returns Simulated data set
#' @export
#'
#' @examples
sim_SB <- function(n, p,
                   rho = 0.5, mu = 0,
                   Dp = 0.5,
                   tau = 0.15, c = 1, sigma_d = sqrt(0.7), tilde_sigma_d = 0.21, lambda = 0.7,
                   V_mu = 0.3, V_tau){
  
  #browser()
  
  # Calculating derived hyperparameters
  l <- sqrt((1 - lambda)/(1 - lambda^p) * lambda^(1:p-1))     # Geometric sequence for decaying effect
  kappa_0 <- l * sqrt(sigma_d^2 - tilde_sigma_d^2)
  kappa_1 <- l * sqrt(sigma_d^2 - tilde_sigma_d^2)
  
  beta0 <- l * sqrt(V_mu)
  beta1 <- l * sqrt(V_tau)
  
  # Simulating covariates
  Sigma_x <- rbind(cbind(diag(p), rho * diag(p)), cbind(rho * diag(p), diag(p)))
  mu_x <- rep(0, 2*p)
  
  X <- mvrnorm(n = n, mu = mu_x, Sigma = Sigma_x)
  X0 <- X[,1:p]
  X1 <- X[,(p+1):(2*p)]
  
  
  # Simulating treatment
  Dp <- 0.5
  D <- rbinom(n, 1, Dp)
  
  
  # Simulating errors
  U <- mvrnorm(n = n, mu = rep(0,2), Sigma = diag(2))
  U0 <- U[,1]
  U1 <- U[,2]
  
  
  # Simulating outcome
  Y0 <- c + X0 %*% beta0 + U0 * c(sqrt(tilde_sigma_d^2 + (X0 %*% kappa_0)^2))
  Y1 <- (c + tau) + X0 %*% beta0 + X1 %*% beta1 + U1 * c(sqrt(tilde_sigma_d^2 + (X1 %*% kappa_1)^2))
  Y <- D * Y1 + (1-D) * Y0
  
  # Return output
  return(data.frame(Y, D, X0, X1))
}

#test_data <- sim_SB(n = 100000, p = 5, V_tau = 0.5)








############################################################
#                 Dukes et al. simulations                 #
############################################################

# Generate uniform random variables:


sim_dukes <- function(n, heterogeneity = c("None", "QuantHMC", "QuantHNMC", "QualHMC", "QualHNMC"),
                      noise_mean = 0, noise_sd = 3){
  
  heterogeneity <- match.arg(heterogeneity)
  
  X1 <- runif(n = n, min = -1, max = 1)
  X2 <- runif(n = n, min = -1, max = 1)
  X3 <- runif(n = n, min = -1, max = 1)
  
  pi_01x <- expit(1/8 * X1 + 1/4*sin(pi*X2))
  A <- rbinom(n = n, size = 1, prob = pi_01x)
  
  hx <- X1 + expit(1/2*(X2 + X3))
  
  if(heterogeneity == "None"){ # No heterogeneity
    gamma <- 3/4
  } else if (heterogeneity == "QuantHMC"){ # Quantitative heterogeneity; monotone CATE
    gamma <- 15*(X3 - 0.5) * (X3 > 0.5)
  } else if (heterogeneity == "QuantHNMC"){ # Quantitative heterogeneity; non-monotone CATE
    gamma <- 3*(1 - X3^2)
  } else if (heterogeneity == "QualHMC"){ # Qualitative heterogeneity; monotone CATE
    gamma <- 3 * sign(X3)*X3^2
  } else if (heterogeneity == "QualHNMC"){ # Qualitative heterogeneity; non-monotone CATE
    gamma <- 3 * cos(3 * pi / 2 * X3)
  }
  
  Y <- hx + A * gamma + rnorm(n = n, mean = noise_mean, sd = noise_sd)
  
  return(data.frame(Y = Y, A = A, X1 = X1, X2 = X2, X3 = X3))
}














  
