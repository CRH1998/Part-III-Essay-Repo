############################################################
#                                                          #
#   The following script includes estimation functions     #
#   for the test of heterogenous effects simulation study  #
#                                                          #
############################################################


############################################################
library(glmnet)
library(SuperLearner)
library(tidyverse)
library(grf)

############################################################


############################################################
#                 Sánchez-Becerra functions                #
############################################################



#' Title
#'
#' @param data The entire data frame
#' @param Y The outcome variable
#' @param D The treatment variable
#' @param xvars Covariates to include; if NULL, all variables except Y and D are used
#' @param k Number of folds
#' @param B Number of repetitions
#' @param alpha glmnet alpha parameter (1 = lasso, 0 = ridge)
#' @param lambda Choice of lambda selection in glmnet: "lambda.1se" or "lambda.min"
#' @param unpenalize_d Whether to leave the treatment main effect unpenalized
#' @param seed Random seed for reproducibility
#'
#' @returns
#' @export
#'
#' @examples
crossfit_cate_lasso <- function(data, Y, D, xvars = NULL,
                                k = 8, B = 1,
                                alpha = 1,
                                lambda = c("lambda.1se", "lambda.min"),
                                unpenalize_d = TRUE,
                                seed = NULL) {
  stopifnot(is.data.frame(data))
  lambda <- match.arg(lambda)
  
  if (!is.null(seed)) set.seed(seed)
  
  # choose covariates automatically if not provided
  if (is.null(xvars)) {
    xvars <- setdiff(names(data), c("Y", "D"))
  }
  
  n <- nrow(data)
  mu0_hat <- rep(NA, n)
  mu1_hat <- rep(NA, n)
  
  # formula: y ~ d * (x1 + x2 + ...)
  rhs <- paste0("D", " * (", paste(xvars, collapse = " + "), ")")
  fml <- as.formula(paste("Y", "~", rhs))
  
  for (b in seq_len(B)) {
    fold <- sample(rep(seq_len(k), length.out = n))
    
    for (i in seq_len(k)) {
      tr_idx <- fold != i
      te_idx <- fold == i
      
      tr <- data[tr_idx, , drop = FALSE]
      te <- data[te_idx, , drop = FALSE]
      
      # Build model matrices with consistent columns via the same terms object
      tt <- terms(fml, data = tr)                 # locks in column structure
      X_tr <- model.matrix(tt, data = tr)
      y_tr <- tr[["Y"]]
      
      # Potential-outcome design matrices (set treatment to 0/1)
      te0 <- te %>% mutate(D = 0)
      te1 <- te %>% mutate(D = 1)
      
      X_te0 <- model.matrix(tt, data = te0)
      X_te1 <- model.matrix(tt, data = te1)
      
      # Optional: no penalization of treatment main effect
      pf <- rep(1, ncol(X_tr))
      if (unpenalize_d) {
        # intercept is already unpenalized by glmnet, but harmless to set:
        pf[colnames(X_tr) == "(Intercept)"] <- 0
        # treatment column is named D
        pf[colnames(X_tr) == "D"] <- 0
      }
      
      fit <- glmnet::cv.glmnet(
        x = X_tr, y = y_tr,
        alpha = alpha,
        penalty.factor = pf
      )
      
      
      mu0_hat[te_idx] <- as.numeric(predict(fit, newx = X_te0, s = lambda))
      mu1_hat[te_idx] <- as.numeric(predict(fit, newx = X_te1, s = lambda))
    }
  }
  
  cate_hat <- mu1_hat - mu0_hat
  ate_hat  <- mean(cate_hat, na.rm = TRUE)
  
  list(
    cate = cate_hat,
    mu0  = mu0_hat,
    mu1  = mu1_hat,
    ate  = ate_hat,
    formula = fml,
    xvars = xvars
  )
}




#' Title
#'
#' @param data The entire data frame
#' @param y The outcome variable
#' @param d The treatment variable
#' @param xvars Covariates to include; if NULL, all variables except Y and D are used
#' @param k Number of folds
#' @param B Number of repetitions
#' @param alpha glmnet alpha parameter (1 = lasso, 0 = ridge)
#' @param lambda Choice of lambda selection in glmnet: "lambda.1se" or "lambda.min"
#' @param unpenalize_d Whether to leave the treatment main effect unpenalized
#' @param seed Random seed for reproducibility
#'
#' @returns
#' @export
#'
#' @examples
crossfit_vcate_lasso <- function(data, y, d, xvars = NULL,
                                k = 8, B = 1,
                                alpha = 1,
                                lambda = c("lambda.1se", "lambda.min"),
                                unpenalize_d = TRUE,
                                seed = NULL) {
  browser()
  stopifnot(is.data.frame(data))
  lambda <- match.arg(lambda)
  
  if (!is.null(seed)) set.seed(seed)
  
  # choose covariates automatically if not provided
  if (is.null(xvars)) {
    xvars <- setdiff(names(data), c("Y", "D"))
  }
  
  n <- nrow(data)
  mu0_hat <- rep(NA, n)
  mu1_hat <- rep(NA, n)
  
  # formula: y ~ d * (x1 + x2 + ...) - here including interactions for heterogeneity
  rhs <- paste0("D", " * (", paste(xvars, collapse = " + "), ")")
  fml <- as.formula(paste("Y", "~", rhs))
  
  for (b in seq_len(B)) {
    fold <- sample(rep(seq_len(k), length.out = n))
    VCATE <- numeric(k)
    
    for (i in seq_len(k)) {
      tr_idx <- fold != i
      te_idx <- fold == i
      
      tr <- data[tr_idx, , drop = FALSE]
      te <- data[te_idx, , drop = FALSE]
      
      # Build model matrices with consistent columns via the same terms object
      tt <- terms(fml, data = tr)                 # locks in column structure, we include interaction effects to estimate heterogeneity
      X_tr <- model.matrix(tt, data = tr)
      y_tr <- tr[["Y"]]
      
      # Potential-outcome design matrices (set treatment to 0/1)
      te0 <- te %>% mutate(D = 0)
      te1 <- te %>% mutate(D = 1)
      
      X_te0 <- model.matrix(tt, data = te0)
      X_te1 <- model.matrix(tt, data = te1)
      
      # Optional: no penalization of treatment main effect
      pf <- rep(1, ncol(X_tr))
      if (unpenalize_d) {
        # intercept is already unpenalized by glmnet, but harmless to set:
        pf[colnames(X_tr) == "(Intercept)"] <- 0
        # treatment column is named D
        pf[colnames(X_tr) == "D"] <- 0
      }
      
      # Fitting the model - this could be modified to include other learners. Maybe use SuperLearner?
      fit <- glmnet::cv.glmnet(
        x = X_tr, y = y_tr,
        alpha = alpha,
        penalty.factor = pf
      )
      
      
      # Estimating mu1_hat and mu0_hat
      mu0_hat_k <- as.numeric(predict(fit, newx = X_te0, s = lambda))
      mu1_hat_k <- as.numeric(predict(fit, newx = X_te1, s = lambda))
      
      # Estimating CATE (tau_hat(X)) and ATE (tau_hat)
      cate_hat_k <- mu1_hat_k - mu0_hat_k
      ate_hat_k <- mean(cate_hat_k)
      
      # Calculating S_hat
      S_hat_k <- cate_hat_k - ate_hat_k
      
      # Estimating M_hat
      pi <- 0.5 # Change this to a function perhaps, if we wish to include the option of estimating propensity
      M_hat_k <- mu0_hat_k + pi * cate_hat_k
      
      # Estimating W_k
      W_hat_k <- cbind(1, M_hat_k, (te$D - pi), (te$D - pi) * S_hat_k)
      
      # Estimating theta_hat
      theta_hat <- solve(crossprod(W_hat_k) / nrow(W_hat_k)) %*% colMeans(W_hat_k * te$Y)
      
      # Estimating V_xnk_hat
      V_xnk_hat <- mean(S_hat_k^2)
      
      # Estimating V_taunk_hat
      V_taunk_hat <- V_xnk_hat * theta_hat[4]^2
      
      VCATE[i] <- V_taunk_hat
      
      mu0_hat[te_idx] <- mu0_hat_k
      mu1_hat[te_idx] <- mu1_hat_k
    }
  }
  
  cate_hat <- mu1_hat - mu0_hat
  ate_hat  <- mean(cate_hat, na.rm = TRUE)
  
  list(
    cate = cate_hat,
    mu0  = mu0_hat,
    mu1  = mu1_hat,
    ate  = ate_hat,
    formula = fml,
    xvars = xvars
  )
}


test_data <- sim_SB(n = 5000, p = 5, V_tau = 0.5, tau = 0.15)
crossfit_cate_lasso(data = test_data, y = test_data$Y, d = test_data$D)
test_data %>% filter(D == 0) %>% summarize(mean(Y))
test_data %>% filter(D == 1) %>% summarize(mean(Y))








