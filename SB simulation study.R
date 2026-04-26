#############################
#     Simulation study      #
#############################



# Testing functions and code
source("/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/Estimation functions.R")
source("/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/Simulations.R")
source("/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/code-JASA/np-qualint-JASA.R")



# Simulate data
# For defining "SL.glmnet_1"
nk <- 250
train_data <- sim_SB(n = nk, J = 2, V_tau = 0.01, tau = 0.15)
test_data <- sim_SB(n = nk, J = 2, V_tau = 0.01, tau = 0.15)
nuisance_param <- nuisance_param_func(train_data = train_data,
                                      eval_data = test_data,
                                      SL.lib = list(create.Learner(
                                        "SL.glmnet", params = list(alpha = 1, nfolds = 20)
                                      )$names),
                                      hal = F,
                                      K = 5)




run_one <- function(n, sim_type = c("SB", "Dukes"), randomized = TRUE, seed = NA,
                    J, V_tau = NA, tau = 0.15, family = gaussian(), K = 2, hal = F,
                    heterogeneity = c("None", "QuantHMC", "QuantHNMC"), heterogeneity_degree = 15, pi_hat = 0.5,
                    analysis_type = c("SB", "Dukes", "both"),
                    SL.lib = list(#create.Learner(
                      #"SL.glmnet", params = list(alpha = 1, nfolds = 10)
                    #)$names, 
                      "SL.glmnet", "SL.xgboost", "SL.ranger")) {
  
  sim_type <- match.arg(sim_type)
  analysis_type <- match.arg(analysis_type)
  
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  
  # Choose simulation type
  if (sim_type == "SB") {
    data <- sim_SB(n = n, J = J, V_tau = V_tau, tau = tau)
  } else if (sim_type == "Dukes") {
    data <- sim_Dukes(n = n, heterogeneity = heterogeneity, 
                      heterogeneity_degree = heterogeneity_degree, randomized = randomized, pi_01x = pi_hat)
    
    V_tau <- 13/768 * heterogeneity_degree^2
  }
  
  ci_low_multifold <- NA
  ci_high_multifold <- NA
  coverage_multifold <- NA
  V_taun <- NA
  beta1n <- NA
  DukesP <- NA
  
  # Choose analysis type
  if (analysis_type == "SB" | analysis_type == "both") {
    SB_CF <- crossfit_pseudoVCATE(data = data,
                                  SL.lib = SL.lib,
                                  randomized = randomized,
                                  SL_cv_folds = K,
                                  family = family, 
                                  hal = hal,
                                  crossfit_folds = 10,
                                  pi_hat = pi_hat,
                                  alpha = 0.025,
                                  Vmax = 0.05,
                                  gridN = 400,
                                  M = 2000,
                                  seed = NA)
    
    ci_low_multifold <- median(SB_CF$ci_low)
    ci_high_multifold <- median(SB_CF$ci_high)
    coverage_multifold <- (ci_low_multifold <= V_tau) & (ci_high_multifold >= V_tau)
    V_taun <- mean(SB_CF$V_taunk_star_hat)
    beta1n <- mean(SB_CF$theta_3)
  }
  
  if (analysis_type == "Dukes" | analysis_type == "both") {
    Y <- data$Y
    A <- data$A
    W <- data[, -c(1, 2)]
    
    DukesP <- QuantTestMulti(y = Y, a = A, w = W, x = W, 
                             randomized = randomized, no.bs = 100,
                             SL.lib = SL.lib, 
                             K = K, 
                             family = family,
                             hal = hal)$p.val
    
    
    
    
  } 
  
  
  return(list(n = n,
              V_tau = V_tau,
              V_taun = V_taun,
              coverage_multifold = coverage_multifold,
              ci_low_multifold = ci_low_multifold,
              ci_high_multifold = ci_high_multifold,
              beta1n = beta1n,
              DukesP = DukesP))
}


# Set parameters for the simulation study
n <- 2500
J <- 4
V_tau <- c(0.01)
heterogeneity_degree <- c(0, 5, 10)

set.seed(16042026)
# Run the simulation study
start.time <- Sys.time()
results <- lapply(heterogeneity_degree, function(heterogeneity_degree) {
  replicate(5, run_one(n = n, sim_type = "Dukes", randomized = TRUE, seed = NA,
                        J = J, V_tau = NA, 
                       heterogeneity = c("QuantHMC"), heterogeneity_degree = heterogeneity_degree,
                        K = 10, family = gaussian(), hal = F,
                        analysis_type = "SB"), simplify = TRUE)
})
end.time <- Sys.time()
print(paste("Simulation time:", end.time - start.time))





# Process results
results_df <- lapply(results, function(sim_result) {
  matrix(unlist(sim_result), ncol = 8, byrow = T)
})
results_df <- data.frame(do.call(rbind, results_df))
colnames(results_df) <- c("n", "V_tau", "V_taun", "coverage_multifold", "ci_low_multifold", "ci_high_multifold", "beta1n", "nonparametricp")


#Save results df as txt
#write.table(results_df, file = "/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/DukesSimulations16042026.txt", row.names = FALSE, col.names = TRUE, sep = "\t")


# Summarize results
summary_results <- results_df %>% 
  group_by(V_tau) %>% 
  summarise(n = mean(n),
            coverage = mean(coverage_multifold),
            power = mean(ci_low_multifold > 0),
            ci_low_multifold = mean(ci_low_multifold),
            ci_high_multifold = mean(ci_high_multifold),
            V_taun = mean(V_taun),
            beta1n = mean(beta1n),
            nonparametriclevel = mean(nonparametricp > 0.05 & V_tau == 0),
            nonparametricpower = mean(nonparametricp < 0.05 & V_tau > 0))


# Plot results
ggplot(summary_results, aes(x = n)) +
  geom_line(aes(y = coverage), color = "blue") +
  geom_line(aes(y = power), color = "red") +
  geom_line(aes(y = nonparametricpower), color = "green") +
  labs(title = "Coverage and Power vs V_tau", y = "Proportion", x = "V_tau") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_line(aes(y = 0.95), color = "black", linetype = "dashed") +
  geom_line(aes(y = 0.05), color = "black", linetype = "dashed")














# Dukes et al
V_tau_values <- c(0, 0.01, 0.05, 0.5)
n <- 2500
J <- 3
data <- sim_SB(n = n, J = J, V_tau = 0.01, tau = 0.15)










