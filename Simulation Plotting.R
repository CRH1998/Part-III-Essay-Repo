# Simulation plotting
library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(patchwork)


#####################################################################################
#                 Different levels of heterogeneity - SB simulations                #
#####################################################################################
SBdiffVCATE <- read_table("Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/Simulation results/simulation_resultsSBsim vs VCATE 02042026.txt")
colnames(SBdiffVCATE) <- c("V_tau", "V_taun", "coverage_multifold", "ci_low_multifold", "ci_high_multifold", "beta1n", "nonparametricp")



sumSBdiffVCATE <- SBdiffVCATE %>% 
  group_by(V_tau) %>% 
  summarise(coverage = mean(coverage_multifold),
            semiparametriclevel = mean(0 < ci_low_multifold & V_tau == 0),
            semiparametricpower = mean(ci_low_multifold > 0),
            nonparametriclevel = mean(nonparametricp < 0.05 & V_tau == 0),
            nonparametricpower = mean(nonparametricp < 0.05),
            ci_low_multifold = median(ci_low_multifold),
            ci_high_multifold = median(ci_high_multifold),
            V_taun = mean(V_taun),
            beta1n = mean(beta1n))


# SB pseudo-VCATE plot
p1 <- ggplot(sumSBdiffVCATE, aes(x = V_tau)) +
  
  # Main curves with legend mapping
  geom_line(aes(y = V_taun), linewidth = 0.6) +
  geom_line(aes(y = ci_low_multifold), linetype = "dashed", linewidth = 0.6) +
  geom_line(aes(y = ci_high_multifold), linetype = "dashed", linewidth = 0.6) +
  
  # Axis labels
  labs(
    title = "Estimate of VCATE with error bounds",
    x = expression(V[tau]),
    y = "Estimated VCATE",
    color = NULL
  ) +
  
  # Axis limits
  coord_cartesian(ylim = c(0, 0.5)) +
  
  # Theme (publication style)
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    legend.text = element_text(size = 10)
  )

# Print plot
print(p1)

# Save to file (PDF recommended for papers)
ggsave("/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/pseudo_VCATE_plot.png", plot = p1, width = 6, height = 4)




# -------- Power plot --------
p_power <- ggplot(sumSBdiffVCATE, aes(x = V_tau)) +
  
  geom_line(aes(y = semiparametricpower, color = "SB"), linewidth = 0.7) +
  geom_line(aes(y = nonparametricpower, color = "Dukes"), linetype = "dashed", linewidth = 0.7) +
  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  
  
  geom_point(aes(y = semiparametricpower, color = "SB"),
             size = 1.8, shape = 1) +
  geom_point(aes(y = nonparametricpower, color = "Dukes"),
             size = 1.8, shape = 2) +
  
  labs(
    title = "Power",
    x = expression(V[tau]),
    y = "Probability",
    color = NULL
  ) +
  
  scale_color_manual(values = c(
    "SB" = "#2e8b57",
    "Dukes" = "#1b3c73"
  )) +
  
  coord_cartesian(ylim = c(0, 1)) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# -------- Level plot --------
p_level <- ggplot(sumSBdiffVCATE, aes(x = V_tau)) +
  
  geom_line(aes(y = semiparametriclevel, color = "SB"), linewidth = 0.7) +
  geom_line(aes(y = nonparametriclevel, color = "Dukes"), linewidth = 0.7, linetype = "dashed") +
  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  
  
  geom_point(aes(y = semiparametriclevel, color = "SB"),
             size = 1.8, shape = 1) +
  geom_point(aes(y = nonparametriclevel, color = "Dukes"),
             size = 1.8, shape = 2) +
  
  labs(
    title = "Level",
    x = expression(V[tau]),
    y = NULL,
    color = NULL
  ) +
  
  scale_color_manual(values = c(
    "SB" = "#2e8b57",
    "Dukes" = "#1b3c73"
  )) +
  
  coord_cartesian(ylim = c(0, 0.075)) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# -------- Combine --------
p_combined <- p_power + p_level

# Print
print(p_combined)

# Save
ggsave("/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/level_power_side_by_side.png", plot = p_combined, width = 8, height = 3)




#####################################################################################
#                 Different levels of n - SB simulations                            #
#####################################################################################

# No heterogeneity
SBdiffn0 <- read_table("Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/Simulation results/simulation_results_diffn04042026_no_heterogeneity.txt")
colnames(SBdiffn0) <- c("n", "V_tau", "V_taun", "coverage_multifold", "ci_low_multifold", "ci_high_multifold", "beta1n", "nonparametricp")
sumSBdiffn0 <- SBdiffn0 %>% 
  group_by(n) %>% 
  summarise(coverage = mean(coverage_multifold),
            semiparametriclevel = mean(0 < ci_low_multifold & V_tau == 0),
            semiparametricpower = mean(ci_low_multifold > 0),
            nonparametriclevel = mean(nonparametricp < 0.05 & V_tau == 0),
            nonparametricpower = mean(nonparametricp < 0.05),
            ci_low_multifold = median(ci_low_multifold),
            ci_high_multifold = median(ci_high_multifold),
            V_taun = mean(V_taun),
            beta1n = mean(beta1n))


p_level0 <- ggplot(sumSBdiffn0, aes(x = n)) +
  
  geom_line(aes(y = semiparametriclevel, color = "SB"), linewidth = 0.7) +
  geom_line(aes(y = nonparametriclevel, color = "Dukes"), linewidth = 0.7, linetype = "dashed") +
  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  
  
  geom_point(aes(y = semiparametriclevel, color = "SB"),
             size = 1.8, shape = 1) +
  geom_point(aes(y = nonparametriclevel, color = "Dukes"),
             size = 1.8, shape = 2) +
  
  labs(
    title = "Level",
    x = expression(n),
    y = NULL,
    color = NULL
  ) +
  
  scale_color_manual(values = c(
    "SB" = "#2e8b57",
    "Dukes" = "#1b3c73"
  )) +
  
  coord_cartesian(ylim = c(0, 0.1)) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )






# 0.01 heterogeneity

SBdiffn500_2500 <- read_table("Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/Simulation results/simulation_results_diffn04042026_0.1heterogeneity.txt")
SBdiffn5000 <- read_table("Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/Simulation results/simulation_results_5000n05042026_0.1heterogeneity.txt")
SBdiffn01 <- rbind(SBdiffn500_2500, SBdiffn5000)
colnames(SBdiffn01) <- c("n", "V_tau", "V_taun", "coverage_multifold", "ci_low_multifold", "ci_high_multifold", "beta1n", "nonparametricp")



sumSBdiffn01 <- SBdiffn01 %>% 
  group_by(n) %>% 
  summarise(coverage = mean(coverage_multifold),
            semiparametriclevel = mean(0 < ci_low_multifold & V_tau == 0),
            semiparametricpower = mean(ci_low_multifold > 0),
            nonparametriclevel = mean(nonparametricp < 0.05 & V_tau == 0),
            nonparametricpower = mean(nonparametricp < 0.05),
            ci_low_multifold = median(ci_low_multifold),
            ci_high_multifold = median(ci_high_multifold),
            V_taun = mean(V_taun),
            beta1n = mean(beta1n))


# -------- Power plot --------
p_power01 <- ggplot(sumSBdiffn01, aes(x = n)) +
  
  geom_line(aes(y = semiparametricpower, color = "SB"), linewidth = 0.7) +
  geom_line(aes(y = nonparametricpower, color = "Dukes"), linetype = "dashed", linewidth = 0.7) +
  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  
  
  geom_point(aes(y = semiparametricpower, color = "SB"),
             size = 1.8, shape = 1) +
  geom_point(aes(y = nonparametricpower, color = "Dukes"),
             size = 1.8, shape = 2) +
  
  labs(
    title = "Power",
    x = expression(n),
    y = "Probability",
    color = NULL
  ) +
  
  scale_color_manual(values = c(
    "SB" = "#2e8b57",
    "Dukes" = "#1b3c73"
  )) +
  
  coord_cartesian(ylim = c(0, 1)) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# -------- Combine --------
p_combined_n <- p_power01 + p_level0

# Print
print(p_combined_n)

ggsave("/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/level_power_side_by_side_n.png", plot = p_combined_n, width = 8, height = 3)




#####################################################################################
#                 Different levels of VCATE - Dukes simulations                     #
#####################################################################################

# DukesdiffVCATE <- read_table("Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/Simulation results/simulation_results_DukesXGBLassoRF0604.txt")
# DukesdiffVCATE$`"nonparametricp"` <- NULL
# colnames(DukesdiffVCATE) <- c("n", "V_tau", "V_taun", "coverage_multifold", "ci_low_multifold", "ci_high_multifold", "beta1n")
# 
# 
# sumDukesdiffVCATE <- DukesdiffVCATE %>% 
#   group_by(V_tau) %>% 
#   summarise(coverage = mean(coverage_multifold),
#             semiparametriclevel = mean(0 < ci_low_multifold & V_tau == 0),
#             semiparametricpower = mean(ci_low_multifold > 0),
#             ci_low_multifold = median(ci_low_multifold),
#             ci_high_multifold = median(ci_high_multifold),
#             V_taun = mean(V_taun),
#             beta1n = mean(beta1n))
# 
# 
# 
# 
# p1_dukes <- ggplot(sumDukesdiffVCATE, aes(x = V_tau)) +
#   
#   # Main curves with legend mapping
#   geom_line(aes(y = V_taun), linewidth = 0.6) +
#   geom_point(aes(y = V_taun), size = 1.8, shape = 1) +
#   geom_line(aes(y = ci_low_multifold), linetype = "dashed", linewidth = 0.6) +
#   geom_line(aes(y = ci_high_multifold), linetype = "dashed", linewidth = 0.6) +
#   
#   # Axis labels
#   labs(
#     title = "Estimate of VCATE with error bounds",
#     x = expression(V[tau]),
#     y = "Estimated VCATE",
#     color = NULL
#   ) +
#   
#   # Axis limits
#   coord_cartesian(ylim = c(0, 5)) +
#   
#   # Theme (publication style)
#   theme_bw(base_size = 12) +
#   theme(
#     legend.position = "right",
#     legend.direction = "vertical",
#     legend.box = "vertical",
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_line(color = "gray90"),
#     legend.text = element_text(size = 10)
#   )
# 
# # Print plot
# print(p1_dukes)
# ggsave("/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/p1_dukes.png", plot = p1_dukes, width = 6, height = 4)
# 
# 
# 
# # -------- Coverage plot --------
# p_coverage_dukes <- ggplot(sumDukesdiffVCATE, aes(x = V_tau)) +
#   
#   geom_line(aes(y = coverage), linewidth = 0.7) +
#   
#   geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray40", linewidth = 0.6) +
#   
#   
#   geom_point(aes(y = coverage),
#              size = 1.8, shape = 1) +
#   
#   labs(
#     title = "Coverage",
#     x = expression(V[tau]),
#     y = "Probability",
#     color = NULL
#   ) +
#   
#   coord_cartesian(ylim = c(0, 1)) +
#   
#   theme_bw(base_size = 12) +
#   theme(
#     legend.position = "none",
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_line(color = "gray90"),
#     plot.title = element_text(face = "bold", hjust = 0.5)
#   )
# 
# 
# # -------- Power plot --------
# p_power_dukes <- ggplot(sumDukesdiffVCATE, aes(x = V_tau)) +
#   
#   geom_line(aes(y = semiparametricpower), linewidth = 0.7) +
#   
#   geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40", linewidth = 0.6) +
#   
#   
#   geom_point(aes(y = semiparametricpower),
#              size = 1.8, shape = 1) +
#   
#   labs(
#     title = "Power",
#     x = expression(V[tau]),
#     y = NULL,
#     color = NULL
#   ) +
#   
#   coord_cartesian(ylim = c(0, 1)) +
#   
#   theme_bw(base_size = 12) +
#   theme(
#     legend.position = "none",
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_line(color = "gray90"),
#     plot.title = element_text(face = "bold", hjust = 0.5)
#   )
# 
# # -------- Level plot --------
# p_level_dukes <- ggplot(sumDukesdiffVCATE, aes(x = V_tau)) +
#   
#   geom_line(aes(y = semiparametriclevel), linewidth = 0.7) +
#   
#   geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40", linewidth = 0.6) +
#   
#   
#   geom_point(aes(y = semiparametriclevel),
#              size = 1.8, shape = 1) +
#   
#   labs(
#     title = "Level",
#     x = expression(V[tau]),
#     y = NULL,
#     color = NULL
#   ) +
#   
#   
#   coord_cartesian(ylim = c(0, 0.075)) +
#   
#   theme_bw(base_size = 12) +
#   theme(
#     legend.position = "right",
#     panel.grid.minor = element_blank(),
#     panel.grid.major = element_line(color = "gray90"),
#     plot.title = element_text(face = "bold", hjust = 0.5)
#   )
# 
# # -------- Combine --------
# p_combined_dukes <- p_coverage_dukes + p_power_dukes + p_level_dukes
# 
# # Print
# print(p_combined_dukes)
# ggsave("/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/p_combined_dukes.png", plot = p_combined_dukes, width = 8, height = 4)
# 







#####################################################################################
#                 Different levels of VCATE - New Dukes simulations                 #
#####################################################################################

DukesdiffVCATE <- read_table("/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/Simulation results/DukesSimulations16042026.txt")
DukesdiffVCATE$`"nonparametricp"` <- NULL
colnames(DukesdiffVCATE) <- c("n", "V_tau", "V_taun", "coverage_multifold", "ci_low_multifold", "ci_high_multifold", "beta1n")


sumDukesdiffVCATE <- DukesdiffVCATE %>% 
  group_by(V_tau) %>% 
  summarise(coverage = mean(coverage_multifold),
            semiparametriclevel = mean(0 < ci_low_multifold & V_tau == 0),
            semiparametricpower = mean(ci_low_multifold > 0),
            ci_low_multifold = median(ci_low_multifold),
            ci_high_multifold = median(ci_high_multifold),
            V_taun = mean(V_taun),
            beta1n = mean(beta1n))


# -------- Coverage plot --------
p_coverage_dukes <- ggplot(sumDukesdiffVCATE, aes(x = V_tau)) +
  
  geom_line(aes(y = coverage), linewidth = 0.7) +
  
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  
  
  geom_point(aes(y = coverage),
             size = 1.8, shape = 1) +
  
  labs(
    title = "Coverage",
    x = expression(V[tau]),
    y = "Probability",
    color = NULL
  ) +
  
  coord_cartesian(ylim = c(0, 1)) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


# -------- Power plot --------
p_power_dukes <- ggplot(sumDukesdiffVCATE, aes(x = V_tau)) +
  
  geom_line(aes(y = semiparametricpower), linewidth = 0.7) +
  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  
  
  geom_point(aes(y = semiparametricpower),
             size = 1.8, shape = 1) +
  
  labs(
    title = "Power",
    x = expression(V[tau]),
    y = NULL,
    color = NULL
  ) +
  
  coord_cartesian(ylim = c(0, 1)) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# -------- Level plot --------
p_level_dukes <- ggplot(sumDukesdiffVCATE, aes(x = V_tau)) +
  
  geom_line(aes(y = semiparametriclevel), linewidth = 0.7) +
  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  
  
  geom_point(aes(y = semiparametriclevel),
             size = 1.8, shape = 1) +
  
  labs(
    title = "Level",
    x = expression(V[tau]),
    y = NULL,
    color = NULL
  ) +
  
  
  coord_cartesian(ylim = c(0, 0.075)) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


# -------- Density plot --------
plot_list <- list()

for (i in 1:nrow(sumDukesdiffVCATE)) {
  V_tau_value <- sumDukesdiffVCATE$V_tau[i]
  V_taun_value <- sumDukesdiffVCATE$V_taun[i]
  
  data_subset <- DukesdiffVCATE %>% filter(V_tau == V_tau_value)
  
  p <- ggplot(data_subset, aes(x = V_taun)) +
    geom_histogram(bins = 30, aes(y = after_stat(density)), fill = "steelblue", color = "black") +
    geom_density(color = "red", size = 1) +
    labs(title = paste("True VCATE is", round(V_tau_value, 3)),
         x = expression(V[taun]),
         y = "Density") +
    theme_minimal() +
    geom_vline(xintercept = V_tau_value, linetype = "dashed", color = "red", linewidth = 0.7)
  
  plot_list[[i]] <- p
}



# -------- Power against local alternatives ------
PALA <- read_table("Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/Simulation results/DukesPowerAgainstLocalAlternatives16042026.txt")
PALA$`"nonparametricp"` <- NULL
colnames(PALA) <- c("n", "V_tau", "V_taun", "coverage_multifold", "ci_low_multifold", "ci_high_multifold", "beta1n")


sumPALA <- PALA %>% 
  group_by(V_tau) %>% 
  summarise(n = mean(n),
            coverage = mean(coverage_multifold),
            semiparametriclevel = mean(0 < ci_low_multifold & V_tau == 0),
            semiparametricpower = mean(ci_low_multifold > 0),
            ci_low_multifold = median(ci_low_multifold),
            ci_high_multifold = median(ci_high_multifold),
            V_taun = mean(V_taun),
            beta1n = mean(beta1n))

# With V_tau on x-axis at botton and n on x-axis at top
p_power_local <- ggplot(sumPALA, aes(x = n)) +
  
  geom_line(aes(y = semiparametricpower), color = "black", linewidth = 0.7) +
  
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  
  
  geom_point(aes(y = semiparametricpower),
             size = 1.8, shape = 1, color = "black") +
  
  labs(
    title = "Power against local alternatives",
    x = "n",
    y = "Probability",
    color = NULL
  ) +
  
  coord_cartesian(ylim = c(0, 0.15)) +
  
  theme_bw(base_size = 20) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# -------- Combine plots --------
dens_p <- (plot_list[[1]] + plot_list[[2]] + plot_list[[3]]) / (plot_list[[4]] + plot_list[[5]] + plot_list[[6]])
par(mfrow = c(1, 4))
p_combined_dukes <- p_coverage_dukes + p_power_dukes + p_level_dukes


# Save the combined plot
ggsave("/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/density_plots_VCATE.png", plot = dens_p, width = 12, height = 7)
ggsave("/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/power_against_local_alternatives.png", plot = p_power_local, width = 8, height = 5)
ggsave("/Users/christianhejstvig-larsen/Desktop/Studie/Master/Cambridge/Essay/R/Part-III-Essay-Repo/p_combined_dukes.png", plot = p_combined_dukes, width = 12, height = 4)
