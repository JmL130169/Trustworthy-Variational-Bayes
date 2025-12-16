# Load packages
suppressPackageStartupMessages({
  
  suppressWarnings(library(MASS))
  suppressWarnings(library(Matrix))
  suppressWarnings(library(clue))
  suppressWarnings(library(stats))
  suppressWarnings(library(boot))
  suppressWarnings(library(gtools))
  suppressWarnings(library(ghyp))
  suppressWarnings(library(ald))
  suppressWarnings(library(statmod))
  suppressWarnings(library(ggplot2))
  suppressWarnings(library(optparse))
  suppressWarnings(library(matrixcalc))
  suppressWarnings(library(data.table))
  suppressWarnings(library(purrr))
  suppressWarnings(library(mclust))
  suppressWarnings(library(mvtnorm))
  suppressWarnings(library(parallel))
  suppressWarnings(library(latex2exp))
  suppressWarnings(library(patchwork))
  
})

# Load functions
source("./TVB_support_gmm.R")

##################################
# Scatter plot for VB estimation #
##################################

data("faithful")
mydata <- faithful
mydata_numeric <- data.frame(lapply(mydata, function(x) as.numeric(as.character(x))))
sum(is.na(mydata_numeric))
mydata <- as.matrix(mydata_numeric)

vb <- vb_gmm_pac(mydata, omega = 1)

mu_1 <- vb$rho[, 1]
mu_2 <- vb$rho[, 2]
Sig_1 <- solve((vb$Phi[, , 1] * vb$nu[1]) * vb$beta[1])
Sig_2 <- solve((vb$Phi[, , 2] * vb$nu[2]) * vb$beta[2])
Z_prob <- vb$r_nk / rowSums(vb$r_nk)
Z_labels <- apply(Z_prob, 1, function(x) which.max(x))
Z <- factor(Z_labels)
df <- data.frame(X1 = faithful[, 1], X2 = faithful[, 2], Cluster = Z)

p1 <- ggplot() +
  geom_point(data = df, aes(x = X1, y = X2, color = Cluster), size = 2) +
  geom_point(aes(x = mu_1[1], y = mu_1[2]), color = 'blue', size = 4, shape = 8) +
  geom_point(aes(x = mu_2[1], y = mu_2[2]), color = 'red', size = 4, shape = 8) +
  labs(x = "eruptions", y = "waiting") +
  theme_minimal() +
  theme(legend.position = "none", aspect.ratio = 1)


#################################
# Density plot (use dictionary) #
#################################

load("~/gmm_realdata.RData")

pi_est <- matrix(0, 2, 2)
lam_est <- matrix(0, 2, 2)
pi_est[1, ] <- gmm_vb_pi$pi_k
lam_est[1, ] <- gmm_vb_pi$lambda
pi_est[2, ] <- gmm_tvb_grid_pi$pi_k
lam_est[2, ] <- gmm_tvb_grid_pi$lambda

# Define model names
model_names <- c("VB", "TVB-grid")

# Plot the posterior of pi
x_vals <- seq(0, 1, length.out = 1000)

df_line <- do.call(rbind, lapply(1:2, function(i) {
  a <- lam_est[i, 1]
  b <- lam_est[i, 2]
  data.frame(
    x = x_vals,
    density = dbeta(x_vals, a, b),
    group = factor(model_names[i])
  )
}))

df_ribbon <- do.call(rbind, lapply(1:2, function(i) {
  a <- lam_est[i, 1]
  b <- lam_est[i, 2]
  lower_ci <- qbeta(0.025, a, b)
  upper_ci <- qbeta(0.975, a, b)
  idx <- x_vals >= lower_ci & x_vals <= upper_ci
  data.frame(
    x = x_vals[idx],
    density = dbeta(x_vals[idx], a, b),
    group = factor(model_names[i])
  )
}))

p2 <- ggplot() +
  geom_ribbon(data = df_ribbon, aes(x = x, ymin = 0, ymax = density, fill = group),
              alpha = 0.1) +
  geom_line(data = df_line, aes(x = x, y = density, color = group), size = 0.5) +
  labs(x = TeX("$\\pi$"),
       y = "Density") +
  scale_color_discrete(name = NULL) +
  scale_fill_discrete(name = NULL) +
  theme_minimal(base_size = 10) +
  coord_cartesian(ylim = c(0, 15)) + 
  theme(legend.position = "right", aspect.ratio = 1)

# Combine plots
plot_tvb <- grid.arrange(p1, p2, ncol = 2)
ggsave("~/plot_tvb.png", plot_tvb, width = 14, height = 6, dpi = 300)