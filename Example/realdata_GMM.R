# Example: Use grid-search TVB on old faithful dataset to build a dictionary

# Import necessary functions
source("~/code/TVB_support_gmm.R")


#############################
#############################
##                         ##
## Part I: Grid Search TVB ##
##                         ##
#############################
#############################


###################
# Make Dictionary #
###################

set.seed(42)

data("faithful")
mydata <- faithful
mydata <- as.matrix(mydata)
target_coverage <- 0.95
omega_grid <- exp(seq(from = log(0.001), to = log(1), length.out = 500))

## Make TVB dictionary
omega_boot <- omega_bootstrap(omega_grid, mydata)

# Save dictionary for future usage
save(omega_boot, file = "~/results/TVB_dict.RData")


##################
# Use Dictionary #
##################

## Case I: calibrated CI of mixing parameter
cover_pi <- numeric(length(omega_grid))

for (k in 1:length(omega_grid)) {
  
    boot_val <- omega_boot[[k]]
    est_cal <- boot_val$vb_results$pi_k[1]
    B <- length(boot_val$vb_boot)
    coverage_count <- 0
    
    for (b in 1:B) {
      
        vb_boot_result <- boot_val$vb_boot[[b]]
        ci <- dir_beta_ci(vb_boot_result$lambda, 1, conf_level = target_coverage)
        if (est_cal < ci[2] && est_cal > ci[1]) {
            coverage_count <- coverage_count + 1
        }
      
    }
    
    cover_boot <- coverage_count / B
    cover_pi[k] <- cover_boot
    
}

diff <- abs(cover_pi - target_coverage)
idx_pi <- which(diff == min(diff))
best_idx_pi <- max(idx_pi)
omega_tvb_grid_pi <- omega_grid[best_idx_pi] # Selected omega for CI of pi
gmm_tvb_grid_pi <- vb_gmm_pac(mydata, m = 2, omega = omega_tvb_grid_pi)
ci_tvb_grid_pi <- dir_beta_ci(gmm_tvb_grid_pi$lambda, 1, conf_level = target_coverage)

## Case II: calibrated CI of mean estimation (\sum\beta_{1})
cover_beta <- numeric(length(omega_grid))

for (k in 1:length(omega_grid)) {
  
    boot_val <- omega_boot[[k]]
    est_cal <- sum(boot_val$vb_results$rho[, 1])
    B <- length(boot_val$vb_boot)
    coverage_count <- 0
    
    for (b in 1:B) {
      
        vb_boot_result <- boot_val$vb_boot[[b]]
        mean_boot <- sum(vb_boot_result$rho[, 1])
        var_boot <- sum(solve((vb_boot_result$Phi[, , 1] * vb_boot_result$nu[1]) * vb_boot_result$beta[1]))
        ci <- normal_ci(mean_boot, var_boot, conf_level = target_coverage)
        if (est_cal < ci[2] && est_cal > ci[1]) {
            coverage_count <- coverage_count + 1
        }
      
    }
    
    cover_boot <- coverage_count / B
    cover_beta[k] <- cover_boot
  
}

diff <- abs(cover_beta - target_coverage)
idx_beta <- which(diff == min(diff))
best_idx_beta <- max(idx_beta)
omega_tvb_grid_beta <- omega_grid[best_idx_beta] # Selected omega for CI of sum of beta_1
gmm_tvb_grid_beta <- vb_gmm_pac(mydata, m = 2, omega = omega_tvb_grid_beta)
mean_tvb_grid <- sum(gmm_tvb_grid_beta$rho[, 1])
var_tvb_grid <- sum(solve((gmm_tvb_grid_beta$Phi[, , 1] * gmm_tvb_grid_beta$nu[1]) * gmm_tvb_grid_beta$beta[1]))
ci_tvb_grid_beta <- normal_ci(mean_tvb_grid, var_tvb_grid, conf_level = target_coverage)


#############################################
#############################################
##                                         ##
## Part II: Sequential TVB and original VB ##
##                                         ##
#############################################
#############################################


### Original VB ###

## Case I: mixing parameter pi
gmm_vb_pi <- vb_gmm_pac(mydata, m = 2)
ci_vb_pi <- dir_beta_ci(gmm_vb_pi$lambda, 1, conf_level = target_coverage)

## Case II: mean estimation (\sum\beta_{1})
gmm_vb_beta <- vb_gmm_pac(mydata, m = 2)
mean_vb_beta <- sum(gmm_vb_beta$rho[, 1])
var_vb_beta <- sum(solve((gmm_vb_beta$Phi[, , 1] * gmm_vb_beta$nu[1]) * gmm_vb_beta$beta[1]))
ci_vb_beta <- normal_ci(mean_vb_beta, var_vb_beta, conf_level = target_coverage)


### Sequential TVB ###

set.seed(42)

## Case I: mixing parameter pi
omega_tvb_seq_pi <- update_omega_pi(mydata, idx = 1)
gmm_tvb_seq_pi <- vb_gmm_pac(mydata, m = 2, omega = omega_tvb_seq_pi$omega)
ci_tvb_seq_pi <- dir_beta_ci(gmm_tvb_seq_pi$lambda, 1, conf_level = target_coverage)

## Case II: mean estimation (\sum\beta_{1})
omega_tvb_seq_beta <- update_omega_beta(mydata, idx = 1)
gmm_tvb_seq_beta <- vb_gmm_pac(mydata, m = 2, omega = omega_tvb_seq_beta$omega)
mean_tvb_seq_beta <- sum(gmm_tvb_seq_beta$rho[, 1])
var_tvb_seq_beta <- sum(solve((gmm_tvb_seq_beta$Phi[, , 1] * gmm_tvb_seq_beta$nu[1]) * gmm_tvb_seq_beta$beta[1]))
ci_tvb_seq_beta <- normal_ci(mean_tvb_seq_beta, var_tvb_seq_beta, conf_level = target_coverage)


# Save the estimation results
save(
  
    ## pi ##
    
    # TVB2
    cover_pi,
    gmm_tvb_grid_pi,
    omega_tvb_grid_pi,
    ci_tvb_grid_pi,
    
    # VB
    gmm_vb_pi,
    ci_vb_pi,
    
    # TVB1
    gmm_tvb_seq_pi,
    omega_tvb_seq_pi,
    ci_tvb_seq_pi,
    
    ## beta ##
    
    # TVB2
    cover_beta,
    gmm_tvb_grid_beta,
    omega_tvb_grid_beta,
    ci_tvb_grid_beta,
    
    # VB
    gmm_vb_beta,
    ci_vb_beta,
    
    # TVB1
    gmm_tvb_seq_beta,
    omega_tvb_seq_beta,
    ci_tvb_seq_beta,
    
    file = "~/results/gmm_realdata.RData"
  
)





