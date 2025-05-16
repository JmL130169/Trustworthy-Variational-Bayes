# Bayesian Gaussian mixture model
# Variational calibration
# Calibration for mixing parameter pi

#############
# Preparing #
#############

# Working directory
setwd("~/VBC/") # Replace with your own file path

# Loading package
library(MASS)
library(Matrix)
library(clue)
library(stats)
library(boot)
library(gtools)
library(ggplot2)
library(optparse)
library(matrixcalc)
library(data.table)
library(purrr)
library(mvtnorm)
library(parallel)

# Option
option_list <- list(
  
    make_option("--name_batch", type = "character", default = NA, action = "store"),
    make_option("--seed", type = "numeric", default = NA, action = "store"),
    make_option("--cluster", type = "numeric", default = NA, action = "store"),
    make_option("--dim", type = "numeric", default = NA, action = "store"),
    make_option("--center", type = "numeric", default = NA, action = "store"),
    make_option("--sample", type = "numeric", default = NA, action = "store"),
    make_option("--prob", type = "numeric", default = NA, action = "store")
  
)

opt <- parse_args(OptionParser(option_list = option_list))

name.batch <- opt$name_batch
seed <- opt$seed
m <- opt$cluster
d <- opt$dim
center <- opt$center
num_sample <- opt$sample
true_pi <- opt$prob

# For code test only
# name.batch <- "sim0"
# seed <- 1
# m <- 2
# d <- 2
# center <- 2
# num_sample <- 500
# true_pi <- 0.65

out.id <- as.numeric(gsub("sim", "", name.batch))

# Make folder
dir.create(file.path("results", name.batch), showWarnings = FALSE, recursive = TRUE)


########################
# Supporting functions #
########################

log_sum_exp <- function(x) {
  
    offset <- max(x)
    s <- log(sum(exp(x - offset))) + offset
    i <- which(!is.finite(s))
    if (length(i) > 0) { s[i] <- offset }
    return(s)
  
}

logB <- function(W, nu){
  
    D <- NCOL(W)
    return(-0.5 * nu * log(det(W)) - 
             (0.5 * nu * D * log(2) + 0.25 * D * (D - 1) * 
                log(pi) + sum(lgamma(0.5 * (nu + 1 - 1:NCOL(W)))))
    )
    
}

logB_uni <- function(W, nu){
  
    D <- NCOL(W)
    return(-0.5 * nu * log(W) - 
             (0.5 * nu * D * log(2) + 0.25 * D * (D - 1) * 
                log(pi) + sum(lgamma(0.5 * (nu + 1 - 1:NCOL(W)))))
    )
  
}

compute_ci <- function(estimate, variance, n, level = 0.95) {
  
    z <- qnorm((1 + level) / 2)
    margin <- z * sqrt(variance)
    ci <- c(estimate - margin, estimate + margin)
    return(ci)
  
}

dir_beta_ci <- function(alpha1, alpha2, conf_level = 0.95) {
  
    # Calculate the lower and upper quantiles
    lower_quantile <- (1 - conf_level) / 2
    upper_quantile <- 1 - lower_quantile
    
    # Calculate the CI using the qbeta function
    ci_lower <- qbeta(lower_quantile, alpha1, alpha2)
    ci_upper <- qbeta(upper_quantile, alpha1, alpha2)
    
    # Return the results
    return(list(
      ci_lower = ci_lower,
      ci_upper = ci_upper
    ))
  
}


##################
# Main functions #
##################

# New version of \alpha-VB update
vb_gmm <- function(data, m = 2, lambda0 = 1/m, 
                   rho0 = c(colMeans(data)), omega = 1, 
                   beta0 = 1, nu0 = NCOL(data) + 2, 
                   Phi0 = cov(data), 
                   max_iter = 200, tol = 1e-4) {
  
    ################################
    # Explanation of each variable #
    ################################
    
    # "lambda": parameter for Dir(\pi)
    # "r_nk": responsibility after E-step
    # "rho_k": mean value of \mu in the Gaussian-Wishart framework
    # "beta_k": variance scale of \mu in the Gaussian-Wishart framework
    # "nu_k": df in Wishart distribution of \Lambda
    # "Phi_k": covariance matrix in Wishart distribution of \Lambda
    
    # "log_Lambda": E[ln|\Lambda_{k}|] in estimating r_nk
    # "log_pi": E[ln\pi_{k}] in estimating r_nk
    
    ##################
    # Initialization #
    ##################
    
    n <- nrow(data)
    d <- ncol(data)
    
    # Initialize parameters
    Phi0_inv <- solve(Phi0)
    L <- rep(-Inf, max_iter) # store ELBO value
    r_nk = log_r_nk = log_rho_nk <- matrix(0, nrow = n, ncol = m)
    x_bar_k <- matrix(0, nrow = d, ncol = m)
    S_k = Phi_k <- array(0, c(d, d, m))
    log_pi = log_Lambda <- rep(0, m)
    
    lambda <- rep(lambda0, m)
    rho_k <- t(kmeans(data, m, nstart = 25)$centers)
    beta_k <- rep(beta0, m)
    nu_k <- rep(nu0, m)
    log_pi  <- digamma(lambda) - digamma(sum(lambda))
    for (k in 1:m) {
        Phi_k[, , k] <- Phi0  # Scale matrix for Wishart
        log_Lambda[k] <- sum(digamma((nu_k[k] + 1 - c(1:d))/2)) + 
          d*log(2) + log(det(Phi_k[, , k]))
    }
    
    
    ##################
    # Main Iteration #
    ##################
    
    for (i in 2:max_iter) {
      
        # Compute responsibilities
        for (k in 1:m) {
            diff <- sweep(data, MARGIN = 2, STATS = rho_k[, k], FUN = "-")
            log_rho_nk[, k] <- omega * log_pi[k] + 0.5 * omega * log_Lambda[k] - 
              0.5 * omega * (d/beta_k[k]) - 
              0.5 * omega * nu_k[k] * diag(diff %*% Phi_k[, , k] %*% t(diff))
            log_rho_nk[, k] <- log_rho_nk[, k] / omega
        }
        
        Z <- apply(log_rho_nk, 1, log_sum_exp)
        log_r_nk <- log_rho_nk - Z
        r_nk <- apply(log_r_nk, 2, exp)
        
        # Update parameters, M-step
        N_k <- colSums(r_nk) + 1e-10
        for (k in 1:m) {
            x_bar_k[, k] <- (r_nk[, k] %*% data) / N_k[k]
            x_cen <- sweep(data, MARGIN = 2, STATS = x_bar_k[, k], FUN = "-")
            S_k[, , k] <- t(x_cen) %*% (x_cen * r_nk[, k]) / N_k[k]
        }
        
        lambda <- lambda0 + omega * N_k
        pi_k <- (lambda0 + omega * N_k) / (m * lambda0 + omega * n)
        beta_k <- beta0 + omega * N_k
        nu_k <- nu0 + omega * N_k # have to +1 ?
        
        for (k in 1:m) {
            rho_k[, k] <- (1/beta_k[k]) * (beta0 * rho0 + omega * N_k[k] * x_bar_k[, k])  
            Phi_k[, , k] <- Phi0_inv + omega * N_k[k] * S_k[, , k] + 
              ((beta0 * N_k[k])/(beta0 + N_k[k])) * 
              tcrossprod((x_bar_k[, k] - rho0))    
            Phi_k[, , k] <- solve(Phi_k[, , k])
        }
        
        log_pi <- digamma(lambda) - digamma(sum(lambda))                      
        for (k in 1:m) {                                              
            log_Lambda[k] <- sum(digamma((nu_k[k] + 1 - 1:d)/2)) + 
              d * log(2) + log(det(Phi_k[, , k])) 
        }
        
        
        ##################
        # Calculate ELBO #
        ##################
        
        lb_px = lb_pml = lb_pml2 = lb_qml <- 0
        for (k in 1:m) {
          
            lb_px <- lb_px + N_k[k] * (log_Lambda[k] - d/beta_k[k] - nu_k[k] * 
                                         matrix.trace(S_k[, , k] %*% Phi_k[, , k]) - 
                                         nu_k[k] * t(x_bar_k[,k] - rho_k[, k]) %*% 
                                         Phi_k[, , k] %*% (x_bar_k[,k] - rho_k[, k]) - d*log(2*pi) ) 
            
            lb_pml <- lb_pml + d * log(beta0/(2 * pi)) + log_Lambda[k] - 
              (d * beta0)/beta_k[k] - beta0 * nu_k[k] * t(rho_k[, k] - rho0) %*% 
              Phi_k[, , k] %*% (rho_k[, k] - rho0)    
            
            lb_pml2 <- lb_pml2 + nu_k[k] * matrix.trace(Phi0_inv %*% Phi_k[, , k]) 
            
            lb_qml <- lb_qml + 0.5 * log_Lambda[k] + 
              0.5 * d * log(beta_k[k]/(2 * pi)) - 
              0.5 * d - (-logB(W = Phi_k[, , k], nu = nu_k[k]) - 
                           0.5 * (nu_k[k] - d - 1) * log_Lambda[k] + 0.5 * nu_k[k] * d)
          
        }
        
        lb_px <- 0.5 * omega * lb_px             
        lb_pml <- 0.5 * lb_pml + m * logB(W = Phi0, nu = nu0) + 
          0.5 * (nu0 - d - 1) * sum(log_Lambda) - 0.5 * lb_pml2
        lb_pz <- sum((omega * r_nk) %*% log_pi)
        lb_qz <- sum(omega * r_nk * log_r_nk)
        lb_pp <- sum((lambda0 - 1) * log_pi) + lgamma(sum(m * lambda0)) -
          m * sum(lgamma(lambda0))
        lb_qp <- sum((lambda - 1) * log_pi) + lgamma(sum(lambda)) - 
          sum(lgamma(lambda))
        
        L[i] <- lb_px + lb_pz + lb_pp + lb_pml - lb_qz - lb_qp - lb_qml
        
        # Check for convergence
        if (i > 2 && abs(L[i] - L[i-1]) < tol) {
            break
        }
      
    }
    
    # Make sure no switch problem
    if (length(pi_k) == 2 & pi_k[1] < 0.5) {
        pi_k <- 1 - pi_k
        lambda <- rev(lambda)
    }
    
    obj <- structure(list(X = data, m = m, n = n, d = d, pi_k = pi_k, 
                          lambda = lambda, r_nk = r_nk,  rho = rho_k, 
                          Phi = Phi_k, beta = beta_k, nu = nu_k, L = L[2:i]), 
                     class = "vb_gmm_pac")
    return(obj)
  
}

# Update omega: Sequential TVB
update_omega_pi <- function(x, initial_alpha = 0, alpha = 0.05, 
                            tolerance = 1e-3, max_iter = 100, 
                            B = 500, target_coverage = 0.95) {
  
    N <- nrow(x)
    half_N <- floor(N/2)
    alpha_param <- initial_alpha
    k <- 0
    omega_all <- c()
    converged <- FALSE
    
    while (k < max_iter && !converged) {
      
        omega <- exp(alpha_param)
        
        train_indices <- sample(1:N, half_N, replace = FALSE)
        test_indices <- setdiff(1:N, train_indices)
        
        x_train <- x[train_indices, ]
        vb_result <- vb_gmm(x_train, omega = omega)
        pi_est <- vb_result$pi_k
        var_est <- vb_result$pi_k * (1 - vb_result$pi_k) / (sum(vb_result$lambda) + 1)
        
        # Bootstrap samples
        x_test <- x[test_indices, ]
        boot_samples <- boot(data = 1:(N - half_N), statistic = function(data, indices) indices, R = B)
        
        # Calculate empirical coverage probability
        coverage_count <- 0
        for (b in 1:B) {
          
            sample_indices <- boot_samples$t[b, ]
            X_sample <- lapply(sample_indices, function(i) x_test[i, ])
            
            # Get posterior distribution for the current bootstrap sample
            posterior <- vb_gmm(do.call(rbind, X_sample), omega = omega)
            
            ci <- dir_beta_ci(posterior$lambda[1], posterior$lambda[2], 
                              conf_level = target_coverage)
            if (pi_est[1] < ci[2] && pi_est[1] > ci[1]) {
                coverage_count <- coverage_count + 1
            }
          
        }
        
        empirical_coverage <- coverage_count / B
        
        if (empirical_coverage > target_coverage | 
            abs(empirical_coverage - target_coverage) < tolerance) {
            converged <- TRUE
        } else {
            # Update omega
            alpha_param <- alpha_param + (empirical_coverage - target_coverage) * (k + 1)^(-0.51)
        }
        
        omega_all[k + 1] <- omega
        k <- k + 1
      
    }
    
    return(list(omega = omega, 
                omega_all = omega_all,
                alpha_param = alpha_param,
                cov_emp = empirical_coverage))
  
}

# Update omega: Grid search TVB (Building dictionary)
boot_tvb <- function(omega, x, B = 500, target_coverage = 0.95) {
  
    N <- nrow(x)
    half_N <- floor(N/2)
    train_indices <- sample(1:N, half_N, replace = FALSE)
    test_indices <- setdiff(1:N, train_indices)
    boot_samples <- boot(data = 1:(N - half_N), 
                         statistic = function(data, indices) indices, R = B)
    x_train <- x[train_indices, ]
    x_test <- x[test_indices, ]
    
    vb_result <- vb_gmm(x_train, omega = omega) # Estimation of training data
    pi_est <- vb_result$pi_k
    coverage_count <- 0
    
    for (b in 1:B) {
      
        sample_indices <- boot_samples$t[b, ]
        X_sample <- lapply(sample_indices, function(i) x_test[i, ])
        vb_boot_result <- vb_gmm(do.call(rbind, X_sample), omega = omega)
        ci <- dir_beta_ci(vb_boot_result$lambda[1], vb_boot_result$lambda[2],
                          conf_level = target_coverage)
        if (pi_est[1] < ci[2] && pi_est[1] > ci[1]) {
            coverage_count <- coverage_count + 1
        }
      
    }
    
    cover_boot <- coverage_count / B
    
    return(cover_boot)
  
}


###################
# Main simulation #
###################

mu1 <- c(0, 0)
mu2 <- c(center, center)
sigma1 <- diag(d)
sigma2 <- diag(d)
target_coverage <- 0.95
omega_grid <- exp(seq(from = log(0.001), to = log(1), length.out = 100))

# Simulate data
seed_use <- seed * 40 + 10000 * out.id
set.seed(seed_use)

component <- rbinom(num_sample, 1, true_pi)
data <- matrix(0, num_sample, d)
for (k in 1:num_sample) {
    if (component[k] == 1) {
        data[k, ] <- rmvnorm(1, mu1, sigma1)
    } else {
        data[k, ] <- rmvnorm(1, mu2, sigma2)
    }
}

# Three strategies to obtain omega 

# S1: Original VB
# No omega need to select!

# S2: Sequential TVB
omega_tvb_seq <- update_omega_pi(data)

# S3: Pre-calculation with random splitting and Bootstrap
emp_cov <- numeric(length(omega_grid))
for (k in 1:length(omega_grid)) {
    emp_cov[k] <- boot_tvb(omega_grid[k], data)
}
diff <- abs(emp_cov - target_coverage)
idx_all <- which(diff == min(diff))
idx_best <- max(idx_all) # Pick the largest omega when some omega's share the coverage
omega_tvb_grid <- omega_grid[idx_best]


# Compare the result
# Original VB
vb_est <- vb_gmm(data)
ci_vb <- dir_beta_ci(vb_est$lambda[1], vb_est$lambda[2])
cov_vb <- ci_vb[1] <= true_pi && true_pi <= ci_vb[2]

# Sequential TVB
vb_tvb_seq <- vb_gmm(data, omega = omega_tvb_seq$omega)
ci_tvb_seq <- dir_beta_ci(vb_tvb_seq$lambda[1], vb_tvb_seq$lambda[2])
cov_tvb_seq <- ci_tvb_seq[1] <= true_pi && true_pi <= ci_tvb_seq[2]

# Grid-search TVB
vb_tvb_grid <- vb_gmm(data, omega = omega_tvb_grid)
ci_tvb_grid <- dir_beta_ci(vb_tvb_grid$lambda[1], vb_tvb_grid$lambda[2])
cov_tvb_grid <- ci_tvb_grid[1] <= true_pi && true_pi <= ci_tvb_grid[2]


# Save the data
save(
  
    data,
    
    vb_est,
    ci_vb,
    cov_vb,
    
    vb_tvb_seq,
    ci_tvb_seq,
    cov_tvb_seq,
    omega_tvb_seq,
    
    vb_tvb_grid,
    ci_tvb_grid,
    cov_tvb_grid,
    omega_tvb_grid,
    
    file = paste0("results/", name.batch, "/", seed_use, ".Rdata")
  
)













