#############
# Preparing #
#############

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
library(datasets)


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

dir_beta_ci <- function(lambda, idx, conf_level = 0.95) {
  
    alpha1 <- lambda[idx]
    alpha2 <- sum(lambda) - alpha1
    
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

# CI for beta
normal_ci <- function(mean, variance, conf_level = 0.95) {
  
    sd <- sqrt(variance)
    
    critical_value <- qnorm((1 + conf_level) / 2)
    
    # Calculate lower and upper bounds
    margin_of_error <- critical_value * sd
    lower_bound <- mean - margin_of_error
    upper_bound <- mean + margin_of_error
    
    # Return the results
    return(list(
      ci_lower = lower_bound,
      ci_upper = upper_bound
    ))
  
}


##################
# Main functions #
##################

vb_gmm_pac <- function(data, m = 2, lambda0 = 1/m, 
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
    sorted_idx <- order(pi_k, decreasing = TRUE)
    pi_k <- pi_k[sorted_idx]
    lambda <- lambda[sorted_idx]
    Phi_k <- Phi_k[, , sorted_idx]
    beta_k <- beta_k[sorted_idx]
    nu_k <- nu_k[sorted_idx]
    rho_k <- rho_k[, sorted_idx]
    r_nk <- r_nk[, sorted_idx]
    
    obj <- structure(list(X = data, m = m, n = n, d = d, pi_k = pi_k, 
                          lambda = lambda, r_nk = r_nk,  rho = rho_k, 
                          Phi = Phi_k, beta = beta_k, nu = nu_k, L = L[2:i]), 
                     class = "vb_gmm_pac")
    return(obj)
  
}

# Bootstrap for a specific omega
boot_gmm <- function(omega, x, m = 2, B = 500) {
    
    set.seed(42)
    N <- nrow(x)
    half_N <- floor(N/2)
    train_indices <- sample(1:N, half_N, replace = FALSE)
    test_indices <- setdiff(1:N, train_indices)
    boot_samples <- boot(data = 1:(N - half_N), 
                         statistic = function(data, indices) indices, R = B)
    x_train <- x[train_indices, ]
    x_test <- x[test_indices, ]
    
    vb_result <- vb_gmm_pac(x_train, m = m, omega = omega) # Estimation of training data
    
    vb_boot <- list()
    for (b in 1:B) {
      
        sample_indices <- boot_samples$t[b, ]
        X_sample <- lapply(sample_indices, function(i) x_test[i, ])
        vb_boot_result <- vb_gmm_pac(do.call(rbind, X_sample), m = m, omega = omega)
        vb_boot[[b]] <- vb_boot_result
      
    }

    return(list(
        vb_results = vb_result,
        vb_boot = vb_boot
    ))
  
}

# Build dictionary
omega_bootstrap <- function(omega_grid, x, m = 2, B = 500) {
  
    omega_boot <- list()
    for (i in 1:length(omega_grid)) {
        omega_boot[[i]] <- boot_gmm(omega_grid[i], x, m = m, B = B)
    }
    
    return(omega_boot)
    
}

# Sequential TVB functions
# For sum(beta)
update_omega_beta <- function(x, initial_alpha = 0, alpha = 0.05, idx = 1,
                              tolerance = 1e-3, max_iter = 100, m = 2,
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
        vb_result <- vb_gmm_pac(x_train, m = m, omega = omega)
        mean_est <- sum(vb_result$rho[, idx])
        
        # Bootstrap samples
        x_test <- x[test_indices, ]
        boot_samples <- boot(data = 1:(N - half_N), statistic = function(data, indices) indices, R = B)
        
        # Calculate empirical coverage probability
        coverage_count <- 0
        for (b in 1:B) {
          
            sample_indices <- boot_samples$t[b, ]
            X_sample <- lapply(sample_indices, function(i) x_test[i, ])
            
            # Get posterior distribution for the current bootstrap sample
            posterior <- vb_gmm_pac(do.call(rbind, X_sample), m = m, omega = omega)
            mean_boot <- sum(posterior$rho[, idx])
            var_boot <- sum(solve((posterior$Phi[, , idx] * posterior$nu[idx]) * posterior$beta[idx]))
            
            ci <- normal_ci(mean_boot, var_boot, conf_level = target_coverage)
            if (mean_est < ci[2] && mean_est > ci[1]) {
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

# For pi
update_omega_pi <- function(x, initial_alpha = 0, alpha = 0.05, idx = 1,
                            tolerance = 1e-3, max_iter = 100, m = 2,
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
        vb_result <- vb_gmm_pac(x_train, m = m, omega = omega)
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
            posterior <- vb_gmm_pac(do.call(rbind, X_sample), m = m, omega = omega)
            
            ci <- dir_beta_ci(posterior$lambda, idx, conf_level = target_coverage)
            if (pi_est[idx] < ci[2] && pi_est[idx] > ci[1]) {
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





