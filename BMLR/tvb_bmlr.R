# Bayesian mixture linear regression
# Variational calibration
# Pre-calculation strategy

#############
# Preparing #
#############

# Working directory
setwd("~/BMLR/")

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

# Option
option_list <- list(
  
    make_option("--name_batch", type = "character", default = NA, action = "store"),
    make_option("--seed", type = "numeric", default = NA, action = "store"),
    make_option("--cluster", type = "numeric", default = NA, action = "store"),
    make_option("--dataset", type = "numeric", default = NA, action = "store"),
    make_option("--dim", type = "numeric", default = NA, action = "store"),
    make_option("--SNR", type = "numeric", default = NA, action = "store"),
    make_option("--mixing", type = "numeric", default = NA, action = "store")
  
)

opt <- parse_args(OptionParser(option_list = option_list))

name.batch <- opt$name_batch
seed <- opt$seed
m <- opt$cluster
N <- opt$dataset
p <- opt$dim
SNR <- opt$SNR
true_pi <- opt$mixing

# For code test only
# name.batch <- "sim0"
# seed <- 1
# m <- 2
# N <- 100
# p <- 10
# SNR <- 1
# true_pi <- 0.65

out.id <- as.numeric(gsub("sim", "", name.batch))

# Make folder
dir.create(file.path("results", name.batch), showWarnings = FALSE, recursive = TRUE)


########################
# Supporting functions #
########################

# Stable numerical performance
log_sum_exp <- function(x) {
  
    # Computes log(sum(exp(x))
    offset <- max(x)
    s <- log(sum(exp(x - offset))) + offset
    i <- which(!is.finite(s))
    if (length(i) > 0) { s[i] <- offset }
    return(s)
  
}

# A general-purpose adder:
add_func <- function(mat_list) {
    Reduce("+", mat_list)
}

# CI for mixing parameter
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

# VB for Gibbs posterior
vb_bmlr <- function(x, y, m = 2, omega = 1, 
                    delta_0 = rep(1/k, k), lambda = 1,
                    alpha_0 = 0.1, beta_0 = 0.1, 
                    max_iter = 300, tol = 1e-5) {
  
    ##################
    # Initialization #
    ##################
    
    # Dimension and necessary values
    N <- length(x)
    p <- ncol(x[[1]])
    xx <- lapply(X = x, FUN = function(x) crossprod(x))
    yy <- unlist(lapply(X = y, FUN = function(y) crossprod(y)))
    xy <- lapply(X = 1:N, FUN = function(i) crossprod(x[[i]], y[[i]]))
    len_y <- unlist(lapply(X = y, FUN = function(y) length(y))) # Obs for each sample
    L <- rep(-Inf, max_iter)
    r_nk = log_r_nk = log_rho_nk <- matrix(0, nrow = N, ncol = m) # Responsibility
    E_ww <- vector("numeric", length = m)
    
    coefficients_list <- lapply(1:N, function(i) {
        fit <- lm(y[[i]] ~ x[[i]] + 0) # Linear regression without intercept
        return(coef(fit))
    })
    
    alpha_k <- rep(alpha_0 + p/2, m)
    W <- do.call(cbind, coefficients_list)
    cl  <- stats::kmeans(t(W), m, nstart = 25)
    m_k <- t(cl$centers) # Mean of each cluster, initialization
    S_k <- array(0, dim = c(p, p, m))
    for (k in 1:m) { 
        S_k[, , k] <- solve(diag(2, p))
    }      
    # Scale of precision matrix
    beta_k <- rep(beta_0, m)   
    # Dirichlet parameter
    delta_k <- delta_0                        
    # Expectation of log Dirichlet  
    e_log_pi <- digamma(delta_k) - digamma(sum(delta_k))     
    mk_Sk <- lapply(X = 1:m, function(k) tcrossprod(m_k[, k]) + S_k[, , k])
    
    
    #############
    # Iteration #
    #############
    for (i in 2:max_iter) {
      
        ##-------------------
        # Variational E-Step
        ##-------------------
        for (k in 1:m) {
            log_rho_nk[, k] <- omega * e_log_pi[k] + omega * lambda * sapply(1:N, function(n) 
              m_k[, k] %*% xy[[n]] - 0.5 * matrix.trace(xx[[n]] %*% mk_Sk[[k]]))
            log_rho_nk[, k] <- log_rho_nk[, k] / omega
        }
        Z <- apply(log_rho_nk, 1, log_sum_exp)
        log_r_nk <- log_rho_nk - Z
        r_nk <- apply(log_r_nk, 2, exp)
        
        ##-------------------
        # Variational M-Step
        ##-------------------
        delta_k <- delta_0 + omega * colSums(r_nk)
        for (k in 1:m) {
          
            w_XX <- lapply(X = 1:N, function(x) xx[[x]] * r_nk[x, k])
            S_k[, , k] <- solve(diag(alpha_k[k]/beta_k[k], p) + omega * lambda * 
                                  add_func(w_XX))
            # Update mean for Gaussian
            w_Xy <- lapply(X = 1:N, function(x) xy[[x]] * r_nk[x, k])
            m_k[, k] <- omega * lambda * S_k[, , k] %*% add_func(w_Xy)
            # Update \beta_k parameter for Gamma
            E_ww[k] <- crossprod(m_k[, k]) + matrix.trace(S_k[, , k])
            beta_k[k]  <- beta_0 + 0.5 * E_ww[k]
          
        }
        # Expected value of mixing proportions
        pi_k <- (delta_0 + omega * colSums(r_nk)) / (m * delta_0 + omega * N)
        # Expectations over \ln\pi
        e_log_pi <- digamma(delta_k) - digamma(sum(delta_k))
        # Compute expectation of E[a]
        E_alpha <- alpha_k / beta_k
        
        ##------------------------
        # Variational lower bound
        ##------------------------
        mk_Sk <- lapply(X = 1:m, function(k) tcrossprod(m_k[, k]) + S_k[, , k])
        
        lb_p_y <- omega * (-0.5 * sum(len_y) * log(2 * pi * (1/lambda)) - 0.5 * lambda * sum(yy) + 
                             sum(sapply(1:m, function(k) lambda * (sum(sapply(1:N, function(n) 
                               r_nk[n, k] * (m_k[, k] %*% xy[[n]] - 
                                               0.5 * matrix.trace(xx[[n]] %*% mk_Sk[[k]]))))))))
        lb_p_w <- sum(-0.5 * p * log(2 * pi) + 0.5 * p * (digamma(alpha_k) - 
                                                            log(beta_k)) - 0.5 * E_alpha * E_ww)
        lb_p_c <- sum((omega * r_nk) %*% e_log_pi)   
        lb_p_pi <- sum((delta_0 - 1) * e_log_pi) + lgamma(sum(delta_0)) - 
          sum(lgamma(delta_0))
        lb_p_tau <- sum(alpha_0 * log(beta_0) + (alpha_0 - 1) * (digamma(alpha_k) - 
                                                                   log(beta_k)) - beta_0 * E_alpha - lgamma(alpha_0))
        lb_q_c <- sum(omega * r_nk * log_r_nk)  
        lb_q_pi <- sum((delta_k - 1) * e_log_pi) + lgamma(sum(delta_k)) - 
          sum(lgamma(delta_k))
        lb_q_w   <- sum(-0.5 * log(sapply(X = 1:m, function(k) det(S_k[, , k]))) - 
                          0.5 * p * (1 + log(2 * pi)))
        lb_q_tau <- sum(-lgamma(alpha_k) + (alpha_k - 1) * digamma(alpha_k) + 
                          log(beta_k) - alpha_k)
        
        # ELBO update
        L[i] <- lb_p_y + lb_p_c + lb_p_pi + lb_p_w + lb_p_tau - lb_q_c - 
          lb_q_pi - lb_q_w - lb_q_tau
        
        # Check for convergence
        if (i > 2 && abs(L[i] - L[i-1]) < tol) {
            break
        }
      
    }
    
    # Make sure no switch problem
    if (length(pi_k) == 2 & pi_k[1] < 0.5) {
        pi_k <- 1 - pi_k
        delta_k <- rev(delta_k)
    }
    
    return_obj <- list(
        m = m_k, 
        S = S_k, 
        delta = delta_k, 
        r_nk = r_nk, 
        lambda = lambda, 
        pi_k = pi_k, 
        beta = beta_k, 
        alpha = alpha_k, 
        L = L[2:i], 
        m = m,
        N = N, 
        dimension = p
    )
    
    return(return_obj)
  
}

# Sequential TVB
update_omega_beta <- function(x, y, initial_alpha = 0, alpha = 0.05, 
                              tolerance = 1e-3, max_iter = 100, 
                              B = 500, target_coverage = 0.95) {
  
    N <- length(x)
    half_N <- floor(N/2)
    alpha_param <- initial_alpha
    k <- 0
    omega_all <- c()
    converged <- FALSE
    
    while (k < max_iter && !converged) {
      
        omega <- exp(alpha_param)
        
        # Split the samples
        train_indices <- sample(1:N, half_N, replace = FALSE)
        test_indices <- setdiff(1:N, train_indices)
        
        # Create training lists
        x_train <- x[train_indices]
        y_train <- y[train_indices]
        
        # Create test lists
        x_test <- x[test_indices]
        y_test <- y[test_indices]
        
        # Train on the first half
        vb_result <- vb_bmlr(x_train, y_train, omega = omega)
        mean_est <- sum(vb_result$m)
        
        # Bootstrap samples from the test set
        test_N <- length(test_indices)
        boot_samples <- boot(data = 1:test_N, statistic = function(data, indices) indices, R = B)
        
        # Calculate empirical coverage probability
        coverage_count <- 0
        for (b in 1:B) {
          
            sample_indices <- boot_samples$t[b, ]
            X_sample <- lapply(sample_indices, function(i) x_test[[i]])
            y_sample <- lapply(sample_indices, function(i) y_test[[i]])
            
            # Get posterior distribution for the current bootstrap sample
            posterior <- vb_bmlr(X_sample, y_sample, omega = omega)
            mean_bs <- sum(posterior$m)
            var_bs <- sum(posterior$S)
            
            ci <- normal_ci(mean_bs, var_bs, conf_level = target_coverage)
            
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

# Grid-search TVB (Building dictionary)
boot_tvb <- function(omega, x, y, B = 500, target_coverage = 0.95) {
  
    N <- length(x)
    half_N <- floor(N/2)
    train_indices <- sample(1:N, half_N, replace = FALSE)
    test_indices <- setdiff(1:N, train_indices)
    test_N <- length(test_indices)
    boot_samples <- boot(data = 1:test_N, statistic = function(data, indices) indices, R = B)
    
    x_train <- x[train_indices]
    y_train <- y[train_indices]
    x_test <- x[test_indices]
    y_test <- y[test_indices]
    
    vb_result <- vb_bmlr(x_train, y_train, omega = omega)
    mean_est <- sum(vb_result$m)
    
    coverage_count <- 0
    
    for (b in 1:B) {
      
        sample_indices <- boot_samples$t[b, ]
        X_sample <- lapply(sample_indices, function(i) x_test[[i]])
        y_sample <- lapply(sample_indices, function(i) y_test[[i]])
        
        # Get posterior distribution for the current bootstrap sample
        posterior <- vb_bmlr(X_sample, y_sample, omega = omega)
        mean_bs <- sum(posterior$m)
        var_bs <- sum(posterior$S)
        
        ci <- normal_ci(mean_bs, var_bs, conf_level = target_coverage)
        
        if (mean_est < ci[2] && mean_est > ci[1]) {
            coverage_count <- coverage_count + 1
        }
      
    }
    
    cover_boot <- coverage_count / B
    
    return(cover_boot)
  
}


###################
# Main simulation #
###################

# True parameters
set.seed(1)
cluster_probs <- c(true_pi, 1 - true_pi)  # Normalize to sum to 1
cluster_assignments <- sample(1:m, N, replace = TRUE, prob = cluster_probs)
a_0 <- 1
b_0 <- 1
tau <- c(1, 1)
cov_mat1 <- matrix(0, nrow = p, ncol = p)
cov_mat2 <- matrix(0, nrow = p, ncol = p)
diag(cov_mat1) <- 1 / tau[1]
diag(cov_mat2) <- 1 / tau[2]

w_1 <- rmvnorm(1, rep(0, p), cov_mat1)
w_2 <- rmvnorm(1, rep(0, p), cov_mat2)

w <- matrix(0, nrow = m, ncol = p)
w[1, ] <- w_1
w[2, ] <- w_2

m_var <- c(1000, 1500, 2000, 2500, 3000) # Different sample size
size <- length(m_var)
target_coverage <- 0.95
omega_grid <- exp(seq(from = log(0.001), to = log(1), length.out = 100))

# Data generation and estimation
cover_vb <- numeric(size)
cover_tvb_grid <- numeric(size)
cover_tvb_seq <- numeric(size)
len_vb <- matrix(0, nrow = size, ncol = 2)
len_tvb_grid <- matrix(0, nrow = size, ncol = 2)
len_tvb_seq <- matrix(0, nrow = size, ncol = 2)
omega_tvb_seq <- numeric(size)
omega_tvb_grid <- numeric(size)
data_x <- vector("list", size)
data_y <- vector("list", size)
est_vb <- vector("list", size)
est_tvb_grid <- vector("list", size)
est_tvb_seq <- vector("list", size)

# Simulate data
seed_use <- seed * 40 + 10000 * out.id
set.seed(seed_use)

for (j in 1:size) {
  
    x <- lapply(1:N, function(i) matrix(rnorm(m_var[j] * p), nrow = m_var[j], ncol = p))
    y <- vector("list", N)
    
    noise_cluster_1 <- rnorm(m_var[j])
    noise_cluster_2 <- rnorm(m_var[j])
    noise <- rbind(noise_cluster_1, noise_cluster_2)
    
    for (l in 1:N) {
        K <- cluster_assignments[l]  
        noise.temp <- (noise[K, ] / norm(as.matrix(noise[K, ]), type = 'F')) *
          norm(x[[l]] %*% w[K, ], type = 'F') / SNR
        y[[l]] <- x[[l]] %*% w[K, ] + noise.temp  # Generate response
    }
    
    data_x[[j]] <- x
    data_y[[j]] <- y
    
    # VB estimation
    est <- vb_bmlr(x, y, omega = 1)
    mean_val <- sum(est$m)
    variance_val <- sum(est$S)
    ci_vb <- normal_ci(mean_val, variance_val, conf_level = 0.95)
    if (ci_vb[1] < sum(w) & ci_vb[2] > sum(w)) {
        cover_vb[j] <- 1
    }
    
    est_vb[[j]] <- est
    len_vb[j, 1] <- as.numeric(ci_vb[1])
    len_vb[j, 2] <- as.numeric(ci_vb[2])
    
    # Grid-search TVB
    emp_cov <- numeric(length(omega_grid))
    for (k in 1:length(omega_grid)) {
        emp_cov[k] <- boot_tvb(omega_grid[k], x, y)
    }
    diff <- abs(emp_cov - target_coverage)
    idx_all <- which(diff == min(diff))
    idx_best <- max(idx_all) # Pick the largest omega when some omega's share the coverage
    omega_select <- omega_grid[idx_best]
    omega_tvb_grid[j] <- omega_select
  
    est <- vb_bmlr(x, y, omega = omega_select)
    mean_val_tvb2 <- sum(est$m)
    variance_val_tvb2 <- sum(est$S)
    ci_tvb_grid <- normal_ci(mean_val_tvb2, variance_val_tvb2, conf_level = 0.95)
    if (ci_tvb_grid[1] < sum(w) & ci_tvb_grid[2] > sum(w)) {
        cover_tvb_grid[j] <- 1
    }
    
    est_tvb_grid[[j]] <- est
    len_tvb_grid[j, 1] <- as.numeric(ci_tvb_grid[1])
    len_tvb_grid[j, 2] <- as.numeric(ci_tvb_grid[2])
    
    # Sequential TVB
    omega_select <- update_omega_beta(x, y)
    omega_tvb_seq[j] <- omega_select$omega
    est <- vb_bmlr(x, y, omega = omega_select$omega)
    mean_val_tvb1 <- sum(est$m)
    variance_val_tvb1 <- sum(est$S)
    ci_tvb_seq <- normal_ci(mean_val_tvb1, variance_val_tvb1, conf_level = 0.95)
    if (ci_tvb_seq[1] < sum(w) & ci_tvb_seq[2] > sum(w)) {
        cover_tvb_seq[j] <- 1
    }
    
    est_tvb_seq[[j]] <- est
    len_tvb_seq[j, 1] <- as.numeric(ci_tvb_seq[1])
    len_tvb_seq[j, 2] <- as.numeric(ci_tvb_seq[2])
  
}

# Save the results
save(
  
    data_x,
    data_y,
    
    est_vb,
    len_vb,
    
    est_tvb_seq,
    len_tvb_seq,
    omega_tvb_seq,
    
    est_tvb_grid,
    len_tvb_grid,
    omega_tvb_grid,
    
    w,
    omega_grid,
    
    file = paste0("results/", name.batch, "/", seed_use, ".Rdata")
  
)











