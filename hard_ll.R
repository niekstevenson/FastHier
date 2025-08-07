loglik_sinmix <- function(Theta, data, model = NULL) {
  X <- data$X
  y <- data$y
  sigma_y <- data$sigma_y

  if (is.null(X) || is.null(y) || is.null(sigma_y)) {
    stop("data must contain X (N x 2), y (length N), and sigma_y (scalar > 0).")
  }

  Theta <- as.matrix(Theta)
  M <- nrow(Theta); D <- ncol(Theta)
  if (D != 5) stop(sprintf("Theta must have 5 columns, got %d", D))
  if (ncol(X) != 2) stop(sprintf("X must have 2 columns, got %d", ncol(X)))
  N <- nrow(X)
  if (length(y) != N) stop(sprintf("length(y)=%d does not match nrow(X)=%d", length(y), N))
  if (!is.finite(sigma_y) || sigma_y <= 0) stop("sigma_y must be a positive finite scalar.")

  # Extract parameters: theta1, theta2, theta3, theta4, theta5
  theta1 <- Theta[, 1]  # exponential amplitude
  theta2 <- Theta[, 2]  # exponential rate
  theta3 <- Theta[, 3]  # sine amplitude  
  theta4 <- Theta[, 4]  # sine frequency
  theta5 <- Theta[, 5]  # sine phase

  # Nonlinear model: y = theta1 * exp(theta2 * x1) + theta3 * sin(theta4 * x2 + theta5) + noise
  x1 <- X[, 1]  # first covariate
  x2 <- X[, 2]  # second covariate
  
  # Compute mean function for each parameter vector (M x N matrix)
  # For each row of Theta (each parameter vector), compute the mean at each observation
  mu <- matrix(NA_real_, M, N)
  for (m in 1:M) {
    exp_term <- theta1[m] * exp(theta2[m] * x1)
    sin_term <- theta3[m] * sin(theta4[m] * x2 + theta5[m])
    mu[m, ] <- exp_term + sin_term
  }

  # Likelihood calculation
  const_per_obs <- -0.5 * log(2 * pi) - log(sigma_y)
  inv_var_half  <- -0.5 / (sigma_y * sigma_y)

  resid <- sweep(mu, 2L, y, `-`)           # M x N
  const_per_obs * N + inv_var_half * rowSums(resid * resid)
}

make_prior_mvn <- function(D = 5, scale = 2, param_names = NULL, Theta_template = NULL) {
  if (!is.null(Theta_template)) {
    D <- ncol(Theta_template)
    if (is.null(param_names)) param_names <- colnames(Theta_template)
  }
  if (is.null(D) || D != 5) stop("This model requires exactly 5 parameters.")
  if (is.null(param_names)) param_names <- c("theta1", "theta2", "theta3", "theta4", "theta5")

  # Tailored priors for each parameter type
  mu <- c(
    theta1 = 1.0,   # exponential amplitude (positive expected)
    theta2 = 0.0,   # exponential rate (can be positive/negative)
    theta3 = 1.0,   # sine amplitude (positive expected)
    theta4 = 1.0,   # sine frequency (positive expected)
    theta5 = 0.0    # sine phase (centered at 0)
  )
  
  # Covariance matrix with mild correlations
  Sigma <- diag(scale^2, 5)
  # Add some correlation between related parameters
  Sigma[1, 3] <- Sigma[3, 1] <- 0.3 * scale^2  # amplitude parameters slightly correlated
  Sigma[2, 4] <- Sigma[4, 2] <- 0.2 * scale^2  # rate and frequency slightly correlated
  
  dimnames(Sigma) <- list(param_names, param_names)

  list(mu = mu, Sigma = Sigma)
}

# --- Synthetic data generator for the 5-parameter model ----
gen_data <- function(N = 200, sigma_y = 0.2, prior_scale = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Generate 2D covariates
  X <- matrix(runif(N * 2, -1, 1), N, 2)  # x1, x2 in [-1, 1]
  
  # True parameter values (moderate values to ensure identifiability)
  theta_true <- c(
    theta1 = 1.5,   # exponential amplitude
    theta2 = 0.8,   # exponential rate  
    theta3 = 1.2,   # sine amplitude
    theta4 = 2.0,   # sine frequency
    theta5 = 0.5    # sine phase
  )
  
  # Generate true mean function
  x1 <- X[, 1]
  x2 <- X[, 2]
  eta <- theta_true[1] * exp(theta_true[2] * x1) + 
         theta_true[3] * sin(theta_true[4] * x2 + theta_true[5])
  
  # Add noise
  y <- eta + rnorm(N, 0, sigma_y)

  # Prior setup
  prior <- make_prior_mvn(D = 5, scale = prior_scale)

  list(N = N, D = 5, X = X, y = y, 
       mu = prior$mu, Sigma = prior$Sigma,
       sigma_y = sigma_y, theta_true = theta_true, eta_true = eta)
}
