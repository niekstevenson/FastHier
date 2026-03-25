normalize_weights <- function(w) {
  w <- pmax(as.numeric(w), 0)
  w / sum(w)
}

weighted_resample_rows <- function(x, w, n_draws = NULL, seed = NULL) {
  x <- as.matrix(x)
  w <- normalize_weights(w)
  if (is.null(n_draws)) {
    n_draws <- nrow(x)
  }
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }
  idx <- sample.int(
    nrow(x),
    size = as.integer(n_draws),
    replace = TRUE,
    prob = w
  )
  x[idx, , drop = FALSE]
}

smc_posteriors <- function(smc, n_draws = NULL, seed = NULL, population_model = NULL) {
  draws <- if (!is.null(smc$Theta)) {
    weighted_resample_rows(smc$Theta, smc$w, n_draws = n_draws, seed = seed)
  } else {
    weighted_resample_rows(smc$theta, smc$w, n_draws = n_draws, seed = seed)
  }

  if (is.null(population_model) || is.null(population_model$theta_to_list)) {
    return(as.data.frame(draws, check.names = FALSE))
  }

  parts <- population_model$theta_to_list(draws)
  lapply(parts, function(x) as.data.frame(as.matrix(x), check.names = FALSE))
}

posterior_intervals <- function(x) {
  x <- as.matrix(x)
  q <- t(apply(x, 2L, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
  data.frame(
    parameter = colnames(x),
    q025 = q[, 1L],
    q500 = q[, 2L],
    q975 = q[, 3L],
    row.names = NULL,
    check.names = FALSE
  )
}

plot_posteriors <- function(x,
                            y = NULL,
                            labels = c("Posterior 1", "Posterior 2"),
                            cols = c("steelblue", "firebrick3"),
                            n_cols = NULL,
                            main_prefix = NULL) {
  x <- as.matrix(x)
  if (!is.null(y)) {
    y <- as.matrix(y)
    if (!is.null(colnames(x)) && !is.null(colnames(y))) {
      y <- y[, colnames(x), drop = FALSE]
    }
  }

  n_params <- ncol(x)
  if (is.null(n_cols)) {
    n_cols <- min(3L, n_params)
  }
  n_rows <- ceiling(n_params / n_cols)

  graphics::par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

  for (j in seq_len(n_params)) {
    dx <- stats::density(x[, j])
    title <- colnames(x)[j]
    if (!is.null(main_prefix)) {
      title <- paste0(main_prefix, ": ", title)
    }

    if (is.null(y)) {
      plot(
        dx,
        lwd = 2,
        col = cols[1L],
        main = title,
        xlab = colnames(x)[j],
        ylab = "Density"
      )
    } else {
      dy <- stats::density(y[, j])
      xlim <- range(c(dx$x, dy$x))
      ylim <- c(0, 1.05 * max(dx$y, dy$y))
      plot(
        dx,
        lwd = 2,
        col = cols[1L],
        xlim = xlim,
        ylim = ylim,
        main = title,
        xlab = colnames(x)[j],
        ylab = "Density"
      )
      lines(dy, lwd = 2, col = cols[2L])
      legend(
        "topright",
        legend = labels,
        col = cols[seq_along(labels)],
        lwd = 2,
        bty = "n",
        cex = 0.85
      )
    }
  }
}
