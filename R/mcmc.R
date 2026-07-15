# =============================================================================
# Generic random-walk Metropolis-Hastings MCMC sampler and diagnostics
# =============================================================================

# Random-walk Metropolis-Hastings MCMC sampler.
# A simple general-purpose sampler. Each iteration proposes a new parameter
# vector by adding independent Gaussian noise to the current state, then
# accepts or rejects via the Metropolis criterion.
# log_post_fn: function taking a named numeric vector and returning a scalar
# log-posterior.
# init: named numeric vector of starting parameter values.
# n_iter: total iterations (including burn-in).
# proposal_sd: named numeric vector of proposal standard deviations, one per
# parameter (tune for 20-40% acceptance rate).
# Returns a list with samples (matrix n_iter x n_params), acceptance_rate,
# and log_post (log-posterior trace).
run_mcmc <- function(log_post_fn, init, n_iter, proposal_sd, ...) {
  n_params <- length(init)
  pnames   <- names(init)
  samples  <- matrix(NA_real_, nrow = n_iter, ncol = n_params,
                     dimnames = list(NULL, pnames))
  log_post_trace <- numeric(n_iter)

  current  <- init
  lp_curr  <- log_post_fn(current, ...)
  accepted <- 0L

  for (i in seq_len(n_iter)) {
    proposal        <- current + stats::rnorm(n_params, 0, proposal_sd)
    names(proposal) <- pnames

    lp_prop <- log_post_fn(proposal, ...)
    log_r   <- lp_prop - lp_curr

    if (is.finite(log_r) && log(stats::runif(1)) < log_r) {
      current  <- proposal
      lp_curr  <- lp_prop
      accepted <- accepted + 1L
    }

    samples[i, ]      <- current
    log_post_trace[i] <- lp_curr
  }

  list(samples         = samples,
       acceptance_rate = accepted / n_iter,
       log_post        = log_post_trace)
}


# Trace plots for MCMC output.
# burnin: iterations to shade (not removed).
plot_mcmc_chains <- function(mcmc_out, burnin = 1000) {
  df           <- as.data.frame(mcmc_out$samples)
  df$iteration <- seq_len(nrow(df))
  param_long   <- tidyr::pivot_longer(df, -"iteration",
                                      names_to  = "parameter",
                                      values_to = "value")

  lp_long <- data.frame(
    iteration = seq_along(mcmc_out$log_post),
    parameter = "log(objective)",
    value     = mcmc_out$log_post
  )

  long <- rbind(lp_long, param_long)
  long$parameter <- factor(long$parameter,
                           levels = c("log(objective)",
                                      setdiff(unique(long$parameter),
                                              "log(objective)")))

  ggplot2::ggplot(long, ggplot2::aes(x = .data$iteration, y = .data$value)) +
    ggplot2::annotate("rect", xmin = 0, xmax = burnin,
                      ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.5) +
    ggplot2::geom_line(linewidth = 0.3, colour = "steelblue") +
    ggplot2::facet_wrap(~parameter, scales = "free_y", ncol = 1) +
    ggplot2::labs(x = "Iteration", y = "Value",
                  caption = paste0("Acceptance rate: ",
                                   round(mcmc_out$acceptance_rate * 100, 1),
                                   "%")) +
    ggplot2::theme_bw()
}
