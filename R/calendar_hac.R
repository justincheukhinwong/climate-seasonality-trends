# Estimate coefficient covariance using Bartlett weights on calendar-day lags.
# `days` must align with the observations used in an unweighted, full-rank lm.
# Missing days contribute no score; no temperatures or residuals are imputed.
calendar_hac <- function(fit, days, lag, adjust = TRUE) {
  if (!inherits(fit, "lm") || inherits(fit, "glm") || inherits(fit, "mlm") ||
      !is.null(fit$weights)) {
    stop("fit must be an unweighted, single-response lm.")
  }
  if (length(lag) != 1L || !is.numeric(lag) || !is.finite(lag) ||
      lag < 0 || lag != floor(lag)) {
    stop("lag must be one nonnegative integer number of calendar days.")
  }
  if (length(adjust) != 1L || !is.logical(adjust) || is.na(adjust)) {
    stop("adjust must be TRUE or FALSE.")
  }

  X <- model.matrix(fit)
  # Stored residuals correspond to fitted rows even when na.exclude was used.
  residual <- fit$residuals
  n <- nrow(X)
  p <- ncol(X)
  if (fit$rank != p || n <= p || length(residual) != n ||
      any(!is.finite(X)) || any(!is.finite(residual))) {
    stop("fit must have full column rank, finite fitted rows, and residual df > 0.")
  }
  if (inherits(days, "Date")) days <- as.numeric(days)
  if (!is.numeric(days) || length(days) != n || any(!is.finite(days)) ||
      any(days != floor(days)) || anyDuplicated(days)) {
    stop("days must contain one unique, finite integer day per fitted observation.")
  }

  # Keep true daily spacing: a gap must not create a spurious lag-one pair.
  calendar_index <- days - min(days) + 1
  calendar_n <- max(calendar_index)
  score <- matrix(0, nrow = calendar_n, ncol = p)
  score[calendar_index, ] <- X * as.numeric(residual)

  # S = sum(s_t s_t') + sum_h [1 - h/(lag+1)] sum_t (s_t s_(t-h)' + transpose).
  meat <- crossprod(score)
  last_lag <- min(lag, calendar_n - 1L)
  if (last_lag > 0L) {
    for (h in seq_len(last_lag)) {
      lag_product <- crossprod(
        score[(h + 1L):calendar_n, , drop = FALSE],
        score[seq_len(calendar_n - h), , drop = FALSE]
      )
      weight <- 1 - h / (lag + 1)
      meat <- meat + weight * (lag_product + t(lag_product))
    }
  }

  # Use observed design rows and observed n for the finite-sample adjustment.
  bread <- chol2inv(chol(crossprod(X)))
  covariance <- bread %*% meat %*% bread
  if (adjust) covariance <- covariance * n / (n - p)
  covariance <- (covariance + t(covariance)) / 2
  dimnames(covariance) <- list(colnames(X), colnames(X))
  covariance
}
