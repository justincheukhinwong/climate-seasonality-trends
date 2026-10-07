# Run from the repository root: Rscript tests/validate_hac.R
source(file.path("R", "calendar_hac.R"))
if (!requireNamespace("sandwich", quietly = TRUE)) {
  stop("Install the validation dependency first: install.packages('sandwich')")
}

assert_close <- function(actual, expected, label, tolerance = 1e-9) {
  error <- max(abs(actual - expected))
  scale <- max(1, max(abs(expected)))
  if (!is.finite(error) || error > tolerance * scale) {
    stop(label, ": maximum absolute difference = ", signif(error, 5))
  }
}

expect_error <- function(expression, label) {
  failed <- inherits(tryCatch(force(expression), error = identity), "error")
  if (!failed) stop(label, ": expected an informative error.")
}

# A heteroskedastic AR process exercises both parts of the HAC estimator.
set.seed(439)
n <- 360L
day <- seq_len(n)
error <- as.numeric(arima.sim(list(ar = 0.65), n = n)) * (1 + day / n)
dense <- data.frame(
  day = day,
  trend = (day - mean(day)) / 365.25,
  annual = sin(2 * pi * day / 365.25),
  temperature = 12 + 0.03 * day + 4 * sin(2 * pi * day / 365.25) + error
)
fit <- lm(temperature ~ trend + annual, data = dense)
for (lag in c(0L, 7L, 30L, 90L)) {
  for (adjust in c(FALSE, TRUE)) {
    expected <- sandwich::NeweyWest(fit, lag = lag, prewhite = FALSE, adjust = adjust)
    actual <- calendar_hac(fit, dense$day, lag = lag, adjust = adjust)
    assert_close(actual, expected, paste("Dense reference", lag, adjust))
  }
}

# Independently assemble all pair weights from calendar distances on gappy data.
gappy <- dense[-c(5L, 11L, 12L, 45:60, 101:135, 250L), ]
gappy_fit <- lm(temperature ~ trend + annual, data = gappy)
X <- model.matrix(gappy_fit)
score <- X * as.numeric(gappy_fit$residuals)
bread <- solve(crossprod(X))
distance <- abs(outer(gappy$day, gappy$day, "-"))
for (lag in c(0L, 7L, 30L, 90L, 500L)) {
  weights <- 1 - distance / (lag + 1)
  weights[weights < 0] <- 0
  expected <- bread %*% crossprod(score, weights %*% score) %*% bread
  expected <- expected * nrow(X) / (nrow(X) - ncol(X))
  actual <- calendar_hac(gappy_fit, gappy$day, lag = lag)
  assert_close(actual, expected, paste("Calendar-pair reference", lag))
  stopifnot(min(eigen(actual, symmetric = TRUE, only.values = TRUE)$values) > -1e-10)
}

# Observation order and calendar origin must not change coefficient covariance.
permutation <- sample(seq_len(nrow(gappy)))
reordered <- gappy[permutation, ]
reordered_fit <- lm(temperature ~ trend + annual, data = reordered)
reference <- calendar_hac(gappy_fit, gappy$day, lag = 30L)
assert_close(calendar_hac(reordered_fit, reordered$day, lag = 30L), reference,
             "Row-order invariance")
assert_close(calendar_hac(gappy_fit, gappy$day + 20000L, lag = 30L), reference,
             "Calendar-origin invariance")
assert_close(calendar_hac(gappy_fit, as.Date("2000-01-01") + gappy$day, lag = 30L),
             reference, "Date input")
assert_close(calendar_hac(gappy_fit, gappy$day, lag = 0L),
             sandwich::vcovHC(gappy_fit, type = "HC1"), "Lag-zero HC1 identity")

# Missing response handling must align days with fitted, rather than input, rows.
with_missing <- dense
with_missing$temperature[c(20L, 80L)] <- NA_real_
missing_fit <- lm(temperature ~ trend + annual, data = with_missing, na.action = na.exclude)
complete <- with_missing[complete.cases(with_missing), ]
complete_fit <- lm(temperature ~ trend + annual, data = complete)
assert_close(calendar_hac(missing_fit, complete$day, lag = 30L),
             calendar_hac(complete_fit, complete$day, lag = 30L), "na.exclude row alignment")

expect_error(calendar_hac(fit, rep(1L, n), 7L), "Duplicate days")
expect_error(calendar_hac(fit, dense$day + 0.5, 7L), "Fractional days")
expect_error(calendar_hac(fit, dense$day[-1L], 7L), "Unaligned day vector")
expect_error(calendar_hac(fit, dense$day, -1L), "Negative lag")
expect_error(calendar_hac(fit, dense$day, 1.5), "Fractional lag")
weighted_fit <- lm(temperature ~ trend + annual, data = dense, weights = rep(2, n))
expect_error(calendar_hac(weighted_fit, dense$day, 7L), "Unsupported weighted fit")
singular_fit <- lm(temperature ~ trend + I(2 * trend), data = dense)
expect_error(calendar_hac(singular_fit, dense$day, 7L), "Rank-deficient fit")

cat("HAC validation passed: dense Newey-West agreement; gappy calendar-pair agreement;\n",
    "positive semidefiniteness, invariances, HC1 identity, missing-row alignment, and input checks.\n",
    sep = "")
