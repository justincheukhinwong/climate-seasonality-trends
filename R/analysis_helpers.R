# Reusable preparation, model summaries, and graphics for the climate analysis.

read_climate_data <- function(path) {
  if (!file.exists(path)) {
    stop("Place uscrn_daily.RData in data/raw/ first. See data/README.md.", call. = FALSE)
  }
  # Load into a private environment so the saved workspace cannot replace analysis objects.
  input <- new.env(parent = baseenv())
  load(path, envir = input)
  if (!exists("dat", input, inherits = FALSE) || !is.data.frame(input$dat)) {
    stop("Expected a data frame named 'dat' inside the RData file.")
  }
  raw <- input$dat
  required <- c("site", "LST_DATE", "T_DAILY_AVG", "LONGITUDE", "LATITUDE")
  if (!all(required %in% names(raw))) stop("The input is missing required columns.")
  d <- data.frame(
    site = as.character(raw$site),
    date = as.Date(as.character(raw$LST_DATE), "%Y%m%d"),
    temperature = as.numeric(raw$T_DAILY_AVG),
    longitude = as.numeric(raw$LONGITUDE), latitude = as.numeric(raw$LATITUDE)
  )
  if (anyNA(d$date) || anyNA(d$site) || any(!nzchar(d$site))) {
    stop("Every record must have a valid station name and calendar date.")
  }
  if (anyDuplicated(d[c("site", "date")])) stop("Duplicate station-date records found.")
  if (any(!is.finite(d$longitude)) || any(!is.finite(d$latitude))) {
    stop("Station coordinates must be finite.")
  }
  # NOAA uses -9999 for this response; other measurement columns are not analyzed here.
  sentinel <- d$temperature == -9999 & !is.na(d$temperature)
  d$temperature[sentinel | !is.finite(d$temperature)] <- NA_real_
  d$state <- substr(d$site, 1, 2)
  d$in_scope <- d$state %in% setdiff(state.abb, c("AK", "HI"))
  d <- d[order(d$site, d$date), ]
  rownames(d) <- NULL

  coverage <- do.call(rbind, lapply(split(d, d$site), function(s) {
    good <- s[!is.na(s$temperature), ]
    valid_span <- if (nrow(good)) as.integer(max(good$date) - min(good$date)) + 1L else NA_integer_
    data.frame(
      site = s$site[1], state = s$state[1], in_scope = s$in_scope[1],
      rows = nrow(s), valid_temperature = nrow(good), missing_temperature = sum(is.na(s$temperature)),
      first_record = min(s$date), last_record = max(s$date),
      first_valid = if (nrow(good)) min(good$date) else as.Date(NA),
      last_valid = if (nrow(good)) max(good$date) else as.Date(NA),
      valid_span_years = (valid_span - 1) / 365.25,
      valid_fraction_within_span = nrow(good) / valid_span,
      missing_calendar_days = as.integer(max(s$date) - min(s$date)) + 1L - nrow(s),
      longitude = median(s$longitude), latitude = median(s$latitude),
      coordinate_versions = nrow(unique(s[c("longitude", "latitude")])),
      longitude_range = diff(range(s$longitude)), latitude_range = diff(range(s$latitude)),
      status = if (!s$in_scope[1]) "outside_contiguous_US" else if (nrow(good) < 100) "fewer_than_100_valid_days" else "included"
    )
  }))
  rownames(coverage) <- NULL
  scope <- d[d$in_scope, ]
  analysis <- scope[!is.na(scope$temperature) & scope$site %in% coverage$site[coverage$status == "included"], ]
  analysis$days <- as.numeric(analysis$date - as.Date("2020-01-01"))
  analysis$years <- analysis$days / 365.25
  analysis <- add_harmonics(analysis)
  audit <- data.frame(
    metric = c("raw_rows", "raw_stations", "temperature_sentinels", "other_missing_response",
               "outside_scope_rows", "outside_scope_stations", "in_scope_rows", "in_scope_missing_response",
               "modeled_rows", "modeled_stations", "duplicate_station_dates"),
    value = c(nrow(raw), length(unique(d$site)), sum(sentinel),
              sum(is.na(d$temperature)) - sum(sentinel), sum(!d$in_scope),
              length(unique(d$site[!d$in_scope])), nrow(scope), sum(is.na(scope$temperature)),
              nrow(analysis), length(unique(analysis$site)), 0)
  )
  list(data = analysis, coverage = coverage, audit = audit)
}

add_harmonics <- function(d) {
  phase <- 2 * pi * d$days / 365.25
  d$sin1 <- sin(phase)
  d$cos1 <- cos(phase)
  d$sin2 <- sin(2 * phase)
  d$cos2 <- cos(2 * phase)
  d
}

fit_stations <- function(d) {
  # All six coefficients vary by station; the same model objects supply every summary.
  lapply(split(d, d$site), function(s) {
    fit <- lm(temperature ~ years + sin1 + cos1 + sin2 + cos2, data = s,
              na.action = na.fail)
    if (fit$rank != 6L) stop("Rank-deficient model at ", s$site[1])
    list(fit = fit, data = s)
  })
}

calendar_residual_cor <- function(errors, days, lag) {
  # Match actual dates rather than treating observations on either side of a gap as adjacent.
  previous <- match(days - lag, days)
  keep <- !is.na(previous)
  if (sum(keep) < 3L) return(NA_real_)
  cor(errors[keep], errors[previous[keep]])
}

summarize_models <- function(models) {
  do.call(rbind, lapply(names(models), function(site) {
    m <- models[[site]]
    fit <- m$fit
    b <- coef(fit)[["years"]]
    iid_se <- sqrt(vcov(fit)["years", "years"])
    crit <- qt(0.975, df.residual(fit))
    se <- vapply(c(7L, 30L, 90L), function(lag) {
      sqrt(calendar_hac(fit, m$data$days, lag)["years", "years"])
    }, numeric(1))
    data.frame(
      site = site, n = nobs(fit), residual_df = df.residual(fit),
      intercept_c = coef(fit)[["(Intercept)"]], trend_c_per_year = b,
      iid_se = iid_se, iid_lower = b - crit * iid_se, iid_upper = b + crit * iid_se,
      hac7_se = se[1], hac30_se = se[2], hac90_se = se[3],
      hac30_lower = b - qnorm(0.975) * se[2], hac30_upper = b + qnorm(0.975) * se[2],
      residual_sd_c = sigma(fit), r_squared = summary(fit)$r.squared,
      residual_cor_day1 = calendar_residual_cor(residuals(fit), m$data$days, 1)
    )
  }))
}

seasonal_features <- function(models) {
  # Hold the trend at its 2020 baseline and vary only seasonal phase over an idealized year.
  grid <- data.frame(days = (0:364) * 365.25 / 365, years = 0)
  grid <- add_harmonics(grid)
  curves <- vapply(models, function(m) as.numeric(predict(m$fit, newdata = grid)), numeric(365))
  features <- data.frame(
    site = names(models), seasonal_range_c = apply(curves, 2, function(x) diff(range(x))),
    coldest_day = apply(curves, 2, which.min), warmest_day = apply(curves, 2, which.max),
    row.names = NULL
  )
  # A July-to-June ordering keeps late December before early January within the same winter.
  features$winter_order <- (features$coldest_day - 183L) %% 365L
  list(features = features, curves = curves)
}

select_cycles <- function(features) {
  index <- c(which.max(features$seasonal_range_c), which.min(features$seasonal_range_c),
             which.min(features$winter_order), which.max(features$winter_order))
  cbind(selection = c("Largest seasonal range", "Smallest seasonal range",
                      "Earliest winter minimum", "Latest winter minimum"), features[index, ])
}

plot_style <- function() {
  par(mar = c(4.2, 4.5, 3.0, 1), mgp = c(2.7, 0.7, 0), las = 1,
      col.axis = "#384759", col.lab = "#24354B", col.main = "#152C45", cex.main = 1.05)
}

plot_exploration <- function(d) {
  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(mfrow = c(1, 2)); plot_style()
  h <- hist(d$temperature, breaks = 80, plot = FALSE)
  plot(h, col = "#387A9C", border = "white", yaxt = "n",
       main = "Observed daily temperatures", xlab = "Daily average temperature (C)", ylab = "Station-days")
  ticks <- pretty(h$counts)
  axis(2, at = ticks, labels = ifelse(ticks == 0, "0", paste0(ticks / 1000, "k")))
  # Average within station first so long records do not receive extra station weight in this display.
  month_day <- format(d$date, "%m-%d")
  keep <- month_day != "02-29"
  station_day <- aggregate(d$temperature[keep], list(site = d$site[keep], month_day = month_day[keep]), mean)
  cycle <- aggregate(x ~ month_day, station_day, mean)
  day <- as.integer(as.Date(paste0("2001-", cycle$month_day)) - as.Date("2001-01-01")) + 1L
  plot(day, cycle$x, type = "l", lwd = 2, col = "#B25E35", xlab = "Day of a non-leap year",
       ylab = "Station-average temperature (C)", main = "Observed seasonal cycle")
}

plot_station_boxes <- function(d) {
  old <- par(no.readonly = TRUE); on.exit(par(old))
  set.seed(123)
  sites <- sort(sample(sort(unique(d$site)), 15))
  x <- d[d$site %in% sites, ]
  x$site <- factor(x$site, levels = sites)
  par(mar = c(9, 4.4, 2.8, 1), las = 2)
  boxplot(temperature ~ site, data = x, outline = FALSE, col = "#A4C9D4",
          names = gsub("_", " ", sites), cex.axis = 0.66, xlab = "", ylab = "Daily average temperature (C)",
          main = "Temperature distributions at 15 reproducibly sampled stations")
}

plot_trends <- function(trends) {
  old <- par(no.readonly = TRUE); on.exit(par(old)); plot_style()
  x <- trends[order(trends$trend_c_per_year), ]
  index <- seq_len(nrow(x))
  limits <- range(x$iid_lower, x$iid_upper, x$hac30_lower, x$hac30_upper)
  plot(index, x$trend_c_per_year, type = "n", ylim = limits,
       xlab = "Station rank by estimated trend", ylab = "Temperature trend (C/year)",
       main = "Same OLS trends, two uncertainty assumptions")
  abline(h = 0, col = "#788696", lty = 2)
  segments(index, x$hac30_lower, index, x$hac30_upper, col = "#B9C9D3", lwd = 2)
  segments(index, x$iid_lower, index, x$iid_upper, col = "#B25E35", lwd = 1)
  points(index, x$trend_c_per_year, pch = 16, cex = 0.55, col = "#153E5A")
  legend("topleft", c("OLS estimate", "95% iid t interval", "95% HAC30 normal interval"),
         pch = c(16, NA, NA), lty = c(NA, 1, 1), lwd = c(NA, 1, 2),
         col = c("#153E5A", "#B25E35", "#B9C9D3"), bty = "n", cex = 0.85)
}

plot_cycles <- function(selected, curves) {
  old <- par(no.readonly = TRUE); on.exit(par(old))
  colors <- c("#285F89", "#B65B35", "#45876F", "#876294")
  par(mfrow = c(2, 2)); plot_style()
  limits <- range(curves[, selected$site, drop = FALSE])
  for (i in seq_len(nrow(selected))) {
    site <- selected$site[i]
    plot(1:365, curves[, site], type = "l", col = colors[i], lwd = 2, ylim = limits,
         xlab = "Day of idealized year", ylab = "Fitted seasonal mean (C)",
         main = paste(selected$selection[i], gsub("_", " ", site), sep = "\n"), cex.main = 0.85)
  }
}

plot_station_map <- function(features, variable, title, legend_title) {
  old <- par(no.readonly = TRUE); on.exit({layout(1); par(old)})
  layout(matrix(c(1, 2), 1), widths = c(5.5, 1))
  par(mar = c(1, 1, 3, 0.5))
  values <- features[[variable]]
  palette <- hcl.colors(80, "YlOrRd", rev = TRUE)
  limits <- range(values)
  index <- 1L + floor(79 * (values - limits[1]) / diff(limits))
  maps::map("state", col = "#C3CBD1", fill = TRUE, bg = "white", border = "white")
  points(features$longitude, features$latitude, pch = 21, cex = 1.0,
         bg = palette[index], col = "#344353", lwd = 0.35)
  title(main = title, col.main = "#152C45", cex.main = 1.05)
  par(mar = c(3, 0.5, 3, 4))
  image(x = 1, y = seq(limits[1], limits[2], length.out = 80),
        z = matrix(seq_len(80), nrow = 1), col = palette, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1, cex.axis = 0.75)
  mtext(legend_title, side = 4, line = 2.7, cex = 0.8)
}

plot_diagnostics <- function(model, site) {
  old <- par(no.readonly = TRUE); on.exit(par(old))
  par(mfrow = c(2, 2)); plot_style()
  fit <- model$fit; residual <- residuals(fit)
  plot(fitted(fit), residual, pch = 16, cex = 0.25, col = "#285F8945",
       xlab = "Fitted temperature (C)", ylab = "Residual (C)", main = "Residuals versus fitted values")
  abline(h = 0, lty = 2, col = "#B25E35")
  qqnorm(residual, pch = 16, cex = 0.25, col = "#387A9C", main = "Normal Q-Q plot")
  qqline(residual, col = "#B25E35", lwd = 2)
  lags <- 1:60
  correlations <- vapply(lags, function(k) calendar_residual_cor(residual, model$data$days, k), numeric(1))
  plot(lags, correlations, type = "h", lwd = 2, col = "#387A9C", ylim = range(c(0, correlations)),
       xlab = "Calendar-day lag", ylab = "Paired residual correlation", main = "Dependence after seasonal adjustment")
  abline(h = 0, col = "#788696")
  month <- factor(format(model$data$date, "%m"), levels = sprintf("%02d", 1:12), labels = month.abb)
  boxplot(residual ~ month, outline = FALSE, col = "#A4C9D4", xlab = "Month", ylab = "Residual (C)",
          main = "Residual spread across seasons", cex.axis = 0.75)
  abline(h = 0, col = "#788696", lty = 2)
}
