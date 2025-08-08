# Function to obtain a crude sequence ratio from exposure and event times
# for each individual
sr_wrapper <- function(id_expo, t_expo, id_event, t_event,
                       window_size = 365, blackout = 0, 
                       obs_period = NULL, 
                       weights = "studypop") {

  # Error handling for observation periods
  # If no observation period is specified, it sets the observation period
  # to the minimum and maximum of the observed event times
  if (is.null(obs_period))
    obs_period <- c(min(t_event), max(t_event))

  expo <- data.table(
    id = id_expo, 
    entry = t_expo, 
    lower = t_expo - window_size, 
    upper = t_expo + window_size
  )
  events <- data.table(id = id_event, t_event = t_event)
  
  d <- events[expo, on = .(id, t_event >= lower, t_event <= upper), nomatch = 0L,
              .(id, t = x.t_event - i.entry, entry = i.entry)]
  d <- subset(d, t != 0)
  d[, first := t > 0]
  
  # Apply blackout period
  d <- subset(d, abs(t) > blackout)
  d <- subset(d, between(entry, obs_period[1] + window_size, 
                                obs_period[2] - window_size))
  
  # Return SR with 95% CI
  first <- sum(d$first)
  last <- nrow(d) - first
  sr <- first / last
  lb_p <- qbeta(0.025, first + 0.5, last + 0.5)
  ub_p <- qbeta(0.975, first + 0.5, last + 0.5)
  
  # Calculate null-effect sequence ratio
  if (weights == "studypop") {
    weight_dates <- d$entry
  } else if (weights == "marginal") {
    weight_dates <- t_expo[t_expo >= (obs_period[1] + window_size) &
                           t_expo <= (obs_period[2] - window_size) &
                           !is.na(t_expo)]
  } else {
    stop("Argument weights must be 'studypop' or 'marginal'")
  }

  null_sr <- calc_null_effect(t_event, weight_dates, obs_period = obs_period, window_size = window_size)
  crude_sr <- c(sr = sr, lb = lb_p / (1 - lb_p), ub = ub_p / (1 - ub_p))
  
  list(
    crude = crude_sr,
    adj   = crude_sr / null_sr,
    null  = null_sr,
    hist  = sr_histogram(d$t, binsize = ifelse(window_size > 90, 30, 7)),
    data  = d
  )
}


# Function to produce a barchart representing the distribution of sequences
# over time
sr_histogram <- function(t, binsize = 30) {
  d <- data.table(t)
  d[, t_cat := fifelse(t > 0, ceiling(t/binsize), floor(t/binsize))]
  
  ggplot(d, aes(t_cat)) +
    geom_bar() +
    geom_vline(xintercept = 0) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_minimal()
}


# Function to calculate null-effect sequence ratios
calc_null_effect <- function(event_times, weight_dates, obs_period, window_size = 365) {
  counts <- data.table(t = event_times)[, .(n = .N), by = t]
  
  study_period <- data.table(t = seq(obs_period[1], obs_period[2], by = 1))
  counts <- counts[study_period, on = "t"]
  counts[is.na(n), n := 0]
  
  counts$before <- frollsum(counts$n, n = window_size + 1, align = "right") - counts$n
  counts$after  <- frollsum(counts$n, n = window_size + 1, align = "left")  - counts$n
  
  weights <- data.table(t = weight_dates)[, .(w = .N), by = t]
  names(weights) <- c("t", "w")
  
  counts[weights, on = "t", nomatch = 0L][, sum(w * after) / sum(w * before)]
}
