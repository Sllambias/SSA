# Basic test of the null effect SR estimation
# Ensures that the trend generated in generate-trend-data.do is removed
# in the estimation of the adjusted SR

library(data.table)
library(ggplot2)
library(extraDistr)

source("ssa-functions.R")

gen_trend_data <- function(nindiv, years) {
  expo_rate <- 1/(50 * 365)
  outc_rate <- 1/(50 * 365)
  max_fu <- 365 * years
  
  data.table(
    id = 1:nindiv,
    outc = ceiling(rgompertz(nindiv, outc_rate, 0.001)),
    expo = ceiling(rexp(nindiv, expo_rate))
  )
}

analyze_trend_data <- function(x, weights) {
  set.seed(x)
  d <- gen_trend_data(1e5, years = 5)
  sr_wrapper(
    d$id, d$expo, d$id, d$outc,
    obs_period = c(0, 365 * 5),
    weights = weights
  )$adj[1]
}

w1 <- sapply(1:250, analyze_trend_data, weights = "studypop")
stopifnot(exp(mean(log(w1))) > 0.99 &
            exp(mean(log(w1))) < 1.01)

w2 <- sapply(1:250, analyze_trend_data, weights = "marginal")
stopifnot(exp(mean(log(w2))) > 0.99 &
            exp(mean(log(w2))) < 1.01)