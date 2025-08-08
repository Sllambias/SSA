# Load functions to perform SSA
library(data.table)
library(ggplot2)
source("ssa-functions.R")


# Set number of individuals to simulate
nindiv <- 1e6

# Simulate data
d <- data.table(
  id = 1:nindiv, 
  t_expo = ceiling(rexp(nindiv, rate = 1/(50  * 365))),
  t_outc = ceiling(rexp(nindiv, rate = 1/(200 * 365)))
)

# Perform sequence symmetry analysis
ssa <- sr_wrapper(
  id_expo  = d$id, 
  t_expo   = d$t_expo, 
  id_event = d$id, 
  t_event  = d$t_outc, 
  obs_period = c(0, 365 * 5)
)

# Show everything that is returned from the function
ssa

# Show histogram
ssa$hist
