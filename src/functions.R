#' -------------------------------------------------------------------------
#' Project:     SRL Growth Modelling 
#'
#' Created by:  Stephen Bradshaw
#' Modified:    23/01/2025
#' Version:     1.0 - Initial Release
#'              
#' Purpose:     Functions for the SRL Growth Modelling project
#' -------------------------------------------------------------------------


###############################################################
########################## FUNCTIONS ##########################
###############################################################

#' test that stan works
#' @return output vector of polygon (assessment area) ID
func_testStanWorks <- function() {
  
  # Prior to the tutorial make sure that the script below runs without error on your R installation.
  # What you need is a working installation of Stan: http://mc-stan.org/ .
  # For installation instructions, see here: 
  # https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
  
  # After installation you should be able to run this script which should output
  # some summary statistics and some pretty plots, :)
  
  # Generating some fake data
  set.seed(123)
  y <- rbinom(30, size = 1, prob = 0.2016)
  
  # Fitting a simple binomial model using Stan
  library(rstan)
  
  model_string <- "
  data {
        int n;
        int y[n];
      }
      parameters {
        real<lower=0, upper=1> theta;
      }
      model {
        y ~ bernoulli(theta);
      }"

stan_samples <- stan(model_code = model_string, data = list(y = y, n = length(y)) )
# stan_samples
traceplot(stan_samples)
plot(stan_samples)

return(stan_samples)
  
}


#' Modification to the standard "in" function
#' Use: x %!in% vector etc
`%!in%` = Negate(`%in%`)



############################################################