### Example adapted from:
# https://gist.github.com/oliviergimenez/68ad17910a62635ff6a062f8ec34292f

# Function 1
gen_data <- function(n = 100, seed = NULL) {
  # n is the sample size
  if(!is.null(seed)) {
    set.seed(seed)} 
  else {
    warning("No seed provided!", call. = FALSE)
    seed <- sample.int(10000, 1)
    set.seed(seed)
    message("Random seed = ", seed, "\n")
  }
  
  x <- sort(rnorm(n)) # covariate values
  int <- 30 # true intercept
  slope <- 10 # true slope
  mu <- int + slope * x # mean
  sigma <- 5 # standard deviation
  y <- rnorm(n, mean = mu, sd = sigma) # response
  data <- data.frame(y = y, x = x)
  
  return(list(data=data,parms=c(int,slope,sigma)))
}

# Function 2
est_post <- function(data, nchains, seed) {
  
  set.seed(seed)
  
  # build models
  m1.jags <- function(){ # model with covariate
    # Likelihood
    for (i in 1:n){
      y[i] ~ dnorm(mu[i], tau)
      mu[i] <- alpha[1] + alpha[2] * x[i]
    }
    # Priors:
    alpha[1] ~ dnorm(0, 0.1) # intercept
    alpha[2] ~ dnorm(0, 0.1) # slope
    sigma ~ dunif(0, 100) # standard deviation
    tau <- 1 / (sigma * sigma) # precision
  }
  
  # initial values for model with covariate
  inits.m1 <- function(){ 
    list(alpha = runif(2), 
         sigma = runif(1))
  }
  
  # parameters to be monitored
  params.m1 <- c("alpha", "sigma") # for model with covariate
  # data
  jags.data.m1 <- list(y = data$y, # for model with covariate
                       x = data$x,
                       n = length(data$y))
  
  # run jags and fit model with covariate
  m1.jags.fit <- jags(data = jags.data.m1, 
                  inits = inits.m1, 
                  parameters.to.save = params.m1, 
                  model.file = m1.jags,
                  n.chains = nchains, 
                  n.iter = 5000, 
                  n.burnin = 1000, 
                  n.thin = 2,
                  progress.bar = "none",
                  quiet = TRUE)
  
  # process output
  samples.m1.coda <- as.mcmc.list(as.mcmc(m1.jags.fit))
  postmeans <- summary(samples.m1.coda[,!grepl("deviance", varnames(samples.m1.coda))])$statistics[,1]
  
  # put results together
  return(list(postmeans=postmeans, jags.fit=m1.jags.fit))
}

# Function 2a
combine_chains <- function(res, n.chains0){
  
  ### combine chains following jags.parallel
  result <- NULL
  model <- NULL
  for (ch in 1:n.chains0) {
    result <- abind(result, res[[ch]]$jagsfit$BUGSoutput$sims.array, 
                    along = 2)
    model[[ch]] <- res[[ch]]$jagsfit$model
  }
  
  result <- R2jags:::as.bugs.array2(result, model.file=res[[1]]$jagsfit$model.file, program = "jags", DIC = TRUE, n.iter = res[[1]]$jagsfit$n.iter, n.burnin = res[[1]]$jagsfit$BUGSoutput$n.burnin, n.thin = res[[1]]$jagsfit$BUGSoutput$n.thin)
  
  out <- list(model = model, BUGSoutput = result, parameters.to.save = res[[1]]$jagsfit$parameters.to.save, model.file = res[[1]]$jagsfit$model.file, n.iter = res[[1]]$jagsfit$n.iter, DIC = TRUE)
  class(out) <- c("rjags.parallel", "rjags")
  
  out.coda <- as.mcmc.list(as.mcmc(out))
  postmeans <- summary(out.coda[,!grepl("deviance", varnames(out.coda))])$statistics[,1]
  
  return(list(postmeans=postmeans))
}

# Function 3
summ_pars <- function(i, est_params, true_params) {
  result <- data.frame(indx = i,
    parameter = c("int", "slope", "sigma"),
    bias = est_params - true_params)
  return(result)
}