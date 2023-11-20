### Example from "Conducting Monte Carlo Simulations in R" by Okan Bulut
# https://doi.org/10.7939/r3-y066-fs55
# https://okanbulut.github.io/simulations_in_r/run.html#parallel-computing

# Function 1
generate_data <- function(nitem, nexaminee, seed = NULL) {
  
  if(!is.null(seed)) {
    set.seed(seed)} 
  else {
    warning("No seed provided!", call. = FALSE)
    seed <- sample.int(10000, 1)
    set.seed(seed)
    message("Random seed = ", seed, "\n")
  }
  
  itempar <- cbind(
    rnorm(nitem, mean = 1.13, sd = 0.25), #a
    rnorm(nitem, mean = 0.21, sd = 0.51), #b
    rnorm(nitem, mean = 0.16, sd = 0.05)) #c
  
  ability <- rnorm(nexaminee, mean = 0, sd = 1)
  
  respdata <- irtoys::sim(ip = itempar, x = ability)
  colnames(respdata) <- paste0("item", 1:nitem)
  
  data <- list(itempar = itempar,
               ability = ability,
               seed = seed,
               respdata = respdata)
  
  return(data)
}

# Function 2
estimate_par <- function(data, guess = -1) {
  # If guessing is fixed
  if(guess >= 0) {
    
    # Model set up
    mod3PL <- mirt::mirt(data, # response data
                         1,    # unidimensional model
                         guess = guess, # fixed guessing
                         verbose = FALSE, # Don't print verbose
                         # Increase the number of EM cycles
                         # Turn off estimation messages
                         technical = list(NCYCLES = 1000,
                                          message = FALSE)) 
  } else {
    mod3PL <- mirt::mirt(data, # response data
                         1,    # unidimensional model
                         itemtype = "3PL", # IRT model
                         verbose = FALSE, # Don't print verbose
                         # Increase the number of EM cycles
                         # Turn off estimation messages
                         technical = list(NCYCLES = 1000,
                                          message = FALSE)) 
  }
  
  # Extract item parameters in typical IRT metric
  itempar_est <- as.data.frame(mirt::coef(mod3PL, IRTpars = TRUE, simplify = TRUE)$item[,1:3])
  return(itempar_est)
}

# Function 3
summarize <- function(est_params, true_params) {
  result <- data.frame(
    parameter = c("a", "b", "c"),
    bias = sapply(1L:3L, function(i) mean((est_params[, i] - true_params[,i]))),
    rmse = sapply(1L:3L, function(i) sqrt(mean((est_params[, i] - true_params[,i])^2))),
    correlation = sapply(1L:3L, function(i) cor(est_params[, i], true_params[,i])))
  return(result)
}