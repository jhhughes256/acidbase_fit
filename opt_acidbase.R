# Set up R environment
# Load library
  library(ggplot2)  # graphical package
  
# Set ggplot2 theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))
  theme_update(plot.title = element_text(hjust = 0.5))
  
# Read in data
  aspirin_data <- read.csv("D:/2019 PostDoc/Des/aspirin_hydrolysis.csv")
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Determine rate constants
# Set constants for reaction
  Kw <- 1E-14
  pKa <- 3.57
  Ka <- 10^-pKa
  err <- 0.1
  
# Set initial parameters (in days^-1)
# Parameters logged to allow easier search for minima
  k_constants <- c(2.5, 0.01, 0.1, 5000)
  init_par <- c(k_constants/24/60, err)
  
  k_pred <- function(par, pH) {
  # Define dependent and independent variables
    x <- pH
    
  # Determine values dependent on independent variable
    H <- 10^-x
    fAH <- H/(H+Ka)
    fA <- Ka/(H+Ka)
    
  # Define parameter values to be used in the equations
    k1H <- par[1]
    k1H20 <- par[2]
    k2H20 <- par[3]
    k2OH <- par[4]

  # Predict dependent value
    k1H*H*fAH + k1H20*fAH + k2H20*fA + k2OH*(Kw/H)*fA
  }
  
  acidbase_fn <- function(par) {
  # Define dependent value
    x <- aspirin_data$pH
    y <- 10^aspirin_data$logk
    
  # Predict dependent value
    yhat <- k_pred(par, x)
    
  # Test likelihood that prediction is correct
    loglik <- dnorm(y, mean = yhat, sd = par[5]*yhat, log = T)
    
  # Return objective function value to be minimised
    return(-1*sum(loglik))
  }
  
# Optimise parameters
  opt_acidbase <- optim(init_par, acidbase_fn, method = "L-BFGS-B",
    lower = c(init_par[1:4]*0.8, 0.001), upper = c(init_par[1:4]*1.2, 1)
    # control = list(
    #   parscale = init_par, fnscale = acidbase_fn(init_par)
    # )
  )
  
# Extract optimised parameters and predict rate constant
  opt_par <- opt_acidbase$par[1:4]*60*24
  pred_data <- data.frame(
    pH = seq(0, 14, by = 0.1)
  )
  
# Predict each reactions contribution to the rate constant
  pred_data$H <- 10^-pred_data$pH
  pred_data$logpred <- log10(k_pred(opt_par, pred_data$pH))
  pred_data$k1H <- log10(with(pred_data, opt_par[1]*H*(H/(H+Ka))))
  pred_data$k1H2O <- log10(with(pred_data, opt_par[2]*(H/(H+Ka))))
  pred_data$k2H2O <- log10(with(pred_data, opt_par[3]*(Ka/(H+Ka))))
  pred_data$k2OH <- log10(with(pred_data, opt_par[4]*(Kw/H)*(Ka/(H+Ka))))
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create plots
# Plot 
  p <- NULL
  p <- ggplot(data = pred_data)
  p <- p + geom_point(aes(x = pH, y = log10(10^(logk)*60*24)), data = aspirin_data, shape = 1)
  p <- p + geom_line(aes(x = pH, y = logpred), size = 1)
  p <- p + geom_line(aes(x = pH, y = k1H), linetype = "dashed")
  p <- p + geom_line(aes(x = pH, y = k1H2O), linetype = "dashed")
  p <- p + geom_line(aes(x = pH, y = k2H2O), linetype = "dashed")
  p <- p + geom_line(aes(x = pH, y = k2OH), linetype = "dashed")
  p <- p + ylab("log k (d^-1)")
  p <- p + coord_cartesian(xlim = c(0.5, 11.5), ylim = c(-4, 2.5))
  p
  