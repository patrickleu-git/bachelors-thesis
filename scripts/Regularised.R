###############################################################################

# Regularised Welfare Maximisation: Function

###############################################################################

# Y: outcome 
# X: matrix of covariates
# D: treatment vector
# S: sensitive attribute vector
# ps: estimated propensity score
# m1: conditional mean of treated individuals
# m0: conditional mean of untreated individuals
# alpha: weight of the male welfare
# tolerance: slackness parameter (discreteness)
# capacity_constraint: maximal number of individuals to be treated
# timelimit: maximum time spent on optimisation in seconds
# lambda: regularisation parameter

optimiseRegularised = function(Y, X, D, S, ps, m1, m0, alpha = mean(1-S), tolerance = 1e-6, 
                               capacity_constraint, timelimit = 5000, lambda = 0.004){
  
  # load libraries
  library(Matrix)
  library(slam)
  library(gurobi)
  
  # compute the doubly robust score (drs)
  G11_hat = (S/mean(S)) * ((D/ps)*(Y-m1)+m1) # treated females
  G01_hat = (S/mean(S)) * (((1-D)/(1-ps))*(Y-m0)+m0) # untreated females
  G10_hat = ((1-S)/mean(1-S)) * ((D/ps)*(Y-m1)+m1) # treated males
  G00_hat = ((1-S)/mean(1-S)) * (((1-D)/(1-ps))*(Y-m0)+m0) # untreated males
  
  GS1 = G11_hat - G01_hat # female
  GS0 = G10_hat - G00_hat # male
  G = alpha*GS0 + (1-alpha)*GS1 # scaled vector
  
  # add intercept for beta0, and normalise
  X = as.matrix(cbind(1, X))
  max_val = max(apply(X, 1, function(x) max(abs(x))))
  XX = X/max_val
  
  # set parameters
  n = nrow(X)
  p = ncol(X)
  
  ## initialise the model
  model = list()
  
  # sense of the optimisation, maximise utility
  model$modelsense = "max"
  
  # objective function contains the welfare & regularisation components
  model$obj = c(G, rep(0, p), rep(-lambda, 2))
  
  # set the linear constraint matrix 
  model$A = rbind(cbind(diag(1, nrow = n), -X, 0, 0), # policies - betas <= 1
                  cbind(diag(1, nrow = n), -X, 0, 0), # policies - betas > 0
                  c(rep(1, n), rep(0, p+2)), # treated <= capacity constraint
                  c(-((1-S)/mean(1-S) - S/mean(S)), rep(0, p), -1, 0), # -w - policies <= 0 (male - female)
                  c(((1-S)/mean(1-S) - S/mean(S)), rep(0, p), 0, -1)) # policies - w <= 0 (male - female)
  
  # the rhs of the constraints, with tolerance (1e-6)
  model$rhs = c(rep(1 - tolerance, n), rep(tolerance, n), capacity_constraint, rep(0, 2))
  
  # set the constraint directions
  model$sense = c(rep("<=", n), rep(">", n), rep("<=", 3))
  
  # define the type of the variables
  model$vtype = c(rep("B", n), rep("C", p+2))
  
  # Put bounds on the parameter space, for z [0,1] and for the betas [-1,1], Inf for the slack variables
  model$ub = c(rep(1, n+p), rep(Inf, 2))
  model$lb = c(rep(0, n), rep(-1, p), rep(0, 2))
  
  # set additional parameters for the optimisation
  params = list(IntFeasTol = 1e-9, FeasibilityTol = 1e-9, TimeLimit = timelimit, # tolerance limits
                BarConvTol = exp(-2), # tolerance on the barrier between primal and dual solution
                Disconnected = 0, # degree of exploitation of (independent) submodels
                Heuristics = 0, # fraction of runtime spent on feasibility heuristics
                NodefileStart = 0.5) # max node memory before writing to hard drive
  
  # solve the model
  result = gurobi(model, params = params)
  
  # extract the values
  beta = result$x[(n+1):(n+p)]
  policies = apply(X, 1, function(x) ifelse(x%*%beta > 0, 1, 0))
  objval = result$objval
  
  return(list(beta = beta, policies = policies, objval = objval, lambda = lambda, 
              result.x = result$x, time = result$runtime, gap = result$mipgap))
  
}



###############################################################################

# Estimation FPT Regularised

###############################################################################

library(slam)
library(Matrix)
library(gurobi)
library(dplyr)

#### Preliminaries 
load("./results/implementation.RData")

Y = df_sample$Y
D = df_sample$D
S = df_sample$S

m1 = df_sample$m1
m0 = df_sample$m0
ps = df_sample$ps

X = df_sample %>% select(C3, X1_D, X2_D, X3_D, X4_D, X5_D)
capacity_constraint = floor(0.33*nrow(X))



#### Regularised Estimation

# regularisation with lambda = 0.004
FPT_regularised_004 = optimiseRegularised(Y, X = cbind(S,X), D, S, ps, m1, m0, alpha = mean(1-S), tolerance = 1e-6,
                                          capacity_constraint, timelimit = 15000, lambda = 0.004)

# regularisation with lambda = 0.010
FPT_regularised_010 = optimiseRegularised(Y, X = cbind(S,X), D, S, ps, m1, m0, alpha = mean(1-S), tolerance = 1e-6,
                                          capacity_constraint, timelimit = 15000, lambda = 0.010)

# regularisation with lambda = 0.100
FPT_regularised_100 = optimiseRegularised(Y, X = cbind(S,X), D, S, ps, m1, m0, alpha = mean(1-S), tolerance = 1e-6,
                                          capacity_constraint, timelimit = 15000, lambda = 0.100)


# save(list = c("FPT_regularised_004", "FPT_regularised_010", "FPT_regularised_100"), file = "./results/fpt_regularised.RData")

# share of treated male / female students with the different approaches
gf = as_tibble(cbind(FPT_regularised_004$policies, FPT_regularised_010$policies, FPT_regularised_100$policies, S))
treatment_dist = gf %>% group_by(S) %>% summarise(D004 = mean(V1), D010 = mean(V2), D100 = mean(V3))
