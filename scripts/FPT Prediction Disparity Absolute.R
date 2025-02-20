###############################################################################

# Fair Policy Targeting (Prediction Disparity Absolute): Function

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
# W_bar: the pareto frontier
# capacity_constraint: maximal number of individuals to be treated
# start: possible starting value for the optimisation
# timelimit: maximum time spent on optimisation in seconds

optimisePredictionDisparityAbs = function(Y, X, D, S, m0, m1, ps, alpha = seq(from = 0.05, to = 0.95, length.out = N), 
                                          tolerance = 1e-6, W_bar, capacity_constraint, start = NA, timelimit = 10000){
  
  # compute the doubly robust score (drs) 
  G11_hat = (S/mean(S)) * ((D/ps)*(Y-m1)+m1) # treated females
  G01_hat = (S/mean(S)) * (((1-D)/(1-ps))*(Y-m0)+m0) # untreated females
  G10_hat = ((1-S)/mean(1-S)) * ((D/ps)*(Y-m1)+m1) # treated males
  G00_hat = ((1-S)/mean(1-S)) * (((1-D)/(1-ps))*(Y-m0)+m0) # untreated males
  
  GS1 = G11_hat - G01_hat # female drs
  GS0 = G10_hat - G00_hat # male drs
  
  # add intercept for beta0, and normalise
  X = as.matrix(cbind(1, X))
  max_val = max(apply(X, 1, function(x) sum(abs(x))))
  XX = X/max_val
  
  # set parameters
  n = nrow(XX)
  p = ncol(XX)
  
  ## initialise the model
  model = list()
  
  # sense of the optimisation, minimise the predictive disparity
  model$modelsense = "min"
  
  # the objective function is the abs prediction disparity, betas and u_j not in objective
  model$obj = c(rep(1, n), rep(0, p+N+2))
  
  # set the linear constraint matrix 
  model$A = rbind(cbind(diag(1, nrow = n), -XX, matrix(0, nrow = n, ncol = N+2)), # policies - betas <= 1
                  cbind(diag(1, nrow = n), -XX, matrix(0, nrow = n, ncol = N+2)), # policies - betas > 0
                  c(rep(1, n), rep(0, p+N+2)), # capacity constraint
                  c(rep(0, n+p), rep(1, N), 0, 0), # constraint on the u_j (C) >= 1
                  c(((1-S)/mean(1-S) - S/mean(S)), rep(0, p+N), -1, 0), # policies - slack_var <= 0
                  c(-((1-S)/mean(1-S) - S/mean(S)), rep(0, p+N), 0, -1)) # -slack_var -policies <= 0
  
  # set the rhs (with tolerance 1e-6)
  model$rhs = c(rep(1-tolerance, n), rep(tolerance, n), capacity_constraint, 1, 0, 0)
  
  # set the constraint directions
  model$sense = c(rep("<=", n), rep(">", n), "<=", ">=", "<=", "<=")
  
  # define the type of the variables
  model$vtype = c(rep("B", n), rep("C", p), rep("B", N), "C", "C")
  
  # put bounds on the parameter space, for z_i, constants, and u_j  [0,1], for the betas [-1,1], Inf for the slack variables
  model$ub =  c(rep(1, n+p+N), Inf, Inf)
  model$lb =  c(rep(0, n), rep(-1, p), rep(0, N+2))
  
  # set up the quadratic constraints (B)
  listnames = c(1:N)
  model$quadcon = sapply(listnames, function(x) NULL)
  
  # initiate a for loop to add the N=sqrt(n) constraints
  for (i in 1:N){
    
    # initiate the Qc matrix
    model$quadcon[[i]]$Qc = matrix(0, nrow = n+p+N+2, ncol = n+p+N+2) # rows and cols must = ncol of A
    # set coefficients in the column of the respective u_j
    model$quadcon[[i]]$Qc[which(S==0), n+p+i] = alpha[i]*GS0[GS0!=0] # male welfare, omit female 0s in GS0
    model$quadcon[[i]]$Qc[which(S==1), n+p+i] = (1-alpha[i])*GS1[GS1!=0] # female welfare, omit male 0s in GS1
    
    # initiate q vector for the linear terms (optimal welfare W_bar)
    model$quadcon[[i]]$q = rep(0, n+p+N+2)
    # insert the respective welfare (- as bring it to lhs)
    model$quadcon[[i]]$q[n+p+i] = -W_bar[i]
    
    # set the rhs
    model$quadcon[[i]]$rhs = -tolerance
    
    # direction of constraint
    model$quadcon[[i]]$sense = ">="
    
  }
  
  if (is.na(start[1]) == F) model$start = start
  
  # set a list of parameters
  params = list(IntFeasTol = 1e-9, FeasibilityTol = 1e-9, TimeLimit = timelimit, # tolerance limits
                BarConvTol = exp(-2), # tolerance on the barrier between primal and dual solution
                Disconnected = 0, # degree of exploitation of (independent) submodels
                Heuristics = 0, # fraction of runtime spent on feasibility heuristics
                NodefileStart = 0.5) # max node memory before writing to hard drive
  
  # solve the model
  result = gurobi(model, params = params)
  
  # extract the values 
  beta = result$x[(n+1):(n+p)] # the betas for deciding on the policy
  u = result$x[((n+p+1)):length(result$x)] # which u_j is equal to 1
  weight = 1-alpha[u == 1] # store the respective optimal weight alpha_j
  policies = apply(X, 1, function(x) ifelse(x%*%beta > 0, 1, 0)) # the estimated policies
  objval = result$objval # the minimised unfairness (prediction disparity)
  
  # return estimated betas, weight of the female group in the optimisation, policies and objective value
  return(list(beta = beta, weight_female = weight, policies = policies, objval = objval, 
              result.x = result$x, time = result$runtime, gap = result$mipgap))
  
}

###############################################################################

# Estimation FPT Prediction Disparity Absolute

###############################################################################

library(slam)
library(Matrix)
library(gurobi)
library(dplyr)

#### Preliminaries 
load("./results/implementation.RData")
load("./results/pareto_frontier.RData")

Y = df_sample$Y
D = df_sample$D
S = df_sample$S

m1 = df_sample$m1
m0 = df_sample$m0
ps = df_sample$ps

X = df_sample %>% select(C3, X1_D, X2_D, X3_D, X4_D, X5_D)
capacity_constraint = floor(0.33*nrow(X))

N = floor(sqrt(nrow(X)))
alpha = seq(from = 0.05, to = 0.95, length.out = N)



#### FPT Prediction Disparity Absolute Estimation

## Build start vector with the help of the pareto frontier values to help find an initial solution in the optimisation
X_start = as.matrix(cbind(1, S, X))

# "list comprehension" to extract the policies from the beta coefficients and the data
policies = t(apply(pareto_beta, 1, function(x) sapply(X_start%*%x, function(y) ifelse(y > 0, 1, 0))))

# compute the objective value, i.e., the unfairness generated by the respective policy
objective = apply(policies, 1, function(x) abs( sum(x*(1-S)/mean(1-S)) - sum(x*S/mean(S)) ))

# find the minimal objective; as there are multiple results, choose the one closest to equal weighting for female and male
minimal_obj = which(objective == min(objective))
minimal_obj = minimal_obj[which.min(abs(minimal_obj - N/2))]

# start value
start_pred_disp_abs = c(policies[minimal_obj,], pareto_beta[minimal_obj,], ifelse(c(1:N) == minimal_obj, 1, 0), NA, NA)

# optimise! after 10'000s gap of 83.6%
FPT_pred_disparity_abs = optimisePredictionDisparityAbs(Y, X = cbind(S, X), D, S, m0, m1, ps, alpha = alpha, tolerance = 1e-6,
                                                        W_bar = pareto_W_bar, capacity_constraint, start = start_pred_disp_abs,
                                                        timelimit = 15000)



# save(list = c("FPT_pred_disparity_abs"), file = "./results/fpt_absolute_disp.RData")
