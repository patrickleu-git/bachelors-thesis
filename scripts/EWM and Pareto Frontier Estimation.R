###############################################################################

# Empirical Welfare Maximisation (EWM): Function

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
# funclass3: policy function class 3 that requires ATE(S==s_1) >= ATE(S==s_2)
# parity_sense: sense of the additional constraint above
# timelimit: maximum time spent on optimisation in seconds

EWM_estimation = function(Y, X, D, S, ps, m1, m0, alpha = mean(1-S),
                          tolerance = 1e-6, capacity_constraint, funclass3 = F, 
                          parity_sense = ">=", timelimit = 5000){
  
  # load libraries
  library(Matrix)
  library(slam)
  library(gurobi)
  
  # compute the doubly robust score (drs)
  # note: one could also add treatment cost here by Y - m_d - treatment.cost
  G11_hat = (S/mean(S)) * ((D/ps)*(Y-m1)+m1) # treated females
  G01_hat = (S/mean(S)) * (((1-D)/(1-ps))*(Y-m0)+m0) # untreated females
  G10_hat = ((1-S)/mean(1-S)) * ((D/ps)*(Y-m1)+m1) # treated males
  G00_hat = ((1-S)/mean(1-S)) * (((1-D)/(1-ps))*(Y-m0)+m0) # untreated males
  
  # compute the difference of the drs as in W_s and full vector G for W_bar
  GS1 = G11_hat - G01_hat # female
  GS0 = G10_hat - G00_hat # male
  G = alpha*GS0 + (1-alpha)*GS1 # EWM: weight corresponding to prevalence in the sample
  
  # include additional column for intercept
  X = as.matrix(cbind(1, X))
  
  # normalise the values by dividing by the max value (if not already normalised)
  max_val = max(apply(X, 1, function(x) max(abs(x))))
  X = as.matrix(X/max_val)
  
  ## Facilitate and speed up computation by "omitting" identical rows
  # find and store all unique rows 
  unique_rows = unique(X)
  
  # store the indexes of the rows which are identical
  index_unique = apply(X, 1, function(y) which(apply(unique_rows, 1, function(x) all(y == x))))
  
  # sum up the drs for individuals which have exactly the same covariate values (i.e. sum up by unique index)
  G_unique = sapply(c(1:nrow(unique_rows)), function(x) sum(G[which(x == index_unique)]))
  GS1_unique = sapply(c(1:nrow(unique_rows)), function(x) sum(GS1[which(x == index_unique)]))
  GS0_unique = sapply(c(1:nrow(unique_rows)), function(x) sum(GS0[which(x == index_unique)]))
  
  # preliminaries
  X_unique = as.matrix(unique_rows) # unique data frame
  n = nrow(X_unique) # number of unique observations
  p = ncol(X_unique) # number of columns (including the one for beta0)
  
  # number of times a certain observation is repeated in the data (required for capacity constraint)
  n_index_unique = sapply(c(1:n), function(x) sum(index_unique == x)) 
  
  
  ## initialise the model
  model  = list()
  
  # sense of optimisation, maximise welfare
  model$modelsense = "max"
  
  # model objective: G_unique as coefficient vector 
  # (betas are not directly in the objective, thus coefficient vector is 0)
  model$obj = c(G_unique, rep(0, p))
  
  # set up the linear constraint matrix (rhs only constants; requirement of gurobi solver)
  A = rbind(cbind(diag(1, nrow = n), -X_unique), # policies - betas <= 1
            cbind(diag(1, nrow = n), -X_unique), # policies - betas > 0
            c(n_index_unique, rep(0, p))) # capacity constraint
  
  # the rhs of the constraints, with tolerance (1e-6) and max treated individuals equal to the capacity constraint
  rhs = c(rep(1 - tolerance, n), rep(tolerance, n), capacity_constraint)
  
  # direction of constraints
  sense = c(rep("<=", n), rep(">", n), "<=")
  
  # additional constraint if we choose policy function class 3
  if (funclass3 == T){
    A = rbind(A, c(GS1_unique - GS0_unique, rep(0, p))) # drs as estimated treatment effect
    rhs = c(rhs, 0) # treatment_effect|S=1 - treatment_effect|S=0 >= 0
    sense = c(sense, parity_sense) # treatment effect of female >= male / female <= male
  }
  
  # combine all parts
  model$A = A
  model$rhs = rhs
  model$sense = sense
  
  # variable types, we have n binary variables (z) and the rest continuous (beta)
  model$vtype= c(rep("B", n), rep("C", p))
  
  # Put bounds on the parameter space, for z [0,1] and for the betas [-1,1] (assumption from paper)
  model$ub= c(rep(1, n), rep(1, p))
  model$lb= c(rep(0, n), rep(-1, p))
  
  # set additional parameters for the optimisation
  params = list(IntFeasTol = 1e-9, FeasibilityTol = 1e-9, TimeLimit = timelimit, # solver tolerance limits
                BarConvTol = exp(-2), # tolerance on the barrier between primal and dual solution
                Disconnected=0, # degree of exploitation of (independent) submodels
                Heuristics=0, # fraction of runtime spent on feasibility heuristics
                NodefileStart = 0.5) # max node memory before writing to hard drive
  
  # start the optimisation
  result= gurobi(model, params = params)
  
  # extract the estimated beta_hats that determine the policy pi
  beta = result$x[(n+1):(n+p)]
  
  # extract the estimated policies pi_hats
  policies = apply(X, 1, function(x) ifelse(x%*%beta > 0, 1, 0))
  
  # extract the welfare
  welfare = result$objval
  
  return(list(welfare = welfare, policies = policies, beta = beta, time = result$runtime, gap = result$mipgap))
  
}



###############################################################################

# Estimation EWM (for comparison with Fair Policy Targeting algorithms)

###############################################################################

#### Prepare data and sample
library(dplyr)
load("./results/implementation.RData")

Y = df_sample$Y
D = df_sample$D
S = df_sample$S

m1 = df_sample$m1
m0 = df_sample$m0
ps = df_sample$ps

X = df_sample %>% select(C3, X1_D, X2_D, X3_D, X4_D, X5_D)
capacity_constraint = floor(0.33*nrow(X))



#### EWM Estimation 

# EWM Policy Class 1 (no constraint)
EWM_P1 = EWM_estimation(Y, X = cbind(S, X), D, S, ps, m1, m0, alpha = mean(1-S),
                        capacity_constraint, tolerance = 1e-6, funclass3 = F, timelimit = 5000)

# EWM Policy Class 2 (fairness through unawareness)
EWM_P2 = EWM_estimation(Y, X, D, S, ps, m1, m0, alpha = mean(1-S), capacity_constraint,
                        tolerance = 1e-6, funclass3 = F, timelimit = 5000)

# EWM Policy Class 3 (additional constraint)
EWM_P3 = EWM_estimation(Y, X, D, S, ps, m1, m0, alpha = mean(1-S), capacity_constraint,
                        tolerance = 1e-6, funclass3 = T, parity_sense = ">=", timelimit = 5000)


# save(list = c("EWM_P1", "EWM_P2", "EWM_P3"), file = "./results/EWM_estimation.RData")



###############################################################################

# Estimation of the Pareto Frontier (for Fair Policy Targeting algorithms)

###############################################################################

# additional set up (discretisation of pareto frontier)
N = floor(sqrt(nrow(X)))
alpha = seq(from = 0.05, to = 0.95, length.out = N)

# load libraries for parallelisation
library(foreach)
library(future)
library(doFuture)
plan(multisession)



# Estimation Pareto Frontier

# use "dofuture" to parallelise computation (individual EWM estimations are independent of each other)
pareto_frontier = foreach(i = alpha, .combine = rbind) %dofuture% {
  
  # EWM estimation for respective alpha
  result = EWM_estimation(Y, X = cbind(S, X), D, S, ps, m1, m0, alpha = i, tolerance = 1e-3,
                          capacity_constraint, funclass3 = F, timelimit = 1000)
  
  # collect results
  c(result[[2]], result[[3]], result[[1]])
  
}

n = nrow(X)
p = ncol(X)

# extract values
pareto_policies = pareto_frontier[,(1:n)]
pareto_beta = pareto_frontier[,(n+1):(n+p+2)]
pareto_W_bar = pareto_frontier[,(n+p+3)]

# save(list = c("pareto_policies", "pareto_beta", "pareto_W_bar"), file = "./results/pareto_frontier.RData")
