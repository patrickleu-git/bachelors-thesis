###############################################################################

# Fair Policy Targeting (Counterfactual Envy): Function

###############################################################################

# Y: outcome 
# X: matrix of covariates
# D: treatment vector
# S: sensitive attribute vector
# ps: estimated propensity score
# m1: conditional mean of treated individuals
# m0: conditional mean of untreated individuals
# m_d_s: conditional mean of individual with D=d, S=s
# alpha: weight of the male welfare
# tolerance: slackness parameter (discreteness)
# W_bar: the pareto frontier
# capacity_constraint: maximal number of individuals to be treated
# start: possible starting value for the optimisation
# timelimit: maximum time spent on optimisation in seconds

optimiseCounterfactualEnvy = function(Y, X, D, S, m0, m1, m11_hat, m01_hat, m10_hat, m00_hat, ps,
                                      alpha = seq(from = 0.05, to = 0.95, length.out = N), tolerance = 1e-6,
                                      W_bar, capacity_constraint, start = NA, timelimit = 10000){
  
  
  # compute the doubly robust score (drs) 
  G11_hat = (S/mean(S)) * ((D/ps)*(Y-m1)+m1) # treated females
  G01_hat = (S/mean(S)) * (((1-D)/(1-ps))*(Y-m0)+m0) # untreated females
  G10_hat = ((1-S)/mean(1-S)) * ((D/ps)*(Y-m1)+m1) # treated males
  G00_hat = ((1-S)/mean(1-S)) * (((1-D)/(1-ps))*(Y-m0)+m0) # untreated males
  
  GS1 = G11_hat - G01_hat # female drs
  GS0 = G10_hat - G00_hat # male drs
  
  # data: intercept, sensitive attribute, covariates (case 2: average score, school rank)
  # note on second column: represents S with counterfactual covariates needed for V_pi(x,s) (x(s), s') 
  X2_female = as.matrix(cbind(1, 1, X[,-1]))
  X2_male = as.matrix(cbind(1, 0, X[,-1]))
  XX = rbind(X2_female, X2_male)
  
  # normalise the data
  max_val = max(apply(XX, 1, function(x) sum(abs(x))))
  XX = XX/max_val
  
  # set parameters
  n = nrow(XX)
  p = ncol(XX)
  
  ## initialise the model
  model = list()
  
  # minimise unfairness (counterfactual envy in this case)
  model$modelsense = "min"
  
  # objective coefficients counterfactual envy; estimator A_hat but with omitted constants!
  A_hat_female = m10_hat*S/mean(S) - m00_hat*S/mean(S) - GS1
  A_hat_male = m11_hat*(1-S)/mean(1-S) - m01_hat*(1-S)/mean(1-S) - GS0
  
  # objective: the coefficients computed above, betas and alphas are not part of objective
  model$obj = c(A_hat_female, A_hat_male, rep(0, p+N))
  
  # set the linear constraint matrix
  model$A = rbind(cbind(diag(1, nrow = n), -XX, matrix(0, nrow = n, ncol = N)), # policies - betas <= 1
                  cbind(diag(1, nrow = n), -XX, matrix(0, nrow = n, ncol = N)), # policies - betas > 0
                  c(S, (1-S), rep(0, p+N)), # capacity constraint (first female then male individuals as in objective) 
                  c(rep(0, n+p), rep(1, N))) # constraint on the u_j (C) >= 1
  
  # set the rhs (with tolerance 1e-6)
  model$rhs = c(rep(1-tolerance, n), rep(tolerance, n), capacity_constraint, 1)
  
  # set the constraint directions
  model$sense = c(rep("<=", n), # policies - betas <= 1
                  rep(">", n), # policies - betas > 0
                  "<=", # capacity constraint 
                  ">=") # constraint on the u_j >=1
  
  # set up the quadratic constraints (B)
  listnames = c(1:N)
  model$quadcon = sapply(listnames, function(x) NULL)
  
  # initiate a for loop to add the N=sqrt(n) constraints
  for (i in 1:N){
    
    # initiate the Qc matrix
    model$quadcon[[i]]$Qc = matrix(0, nrow = n+p+N, ncol = n+p+N)
    # set coefficients in the column of the respective u_j
    model$quadcon[[i]]$Qc[1:(n/2), n+p+i] = (1-alpha[i])*GS1 # female weighted with alpha
    model$quadcon[[i]]$Qc[(n/2 + 1):n, n+p+i] = alpha[i]*GS0 # male weighted with 1-alpha
    
    # initiate q vector for the linear terms (optimal welfare W_bar)
    model$quadcon[[i]]$q = rep(0, n+p+N)
    # insert the respective welfare (- as bring to lhs)
    model$quadcon[[i]]$q[n+p+i] = -W_bar[i]
    
    # set the rhs
    model$quadcon[[i]]$rhs = -tolerance
    
    # direction of constraint
    model$quadcon[[i]]$sense = ">="
    
  }
  
  # set the variable type
  model$vtype = c(rep("B", n), rep("C", p), rep("B", N))
  
  # set bounds: for z_i and u_j [0,1], for betas [-1,1]
  model$ub = rep(1, n+p+N)
  model$lb = c(rep(0, n), rep(-1, p), rep(0, N))
  
  if (is.na(start[1]) == F) model$start = start
  
  # set a list of parameters
  params = list(IntFeasTol = 1e-9, FeasibilityTol = 1e-9, TimeLimit = timelimit, # tolerance limits
                BarConvTol = exp(-2), # tolerance on the barrier between primal and dual solution
                Disconnected=0, # degree of exploitation of (independent) submodels
                Heuristics=0, # fraction of runtime spent on feasibility heuristics
                NodefileStart = 0.5) # max node memory before writing to hard drive
  
  # results
  result = gurobi(model, params = params)
  
  # extract betas, weight, policies, objective value
  beta = result$x[(n+1):(n+p)] # the betas for deciding on the policy
  u = result$x[(n+p+2):length(result$x)] # which u_j is equal to 1
  weight = (1-alpha[u == 1]) # store the respective optimal weight alpha_j
  policies = apply(cbind(1, as.matrix(X)), 1, function(x) ifelse(x%*%beta > 0, 1, 0)) # the estimated policies
  objval = result$objval # the minimised unfairness (counterfactual envy)
  
  # return estimated betas, weight of the female group in the optimisation, policies and objective value
  return(list(beta = beta, weight_female = weight, policies = policies, objval = objval, 
              result.x = result$x, time = result$runtime, gap = result$mipgap))
  
}

###############################################################################

# Estimation FPT Counterfactual Envy

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
m11_hat = df_sample$m11_hat
m10_hat = df_sample$m10_hat
m01_hat = df_sample$m01_hat
m00_hat = df_sample$m00_hat
ps = df_sample$ps

X = df_sample %>% select(C3, X1_D, X2_D, X3_D, X4_D, X5_D)
capacity_constraint = floor(0.33*nrow(X))

N = floor(sqrt(nrow(X)))
alpha = seq(from = 0.05, to = 0.95, length.out = N)

G11_hat = (S/mean(S)) * ((D/ps)*(Y-m1)+m1) # treated females
G01_hat = (S/mean(S)) * (((1-D)/(1-ps))*(Y-m0)+m0) # untreated females
G10_hat = ((1-S)/mean(1-S)) * ((D/ps)*(Y-m1)+m1) # treated males
G00_hat = ((1-S)/mean(1-S)) * (((1-D)/(1-ps))*(Y-m0)+m0) # untreated males

GS1 = G11_hat - G01_hat # female drs
GS0 = G10_hat - G00_hat # male drs


#### FPT Counterfactual Envy Estimation

## Build start vector with the help of the pareto frontier values to help find an initial solution in the optimisation
Xfem = as.matrix(cbind(1, 1, X))
Xmale = as.matrix(cbind(1, 0, X))

# compute the policies that female / male individuals would be prescribed
fem_policies = t(apply(pareto_beta, 1, function(x) sapply(Xfem%*%x, function(y) ifelse(y > 0, 1, 0))))
male_policies = t(apply(pareto_beta, 1, function(x) sapply(Xmale%*%x, function(y) ifelse(y > 0, 1, 0))))

# compute the welfare female / male individuals would achieve with their allocated policy
fem_welfare = apply(fem_policies, 1, function(x) sum(GS1*x))
male_welfare = apply(male_policies, 1, function(x) sum(GS0*x))

# objective for the start vector 
objective_fem = apply(fem_policies, 1, function(x) sum(m10_hat*x*S)/mean(S) +  sum(m00_hat*(1 - x)*S)/mean(S)) -  male_welfare
objective_male = apply(male_policies, 1, function(x) sum(m11_hat*x*(1 - S))/(1 - mean(S)) +  sum(m01_hat*(1 - x)*(1 - S))/(1 - mean(S))) - fem_welfare
objective_start = objective_fem + objective_male

# compute the least unfair objective combination, as there might be multiple take the one closest to equal weighting
least_unfair = which(objective_start == min(objective_start))
least_unfair = least_unfair[which.min(abs(least_unfair - N/2))]

# combine elements to the start vector
start_envy = c(fem_policies[least_unfair,], male_policies[least_unfair,], pareto_beta[least_unfair,], ifelse(c(1:N) == least_unfair, 1, 0))

# optimise!
FPT_counterfactual_envy = optimiseCounterfactualEnvy(Y, X = cbind(S, X), D, S, m0, m1, m11_hat, m01_hat, m10_hat, m00_hat, 
                                                     ps, alpha = alpha, tolerance = 1e-6, capacity_constraint, 
                                                     W_bar = pareto_W_bar, start = start_envy, timelimit = 15000)


# save(list = c("FPT_counterfactual_envy"), file = "./results/fpt_counterfactual_envy.RData")
