###############################################################################

# Estimation of the Pareto Frontier for the Frontier Plots

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

N = 100
alpha = seq(from = 0.05, to = 0.95, length.out = N)


# load libraries
library(slam)
library(gurobi)
library(foreach)
library(future)
library(doFuture)
plan(multisession)

# policy 1 and female >= male
frontier_P1_fem = foreach(i = alpha, .combine = rbind) %dofuture% {
  
  result = EWM_estimation(Y, X = cbind(S,X), D, S, ps, m1, m0, alpha = i, tolerance = 1e-3,
                          capacity_constraint, funclass3 = T, parity_sense = ">=", timelimit = 500)
  
  c(result[[2]], result[[3]], result[[1]])
}

# policy 1 and female <= male
frontier_P1_male = foreach(i = alpha, .combine = rbind) %dofuture% {
  
  result = EWM_estimation(Y, X = cbind(S,X), D, S, ps, m1, m0, alpha = i, tolerance = 1e-3,
                          capacity_constraint, funclass3 = T, parity_sense = "<=", timelimit = 500)
  
  c(result[[2]], result[[3]], result[[1]])
}

# policy 2 (not using S in decision rule) and female >= male
frontier_P2_fem = foreach(i = alpha, .combine = rbind) %dofuture% {
  
  result = EWM_estimation(Y, X, D, S, ps, m1, m0, alpha = i, tolerance = 1e-3,
                          capacity_constraint, funclass3 = T, parity_sense = ">=", timelimit = 500)
  
  c(result[[2]], result[[3]], result[[1]])
}

# policy 2 (not using S in decision rule) and female <= male
frontier_P2_male = foreach(i = alpha, .combine = rbind) %dofuture% {
  
  result = EWM_estimation(Y, X, D, S, ps, m1, m0, alpha = i, tolerance = 1e-3,
                          capacity_constraint, funclass3 = T, parity_sense = "<=", timelimit = 500)
  
  c(result[[2]], result[[3]], result[[1]])
}

# save(list = c("frontier_P1_fem", "frontier_P1_male", "frontier_P2_fem", "frontier_P2_male"), file = "./results/frontier_plot.RData")

###############################################################################

# Pareto Frontier Plots

###############################################################################

# load library for plot
library(ggplot2)
library(ggthemes)
library(extrafont)
font_import()

# function for the comparison of the welfares and detection of pareto dominant allocations;
# corresponding largely to Viviano & Bradic's (2024) function in supplementary materials
pareto_dominance = function(female_function, male_function){
  
  # combine the solutions into matrix
  combined = rbind(female_function, male_function)
  
  # dummy vector indicating whether an allocation is not dominated
  dominated = rep(0, nrow(combined))
  
  # initiate for loops, compare each allocation i against all others j
  for(i in 1:nrow(combined)){
    indicator = 0
    for(j in 1:nrow(combined)){
      if(i != j){
        # dominated is 1 if i < j for both, female and male individuals, otherwise 0
        indicator = max(combined[i,1] < combined[j,1] & combined[i,2] < combined[j,2], indicator)
      }
    }
    # vector position of dominated allocations
    dominated[i] = indicator
  }
  
  # return only non dominated allocations
  return(combined[dominated == 0,])
}


# difference of the doubly robust scores for female and male (welfare improvement)
# (multiply with policy if pi = 1; otherwise only baseline welfare)
Wfem = S/mean(S) * (((D/ps)*(Y-m1)+m1) - ((1-D)/(1-ps)*(Y-m0)+m0))
Wmale = (1-S)/mean(1-S) * (((D/ps)*(Y-m1)+m1) - ((1-D)/(1-ps)*(Y-m0)+m0))

# baseline effect, i.e. doubly robust score of untreated individuals 
# (multiply with S=s/mean(S=s) )
baseline_fem = S/mean(S) * ((1-D)/(1-ps)*(Y-m0)+m0)
baseline_male = (1-S)/mean(1-S) * ((1-D)/(1-ps)*(Y-m0)+m0)

# covariates as matrix; dimensions of the matrix
X = as.matrix(df_sample %>% select(C3, X1_D, X2_D, X3_D, X4_D, X5_D))
n = nrow(X)
p = ncol(X)+1

# data frame with the covariates for female / male
df_fem = as.matrix(cbind(1,1,X))
df_male = as.matrix(cbind(1,0,X))

# load the results of the frontier computations, store the betas
load("./results/frontier_plot.RData")
beta_P1_fem = frontier_P1_fem[, (n+1):(n+p+1)] # p+1 due to added S
beta_P1_male = frontier_P1_male[, (n+1):(n+p+1)] # p+1 due to added S
beta_P2_fem = frontier_P2_fem[, (n+1):(n+p)]
beta_P2_male = frontier_P2_male[, (n+1):(n+p)]


# Compute the welfare

# Policy function class 1

# with female >= male
welfare_fem1 = apply(beta_P1_fem, 1, function(x) mean(sapply(df_fem%*%x, function(y) ifelse(y >= 0, 1, 0))*Wfem + baseline_fem))
welfare_male1 = apply(beta_P1_fem, 1, function(x) mean(sapply(df_male%*%x, function(y) ifelse(y >= 0, 1, 0))*Wmale + baseline_male))
welfare_fem_function = cbind(welfare_fem1, welfare_male1)

# with female <= male
welfare_fem2 = apply(beta_P1_male, 1, function(x) mean(sapply(df_fem%*%x, function(y) ifelse(y >= 0, 1, 0))*Wfem + baseline_fem))
welfare_male2 = apply(beta_P1_male, 1, function(x) mean(sapply(df_male%*%x, function(y) ifelse(y >= 0, 1, 0))*Wmale + baseline_male))
welfare_male_function = cbind(welfare_fem2, welfare_male2)

# compute pareto frontier
welfare1 = pareto_dominance(welfare_fem_function, welfare_male_function)
# add minimum values for nicer plot
welfare1 = rbind(c(0.4, max(welfare1[,2])), welfare1, c(max(welfare1[,1]), 0.45))


# Policy function class 2

# with female >= male
welfare_fem1 = apply(beta_P2_fem, 1, function(x) mean(sapply(df_fem[,-2]%*%x, function(y) ifelse(y >= 0, 1, 0))*Wfem + baseline_fem))
welfare_male1 = apply(beta_P2_fem, 1, function(x) mean(sapply(df_male[,-2]%*%x, function(y) ifelse(y >= 0, 1, 0))*Wmale + baseline_male))
welfare_fem_function = cbind(welfare_fem1, welfare_male1)

# with female <= male
welfare_fem2 = apply(beta_P2_male, 1, function(x) mean(sapply(df_fem[,-2]%*%x, function(y) ifelse(y >= 0, 1, 0))*Wfem + baseline_fem))
welfare_male2 = apply(beta_P2_male, 1, function(x) mean(sapply(df_male[,-2]%*%x, function(y) ifelse(y >= 0, 1, 0))*Wmale + baseline_male))
welfare_male_function = cbind(welfare_fem2, welfare_male2)

# compute pareto frontier
welfare2 = pareto_dominance(welfare_fem_function, welfare_male_function)
# add minimum values for nicer plot
welfare2 = rbind(c(0.4, max(welfare2[,2])), welfare2, c(max(welfare2[,1]), 0.45))


# Policy function class 3

# only female >= male
welfare3 = cbind(welfare_fem1, welfare_male1)
# add minimum values for nicer plot
welfare3 = rbind(c(0.4, max(welfare3[,2])), welfare3, c(max(welfare3[,1]), 0.45))


# name the allocations after the policy function class
policy_class = c(rep("Class 1", dim(welfare1)[1]), rep("Class 2", dim(welfare2)[1]), rep("Class 3", dim(welfare3)[1]))

# create tibble
df_pareto = as.data.frame(cbind(rbind(welfare1, welfare2, welfare3), policy_class))

# update column names
names(df_pareto) = c('Wfemale', 'Wmale', 'Class')

# coerce columns from character to numeric
df_pareto = transform(df_pareto, Wfemale = as.numeric(Wfemale), Wmale = as.numeric(Wmale))


# create plot
ggplot(data = df_pareto, aes(x = Wfemale, y = Wmale)) +
  geom_line(aes(color = Class), linetype = 1, linewidth = 0.7, show.legend = F) +
  geom_ribbon(aes(ymin = min(Wmale), ymax = Wmale, fill = Class)) +
  scale_color_manual(values=c("#69b3a2", "#404080", "darkgrey")) +
  scale_fill_manual(values = alpha(c("#69b3a2", "#404080", "darkgrey"), 0.4),
                    labels = c("Policy Class 1", "Policy Class 2", "Policy Class 3"),
                    name = "Policy Function Class") +
  xlab("Welfare female individuals") + 
  ylab("Welfare male individuals") +
  theme_hc() +
  theme(text = element_text(family = "Times New Roman"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12))

ggsave("frontier_plot.png", plot = last_plot(), path = "./results/", height = 7, width = 12)
