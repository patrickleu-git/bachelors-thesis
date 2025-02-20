###############################################################################

# Welfare Comparison: Empirical Application

###############################################################################

load("./results/implementation.RData")

Y = df_sample$Y
D = df_sample$D
S = df_sample$S

m1 = df_sample$m1
m0 = df_sample$m0
ps = df_sample$ps

# difference of the doubly robust scores for female and male (welfare improvement)
# (multiply with policy if pi = 1; otherwise only baseline welfare)
Wfem = S/mean(S) * (((D/ps)*(Y-m1)+m1) - ((1-D)/(1-ps)*(Y-m0)+m0))
Wmale = (1-S)/mean(1-S) * (((D/ps)*(Y-m1)+m1) - ((1-D)/(1-ps)*(Y-m0)+m0))

# baseline effect, i.e. doubly robust score of untreated individuals 
# (multiply with S=s/mean(S=s) )
baseline_fem = S/mean(S) * ((1-D)/(1-ps)*(Y-m0)+m0)
baseline_male = (1-S)/mean(1-S) * ((1-D)/(1-ps)*(Y-m0)+m0)

# covariates as matrix
X = as.matrix(df_sample %>% select(C3, X1_D, X2_D, X3_D, X4_D, X5_D))

# data frame with the covariates for female / male
df_fem = as.matrix(cbind(1,1,X))
df_male = as.matrix(cbind(1,0,X))

###############################################################################
# Envy
load("./results/fpt_counterfactual_envy.RData")

# betas for policies and weight of female group
beta_envy = FPT_counterfactual_envy$beta
weight_envy = FPT_counterfactual_envy$weight_female
# welfare under the envy measure for female and male individuals
# note: matrix multiplication to extract policy, multiply with welfare components 
welfare_fem_envy = mean(sapply(df_fem%*%beta_envy, function(y) ifelse(y >= 0, 1, 0))*Wfem + baseline_fem)
welfare_male_envy = mean(sapply(df_male%*%beta_envy, function(y) ifelse(y >= 0, 1, 0))*Wmale + baseline_male)


###############################################################################
# Prediction Disparity
load("./results/fpt_disp.RData")

# betas for policies and weight of female group
beta_pred_disp = FPT_pred_disparity$beta
weight_pred_disp = FPT_pred_disparity$weight_female
# welfare under the envy measure for female and male individuals
# note: matrix multiplication to extract policy, multiply with welfare components 
welfare_fem_pred_disp = mean(sapply(df_fem%*%beta_pred_disp, function(y) ifelse(y >= 0, 1, 0))*Wfem + baseline_fem)
welfare_male_pred_disp = mean(sapply(df_male%*%beta_pred_disp, function(y) ifelse(y >= 0, 1, 0))*Wmale + baseline_male)


###############################################################################
# Prediction Disparity Absolute
load("./results/fpt_absolute_disp.RData")

# betas for policies and weight of female group
beta_pred_disp_abs = FPT_pred_disparity_abs$beta
weight_pred_disp_abs = FPT_pred_disparity_abs$weight_female
# welfare under the envy measure for female and male individuals
# note: matrix multiplication to extract policy, multiply with welfare components 
welfare_fem_pred_disp_abs = mean(sapply(df_fem%*%beta_pred_disp_abs, function(y) ifelse(y >= 0, 1, 0))*Wfem + baseline_fem)
welfare_male_pred_disp_abs = mean(sapply(df_male%*%beta_pred_disp_abs, function(y) ifelse(y >= 0, 1, 0))*Wmale + baseline_male)


###############################################################################
# Regularisation Approach
load("./results/fpt_regularised.RData")

# betas for policies 
beta_004 = FPT_regularised_004$beta
beta_010 = FPT_regularised_010$beta
beta_100 = FPT_regularised_100$beta

# welfare under the envy measure for female and male individuals
# note: matrix multiplication to extract policy, multiply with welfare components 
welfare_fem_regularised_004 = mean(sapply(df_fem%*%beta_004, function(y) ifelse(y >= 0, 1, 0))*Wfem + baseline_fem)
welfare_male_regularised_004 = mean(sapply(df_male%*%beta_004, function(y) ifelse(y >= 0, 1, 0))*Wmale + baseline_male)

welfare_fem_regularised_010 = mean(sapply(df_fem%*%beta_010, function(y) ifelse(y >= 0, 1, 0))*Wfem + baseline_fem)
welfare_male_regularised_010 = mean(sapply(df_male%*%beta_010, function(y) ifelse(y >= 0, 1, 0))*Wmale + baseline_male)

welfare_fem_regularised_100 = mean(sapply(df_fem%*%beta_100, function(y) ifelse(y >= 0, 1, 0))*Wfem + baseline_fem)
welfare_male_regularised_100 = mean(sapply(df_male%*%beta_100, function(y) ifelse(y >= 0, 1, 0))*Wmale + baseline_male)


###############################################################################
# EWM
load("./results/EWM_estimation.RData")

# betas for policies
beta_P1 = EWM_P1$beta
beta_P2 = EWM_P2$beta
beta_P3 = EWM_P3$beta

# welfare under the envy measure for female and male individuals
# note: matrix multiplication to extract policy, multiply with welfare components 
welfare_fem_EWM_P1 = mean(sapply(df_fem%*%beta_P1, function(y) ifelse(y >= 0, 1, 0))*Wfem + baseline_fem)
welfare_male_EWM_P1 = mean(sapply(df_male%*%beta_P1, function(y) ifelse(y >= 0, 1, 0))*Wmale + baseline_male)

welfare_fem_EWM_P2 = mean(sapply(cbind(1,X)%*%beta_P2, function(y) ifelse(y >= 0, 1, 0))*Wfem + baseline_fem)
welfare_male_EWM_P2 = mean(sapply(cbind(1,X)%*%beta_P2, function(y) ifelse(y >= 0, 1, 0))*Wmale + baseline_male)

welfare_fem_EWM_P3 = mean(sapply(cbind(1,X)%*%beta_P3, function(y) ifelse(y >= 0, 1, 0))*Wfem + baseline_fem)
welfare_male_EWM_P3 = mean(sapply(cbind(1,X)%*%beta_P3, function(y) ifelse(y >= 0, 1, 0))*Wmale + baseline_male)


###############################################################################

## Combine all elements

# combine the welfares in a vector
welfares_fem = round(c(welfare_fem_envy, welfare_fem_pred_disp, welfare_fem_pred_disp_abs, 
                 welfare_fem_EWM_P1, welfare_fem_EWM_P2, welfare_fem_EWM_P3,
                 welfare_fem_regularised_004, welfare_fem_regularised_010, welfare_fem_regularised_100),3)

welfares_male = round(c(welfare_male_envy, welfare_male_pred_disp, welfare_male_pred_disp_abs, 
                  welfare_male_EWM_P1, welfare_male_EWM_P2, welfare_male_EWM_P3,
                  welfare_male_regularised_004, welfare_male_regularised_010, welfare_male_regularised_100),3)


# the importance weights (weight for EWM is mean(S) by definition); lambda for regularised approach
weights = round(c(weight_envy, weight_pred_disp, weight_pred_disp_abs, 
                  rep(mean(S), 3), 0.004, 0.010, 0.100),3)

times = round(c(FPT_counterfactual_envy$time, FPT_pred_disparity$time, FPT_pred_disparity_abs$time,
                EWM_P1$time, EWM_P2$time, EWM_P3$time, 
                FPT_regularised_004$time, FPT_regularised_010$time, FPT_regularised_100$time),1)

gaps = round(c(FPT_counterfactual_envy$gap, FPT_pred_disparity$gap, FPT_pred_disparity_abs$gap,
               EWM_P1$gap, EWM_P2$gap, EWM_P3$gap, 
               FPT_regularised_004$gap, FPT_regularised_010$gap, FPT_regularised_100$gap)*100, 3)

## Construct table 
final_table = cbind(welfares_fem, welfares_male, weights, times, gaps)

colnames(final_table) = c("Welfare Female", "Welfare Male", "Importance Weight", "Time (s)", "Gap (%)")

rownames(final_table) = c("Counterfactual Envy", "Prediction Disparity", "Prediction Disparity Abs", 
                          paste("Welfare Max.", c(1:3)), paste("Regularised", c(1:3)))

# xtable(final_table)






