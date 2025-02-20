###########################################################

# Implementation: National Study of Learning Mindset

###########################################################

# load libraries
library(tidyverse)
library(data.table)
library(mltools)
library(caret)
library(glmnet)
library(extrafont)
font_import()

#### Data preparation
data = read_csv("./data/synthetic_data.csv")

# prepare data 
data = data %>% mutate(S = ifelse(data$C2 == 2, 1, 0), # recode S (0,1)
                       D = Z,                          # treatment vector D
                       Y = 1 / (1+exp(-data$Y)))       # normalise Y (0,1)

# one hot encoding for categorical data
data$C1 = as_factor(data$C1)
data$XC = as_factor(data$XC)
data = one_hot(as.data.table(data))



#### Propensity score estimation

# function for the propensity score estimation
# df: data frame for the estimation
# D: the treatment vector
# K: the number of folds for cross-fitting
estimate_ps = function(df, D, K = 5) {
  
  # load required libraries
  library(caret)
  library(glmnet)
  
  # create the folds for the cross-fitting
  set.seed(777)
  folds = createFolds(D, k = K)
  
  # storage for the predictions
  ps = rep(0, nrow(df))
  
  # initiate for loop
  for (i in 1:K){
    
    # split data in training and test set
    test_indices = unlist(folds[i])
    train_indices = unlist(folds[-i])
    
    # remove D as the first column of df
    x_train = as.matrix(df[train_indices, -1])
    y_train = D[train_indices]
    x_test = as.matrix(df[test_indices, -1])
    
    # cross validation for lambda
    ps_reg = cv.glmnet(x_train, as.factor(y_train), family = "binomial")
    
    # make predictions using the test set and optimal lambda
    ps[test_indices] = predict(ps_reg, s = ps_reg$lambda.min, newx = x_test, type = "response")
    
  }
  
  # return the ps
  return(ps)
  
}

# prepare the df for the propensity score estimation
df_ps = data %>% select(D, S, everything(), -Y, -C2, -schoolid, -Z)

# estimate and add the propensity score to df_ps
df_ps$ps = estimate_ps(df = df_ps, D = df_ps$D, K = 5)

# plot to check overlap between treated and untreated (Positivity Check)
ggplot(df_ps) +
  geom_histogram(aes(x = ps, fill = "darkgrey"), color = "white", alpha = 0.6, position = 'identity', show.legend = T) +
  geom_histogram(aes(x = ps, fill = as.factor(D)), colour = "white", alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("#69b3a2", "#404080", "darkgrey"), 
                    labels = c("untreated", "treated", "cumulative"), name = "Treatment Status")+
  labs(x = "Propensity Score", y = "Count") +
  theme_hc() +
  theme(text = element_text(family = "Times New Roman"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14))

ggsave("positivity_check.png", plot = last_plot(), path = "./results/", height = 7, width = 12)

#### Conditional mean estimation

# function for the conditional mean estimation
# Y: outcome variable
# S: sensitive attribute
# df: data frame for the estimation
# K: number of folds for crossfitting
estimate_conditional_mean = function(Y, S, df_mean_reg, K = 5){
  
  # create folds for cross-fitting
  set.seed(777)
  folds = createFolds(Y, k = K)
  
  # create storage for estimations m_d_s
  m11_hat = m10_hat = m01_hat = m00_hat = rep(0, length(Y))
  
  # initiate the for loop over the folds
  for (i in 1:K){
    
    # split data in training and testing set
    train_indices = unlist(folds[-i])
    test_indices = unlist(folds[i])
    
    # main regression for the prediction
    mean_reg = cv.glmnet(as.matrix(df_mean_reg[train_indices, ]), Y[train_indices], family = "gaussian")
    lambda = mean_reg$lambda.min
    
    # m11_hat: test data (assuming D = 1, S = 1) and prediction
    df_m11 = cbind(1, 1, df_mean_reg[,-c(1,2)])
    m11_hat[test_indices] = predict(mean_reg, s = lambda,  newx = as.matrix(df_m11)[test_indices, ], type = "response")
    
    # m01_hat: test data (assuming D = 0, S = 1) and prediction
    df_m01 = cbind(0, 1, df_mean_reg[,-c(1,2)])
    m01_hat[test_indices] = predict(mean_reg, s = lambda, newx = as.matrix(df_m01)[test_indices, ], type = "response")
    
    # m10_hat: test data (assuming D = 1, S = 1) and prediction
    df_m10 = cbind(1, 0, df_mean_reg[,-c(1,2)])
    m10_hat[test_indices] = predict(mean_reg, s = lambda, newx = as.matrix(df_m10)[test_indices, ], type = "response")
    
    # m00_hat: test data (assuming D = 0, S = 0) and prediction
    df_m00 = cbind(0, 0, df_mean_reg[,-c(1,2)])
    m00_hat[test_indices] = predict(mean_reg, s = lambda, newx = as.matrix(df_m00)[test_indices, ], type = "response")
    
  }
  
  # compute the conditional means of treated and untreated individuals
  m1 = S*m11_hat + (1-S)*m10_hat
  m0 = S*m01_hat + (1-S)*m00_hat
  
  return(list(m1 = m1, m0 = m0, m11_hat = m11_hat, m01_hat = m01_hat, m10_hat = m10_hat, m00_hat = m00_hat))
  
}

# prepare the df for the conditional mean estimation
df_mean_reg = data %>% select(D, S, everything(), -Y, -C2, -Z)

# estimate the conditional means
elements = estimate_conditional_mean(Y = data$Y, S = data$S, df_mean_reg = df_mean_reg, K = 5)



#### Prepare final dataset to be used for algorithms
df = data %>% select(Y, D, S, C3, paste0("X", c(1:5))) %>% 
  mutate(X1_D = ifelse(X1 > mean(X1), 1, 0),
         X2_D = ifelse(X2 > mean(X2), 1, 0),
         X3_D = ifelse(X3 > mean(X3), 1, 0),
         X4_D = ifelse(X4 > mean(X4), 1, 0),
         X5_D = ifelse(X5 > mean(X5), 1, 0),
         ps = df_ps$ps, m1 = elements$m1, m0 = elements$m0,
         m11_hat = elements$m11_hat,
         m10_hat = elements$m10_hat,
         m01_hat = elements$m01_hat,
         m00_hat = elements$m00_hat)

# take a random sample of n = 500 observations for the optimisation
set.seed(000)
df_sample = slice_sample(.data = df, n = 500)


# comparison between sample and full dataset for important variables
df_sample_comparison = df_sample %>% select(D, S, ps, m1, m0)
data_comparison = data %>% select(D, S) %>% 
  mutate(ps = df_ps$ps, m1 = elements$m1, m0 = elements$m0)

# xtable(summary(df_sample_comparison))
# xtable(summary(data_comparison))
# save(list = c("df_sample"), file = "./results/implementation.RData")
