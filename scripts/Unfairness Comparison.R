###############################################################################

# Unfairness Plots

###############################################################################

# load required libraries
library(tidyverse)
library(ggthemes)
library(extrafont)
font_import()

# load the required results
load("./results/implementation.RData")
load("./results/fpt_counterfactual_envy.RData")
load("./results/fpt_disp.RData")
load("./results/fpt_absolute_disp.RData")
load("./results/fpt_regularised.RData")
load("./results/EWM_estimation.RData")

# some preliminary parametere definitions
Y = df_sample$Y
D = df_sample$D
S = df_sample$S

# Unfairness comparison with V = prediction disparity

# V with the prediciton disparity objective
V_pred_disp = (sum((1-S)*FPT_pred_disparity$policies)/sum(1-S) 
               - sum(S*FPT_pred_disparity$policies)/sum(S))

# V with the EWM method
V_EWM1 = (sum((1-S)*EWM_P1$policies)/sum(1-S) 
          - sum(S*EWM_P1$policies)/sum(S))

V_EWM2 = (sum((1-S)*EWM_P2$policies)/sum(1-S) 
          - sum(S*EWM_P2$policies)/sum(S))

V_EWM3 = (sum((1-S)*EWM_P3$policies)/sum(1-S) 
          - sum(S*EWM_P3$policies)/sum(S))

# V with the regularised approach  
V_regularised004 = (sum((1-S)*FPT_regularised_004$policies)/sum(1-S) 
                    - sum(S*FPT_regularised_004$policies)/sum(S))

V_regularised010 = (sum((1-S)*FPT_regularised_010$policies)/sum(1-S) 
                    - sum(S*FPT_regularised_010$policies)/sum(S))

V_regularised100 = (sum((1-S)*FPT_regularised_100$policies)/sum(1-S) 
                    - sum(S*FPT_regularised_100$policies)/sum(S))


# Unfairness with V = prediction disparity absolute

# Prediction disparity absolute comparison
V_pred_disp_abs = abs( sum((1-S)*FPT_pred_disparity_abs$policies)/sum(1-S) 
                       - sum(S*FPT_pred_disparity_abs$policies)/sum(S) )

# V with the EWM method
V_EWM1_abs = abs( sum((1-S)*EWM_P1$policies)/sum(1-S) 
                  - sum(S*EWM_P1$policies)/sum(S) )

V_EWM2_abs = abs( sum((1-S)*EWM_P2$policies)/sum(1-S) 
                  - sum(S*EWM_P2$policies)/sum(S) )

V_EWM3_abs = abs( sum((1-S)*EWM_P3$policies)/sum(1-S) 
                  - sum(S*EWM_P3$policies)/sum(S) )

# V with the regularised approach  
V_regularised004_abs = abs( sum((1-S)*FPT_regularised_004$policies)/sum(1-S) 
                            - sum(S*FPT_regularised_004$policies)/sum(S) )

V_regularised010_abs = abs( sum((1-S)*FPT_regularised_010$policies)/sum(1-S) 
                            - sum(S*FPT_regularised_010$policies)/sum(S) )

V_regularised100_abs = abs( sum((1-S)*FPT_regularised_100$policies)/sum(1-S) 
                            - sum(S*FPT_regularised_100$policies)/sum(S) )


# combine into dataframe
V = c(V_pred_disp, V_EWM1, V_EWM2, V_EWM3, V_regularised004, V_regularised010, 
      V_regularised100, V_pred_disp_abs, V_EWM1_abs, V_EWM2_abs, V_EWM3_abs, 
      V_regularised004_abs, V_regularised010_abs, V_regularised100_abs)

names(V) = rep(c("Prediction Disparity", paste("Welfare Max.", c(1:3)), 
             paste("Regularised", c(1:3))),2)

df_V = enframe(V)

df_V$metric = as.factor(c(rep("Prediction Disparity", 7), rep("Prediction Disparity Absolute", 7)))

df_V = df_V %>% 
  mutate(name = factor(name, levels = c("Prediction Disparity", paste("Welfare Max.", c(1:3)), 
                                        paste("Regularised", c(1:3)))))

# create the plot
ggplot(df_V, aes(x = metric, y = value, fill = name, colour = name)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.6) +
  scale_y_continuous(breaks = c(0.15, 0, -0.15, -0.30, -0.45, -0.6))+
  scale_fill_manual(values=c("#507872", "#69b3a2", "#acc6aa", "darkgrey", 
                             "#595959", "#6a5acd", "#404080"), name = "Method")+
  scale_colour_manual(values=alpha(c("#507872", "#69b3a2", "#acc6aa", "darkgrey", 
                               "#595959", "#6a5acd", "#404080"), 0.6), guide = "none")+
  xlab("") +
  ylab("Unfairness Level")+
  theme_hc()+
  theme(text = element_text(family = "Times New Roman"),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12))


ggsave("unfairness.png", plot = last_plot(), path = "./results/", height = 7, width = 12)
























