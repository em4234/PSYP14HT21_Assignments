#####################################################
##### ASSIGMENT PART 1: HIERARCHICAL REGRESSION #####
#####################################################

# Package loading 

library(tidyverse)
library(gsheet)
library(psych)	
library(gridExtra)
library(dplyr)
library(ggplot2)
library(car) 
library(lmtest) 
library(sandwich) 
library(boot) 
library(lmboot)
library(cAIC4)
library(r2glmm)
library(lme4)
library(lmerTest)
library(MuMIn)

# Costume function to obtain regression coefficients source: https://www.statmethods.net/advstats/bootstrapping.html

bs_to_boot <- function(model, data, indices) {
  d <- data[indices, ] # allows boot to select sample
  fit <- lm(formula(model), data = d)
  return(coef(fit))
}
# function to obtain adjusted R^2 source: https://www.statmethods.net/advstats/bootstrapping.html
# (partially modified)

adjR2_to_boot <- function(model, data, indices) {
  d <- data[indices, ] # allows boot to select sample
  fit <- lm(formula(model), data = d)
  return(summary(fit)$adj.r.squared)
}
# Computing the booststrap BCa (bias-corrected and
# accelerated) bootstrap confidence intervals by Elfron
# (1987) This is useful if there is bias or skew in the
# residuals.

confint.boot <- function(model, data = NULL, R = 1000) {
  if (is.null(data)) {
    data = eval(parse(text = as.character(model$call[3])))
  }
  boot.ci_output_table = as.data.frame(matrix(NA, nrow = length(coef(model)),
                                              ncol = 2))
  row.names(boot.ci_output_table) = names(coef(model))
  names(boot.ci_output_table) = c("boot 2.5 %", "boot 97.5 %")
  2
  results.boot = results <- boot(data = data, statistic = bs_to_boot,
                                 R = 1000, model = model)
  for (i in 1:length(coef(model))) {
    boot.ci_output_table[i, ] = unlist(unlist(boot.ci(results.boot,
                                                      type = "bca", index = i))[c("bca4", "bca5")])
  }
  return(boot.ci_output_table)
}
# Computing the booststrapped confidence interval for a
# linear model using wild bottstrapping as descibed by Wu
# (1986) <doi:10.1214/aos/1176350142> requires the lmboot
# pakcage

wild.boot.confint <- function(model, data = NULL, B = 1000) {
  if (is.null(data)) {
    data = eval(parse(text = as.character(model$call[3])))
  }
  wild_boot_estimates = wild.boot(formula(model), data = data,
                                  B = B)
  result = t(apply(wild_boot_estimates[[1]], 2, function(x) quantile(x,
                                                                     probs = c(0.025, 0.975))))
  return(result)
}

# Collection of dataset   

raw_data <- read.csv("https://raw.githubusercontent.com/kekecsz/PSYP14-Advanced-Scientific-Methods/main/Home%20assignment/home_sample_1.csv")
view(raw_data)
raw_data

# Descriptive analysis of data 

describe(raw_data)
dim(raw_data)
str(raw_data)
names(raw_data)

row.names(raw_data)
nrow(raw_data)
ncol(raw_data)
head(raw_data)
tail(raw_data)

raw_data %>% 
  summary()

# Checking data for coding errors and outliers 
ggplot(raw_data, aes(y=pain))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(raw_data, aes(y=age))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(raw_data, aes(y=sex))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(raw_data, aes(y=STAI_trait))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(raw_data, aes(y=pain_cat))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(raw_data, aes(y=mindfulness))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(raw_data, aes(y=cortisol_serum))+stat_boxplot(geom="errorbar")+geom_boxplot()

# Identiying and correcting coding errors   
raw_data$pain[raw_data$pain>10]
raw_data$pain[raw_data$pain==55]<-5
view(raw_data)

raw_data$STAI_trait[raw_data$STAI_trait<20]
raw_data$STAI_trait[raw_data$STAI_trait>80]
raw_data$STAI_trait[raw_data$STAI_trait==4.2]<-42
view(raw_data)

data1 <- raw_data
view(data1)

# Creating model 1 and 2 

model_1 <- lm(pain ~ sex + age, data = data1)

# Checking model 2 for outliers by using Cook's distance

model_2 <- lm(pain ~ sex + age + STAI_trait + pain_cat + mindfulness + cortisol_serum, data = data1)

model_2 %>%
  plot(which = 5)

model_2 %>%
  plot(which = 4)

# Criteria: Cooks distance > 4/N(150)=0.026
data1 <- data1 %>% 
  slice(-47, -74, -86)

# Corrected model 1 and 2 

model_1_corrected <- lm(pain ~ sex + age, data = data1)
model_2_corrected <- lm(pain ~ sex + age + STAI_trait + pain_cat + mindfulness + cortisol_serum, data = data1)
model_1_corrected
model_2_corrected

# Checking model_2_corrected to see if the assumptions of normality, linearity, homoscedasticty and that there is no excess multicollinearity

# 1. Assumption of normality 

  # QQ-plot

model_2_corrected %>%
  plot(which = 2)

  # Histogram 

residuals_model_2_corrected = enframe(residuals(model_2_corrected))
residuals_model_2_corrected %>%
  ggplot() + aes(x = value) + geom_histogram()

  # Skew and kurtosis

describe(residuals(model_2_corrected))

# Skew and kurtosis inside the range of -1 to 1 (skew=-0.14, kurtosis=0)
# The assumption of normality not violated 

# 2. Assumption of linearity

model_2_corrected %>%
  residualPlots()

# No significance in the results, linearity-assumption not violated 

# 3. Assumption of homoscedasticty

model_2_corrected %>%
  plot(which = 3)

bptest(model_2_corrected)
ncvTest(model_2_corrected)

  # No significance in Breush-Pagan nor NCV test - assumption not violated

# 4. 

model_2_corrected %>%
  vif()

  # VIF-value not above 3 

# Adjusted R Squared of model 1 and model 2

summary(model_1_corrected)$adj.r.squared
summary(model_2_corrected)$adj.r.squared

# Compare of model fit 

AIC(model_1_corrected)
AIC(model_2_corrected)

# Summary
summary(model_1_corrected)
summary(model_2_corrected)

coef(model_1_corrected)
coef(model_2_corrected)

confint(model_1_corrected)
confint(model_2_corrected)

# Costum function for model regression coefficients  
coef_table = function(model) {
  require(lm.beta)
  mod_sum = summary(model)
  mod_sum_p_values = as.character(round(mod_sum$coefficients[,
                                                             4], 3))
  mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values !=
                     "1"] = substr(mod_sum_p_values[mod_sum_p_values != "0" &
                                                      mod_sum_p_values != "1"], 2, nchar(mod_sum_p_values[mod_sum_p_values !=
                                                                                                            "0" & mod_sum_p_values != "1"]))
  mod_sum_p_values[mod_sum_p_values == "0"] = "<.001"
  mod_sum_table = cbind(as.data.frame(round(cbind(coef(model),
                                                  confint(model), c(0, lm.beta(model)$standardized.coefficients[c(2:length(model$coefficients))])),
                                            2)), mod_sum_p_values)
  names(mod_sum_table) = c("b", "95%CI lb", "95%CI ub", "Std.Beta",
                           "p-value")
  mod_sum_table["(Intercept)", "Std.Beta"] = "0"
  return(mod_sum_table)
}

# Table models regression coefficients 

coef_table(model_1_corrected)
coef_table(model_2_corrected)

# Anova to compare the models in residual errors and degrees of freedom

anova(model_1_corrected, model_2_corrected)
AIC(model_1_corrected, model_2_corrected)


##################################################
##### ASSIGNMENT PART 2: BACKWARD REGRESSION #####
##################################################


# Creating new model with backward regression
theory_based_model <- lm(pain ~ sex + age + STAI_trait + pain_cat + cortisol_serum + mindfulness + weight + IQ + household_income, data=data1)
model_3_backward <- step(model_3, direction = "backward")
backward_model <- lm(pain ~ age + pain_cat + cortisol_serum + mindfulness, data=data1)

# Summary of data 
summary(theory_based_model)
summary(backward_model)

coef_table(theory_based_model)
coef_table(backward_model)

# AIC, R squared and ANOVA
AIC(theory_based_model)
AIC(backward_model)

summary(theory_based_model)$adj.r.squared
summary(backward_model)$adj.r.squared

anova(theory_based_model, backward_model)

# Prediction performance 

# Calculation of predicted values 
pred_test_theory_1 <- predict(theory_based_model, data1)
pred_test_backward_1 <- predict(backward_model, data1)

pred_test_theory_1
pred_test_backward_1

summary(pred_test_theory_1)
summary(pred_test_backward_1)

# Sums of squares residuals 
RSS_test_theory_1 = sum((data1[, "pain"] - pred_test_theory_1)^2)
RSS_test_backward_1 = sum((data1[, "pain"] - pred_test_backward_1)^2)

RSS_test_theory_1
RSS_test_backward_1

# Testing of model on new data 
data2 <- read.csv("https://raw.githubusercontent.com/kekecsz/PSYP14-Advanced-Scientific-Methods/main/Home%20assignment/home_sample_2.csv")
view(data2)
data2

# Checking dataset for coding errors and outliers 
ggplot(data2, aes(y=pain))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(data2, aes(y=age))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(data2, aes(y=sex))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(data2, aes(y=STAI_trait))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(data2, aes(y=pain_cat))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(data2, aes(y=mindfulness))+stat_boxplot(geom="errorbar")+geom_boxplot()
ggplot(data2, aes(y=cortisol_serum))+stat_boxplot(geom="errorbar")+geom_boxplot()

# Calculation of predicted values 
pred_test_theory_2 <- predict(theory_based_model, data2)
pred_test_backward_2 <- predict(backward_model, data2)

pred_test_theory_2
pred_test_backward_2

summary(pred_test_theory_2)
summary(pred_test_backward_2)

# Sums of squares residuals 
RSS_test_theory_2 = sum((data2[, "pain"] - pred_test_theory_2)^2)
RSS_test_backward_2 = sum((data2[, "pain"] - pred_test_backward_2)^2)

RSS_test_theory_2
RSS_test_backward_2

# Backwards model more errors than theory-based model 


################################################
##### ASSIGMENT PART 3: LINEAR MIXED MODEL #####
################################################


# Collection of data 
data3 <- read.csv("https://raw.githubusercontent.com/kekecsz/PSYP14-Advanced-Scientific-Methods/main/Home%20assignment/home_sample_3.csv")
data4 <- read.csv("https://raw.githubusercontent.com/kekecsz/PSYP14-Advanced-Scientific-Methods/main/Home%20assignment/home_sample_4.csv")
view(data3)
view(data4)

# Costum function 
stdCoef.merMod <- function(object) {
  sdy <- sd(getME(object, "y"))
  sdx <- apply(getME(object, "X"), 2, sd)
  sc <- fixef(object) * sdx/sdy
  se.fixef <- coef(summary(object))[, "Std. Error"]
  se <- se.fixef * sdx/sdy
  return(data.frame(stdcoef = sc, stdse = se))
}

# Converting hospital into factor 
data3 = data3 %>%
  mutate(hospital = factor(hospital))

# Check dataset 
----

# Assigning data to two different dataframes
data_pain_int <- data3
view(data_pain_int)

data_pain_slope <- data3
view(data_pain_slope)

# Data visualization 
data3 %>%
  ggplot() + aes(y = pain, x = sex) + geom_point(aes(color = hospital), size = 4) + geom_smooth(method = "lm", se = F)

data3 %>%
  ggplot() + aes(y = pain, x = age) + geom_point(aes(color = hospital), size = 4) + geom_smooth(method = "lm", se = F)

data3 %>%
  ggplot() + aes(y = pain, x = STAI_trait) + geom_point(aes(color = hospital), size = 4) + geom_smooth(method = "lm", se = F)

data3 %>%
  ggplot() + aes(y = pain, x = pain_cat) + geom_point(aes(color = hospital), size = 4) + geom_smooth(method = "lm", se = F)

data3 %>%
  ggplot() + aes(y = pain, x = cortisol_serum) + geom_point(aes(color = hospital), size = 4) + geom_smooth(method = "lm", se = F)

data3 %>%
  ggplot() + aes(y = pain, x = mindfulness) + geom_point(aes(color = hospital), size = 4) + geom_smooth(method = "lm", se = F)

# Assigning data to two different dataframes
data_pain_int <- data3
view(data_pain_int)

data_pain_slope <- data3
view(data_pain_slope)

# 1. Fixed effect model 
model_fixed_int = lm(pain ~ sex + age + STAI_trait + pain_cat + cortisol_serum + mindfulness + hospital, data = data_int_pred)
model_fixed_int

# 2. Random intercept model 
model_rnd_int = lmer(pain ~ sex + age + STAI_trait + pain_cat + cortisol_serum + mindfulness + (1|hospital), data = data_pain_int)
model_rnd_int

data_int_pred <- data_int %>%
  mutate(pred_int = predict(model_3_rnd_int))

data_int_pred

# Plotting data_int_pred
plot_1 <- data_int_pred %>%
  ggplot() + aes(y = pain, x = sex, group = hospital) +
  geom_point(aes(color = hospital), size = 4) + geom_line(color = "red", aes(y = pred_int, x = age)) + facet_wrap(~hospital, ncol = 2)

plot_2 <- data_int_pred %>%
  ggplot() + aes(y = pain, x = age, group = hospital) +
  geom_point(aes(color = hospital), size = 4) + geom_line(color = "red", aes(y = pred_int, x = age)) + facet_wrap(~hospital, ncol = 2)

plot_3 <- data_int_pred %>%
  ggplot() + aes(y = pain, x = STAI_trait, group = hospital) +
  geom_point(aes(color = hospital), size = 4) + geom_line(color = "red", aes(y = pred_int, x = age)) + facet_wrap(~hospital, ncol = 2)

plot_4 <- data_int_pred %>%
  ggplot() + aes(y = pain, x = pain_cat, group = hospital) +
  geom_point(aes(color = hospital), size = 4) + geom_line(color = "red", aes(y = pred_int, x = age)) + facet_wrap(~hospital, ncol = 2)

plot_5 <- data_int_pred %>%
  ggplot() + aes(y = pain, x = cortisol_serum, group = hospital) +
  geom_point(aes(color = hospital), size = 4) + geom_line(color = "red", aes(y = pred_int, x = age)) + facet_wrap(~hospital, ncol = 2)

plot_6 <- data_int_pred %>%
  ggplot() + aes(y = pain, x = mindfulness, group = hospital) +
  geom_point(aes(color = hospital), size = 4) + geom_line(color = "red", aes(y = pred_int, x = age)) + facet_wrap(~hospital, ncol = 2)

grid.arrange(plot_1, plot_2, nrow = 1)
grid.arrange(plot_3, plot_4, nrow = 1)
grid.arrange(plot_5, plot_6, nrow = 1)

# Exploring clustering in data 
plot_1_int <- data_int_pred %>%
  ggplot() + aes(y = pain, x = sex, color = hospital) +
  geom_point(size = 4) + geom_smooth(method = "lm", se = F,
                                     fullrange = TRUE)

plot_2_int <- data_int_pred %>%
  ggplot() + aes(y = pain, x = age, color = hospital) +
  geom_point(size = 4) + geom_smooth(method = "lm", se = F,
                                     fullrange = TRUE)

plot_3_int <- data_int_pred %>%
  ggplot() + aes(y = pain, x = STAI_trait, color = hospital) +
  geom_point(size = 4) + geom_smooth(method = "lm", se = F,
                                     fullrange = TRUE)

plot_4_int <- data_int_pred %>%
  ggplot() + aes(y = pain, x = pain_cat, color = hospital) +
  geom_point(size = 4) + geom_smooth(method = "lm", se = F,
                                     fullrange = TRUE)

plot_5_int <- data_int_pred %>%
  ggplot() + aes(y = pain, x = cortisol_serum, color = hospital) +
  geom_point(size = 4) + geom_smooth(method = "lm", se = F,
                                     fullrange = TRUE)

plot_6_int <- data_int_pred %>%
  ggplot() + aes(y = pain, x = mindfulness, color = hospital) +
  geom_point(size = 4) + geom_smooth(method = "lm", se = F,
                                     fullrange = TRUE)

grid.arrange(plot_1_int, plot_2_int, nrow = 1)
grid.arrange(plot_3_int, plot_4_int, nrow = 1)
grid.arrange(plot_5_int, plot_6_int, nrow = 1)

# Plotting random intercept
plot_1_int + xlim(-1, 50) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
plot_2_int + xlim(-1, 50) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
plot_3_int + xlim(-1, 50) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
plot_4_int + xlim(-1, 50) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
plot_5_int + xlim(-1, 10) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
plot_6_int + xlim(-1, 10) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

# Plotting slope 
slope_plot_1 = data_slope %>%
  ggplot() + aes(y = pain, x = sex, color = hospital) +
  geom_point(size = 4) + geom_smooth(method = "lm", se = F,
                                     fullrange = TRUE) + xlim(-1, 50) + geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

slope_plot_2 = data_slope %>%
  ggplot() + aes(y = pain, x = age, color = hospital) +
  geom_point(size = 4) + geom_smooth(method = "lm", se = F,
                                     fullrange = TRUE) + xlim(-1, 50) + geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

slope_plot_3 = data_slope %>%
  ggplot() + aes(y = pain, x = STAI_trait, color = hospital) +
  geom_point(size = 4) + geom_smooth(method = "lm", se = F,
                                     fullrange = TRUE) + xlim(-1, 50) + geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

slope_plot_4 = data_slope %>%
  ggplot() + aes(y = pain, x = pain_cat, color = hospital) +
  geom_point(size = 4) + geom_smooth(method = "lm", se = F,
                                     fullrange = TRUE) + xlim(-1, 50) + geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

slope_plot_5 = data_slope %>%
  ggplot() + aes(y = pain, x = cortisol_serum, color = hospital) +
  geom_point(size = 4) + geom_smooth(method = "lm", se = F,
                                     fullrange = TRUE) + xlim(-1, 10) + geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

slope_plot_6 = data_slope %>%
  ggplot() + aes(y = pain, x = mindfulness, color = hospital) +
  geom_point(size = 4) + geom_smooth(method = "lm", se = F,
                                     fullrange = TRUE) + xlim(-1, 10) + geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)

grid.arrange(slope_plot_1, slope_plot_2, nrow = 1)
grid.arrange(slope_plot_3, slope_plot_4, nrow = 1)
grid.arrange(slope_plot_5, slope_plot_6, nrow = 1)

# 3. Random slope model 
mod_rnd_slope = lmer(pain ~ sex + age + STAI_trait + pain_cat + cortisol_serum + mindfulness + (sex | hospital),
                     data = data_slope)
mod_rnd_slope

# Calculating residuals 
sum(residuals(model_3)^2)
sum(residuals(model_3_rnd_int)^2)

# Statistics of models 
coef_table(model_3)
stdCoef.merMod(model_3_rnd_int)

# cAIC comparison of models 
cAIC(model_3)$caic
cAIC(model_3_rnd_int)$caic

# ANOVA 
anova(model_3, model_3_rnd_int)

