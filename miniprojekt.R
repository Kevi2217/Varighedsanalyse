library(dplyr)
library(tidyr)
library(lubridate)
library(openxlsx)
library(ggplot2)
library(leaflet)
library(ggmap)
library(mapview)
library(sp)
library(mapdeck)
library(RColorBrewer)
library(ggcorrplot)
library(gridExtra)
library(survival)
library(survminer)
library(forestmodel)

load("~/Desktop/varighedsanalyse miniprojekt/melanoma_data.RData")

##################################### OPG 2 ####################################
# Function to plot continuous variables (boxplot & histogram)
plot_continuous <- function(data, col_name) {
  col_data <- data[[col_name]]
  
  # Create a combined boxplot and histogram layout
  p1 <- ggplot(data, aes_string(x = col_name)) + 
    geom_boxplot() + 
    ggtitle(paste("Boxplot of", col_name))
  
  p2 <- ggplot(data, aes_string(x = col_name)) + 
    geom_histogram(binwidth = (max(melanoma30[[col_name]], na.rm = TRUE) - min(melanoma30[[col_name]], na.rm = TRUE)) / 30,
                   fill = "blue", alpha = 0.7) + 
    ggtitle(paste("Histogram of", col_name))
  
  # Arrange plots side by side
  gridExtra::grid.arrange(p1, p2, ncol = 2)
}

# Function to plot categorical variables (bar plot)
plot_categorical <- function(data, col_name) {
  ggplot(data, aes_string(x = col_name)) + 
    geom_bar(fill = "blue", alpha = 0.7) + 
    ggtitle(paste("Bar plot of", col_name)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for readability
}

# Main function to loop through specified column indices and make appropriate plots
plot_outliers_cont <- function(data, col_indices) {
  for (col_index in col_indices) {
    col_name <- colnames(data)[col_index]
      print(paste("Continuous column:", col_name))
      plot_continuous(data, col_name)
    }
  }


plot_outliers_cat <- function(data, col_indices) {
  for (col_index in col_indices) {
    col_name <- colnames(data)[col_index]
    print(paste("Categorical column:", col_name))
    plot_categorical(data, col_name)    
  }
}

plot_outliers_cont(melanoma30, c(8, 10, 11))
# plot_categorical(melanoma30, "status")
# plot_categorical(melanoma30, "dead")
plot_categorical(melanoma30, "ici")
plot_categorical(melanoma30, "epicell")
plot_categorical(melanoma30, "ulceration")
plot_categorical(melanoma30, "sex")
plot_categorical(melanoma30, "invas2")




##################################### OPG 3 ####################################

# Creating categories for thickness to turn it into a categorical value
melanoma30 <- melanoma30 %>% 
  dplyr::mutate(thick_class = ifelse(thickness <= 97, "very small",
                              ifelse(thickness > 97 & thickness <= 194, "small",
                              ifelse(thickness > 194 & thickness <= 356, "medium",
                              ifelse(thickness > 356, "large",
                                     "unknown"))))) %>% 
  dplyr::mutate(dead = ifelse(dead == "dead", 1, 0))

# Create a Surv-object for survival analysis
df_surv <- Surv(melanoma30$time, melanoma30$dead)

# Estimate Kaplan-Meier survival curves for the thickness categories
kaplan_fit <- survfit(df_surv ~ thick_class, data = melanoma30)


# Plot the survival curves
ggsurvplot(kaplan_fit, data = melanoma30, pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, legend.title = "Tumor Thickness",
           title = "Kaplan-Meier Curves for Tumor Thickness Categories",
           xlab = "Time", ylab = "Survival Probability")

# Perform a log-rank test to compare survival distributions
survdiff(df_surv ~ thick_class, data = melanoma30)
# log_rank_test suggests that there is a significant difference of the curves 

##################################### OPG 4 ####################################

cox_model <- coxph(df_surv ~ ici + epicell + ulceration + thickness +
                     sex + age + invas2,
                   data = melanoma30)
summary(cox_model)

cox_model_mr_thickness <- coxph(df_surv ~ ici + epicell + ulceration +
                     sex + age + invas2,
                   data = melanoma30)

cox_model_mr_age <- coxph(df_surv ~ ici + epicell + ulceration + thickness +
                                  sex + invas2,
                                data = melanoma30)
# Get martingale residuals
martingale_res_thickness <- residuals(cox_model_mr_thickness, type = "martingale")
martingale_res_age <- residuals(cox_model_mr_age, type = "martingale")

# Plot martingale residuals vs. age
ggplot(melanoma30, aes(x = age, y = martingale_res_age)) +
  geom_point() +
  geom_smooth(method = "loess", color = "red") +
  labs(title = "Martingale Residuals vs. Age")

# Plot martingale residuals vs. tumor thickness
ggplot(melanoma30, aes(x = thickness, y = martingale_res_thickness)) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue") +
  labs(title = "Martingale Residuals vs. Tumor Thickness")

# Plot of log(thickness) to martingale_res
ggplot(melanoma30, aes(x = log(thickness), y = martingale_res_thickness)) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue") +
  labs(title = "Martingale Residuals vs. log(Tumor Thickness)")

# Model 1: Log transformation for tumor thickness
cox_model_log <- coxph(df_surv ~ ici + epicell + ulceration + log(thickness) +
                         sex + age + invas2,
                       data = melanoma30)

# Model 2: Categorize age into quartiles
# melanoma30$age_cat <- cut(melanoma30$age, breaks = quantile(melanoma30$age,
#                                                             probs = seq(0, 1, 0.25), na.rm = TRUE),
#                   include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))
# 
# cox_model_age_cat <- coxph(df_surv ~ ici + epicell + ulceration + thickness +
#                              sex + age_cat + invas2, 
#                            data = melanoma30)

# According to AIC the cox_model_log is the best


AIC(cox_model, cox_model_log)


##################################### OPG 5 ####################################
# Checking Proportional Hazards Assumption
zph_test <- cox.zph(cox_model_log)
print(zph_test)
plot(zph_test)



# We check for any outliers in our fitted model so far cox_model_log2
# We find that none of the observations have a larger deviance than 3
deviance_res <- residuals(cox_model_log, type = "deviance")
which(abs(deviance_res) > 3)

plot(melanoma30$id, deviance_res)

melanoma30_reduced <- melanoma30[-2,]

df_surv_reduced <- Surv(melanoma30_reduced$time, melanoma30_reduced$dead)

cox_model_log_wo_out <- coxph(df_surv_reduced ~ ici + epicell + ulceration + log(thickness) +
                         sex + age + invas2,
                       data = melanoma30_reduced)

##################################### OPG 6 ####################################
# We test in three diferent ways to see if the null hypothesis of no effect of tumor thickness is true
# Creating the model without thickness
cox_wo_out_and_thickness_model_log <- coxph(df_surv_reduced ~ ici + epicell + ulceration +
                          sex + age + invas2,
                        data = melanoma30_reduced)
# a)
summary(cox_model_log_wo_out)

# b)
anova(cox_wo_out_and_thickness_model_log, cox_model_log_wo_out)

##################################### OPG 7 ####################################

new_patient <- data.frame(
  age = 57,
  ici = 2,
  epicell = "yes",
  sex = "male",
  ulceration = "yes",
  thickness = 750,
  invas2 = "Clark I-III"
)

# Estimate the baseline survival function using the Cox model
new_patient_fit <- survfit(cox_model_log_wo_out, newdata = new_patient)


# Predict the survival probability at time 8 (difference between 57 and 65 years)
surv_prob <- summary(new_patient_fit, times = 2920)

print(surv_prob$surv)

##################################### OPG 8 ####################################
# We compare the standard cox_model to the exponential and Weibull models (parametric models)
cox_base_hazard <- basehaz(cox_model)

# Plot the cumulative baseline hazard
plot(cox_base_hazard$time, cox_base_hazard$hazard,
     type = "l", xlab = "Time", ylab = "Cumulative Hazard", 
     main = "Cumulative Baseline Hazard from Cox Model")

# Check for Weibull (log-log plot) If straight line then Weibull could be good fit for data
plot(log(cox_base_hazard$time), log(cox_base_hazard$hazard),
     type = "l", xlab = "log(Time)", ylab = "log(Cumulative Hazard)",
     main = "Log-Log Plot (Weibull Check)")

# We fit the exponential and Weibull to compare their AIC's to the cox_model
exp_model <- survreg(df_surv ~ age + ici + epicell + ulceration +
                       sex + thickness + invas2,
                     data = melanoma30, dist = "exponential")
summary(exp_model)

weibull_model <- survreg(df_surv ~ age + ici + epicell + ulceration +
                           sex + thickness + invas2,
                         data = melanoma30, dist = "weibull")
summary(weibull_model)


# AIC below shows that the cox_model is the best fit for the melanoma30 data
AIC(cox_model, exp_model, weibull_model)





