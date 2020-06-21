library(readr)
library(pscl)
library(car)
library(ggplot2)
library(dplyr)
library(tidyr)
library(GGally)
library(lmtest)
library(pROC)

# Replace with your directory
setwd("../")

attach(mtcars)

loadData = function(path){
  d = read_csv(path)
  return(d)
}

# Functions for IPW
sweep = function(range, data, transit_col, fm){
  output = vector()
  max_weights = vector()
  for (val in range){
    data = data %>% mutate(treatment=(pmin(transit_col)<=val))
    glm.1 = glm(formula=as.formula(fm), family="binomial",data=data)
    pred.used = predict(glm.1, type = "response")
    weights.used = data$treatment/pred.used + (1 - data$treatment)/(1 - pred.used)
    trt_eff = mean(data$average_ridership * data$treatment/pred.used) - mean(data$average_ridership * (1 - data$treatment)/(1 - pred.used))
    output = append(output, trt_eff)
    max_weights = append(max_weights, max(weights.used))
  }
  return(list("output"=output,"max.weights"=max_weights))
}

sweep.boot = function(range, data, transit_col, fm){
  n.boot = 1000
  output = matrix(NA, nrow = length(range), ncol = n.boot)
  for(j in 1:n.boot){
    index = sample(1:nrow(data), nrow(data), replace = TRUE)
    data.boot = data[index, ]
    for (i in 1:length(range)){
      val = range[i]
      data.boot = data.boot %>% mutate(treatment=(pmin(transit_col)<=val))
      glm.1 = glm(formula = as.formula(fm), family="binomial",data=data.boot)
      pred.used = predict(glm.1, type = "response")
      weights.used = data.boot$treatment/pred.used + (1 - data.boot$treatment)/(1 - pred.used)
      trt_eff = mean(data.boot$average_ridership * data.boot$treatment/pred.used) - mean(data.boot$average_ridership * (1 - data.boot$treatment)/(1 - pred.used))
      output[i, j] = trt_eff
    }
  }
  return(output)
}

# load in the imputed data
d_data = loadData("Data/dc_summer_bike_station_data.csv")
colnames(d_data)

sub_data = d_data %>% select(average_ridership,percent_minority,percent_employed,total_pop,income_per_cap,
                             percent_poverty,percent_transit,percent_walk,total_num_jobs,average_pop,average_num_jobs,
                             mean_commute,bus_dist,rail_dist,professional,service,
                             office,construction,production)

sum(!complete.cases(sub_data))

# Plot histograms for all the variables
sub_data %>% gather() %>% ggplot(aes(value)) + facet_wrap(~ key, scales = "free") + geom_density()
# average_num_jobs, construction, bus_dist, average_ridership, production and total_num_jobs is highly skewed to the right
# average_pop, percent_walk, percent_poverty is skewed to the right
# professional is slightly skewed to the left

# Plot logged distributions
is.finite.data.frame <- function(obj){
  sapply(obj,FUN = function(x) sum(is.finite(x)))
}
full_logs = log(sub_data)
is.finite.data.frame(full_logs)
full_logs = full_logs[is.finite(rowSums(full_logs)),] 
nrow(sub_data) - nrow(full_logs) # Note removing 161 rows to keep finite
full_logs %>% gather() %>% ggplot(aes(value)) + facet_wrap(~ key, scales = "free") + geom_density()
# Distributions improve for average_num_jobs, average_pop, bus_dist, average_ridership, bus_dist, construction,
# percent_walk, total_num_jobs, percent_poverty
# average_ridership is highly skewed -> take log
sub_data =sub_data %>% mutate(log_average_ridership=log(average_ridership+1))
colnames(sub_data)

# Plot relationship between variables and average_ridership
sub_data %>%
  gather(-log_average_ridership, key = "var", value = "value") %>%
  ggplot(aes(x = value, y = log_average_ridership)) +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
# with average_ridership: everything looks pretty much like a horizontal line at y=0 because average_ridership is so skewed
# with log_average_ridership: percent_walk looks non-linear

# Plot relationship between logged variables and average_ridership
full_logs %>%
  gather(-average_ridership, key = "var", value = "value") %>%
  ggplot(aes(x = value, y = average_ridership)) +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
# percent_walk, income_per_cap, percent_employed seem to have a slightly positive relationship
# construction, mean_commute, office, percent_minority seem to have a slightly negative relationship 

# Correlation matrix
ggcorr(sub_data,label=TRUE)
# with log_average_ridership: rail_dist (-0.5), mean_commute (-0.5), percent_walk (0.6)
# average_num_jobs (0), percent_poverty (0)

# Correlation matrix of logged variables
ggcorr(full_logs,label=TRUE)
# percent_walk (0.6), mean_commute (-0.5)

##################################### RAIL ###############################################
# Rail distance - start by building logistic regression for 50th percentile
quant_50 = quantile(sub_data$rail_dist, 0.50)
quant_50
sub_data = sub_data %>% mutate(treatment=(pmin(rail_dist)<=quant_50))

sub_data %>%
  gather(-treatment, key = "var", value = "value") %>%
  ggplot(aes(x = treatment, y = value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
# Largest differences appear in percent_transit, total_num_jobs, log_average_ridership, percent_walk, percent_employed, rail_dist
# Smallest differences appear in bus_dist, office, percent_minority, professional, construction, service, total_pop, production, percent_poverty

# BASELINE MODEL
glm.1 = glm(treatment~.-rail_dist-bus_dist-average_ridership-log_average_ridership, family="binomial",data=sub_data)
summary(glm.1) # coeffs for percent_minority, total_pop, percent_transit, percent_walk, total_num_jobs are significant, 
               # AIC: 548.79
glm.1.probs = predict(glm.1,type="response")
glm.1.pred = (glm.1.probs > 0.5)
table(glm.1.pred, sub_data$treatment) # pretty good
vif(glm.1) # vars with value over 5: professional, service, office, construction, production
#           remove: service, office, construction, production

# Subset of variables
glm.2 = glm(treatment~professional+average_pop+percent_minority+total_num_jobs+income_per_cap+average_num_jobs+total_pop+percent_employed+percent_poverty+percent_transit+percent_walk+mean_commute, family="binomial",data=sub_data)
summary(glm.2) # professional, percent_minority, total_pop, total_num_jobs, percent_transit, and percent_walk are significant
               # AIC: 549.62
vif(glm.2) # all below 5
glm.2.probs = predict(glm.2,type="response")
glm.2.pred = (glm.2.probs > 0.5)
table(glm.2.pred, sub_data$treatment) # pretty accurate
lrtest(glm.1,glm.2) # borderline p-value (0.065)
plot(glm.2) 
roc.2 = roc(sub_data$treatment,glm.2.probs)
auc(roc.2) # Area under the curve: 0.8187, not bad (maybe could be better)
ggroc(data=roc.2)

# Testing logs (logged average_num_jobs, total_num_jobs, and average_pop based on EDA above)
# After some trial and error, I also removed total_num_jobs because it had a high VIF with average_num_jobs
# and resulted in a better model when removed.
glm.3 = glm(treatment~log(average_pop)+total_pop+log(average_num_jobs)+income_per_cap+percent_employed+percent_poverty+percent_transit+percent_walk+mean_commute+professional, family="binomial",data=sub_data)
summary(glm.3) # log(average_pop), log(average_num_jobs), percent_transit, and percent_walk are significant
               # AIC: 551.08
vif(glm.3) # all below 5
glm.3.probs = predict(glm.3,type="response")
glm.3.pred = (glm.3.probs > 0.5)
table(glm.3.pred, sub_data$treatment) # pretty accurate
plot(glm.3) # there doesn't seem to be any high leverage points
roc.3 = roc(sub_data$treatment,glm.3.probs)
auc(roc.3) # Area under the curve: 0.8211, a improvement
ggroc(data=roc.3)
# will keep the logs

# Testing polynomial fit for most significant variables percent_transit, percent_walk
# I experimented with different combinations of poly fit for those two variables
glm.4 = glm(treatment~poly(percent_transit,degree=2)+percent_walk+log(average_pop)+total_pop+log(average_num_jobs)+income_per_cap+percent_employed+percent_poverty+mean_commute+professional, family="binomial",data=sub_data)
summary(glm.4) # log(average_pop), log(average_num_jobs), percent_transit (both degrees), and percent_walk are significant
               # AIC: 547.17
vif(glm.4) # all below 5
glm.4.probs = predict(glm.4,type="response")
glm.4.pred = (glm.4.probs > 0.5)
table(glm.4.pred, sub_data$treatment) # pretty accurate
plot(glm.4) # there doesn't seem to be any high leverage points
roc.4 = roc(sub_data$treatment,glm.4.probs)
auc(roc.4) # Area under the curve: 0.825, slight improvement
ggroc(data=roc.4)
# I will keep the poly fit

# IPW resuslts look the same with and without the poly fit, so I will keep it in the model
formula_rail = "treatment~poly(percent_transit,degree=2)+percent_walk+log(average_pop)+total_pop+log(average_num_jobs)+income_per_cap+percent_employed+percent_poverty+mean_commute+professional"
range.rail = quantile(sub_data$rail_dist, probs = seq(0.25, 0.75, by= 0.01)) # decile
output.rail = sweep.boot(range.rail, sub_data, sub_data$rail_dist, formula_rail)
results.rail = sweep(range.rail, sub_data, sub_data$rail_dist, formula_rail)
output.rail.mean = results.rail$output
rail.max.weights = results.rail$max.weights
upper.rail.ci = output.rail.mean + 1.96 * apply(output.rail,1, sd)
lower.rail.ci = output.rail.mean - 1.96 * apply(output.rail,1, sd)

par(mfrow=c(1,2), cex=1.5, lwd=1, bg=NA)
plot(range.rail, output.rail.mean, type="l", xlab='Percentile (km)', ylab='Treatment effect (trips)', ylim = c(min(lower.rail.ci), max(upper.rail.ci)), lwd=2)
lines(range.rail, upper.rail.ci, lwd=2, lty=3)
lines(range.rail, lower.rail.ci, lwd=2, lty=3)
plot(range.rail,rail.max.weights, type="l", xlab='Percentile (km)', ylab = 'Max IPW Weight', lwd=2)
par(mfrow=c(1,1), cex=1.5, lwd=1, bg=NA)

#write.csv(range.rail, 'Data/dc_rail_IPW_range.csv')
#write.csv(output.rail, 'Data/dc_rail_IPW_sweep.csv')
#write.csv(results.rail, 'Data/dc_rail_IPW_res.csv')
#write.csv(upper.rail.ci, 'Data/dc_IPW_upper.csv')
#write.csv(lower.rail.ci, 'Data/dc_rail_IPW_lower.csv')
