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

d_data = loadData("Data/philly_summer_bike_station_data.csv")
colnames(d_data)

# Get all the variables that we are intrested in seeing the distributions of a potentially include
# in our models
sub_data = d_data %>% select(average_ridership,percent_minority,percent_employed,total_pop,income_per_cap,
                             percent_poverty,percent_transit,percent_walk,total_num_jobs,average_pop,average_num_jobs,
                             mean_commute,bus_dist,rail_dist,trolley_dist,high_speed_dist,professional,service,
                             office,construction,production)

sum(!complete.cases(sub_data)) # all samples are complete
# Plot histograms for all the variables
sub_data %>% gather() %>% ggplot(aes(value)) + facet_wrap(~ key, scales = "free") + geom_density()
# average_num_jobs and total_num_jobs is highly skewed to the right
# average_pop, bus_dist, construction is skewed to the right

# Plot logged distributions
is.finite.data.frame <- function(obj){
  sapply(obj,FUN = function(x) sum(is.finite(x)))
}
full_logs = log(sub_data)
is.finite.data.frame(full_logs)
full_logs = full_logs[is.finite(rowSums(full_logs)),] 
nrow(sub_data) - nrow(full_logs) # Note removing 19 rows to keep finite
full_logs %>% gather() %>% ggplot(aes(value)) + facet_wrap(~ key, scales = "free") + geom_density()
# Distributions improve for average_num_jobs, average_pop, bus_dist, construction, production, percent_poverty

hist(sub_data$average_ridership) # slightly skewed - might not need to take log

# Plot relationship between variables and average_ridership
sub_data %>%
  gather(-average_ridership, key = "var", value = "value") %>%
  ggplot(aes(x = value, y = average_ridership)) +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
# Not many super obvious associations
# percent_employed, percent_walk, and income_per_cap seem to have a positive relationship with average_ridership
# professional and mean_commute could have a non linear relationship with average_ridership
# percent_minority, and service seem to have a negative relationship with average_ridership

# Plot relationship between logged variables and average_ridership
full_logs %>%
  gather(-average_ridership, key = "var", value = "value") %>%
  ggplot(aes(x = value, y = average_ridership)) +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
# Relationship between logged average_ridership and logged professional, percent_employed, production, mean_commute, 
# percent_transit, service, percent_walk appears stronger 

# Correlation matrix
ggcorr(sub_data,label=TRUE)
# 0.6 between average_ridership and professional
# 0.5 between average_ridership and income_per_cap
# -0.6 between average_ridership and service
# -0.5 between average_ridership and production, rail_dist, percent_minority
# 0 between average_ridership and trolley_dist
# 0.1 between average_ridership and total_pop and total_num_jobs ---> use average_pop and average_num_jobs
# -0.2 between average_ridership and bus_dist

# Correlation matrix of logged variables
ggcorr(full_logs,label=TRUE)
# 0.6 between average_ridership and professional
# 0.5 between average_ridership and percent_walk, income_per_cap
# -0.6 between average_ridership and service
# -0.5 between average_ridership and production, percent_transit, percent_minority
# 0.1 between average_ridership and total_pop, trolley_dist
# -0.2 between average_ridership and bus_dist, high_speed_dist

##################################### RAIL ###############################################

# Rail distance - start by building logistic regression for 50th percentile
quant_50 = quantile(sub_data$rail_dist, 0.50)
sub_data = sub_data %>% mutate(treatment=(pmin(rail_dist)<=quant_50))

sub_data %>%
  gather(-treatment, key = "var", value = "value") %>%
  ggplot(aes(x = treatment, y = value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
# Largest differences appear in average_ridership, income_per_cap, percent_walk, professional, service
# Smallest differences appear in bus_dist, percent_poverty, total_pop, trolley_dist

# BASELINE MODEL
glm.1 = glm(treatment~.-rail_dist-trolley_dist-bus_dist-high_speed_dist-average_ridership, family="binomial",data=sub_data)
summary(glm.1) # only coefficient for percent_walk is significant, AIC: 101.52
glm.1.probs = predict(glm.1,type="response")
glm.1.pred = (glm.1.probs > 0.5)
table(glm.1.pred, sub_data$treatment)
vif(glm.1) # percent_minority, percent_employed, income_per_cap, total_num_jobs, professional, service, office
           # construction, production have values above 5 --> remove total_num_jobs, serivce, office, construction, production, income_per_cap

# Subset of variables
glm.2 = glm(treatment~average_pop+average_num_jobs+total_pop+percent_minority+percent_employed+percent_poverty+percent_transit+percent_walk+mean_commute+professional, family="binomial",data=sub_data)
summary(glm.2) # average_num_jobs, total_pop, and percent_walk are significant, AIC:96.709
vif(glm.2) # professional still has value above 6, so remove
glm.2.probs = predict(glm.2,type="response")
glm.2.pred = (glm.2.probs > 0.5)
table(glm.2.pred, sub_data$treatment)
lrtest(glm.1,glm.2) # high p-value so fail to reject the null hypothesis. larger model is not necessarily better

# Model after removing professional
glm.3 = glm(treatment~average_pop+average_num_jobs+total_pop+percent_minority+percent_employed+percent_poverty+percent_transit+percent_walk+mean_commute, family="binomial",data=sub_data)
summary(glm.3) # average_num_jobs, total_pop, percent_employed, and percent_walk are significant, AIC:97.17
vif(glm.3) # they are all below 5
glm.3.probs = predict(glm.3,type="response")
glm.3.pred = (glm.3.probs > 0.5)
table(glm.3.pred, sub_data$treatment)
lrtest(glm.2,glm.3) # high p-value so fail to reject the null hypothesis. larger model is not necessarily better
plot(glm.3) # residuals are mainly centered around 0 there seems to be a bit of a > trend which may indicate a polynomial fit
roc.3 = roc(sub_data$treatment,glm.3.probs)
auc(roc.3) # Area under the curve: 0.9451 which is really good
ggroc(data=roc.3)

# Testing logs (logged average_num_jobs and average_pop based on EDA above)
glm.4 = glm(treatment~log(average_pop)+log(average_num_jobs)+total_pop+percent_minority+percent_employed+percent_poverty+percent_transit+percent_walk+mean_commute, family="binomial",data=sub_data)
summary(glm.4) # average_num_jobs, total_pop, percent_employed, and percent_walk are significant, AIC: 95.277 (significance increased for average_num_jobs and percent_walk)
vif(glm.4) # they are all below 5
glm.4.probs = predict(glm.4,type="response")
glm.4.pred = (glm.4.probs > 0.5)
table(glm.4.pred, sub_data$treatment)
plot(glm.4) # residuals are mainly centered around 0, 111 looks like an outlier
roc.4 = roc(sub_data$treatment,glm.4.probs)
auc(roc.4) # Area under the curve: 0.9415
ggroc(data=roc.4)
# Log didn't help all that much

# Testing polynomial fit for most significant variables average_num_jobs, total_pop, and percent_walk
glm.5 = glm(treatment~average_pop+poly(average_num_jobs,degree=2)+poly(total_pop,degree=2)+percent_minority+percent_employed+percent_poverty+percent_transit+poly(percent_walk,degree=2)+mean_commute, family="binomial",data=sub_data)
summary(glm.5) # total_pop and percent_walk are significant, AIC: 99.126
vif(glm.5) # they are all below 5
glm.5.probs = predict(glm.5,type="response")
glm.5.pred = (glm.5.probs > 0.5)
table(glm.5.pred, sub_data$treatment)
plot(glm.5) # residuals are mainly centered around 0 trend seems to still be there
roc.5 = roc(sub_data$treatment,glm.5.probs)
auc(roc.5) # Area under the curve: 0.9406
ggroc(data=roc.5)
# Polynomials didn't help all that much

formula_rail = "treatment~average_pop+average_num_jobs+total_pop+percent_minority+percent_employed+percent_poverty+percent_transit+percent_walk+mean_commute"
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

#write.csv(range.rail, 'Data/philly_rail_IPW_range.csv')
#write.csv(output.rail, 'Data/philly_rail_IPW_sweep.csv')
#write.csv(results.rail, 'Data/philly_rail_IPW_res.csv')
#write.csv(upper.rail.ci, 'Data/philly_rail_IPW_upper.csv')
#write.csv(lower.rail.ci, 'Data/philly_rail_IPW_lower.csv')

##################################### TROLLEY ###############################################
# Trolley distance - start by building logistic regression for 50th percentile
quant_50 = quantile(sub_data$trolley_dist, 0.50)
sub_data = sub_data %>% mutate(treatment=(pmin(trolley_dist)<=quant_50))

sub_data %>%
  gather(-treatment, key = "var", value = "value") %>%
  ggplot(aes(x = treatment, y = value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
# Largest differences appear in percent_walk, percent_employed, income_per_cap
# Smallest differences appear in average_num_jobs, average_pop, construction (most variables have similar means between treatment and control groups)

# BASELINE MODEL
glm.1 = glm(treatment~.-rail_dist-trolley_dist-bus_dist-high_speed_dist-average_ridership, family="binomial",data=sub_data)
summary(glm.1) # only coefficient for percent_poverty and average_pop are significant, AIC: 157.18
glm.1.probs = predict(glm.1,type="response")
glm.1.pred = (glm.1.probs > 0.5)
table(glm.1.pred, sub_data$treatment)
vif(glm.1) # percent_minority, percent_employed, income_per_cap, total_num_jobs, professional, service, office
# construction, production have values above 5 --> remove service, office, construction, production, income_per_cap, professional, percent_minority

# Subset of variables
glm.2 = glm(treatment~average_pop+average_num_jobs+total_pop+total_num_jobs+percent_employed+percent_poverty+percent_transit+percent_walk+mean_commute, family="binomial",data=sub_data)
summary(glm.2) # average_pop, total_pop, and percent_poverty are significant, AIC: 150.31
vif(glm.2) # all below 5
glm.2.probs = predict(glm.2,type="response")
glm.2.pred = (glm.2.probs > 0.5)
table(glm.2.pred, sub_data$treatment)
lrtest(glm.1,glm.2) # high p-value so fail to reject the null hypothesis. larger model is not necessarily better
plot(glm.2) # residuals are mainly centered around 0
roc.2 = roc(sub_data$treatment,glm.2.probs)
auc(roc.2) # Area under the curve: 0.8096 which is really good
ggroc(data=roc.2)

# Testing logs (logged average_num_jobs total_num_jobs and average_pop based on EDA above)
glm.3 = glm(treatment~log(average_pop)+log(average_num_jobs)+total_pop+log(total_num_jobs)+percent_employed+percent_poverty+percent_transit+percent_walk+mean_commute, family="binomial",data=sub_data)
summary(glm.3) # average_pop total_pop and total_num_jobs are significant, AIC: 153.7
vif(glm.3) # log(average_num_jobs) and log(total_num_jobs) both have values much larger than 5 --> remove log(total_num_jobs)
# After removing total_num_jobs
glm.4 = glm(treatment~log(average_pop)+log(average_num_jobs)+total_pop+percent_employed+percent_poverty+percent_transit+percent_walk+mean_commute, family="binomial",data=sub_data)
summary(glm.4) # average_pop and average_num_jobs are significant, AIC: 155.26
vif(glm.4) # values all below 5
glm.4.probs = predict(glm.4,type="response")
glm.4.pred = (glm.4.probs > 0.5)
table(glm.4.pred, sub_data$treatment)
plot(glm.4) # residuals are mainly centered around 0
roc.4 = roc(sub_data$treatment,glm.4.probs)
auc(roc.4) # Area under the curve: 0.7921
ggroc(data=roc.4)
# Log didn't help all that much

# Testing polynomial fit for most significant variables average_pop, total_pop, and percent_poverty
glm.5 = glm(treatment~poly(average_pop,degree=2)+average_num_jobs+poly(total_pop,degree=2)+total_num_jobs+percent_employed+poly(percent_poverty,degree=2)+percent_transit+percent_walk+mean_commute, family="binomial",data=sub_data)
summary(glm.5) # average_pop, total_pop and percent_poverty are significant, AIC: 154.56
vif(glm.5) # they are all below 5
glm.5.probs = predict(glm.5,type="response")
glm.5.pred = (glm.5.probs > 0.5)
table(glm.5.pred, sub_data$treatment) # there seems to be a lot of false positives and false negatives
plot(glm.5) # residuals don't look super good
roc.5 = roc(sub_data$treatment,glm.5.probs)
auc(roc.5) # Area under the curve: 0.6764
ggroc(data=roc.5)

# Polynomial fit for average_pop and total_pop
glm.6 = glm(treatment~poly(average_pop,degree=2)+average_num_jobs+poly(total_pop,degree=2)+total_num_jobs+percent_employed+percent_poverty+percent_transit+percent_walk+mean_commute, family="binomial",data=sub_data)
summary(glm.6) # percent_poverty and average_pop are signifciant, AIC: 153.7
vif(glm.6) # values all below 5
glm.6.probs = predict(glm.6,type="response")
glm.6.pred = (glm.6.probs > 0.5)
table(glm.6.pred, sub_data$treatment)
plot(glm.6)
roc.6 = roc(sub_data$treatment,glm.6.probs)
auc(roc.6) # Area under the curve: 0.821, which is a slight improvement
ggroc(data=roc.6)

formula_trolley = "treatment~poly(average_pop,degree=2)+average_num_jobs+poly(total_pop,degree=2)+total_num_jobs+percent_employed+percent_poverty+percent_transit+percent_walk+mean_commute"
range.trolley = quantile(sub_data$trolley_dist, probs = seq(0.25, 0.75, by= 0.01)) # decile
output.trolley = sweep.boot(range.trolley, sub_data, sub_data$trolley_dist, formula_trolley)
results.trolley = sweep(range.trolley, sub_data, sub_data$trolley_dist, formula_trolley)
output.trolley.mean = results.trolley$output
trolley.max.weights = results.trolley$max.weights
upper.trolley.ci = output.trolley.mean + 1.96 * apply(output.trolley,1, sd)
lower.trolley.ci = output.trolley.mean - 1.96 * apply(output.trolley,1, sd)

par(mfrow=c(1,2), cex=1.5, lwd=1, bg=NA)
plot(range.trolley, output.trolley.mean, type="l", xlab='Percentile (km)', ylab='Treatment effect (trips)', ylim = c(min(lower.trolley.ci), max(upper.trolley.ci)), lwd=2)
lines(range.trolley, upper.trolley.ci, lwd=2, lty=3)
lines(range.trolley, lower.trolley.ci, lwd=2, lty=3)
plot(range.trolley,trolley.max.weights, type="l", xlab='Percentile (km)', ylab = 'Max IPW Weight', lwd=2)

#write.csv(range.trolley, 'Data/philly_trolley_IPW_range.csv')
#write.csv(output.trolley, 'Data/philly_trolley_IPW_sweep.csv')
#write.csv(results.trolley, 'Data/philly_trolley_IPW_res.csv')
#write.csv(upper.trolley.ci, 'Data/philly_trolley_IPW_upper.csv')
#write.csv(lower.trolley.ci, 'Data/philly_trolley_IPW_lower.csv')

par(mfrow=c(1,1), cex=1.5, lwd=1, bg=NA)
##################################### HIGH SPEED ###############################################
# High speed distance - start by building logistic regression for 50th percentile
quant_50 = quantile(sub_data$high_speed_dist, 0.50)
sub_data = sub_data %>% mutate(treatment=(pmin(high_speed_dist)<=quant_50))

sub_data %>%
  gather(-treatment, key = "var", value = "value") %>%
  ggplot(aes(x = treatment, y = value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
# Largest differences appear in professional, percent_employed, income_per_cap
# Smallest differences appear in average_num_jobs, office, mean_commute, construction, total_pop

# BASELINE MODEL
glm.1 = glm(treatment~.-rail_dist-trolley_dist-bus_dist-high_speed_dist-average_ridership, family="binomial",data=sub_data)
summary(glm.1) # only coefficient for percent_employed, percent_poverty, average_num_jobs, mean_commute are significant, AIC: 160.34
glm.1.probs = predict(glm.1,type="response")
glm.1.pred = (glm.1.probs > 0.5)
table(glm.1.pred, sub_data$treatment)
vif(glm.1) # percent_minority, percent_employed, income_per_cap, percent_poverty, percent_walk, total_num_jobs
           # average_num_jobs, professional, service, office, construction, production have values above 5
           # remove professional, service, office, construction, production, total_num_jobs, percent_walk, percent_poverty
           # percent_minority

# Subset of variables (deleted different combinations of variables until inflation factors looked good and auc was highest)
glm.2 = glm(treatment ~ average_pop + average_num_jobs + total_pop + percent_employed + percent_transit + mean_commute, family="binomial",data=sub_data)
summary(glm.2) # average_num_jobs, percent_employed, percent_transit, mean_commute are significant, AIC: 156.11
vif(glm.2) # all below 3
glm.2.probs = predict(glm.2,type="response")
glm.2.pred = (glm.2.probs > 0.5)
table(glm.2.pred, sub_data$treatment)
lrtest(glm.1,glm.2) # high p-value so fail to reject the null hypothesis. larger model is not necessarily better
plot(glm.2) # residuals look pretty good
roc.2 = roc(sub_data$treatment,glm.2.probs)
auc(roc.2) # Area under the curve: 0.7483 which is pretty good
ggroc(data=roc.2)

# Testing logs (logged average_num_jobs and average_pop based on EDA above)
glm.3 = glm(treatment ~ log(average_pop) + log(average_num_jobs) + total_pop + percent_employed + percent_transit + mean_commute, family="binomial",data=sub_data)
summary(glm.3) # average_num_jobs, percent_employed, percent_transit, mean_commute are significant, AIC: 155.03
vif(glm.3) # values all below 3
glm.3.probs = predict(glm.3,type="response")
glm.3.pred = (glm.3.probs > 0.5)
table(glm.3.pred, sub_data$treatment)
plot(glm.3) # residuals are mainly centered around 0 123 may be an outlier worth investigating
roc.3 = roc(sub_data$treatment,glm.3.probs)
auc(roc.3) # Area under the curve: 0.7761 which is pretty good
ggroc(data=roc.3)
# Logs appear to be an improvement so will keep in model

# Testing polynomial fit for most significant variables percent_employed
glm.4 = glm(treatment ~ log(average_pop) + log(average_num_jobs) + total_pop + poly(percent_employed,degree=2) + percent_transit + mean_commute, family="binomial",data=sub_data)
summary(glm.4) # average_pop, total_pop and percent_poverty are significant, AIC: 156.9
vif(glm.4) # they are all below 3
glm.4.probs = predict(glm.4,type="response")
glm.4.pred = (glm.4.probs > 0.5)
table(glm.4.pred, sub_data$treatment)
plot(glm.4) # residual plot looks the same
roc.4 = roc(sub_data$treatment,glm.4.probs)
auc(roc.4) # Area under the curve: 0.7725
ggroc(data=roc.4) 
# Polynomials didn't seem to help much

formula_high = "treatment ~ log(average_pop) + log(average_num_jobs) + total_pop + percent_employed + percent_transit + mean_commute"
range.high = quantile(sub_data$high_speed_dist, probs = seq(0.25, 0.75, by= 0.01)) # decile
output.high = sweep.boot(range.high, sub_data, sub_data$high_speed_dist, formula_high)
results.high = sweep(range.high, sub_data, sub_data$high_speed_dist, formula_high)
output.high.mean = results.high$output
high.max.weights = results.high$max.weights
upper.high.ci = output.high.mean + 1.96 * apply(output.high,1, sd)
lower.high.ci = output.high.mean - 1.96 * apply(output.high,1, sd)

par(mfrow=c(1,2), cex=1.5, lwd=1, bg=NA)
plot(range.high, output.high.mean, type="l", xlab='Percentile (km)', ylab='Treatment effect (trips)', ylim = c(min(lower.high.ci), max(upper.high.ci)), lwd=2)
lines(range.high, upper.high.ci, lwd=2, lty=3)
lines(range.high, lower.high.ci, lwd=2, lty=3)
plot(range.high, high.max.weights, type="l", xlab='Percentile (km)', ylab = 'Max IPW Weight', lwd=2)

#write.csv(range.high, 'Data/philly_high_IPW_range.csv')
#write.csv(output.high, 'Data/philly_high_IPW_sweep.csv')
#write.csv(results.high, 'Data/philly_high_IPW_res.csv')
#write.csv(upper.high.ci, 'Data/philly_high_IPW_upper.csv')
#write.csv(lower.high.ci, 'Data/philly_high_IPW_lower.csv')

par(mfrow=c(1,1), cex=1.5, lwd=1, bg=NA)

##################################### COMBO TRANSPORTATION ###############################################

sub_data = sub_data %>% mutate(treatment=(pmin(high_speed_dist,trolley_dist,rail_dist)<=0.402)) 
# not sure if we should use 0.402 because the 25th percentile for rail is around 0.7 km
summary(sub_data$treatment) # Not bad breakdown

sub_data %>%
  gather(-treatment, key = "var", value = "value") %>%
  ggplot(aes(x = treatment, y = value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
# Largest differences appear in percent_employed, percent_poverty, income_per_cap
# Smallest differences appear in average_num_jobs, average_pop, average_ridership (many variables don't have much diff)

# BASELINE MODEL
glm.1 = glm(treatment~.-rail_dist-trolley_dist-bus_dist-high_speed_dist-average_ridership, family="binomial",data=sub_data)
summary(glm.1) # only coefficient for percent_poverty, percent_walk, total_num_jobs are significant, AIC: 154.29
glm.1.probs = predict(glm.1,type="response")
glm.1.pred = (glm.1.probs > 0.5)
table(glm.1.pred, sub_data$treatment)
vif(glm.1) # vars with value over 5: percent_minority, percent_employed, income_per_cap, percent_walk, total_num_jobs
           # professional, service, office, construction, production
           # remove: professional, service, office, construction, production, total_num_jobs, income_per_cap, percent_minority

# Subset of variables
glm.2 = glm(treatment~average_pop+average_num_jobs+total_pop+percent_employed+percent_poverty+percent_transit+percent_walk+mean_commute, family="binomial",data=sub_data)
summary(glm.2) # average_num_jobs, percent_poverty, mean_commute are significant, AIC: 145.95
vif(glm.2) # everything is below 5
glm.2.probs = predict(glm.2,type="response")
glm.2.pred = (glm.2.probs > 0.5)
table(glm.2.pred, sub_data$treatment)
lrtest(glm.1,glm.2) # high p-value so fail to reject the null hypothesis. larger model is not necessarily better
plot(glm.2) # residuals look pretty good
roc.2 = roc(sub_data$treatment,glm.2.probs)
auc(roc.2) # Area under the curve: 0.8168 which is really good
ggroc(data=roc.2)

# Testing logs (logged average_num_jobs and average_pop based on EDA above)
glm.3 = glm(treatment~log(average_pop)+log(average_num_jobs)+total_pop+percent_employed+percent_poverty+percent_transit+percent_walk+mean_commute, family="binomial",data=sub_data)
summary(glm.3) # average_num_jobs, percent_poverty, percent_walk are significant, AIC: 142.7
vif(glm.3) # values all below 3
glm.3.probs = predict(glm.3,type="response")
glm.3.pred = (glm.3.probs > 0.5)
table(glm.3.pred, sub_data$treatment)
plot(glm.3)
roc.3 = roc(sub_data$treatment,glm.3.probs)
auc(roc.3) # Area under the curve: 0.8284
ggroc(data=roc.3)
# Logs appear to be an improvement so will keep in model

# Testing polynomial fit for most significant variables percent_poverty
glm.4 = glm(treatment~log(average_pop)+log(average_num_jobs)+total_pop+percent_employed+poly(percent_poverty,degree=2)+percent_transit+percent_walk+mean_commute, family="binomial",data=sub_data)
summary(glm.4) # average_num_jobs, percent_poverty, percent_walk are significant, AIC: 144.7
vif(glm.4) # they are all below 3
glm.4.probs = predict(glm.4,type="response")
glm.4.pred = (glm.4.probs > 0.5)
table(glm.4.pred, sub_data$treatment)
plot(glm.4) # residual plot looks the same
roc.4 = roc(sub_data$treatment,glm.4.probs)
auc(roc.4) # Area under the curve: 0.8287
ggroc(data=roc.4) 
# AUC slightly improve so will keep poly in model

formula_trans = "treatment~log(average_pop)+log(average_num_jobs)+total_pop+percent_employed+poly(percent_poverty,degree=2)+percent_transit+percent_walk+mean_commute"
range.trans = seq(0.2, 0.804, by= 0.01) 
output.trans = sweep.boot(range.trans, sub_data, pmin(sub_data$trolley_dist,sub_data$high_speed_dist,sub_data$rail_dist),formula_trans)
results.trans = sweep(range.trans, sub_data, pmin(sub_data$trolley_dist,sub_data$high_speed_dist,sub_data$rail_dist),formula_trans)
output.trans.mean = results.trans$output
trans.max.weights = results.trans$max.weights
upper.trans.ci = output.trans.mean + 1.96 * apply(output.trans,1, sd)
lower.trans.ci = output.trans.mean - 1.96 * apply(output.trans,1, sd)

par(mfrow=c(1,2), cex=1.5, lwd=1, bg=NA)
plot(range.trans, output.trans.mean, type="l", xlab='Percentile (km)', ylab='Treatment effect (trips)', ylim = c(min(lower.trans.ci), max(upper.trans.ci)), lwd=2)
lines(range.trans, upper.trans.ci, lwd=2, lty=3)
lines(range.trans, lower.trans.ci, lwd=2, lty=3)
plot(range.trans,trans.max.weights, type="l", xlab='Percentile (km)', ylab = 'Max IPW Weight', lwd=2)

#write.csv(range.trans, 'Data/philly_combo_IPW_range.csv')
#write.csv(output.trans, 'Data/philly_combo_IPW_sweep.csv')
#write.csv(results.trans, 'Data/philly_combo_IPW_res.csv')
#write.csv(upper.trans.ci, 'Data/philly_combo_IPW_upper.csv')
#write.csv(lower.trans.ci, 'Data/philly_combo_IPW_lower.csv')
