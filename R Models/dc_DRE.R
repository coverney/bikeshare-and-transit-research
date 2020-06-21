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

d_data = loadData("Data/dc_summer_bike_station_data.csv")
colnames(d_data)

# Get all the variables that we are intrested in seeing the distributions of a potentially include
# in our models
sub_data = d_data %>% select(average_ridership,percent_minority,percent_employed,total_pop,income_per_cap,
                             percent_poverty,percent_transit,percent_walk,total_num_jobs,average_pop,average_num_jobs,
                             mean_commute,bus_dist,rail_dist,professional,service,
                             office,construction,production)

sum(!complete.cases(sub_data)) # should be 0

##################################### Rail ###############################################
sub_data = sub_data %>% mutate(treatment=(pmin(rail_dist)<=0.402)) 
summary(sub_data$treatment) # Not the best breakdown: 335 false and 166 true

########################################## TREATED GROUP ################################################
data_treated = sub_data[(sub_data$treatment == TRUE),]

# Build correlation matrix
ggcorr(data_treated,label=TRUE) # percent_walk (0.4), mean_commute (-0.4)

data_treated %>% gather() %>% ggplot(aes(value)) + facet_wrap(~ key, scales = "free") + geom_density()
# production, total_num_jobs, bus_dist, construction skewed to right; professional skewed to left

hist(data_treated$average_ridership) # skewed 
hist(log(data_treated$average_ridership)) # more normal - will prob log

data_treated =data_treated %>% mutate(log_average_ridership=log(average_ridership+1))
ggcorr(data_treated,label=TRUE) # stronger relationships: mean_commute (-0.6), percent_walk (0.6)

# I will use log_average_ridership as my dependent variable in the linear regression model

data_treated %>%
  gather(-average_ridership, key = "var", value = "value") %>%
  ggplot(aes(x = value, y = log(average_ridership))) +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
# percent_walk and total_num_jobs seem to have a non linear relationship with log_average_ridership
# mean_commute has negative relationship
# income has positive relationship

# Build model for treated group, start with baseline model
lm.1 = lm(log_average_ridership~.-average_ridership-rail_dist-bus_dist-treatment, data=data_treated)
summary(lm.1) # only coefficients for percent_minority, percent_employed, percent_poverty, average_pop
              # average_num_jobs (really small coef), mean_commute are sig. 
              # R-squared:0.6368 
vif(lm.1) # vars with value over 5: total_num_jobs, average_num_jobs, professional, service, office, construction, production
          # remove: total_num_jobs (average_num_jobs is sig), service, professional, construction, production
# The reason why I kept office instead of professional was because when I was testing out the subset of variables
# I tried each of the types of jobs and office resulted in the model with the highest adjusted r^2 and p-value

# Subset of variables
lm.2 = lm(log_average_ridership~office+average_num_jobs+percent_walk+income_per_cap+percent_poverty+average_pop+total_pop+percent_minority+percent_transit+mean_commute+percent_employed, data=data_treated)
summary(lm.2) # office, percent_walk, percent_poverty, average_num_jobs, percent_minority, percent_transit, mean_commute, percent_employed are significant, 
              # Adjusted R-squared:  0.6257 
vif(lm.2) # all below 5
anova(lm.1,lm.2) # borderline p-value (0.09)

# Try logging average_num_jobs
lm.3 = lm(log_average_ridership~office+log(average_num_jobs)+percent_walk+income_per_cap+percent_poverty+average_pop+total_pop+percent_minority+percent_transit+mean_commute+percent_employed, data=data_treated)
summary(lm.3) # office, percent_poverty, log(average_num_jobs), percent_minority, percent_transit, mean_commute, percent_employed are significant, 
              # Adjusted R-squared: 0.6066 
vif(lm.3) # all below 5
# Logging didn't improve model

plot(lm.2) # Pretty good
# Residuals for variables
plot(data_treated$office, residuals(lm.2)) # no trend
plot(data_treated$average_pop, residuals(lm.2)) # no trend
plot(data_treated$total_pop, residuals(lm.2)) # no trend
plot(data_treated$income_per_cap, residuals(lm.2)) # no trend
plot(data_treated$percent_walk, residuals(lm.2)) # no trend, might be a little >
plot(data_treated$percent_minority, residuals(lm.2)) # no trend
plot(data_treated$mean_commute, residuals(lm.2)) # no trend
plot(data_treated$percent_transit, residuals(lm.2)) # no trend
plot(data_treated$percent_employed, residuals(lm.2)) # no trend
plot(data_treated$percent_poverty, residuals(lm.2)) # no trend
plot(data_treated$average_num_jobs, residuals(lm.2)) # no trend

# Logging percent_transit didn't improve the model

# Try poly fit for percent_walk
lm.4 = lm(log_average_ridership~office+poly(percent_walk, degree=2)+income_per_cap+percent_poverty+average_pop+total_pop+average_num_jobs+percent_minority+percent_transit+mean_commute+percent_employed, data=data_treated)
summary(lm.4) # office, poly(percent_walk, degree = 2)1, percent_poverty, average_num_jobs, percent_minority, mean_commute, percent_employed are significant, 
              # Adjusted R-squared: 0.6244
vif(lm.4) # percent_walk has value of 7
# model doesn't improve with poly fit, so won't use

formula_treated = "log_average_ridership~office+percent_walk+income_per_cap+percent_poverty+average_pop+total_pop+average_num_jobs+percent_minority+percent_transit+mean_commute+percent_employed"

########################################## CONTROL GROUP ################################################

data_control = sub_data[(sub_data$treatment == FALSE),]

# Build correlation matrix
ggcorr(data_control,label=TRUE) # nothing is correlated (highest value is +- 0.2)

data_control %>% gather() %>% ggplot(aes(value)) + facet_wrap(~ key, scales = "free") + geom_density()
# A few heavy tailed distributions - average_num_jobs, percent_walk, total_num_jobs, production, percent_poverty
hist(data_control$average_ridership) # very skewed
hist(log(data_control$average_ridership)) # less skewed, will use log in model

data_control =data_control %>% mutate(log_average_ridership=log(average_ridership+1))
ggcorr(data_control,label=TRUE) # much stronger relationships, rail_dist (-0.5), percent_walk (0.5), mean_commute (-0.5)

data_control %>%
  gather(-average_ridership, key = "var", value = "value") %>%
  ggplot(aes(x = value, y = log(average_ridership))) +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
# there might be a non linear relationship with percent_walk, income looks slightly positive


# Build model for control group, baseline model
lm.1 = lm(log_average_ridership~.-average_ridership-rail_dist-bus_dist-treatment, data=data_control)
summary(lm.1) # percent_minority, percent_employed, total_pop, income_per_cap, percent_poverty, percent_transit, percent_walk, total_num_jobs, 
              # average_num_jobs, mean_commute, service, office are sig. 
              # Adjusted R-squared:  0.581
vif(lm.1) # vars with value over 5: professional, service
          # remove professional, which I got from trying to remove different combinations of occupation variables

# Subset of variables
lm.2 = lm(log_average_ridership~service+construction+production+total_num_jobs+office+percent_walk+income_per_cap+percent_poverty+average_pop+total_pop+average_num_jobs+percent_minority+percent_transit+mean_commute+percent_employed, data=data_control)
summary(lm.2) # service, construction, production, total_num_jobs, office, percent_walk, income_per_cap
              # percent_poverty, total_pop, average_num_jobs, percent_transit, percent_minority, mean_commute, percent_employed are significant
              # Adjusted R-squared: 0.5823 
vif(lm.2) # all below 5
anova(lm.1,lm.2) # high p-value

plot(lm.2) # Pretty good
# Residuals for variables
plot(data_control$office, residuals(lm.2)) # no trend
plot(data_control$service, residuals(lm.2)) # no trend
plot(data_control$construction, residuals(lm.2)) # no trend
plot(data_control$production, residuals(lm.2)) # no trend
plot(data_control$total_num_jobs, residuals(lm.2)) # no trend
plot(data_control$average_pop, residuals(lm.2)) # no trend
plot(data_control$total_pop, residuals(lm.2)) # no trend
plot(data_control$income_per_cap, residuals(lm.2)) # no trend
plot(data_control$percent_walk, residuals(lm.2)) # no trend
plot(data_control$percent_minority, residuals(lm.2)) # no trend
plot(data_control$mean_commute, residuals(lm.2)) # no trend
plot(data_control$percent_transit, residuals(lm.2)) # no trend
plot(data_control$percent_employed, residuals(lm.2)) # no trend
plot(data_control$percent_poverty, residuals(lm.2)) # no trend
plot(data_control$average_num_jobs, residuals(lm.2)) # no trend

formula_control = "log_average_ridership~service+construction+production+total_num_jobs+office+percent_walk+income_per_cap+percent_poverty+average_pop+total_pop+average_num_jobs+percent_minority+percent_transit+mean_commute+percent_employed"

formula_prob = "treatment~poly(percent_transit,degree=2)+percent_walk+log(average_pop)+total_pop+log(average_num_jobs)+income_per_cap+percent_employed+percent_poverty+mean_commute+professional"
dre = function(data, transit_col, val, formula_prob, formula_control, formula_treated){
  data = data %>% mutate(treatment=(pmin(transit_col)<=val))
  glm.1 = glm(formula=as.formula(formula_prob), family="binomial",data=data)
  lm.control = lm(formula=as.formula(formula_control),data=data)
  lm.treated = lm(formula=as.formula(formula_treated),data=data)
  pred.probs = predict(glm.1, type = "response")
  weights.used = data$treatment/pred.probs + (1 - data$treatment)/(1 - pred.probs)
  pred.treated = exp(predict(lm.treated, type="response"))-1
  pred.control = exp(predict(lm.control, type="response"))-1
  #pred.treated = predict(lm.treated, type="response") 
  #pred.control = predict(lm.control, type="response")
  treat.est = pred.treated - (pred.treated-data$average_ridership)*(data$treatment/pred.probs)
  control.est = pred.control - (pred.control-data$average_ridership)*((1-data$treatment)/(1-pred.probs))
  trt_eff = mean(treat.est) - mean(control.est)
  #return(weights.used) # return the weights
  return(list("trt.eff"=trt_eff,"max.weight"=max(weights.used)))
}

sweep = function(range, data, transit_col, formula_prob, formula_control, formula_treated){
  output = vector()
  max_weights = vector()
  for (val in range){
    results = dre(data, transit_col, val, formula_prob, formula_control, formula_treated)
    output = append(output, results$trt.eff)
    max_weights = append(max_weights, results$max.weight)
  }
  return(list("output"=output,"max.weights"=max_weights))
}

sweep.boot = function(range, data, transit_col, formula_prob, formula_control, formula_treated){
  n.boot = 1000
  output = matrix(NA, nrow = length(range), ncol = n.boot)
  for(j in 1:n.boot){
    index = sample(1:nrow(data), nrow(data), replace = TRUE)
    data.boot = data[index, ]
    results = sweep(range, data.boot, transit_col, formula_prob, formula_control, formula_treated)
    output[, j] = results$output
  }
  return(output)
}

sub_data = sub_data %>% mutate(log_average_ridership = log(average_ridership+1))
dre(sub_data, pmin(sub_data$rail_dist), 0.402, formula_prob, formula_control, formula_treated)
# $`trt.eff`
# [1] 11.67849
# $max.weight
# [1] 36.88149

range.trans = seq(0.201, 0.804, by= 0.01) 
output.trans = sweep.boot(range.trans, sub_data, pmin(sub_data$rail_dist),formula_prob, formula_control, formula_treated)
results.trans = sweep(range.trans, sub_data, pmin(sub_data$rail_dist),formula_prob, formula_control, formula_treated)
output.trans.mean = results.trans$output
trans.max.weights = results.trans$max.weights
upper.trans.ci = output.trans.mean + 1.96 * apply(output.trans,1, sd)
lower.trans.ci = output.trans.mean - 1.96 * apply(output.trans,1, sd)

par(mfrow=c(1,2), cex=1.5, lwd=1, bg=NA)
plot(range.rail, output.rail.mean, type="l", xlab='Percentile (km)', ylab='Treatment effect (trips)', ylim = c(min(lower.rail.ci), max(upper.rail.ci)), lwd=2)
lines(range.rail, upper.rail.ci, lwd=2, lty=3)
lines(range.rail, lower.rail.ci, lwd=2, lty=3)
plot(range.rail,rail.max.weights, type="l", xlab='Percentile (km)', ylab = 'Max IPW Weight', lwd=2)
par(mfrow=c(1,1), cex=1.5, lwd=1, bg=NA)

#write.csv(range.trans, 'Data/dc_rail_DRE_range.csv')
#write.csv(output.trans, 'Data/dc_rail_DRE_sweep.csv')
#write.csv(results.trans, 'Data/dc_rail_DRE_res.csv')
#write.csv(upper.trans.ci, 'Data/dc_rail_DRE_upper.csv')
#write.csv(lower.trans.ci, 'Data/dc_rail_DRE_lower.csv')

sweep_W = function(range, data, transit_col, formula_prob, formula_control, formula_treated){
  weights = matrix(NA, nrow = nrow(data), ncol = length(range))
  for (j in 1:length(range)){
    val = range[j]
    weight = dre(data, transit_col, val, formula_prob, formula_control, formula_treated)
    weights[, j] = weight
  }
  return(weights)
}
range.trans = seq(0.201, 0.804, by= 0.01) 
weights = sweep_W(range.trans, sub_data, pmin(sub_data$rail_dist),formula_prob, formula_control, formula_treated)
# write.csv(weights, 'Data/DRE/dc_rail_DRE_weights.csv')
