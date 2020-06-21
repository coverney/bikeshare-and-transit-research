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

d_data = loadData("Data/philly_summer_bike_station_data.csv")
colnames(d_data)

# Get all the variables that we are intrested in seeing the distributions of a potentially include
# in our models
sub_data = d_data %>% select(average_ridership,percent_minority,percent_employed,total_pop,income_per_cap,
                             percent_poverty,percent_transit,percent_walk,total_num_jobs,average_pop,average_num_jobs,
                             mean_commute,bus_dist,rail_dist,trolley_dist,high_speed_dist,professional,service,
                             office,construction,production)

sum(!complete.cases(sub_data)) # all samples are complete

##################################### COMBO TRANSPORTATION ###############################################

sub_data = sub_data %>% mutate(treatment=(pmin(high_speed_dist,trolley_dist,rail_dist)<=0.402)) 
summary(sub_data$treatment) # Not bad breakdown: 54 false and 70 true

########################################## TREATED GROUP ################################################
data_treated = sub_data[(sub_data$treatment == TRUE),]

# Build correlation matrix
ggcorr(data_treated,label=TRUE)
# Interesting things to note - average ridership has some relationship with income(0.5), num_jobs(0.5), rail_dist(-0.6)
#                              professional(0.6), service(-0.5), production(-0.5)

data_treated %>% gather() %>% ggplot(aes(value)) + facet_wrap(~ key, scales = "free") + geom_density()
# A few heavy tailed distributions - average_number_jobs, average_pop, average ridership, bus_dist, total_num_jobs, trolley_dist
hist(data_treated$average_ridership) # slightly skewed 
hist(log(data_treated$average_ridership)) # slightly more normal - may or may not use transformation

data_treated =data_treated %>% mutate(log_average_ridership=log(average_ridership+1))
ggcorr(data_treated,label=TRUE) # has stronger relationships: production(-0.6), service(-0.6), professional(0.7),
#                                                             rail_dist(-0.6), num_jobs(0.5), percent_walk(0.5),
#                                                             income_per_cap(0.5), percent_minority(-0.5)
# I will use log_average_ridership as my dependent variable in the linear regression model

data_treated %>%
  gather(-average_ridership, key = "var", value = "value") %>%
  ggplot(aes(x = value, y = log(average_ridership))) +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
# Number of jobs might also be logged and maybe average_pop and bus_dist

# Build model for treated group, start with baseline model
lm.1 = lm(log_average_ridership~.-average_ridership-rail_dist-trolley_dist-bus_dist-high_speed_dist-treatment, data=data_treated)
summary(lm.1) # only coefficients for percent_employed, percent_walk, and total_num_jobs(coeff really small) are sig. R-squared:0.5428 
vif(lm.1) # vars with value over 5: percent_minority, percent_employed, income_per_cap, percent_poverty, percent_walk,
#                                   total_num_jobs, average_num_jobs, professional, service, office, construction, production
#           remove: service, office, construction, production, average_num_jobs, percent_poverty, income_per_cap, percent_walk   
# I kept professional because it had a pretty strong relationship with log_average_ridership from EDA

# Subset of variables
lm.2 = lm(log_average_ridership~professional+average_pop+total_pop+total_num_jobs+percent_minority+percent_transit+mean_commute+percent_employed, data=data_treated)
summary(lm.2) # professional, average_pop, total_num_jobs, percent_employed are significant, Adjusted R-squared:  0.5413 
vif(lm.2) # all below 5
anova(lm.1,lm.2) # p-value 0.4289, high

# log number of jobs
lm.3 = lm(log_average_ridership~professional+average_pop+total_pop+log(total_num_jobs)+percent_minority+percent_transit+mean_commute+percent_employed, data=data_treated)
summary(lm.3) # professional and log(total_num_jobs) are significant, Adjusted R-squared:  0.4387 
vif(lm.3) # all below 5
# logging didn't improve model

plot(lm.2) # Pretty good
# Residuals for variables
plot(data_treated$professional, residuals(lm.2)) # no trend
plot(data_treated$average_pop, residuals(lm.2)) # only small evidence of > trend (considering how only one point forms the tip)
plot(data_treated$total_pop, residuals(lm.2)) # no trend
plot(data_treated$total_num_jobs, residuals(lm.2)) # no trend
plot(data_treated$percent_minority, residuals(lm.2)) # no trend
plot(data_treated$mean_commute, residuals(lm.2)) # no trend
plot(data_treated$percent_transit, residuals(lm.2)) # no trend
plot(data_treated$percent_employed, residuals(lm.2)) # no trend

formula_treated = "log_average_ridership~professional+average_pop+total_pop+total_num_jobs+percent_minority+percent_transit+mean_commute+percent_employed"

########################################## CONTROL GROUP ################################################

data_control = sub_data[(sub_data$treatment == FALSE),]

# Build correlation matrix
ggcorr(data_control,label=TRUE)
# interestingly, there seems to be stronger relationships with average_ridership: percent_minority(-0.6), percent_employed(0.7),
#               income_per_cap(0.7), percent_transit(-0.5), professional(0.7), service(-0.6), production(-0.6)
# --> could the cutoff we use be too low?

data_control %>% gather() %>% ggplot(aes(value)) + facet_wrap(~ key, scales = "free") + geom_density()
# A few heavy tailed distributions - average_num_jobs, construction, total_num_jobs
hist(data_control$average_ridership) # doesn't look too skewed 
hist(log(data_control$average_ridership)) # looks more skewed --> probably won't transform

data_control %>%
  gather(-average_ridership, key = "var", value = "value") %>%
  ggplot(aes(x = value, y = average_ridership)) +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
# might log number of jobs

# Build model for control group, baseline model
lm.1 = lm(average_ridership~.-rail_dist-trolley_dist-bus_dist-high_speed_dist-treatment, data=data_control)
summary(lm.1) # Only total_num_jobs is sig, Adjusted R-squared:  0.4227
vif(lm.1) # vars with value over 5: percent_minority, percent_employed, income_per_cap, percent_transit, total_num_jobs,
#                                   mean_commute, professional, service, office, construction, production
#           remove: professional, service, office, construction, production, income_per_cap, average_num_jobs, average_pop

# Subset of variables
lm.2 = lm(average_ridership~percent_employed+percent_walk+total_num_jobs+total_pop+percent_poverty+percent_transit+mean_commute+percent_minority, data=data_control)
summary(lm.2) # percent_employed is significant, Adjusted R-squared:  0.4733
vif(lm.2) # all below 6 (only percent_minority is slightly below 5, removing it decreases r-squared)
anova(lm.1,lm.2) # p-value is high
# logging total_num_jobs doesn't improve the model

plot(lm.2) # Pretty good
# Residuals for variables - no trends
plot(data_control$total_pop, residuals(lm.2)) # no trend 
plot(data_control$total_num_jobs, residuals(lm.2)) # no trend with some outliers
plot(data_control$percent_minority, residuals(lm.2)) # no trend
plot(data_control$percent_poverty, residuals(lm.2)) # no trend
plot(data_control$percent_employed, residuals(lm.2)) # no trend
plot(data_control$percent_walk, residuals(lm.2)) # no trend
plot(data_control$mean_commute, residuals(lm.2)) # no trend
plot(data_control$percent_transit, residuals(lm.2)) # no trend

formula_control = "average_ridership~percent_employed+percent_walk+total_num_jobs+total_pop+percent_poverty+percent_transit+mean_commute+percent_minority"

########################################## TREATED GROUP ################################################
formula_prob = "treatment~log(average_pop)+log(average_num_jobs)+total_pop+percent_employed+poly(percent_poverty,degree=2)+percent_transit+percent_walk+mean_commute"

dre = function(data, transit_col, val, formula_prob, formula_control, formula_treated){
  data = data %>% mutate(treatment=(pmin(transit_col)<=val))
  glm.1 = glm(formula=as.formula(formula_prob), family="binomial",data=data)
  lm.control = lm(formula=as.formula(formula_control),data=data)
  lm.treated = lm(formula=as.formula(formula_treated),data=data)
  pred.probs = predict(glm.1, type = "response")
  weights.used = data$treatment/pred.probs + (1 - data$treatment)/(1 - pred.probs)
  pred.treated = exp(predict(lm.treated, type="response"))-1 # Take exp here??
  pred.control = predict(lm.control, type="response")
  #pred.treated = predict(lm.treated, type="response") # Take exp here??
  #pred.control = predict(lm.control, type="response")
  treat.est = pred.treated - (pred.treated-data$average_ridership)*(data$treatment/pred.probs)
  control.est = pred.control - (pred.control-data$average_ridership)*((1-data$treatment)/(1-pred.probs))
  trt_eff = mean(treat.est) - mean(control.est)
  #return(weights.used) # return the weights
  return(list("trt.eff"=trt_eff,"max.weight"=max(weights.used)))
}

sub_data = sub_data %>% mutate(log_average_ridership = log(average_ridership+1))
dre(sub_data, pmin(sub_data$trolley_dist,sub_data$high_speed_dist,sub_data$rail_dist), 0.402, formula_prob, formula_control, formula_treated)

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

range.trans = seq(0.201, 0.804, by= 0.01) 
output.trans = sweep.boot(range.trans, sub_data, pmin(sub_data$trolley_dist,sub_data$high_speed_dist,sub_data$rail_dist),formula_prob, formula_control, formula_treated)
results.trans = sweep(range.trans, sub_data, pmin(sub_data$trolley_dist,sub_data$high_speed_dist,sub_data$rail_dist),formula_prob, formula_control, formula_treated)
output.trans.mean = results.trans$output
trans.max.weights = results.trans$max.weights
upper.trans.ci = output.trans.mean + 1.96 * apply(output.trans,1, sd)
lower.trans.ci = output.trans.mean - 1.96 * apply(output.trans,1, sd)

par(mfrow=c(1,2), cex=1.5, lwd=1, bg=NA)
plot(range.trans, output.trans.mean, type="l", xlab='Percentile (km)', ylab='Treatment effect (trips)', ylim = c(min(lower.trans.ci), max(upper.trans.ci)), lwd=2)
lines(range.trans, upper.trans.ci, lwd=2, lty=3)
lines(range.trans, lower.trans.ci, lwd=2, lty=3)
plot(range.trans,trans.max.weights, type="l", xlab='Percentile (km)', ylab = 'Max IPW Weight', lwd=2)
par(mfrow=c(1,1), cex=1.5, lwd=1, bg=NA)

#write.csv(range.trans, 'Data/philly_combo_DRE_range.csv')
#write.csv(output.trans, 'Data/philly_combo_DRE_sweep.csv')
#write.csv(results.trans, 'Data/philly_combo_DRE_res.csv')
#write.csv(upper.trans.ci, 'Data/philly_combo_DRE_upper.csv')
#write.csv(lower.trans.ci, 'Data/philly_combo_DRE_lower.csv')

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
weights = sweep_W(range.trans, sub_data, pmin(sub_data$trolley_dist,sub_data$high_speed_dist,sub_data$rail_dist),formula_prob, formula_control, formula_treated)
#write.csv(weights, 'Data/DRE/philly_combo_DRE_weights.csv')
