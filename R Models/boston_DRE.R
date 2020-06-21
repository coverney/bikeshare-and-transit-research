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

d_data = loadData("Data/boston_summer_bike_station_data.csv")
colnames(d_data)

sub_data = d_data %>% select(percent_minority,percent_employed,average_ridership,
                             total_pop,income_per_cap,percent_poverty,total_num_jobs,
                             percent_transit,mean_commute,average_pop,average_num_jobs,
                             tram_dist,bus_dist,rail_dist,subway_dist,
                             ferry_dist,professional,service,office,construction,
                             production,walk)

sub_data[!complete.cases(sub_data),"percent_employed"] = 0 # population = 0 so set percent_employed=0
sum(!complete.cases(sub_data))

##################################### COMBO TRANSPORTATION ###############################################

sub_data = sub_data %>% mutate(treatment=(pmin(tram_dist,subway_dist,rail_dist)<=0.402))
summary(sub_data$treatment) # Not bad breakdown

########### TREATED GROUP ########################################
data_treated = sub_data[(sub_data$treatment == TRUE),]

# Build correlation matrix
ggcorr(data_treated,label=TRUE)
# Interesting things to note - some relationship with walk, mean_commute, type of jobs

data_treated %>% gather() %>% ggplot(aes(value)) + facet_wrap(~ key, scales = "free") + geom_density()
# A few heavy tailed distributions - average_number_jobs, average_pop, average ridership

hist(data_treated$average_ridership) # slightly skewed 
hist(log(data_treated$average_ridership)) # more normal - likely will use transformation

data_treated =data_treated %>% mutate(log_average_ridership=log(average_ridership+1))
ggcorr(data_treated,label=TRUE) # stronger relationships

data_treated %>%
  gather(-average_ridership, key = "var", value = "value") %>%
  ggplot(aes(x = value, y = log(average_ridership))) +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
# Jobs show also might take log - population not so much

# Build model for treated group
lm.1 = lm(log_average_ridership~.-average_ridership-rail_dist-subway_dist-ferry_dist-tram_dist-bus_dist-treatment, data=data_treated)
summary(lm.1)
vif(lm.1)

lm.2 = lm(log_average_ridership~percent_minority+professional+office+average_pop+average_num_jobs+walk+mean_commute+percent_transit, data=data_treated)
summary(lm.2)
vif(lm.2)
anova(lm.1,lm.2)

lm.3 = lm(log_average_ridership~percent_minority+professional+office+average_pop+log(average_num_jobs)+walk+mean_commute+percent_transit, data=data_treated)
summary(lm.3)
vif(lm.3)
# No evidence for polynomials or interactions for significant variables

plot(lm.3) # Pretty good
# Residuals for variables - no trends
plot(data_treated$percent_minority, residuals(lm.3))
plot(data_treated$professional, residuals(lm.3))
plot(data_treated$office, residuals(lm.3))
plot(data_treated$average_pop, residuals(lm.3))
plot(log(data_treated$average_num_jobs), residuals(lm.3))
plot(data_treated$walk, residuals(lm.3))
plot(data_treated$mean_commute, residuals(lm.3))
plot(data_treated$percent_transit, residuals(lm.3))

formula_treated = "log_average_ridership~percent_minority+professional+office+average_pop+log(average_num_jobs)+walk+mean_commute+percent_transit"

########### CONTROL GROUP ########################################

data_control = sub_data[(sub_data$treatment == FALSE),]

# Build correlation matrix
ggcorr(data_control,label=TRUE)
# Interesting things to note - somewhat weaker relationships in general

data_control %>% gather() %>% ggplot(aes(value)) + facet_wrap(~ key, scales = "free") + geom_density()
# Same heavy tailed distributions

hist(data_control$average_ridership) # more higly skewed than treatment
hist(log(data_control$average_ridership)) # more normal - likely will use transformation

data_control = data_control %>% mutate(log_average_ridership=log(average_ridership+1))
ggcorr(data_control,label=TRUE) # stronger relationships

data_control %>%
  gather(-log_average_ridership, key = "var", value = "value") %>%
  ggplot(aes(x = value, y = log_average_ridership)) +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
# Maybe take log of jobs - not as much evidence

# Build model for control group
lm.1 = lm(log_average_ridership~.-average_ridership-rail_dist-subway_dist-ferry_dist-tram_dist-bus_dist-treatment, data=data_control)
summary(lm.1) # Many more significant variables

lm.2 = lm(log_average_ridership~percent_employed+percent_minority+total_num_jobs+percent_transit+mean_commute+walk+construction+service+office+professional, data=data_control)
summary(lm.2)
# Could only remove population and income
anova(lm.1,lm.2)

# Tested log of jobs - not helpful

plot(lm.2) # Pretty good
# Residuals for variables - no trends
plot(data_control$percent_minority, residuals(lm.2))
plot(data_control$professional, residuals(lm.2))
plot(data_control$office, residuals(lm.2))
plot(data_control$service, residuals(lm.2))
plot(data_control$construction, residuals(lm.2))
plot(data_control$total_num_jobs, residuals(lm.2))
plot(data_control$percent_employed, residuals(lm.2))
plot(data_control$walk, residuals(lm.2))
plot(data_control$mean_commute, residuals(lm.2))
plot(data_control$percent_transit, residuals(lm.2))

lm.4 = lm(log_average_ridership~percent_employed+poly(percent_minority,degree=2)+total_num_jobs+percent_transit+mean_commute+poly(walk,degree=2)+professional+office+production+service+construction, data=data_control)
summary(lm.4)
anova(lm.2,lm.4)

lm.5 = lm(log_average_ridership~percent_employed+poly(percent_minority,degree=2)+total_num_jobs+percent_transit+mean_commute+poly(walk,degree=2)+professional+office+production+service+construction+total_num_jobs*percent_transit, data=data_control)
summary(lm.5)
anova(lm.4,lm.5)

formula_control = "log_average_ridership~percent_employed+poly(percent_minority,degree=2)+total_num_jobs+percent_transit+mean_commute+poly(walk,degree=2)+professional+office+production+service+construction+total_num_jobs*percent_transit"

formula_prob = "treatment~poly(walk,degree=2)+poly(percent_transit,degree=2)+mean_commute+income_per_cap+average_pop+average_num_jobs"

dre = function(data, transit_col, val, formula_prob, formula_control, formula_treated){
  data = data %>% mutate(treatment=(pmin(transit_col)<=val))
  glm.1 = glm(formula=as.formula(formula_prob), family="binomial",data=data)
  lm.control = lm(formula=as.formula(formula_control),data=data)
  lm.treated = lm(formula=as.formula(formula_treated),data=data)
  pred.probs = predict(glm.1, type = "response")
  weights.used = data$treatment/pred.probs + (1 - data$treatment)/(1 - pred.probs)
  pred.treated = exp(predict(lm.treated, type="response"))-1 # Take exp here??
  pred.control = exp(predict(lm.control, type="response"))-1
  #pred.treated = predict(lm.treated, type="response") # Take exp here??
  #pred.control = predict(lm.control, type="response")
  treat.est = pred.treated - (pred.treated-data$average_ridership)*(data$treatment/pred.probs)
  control.est = pred.control - (pred.control-data$average_ridership)*((1-data$treatment)/(1-pred.probs))
  trt_eff = mean(treat.est) - mean(control.est)
  #return(weights.used) # return the weights
  return(list("trt.eff"=trt_eff,"max.weight"=max(weights.used)))
}

sub_data = sub_data %>% mutate(log_average_ridership = log(average_ridership+1))
dre(sub_data, pmin(sub_data$tram_dist,sub_data$subway_dist,sub_data$rail_dist), 0.408, formula_prob, formula_control, formula_treated)

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
output.trans = sweep.boot(range.trans, sub_data, pmin(sub_data$tram_dist,sub_data$subway_dist,sub_data$rail_dist),formula_prob, formula_control, formula_treated)
results.trans = sweep(range.trans, sub_data, pmin(sub_data$tram_dist,sub_data$subway_dist,sub_data$rail_dist),formula_prob, formula_control, formula_treated)
output.trans.mean = results.trans$output
trans.max.weights = results.trans$max.weights
upper.trans.ci = output.trans.mean + 1.96 * apply(output.trans,1, sd)
lower.trans.ci = output.trans.mean - 1.96 * apply(output.trans,1, sd)

par(mfrow=c(1,2), cex=1.5, lwd=1, bg=NA)
plot(range.trans, output.trans.mean, type="l", xlab='Percentile (km)', ylab='Treatment effect (trips)', ylim = c(min(lower.trans.ci), max(upper.trans.ci)), lwd=2)
lines(range.trans, upper.trans.ci, lwd=2, lty=3)
lines(range.trans, lower.trans.ci, lwd=2, lty=3)
plot(range.trans,trans.max.weights, type="l", xlab='Percentile (km)', ylab = 'Max IPW Weight', lwd=2)

#write.csv(range.trans, 'Data/DRE/boston_combo_DRE_range.csv')
#write.csv(output.trans, 'Data/DRE/boston_combo_DRE_sweep.csv')
#write.csv(results.trans, 'Data/DRE/boston_combo_DRE_res.csv')
#write.csv(upper.trans.ci, 'Data/DRE/boston_combo_DRE_upper.csv')
#write.csv(lower.trans.ci, 'Data/DRE/boston_combo_DRE_lower.csv')

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
weights = sweep_W(range.trans, sub_data, pmin(sub_data$tram_dist,sub_data$subway_dist,sub_data$rail_dist),formula_prob, formula_control, formula_treated)
# write.csv(weights, 'Data/DRE/boston_combo_DRE_weights.csv')

