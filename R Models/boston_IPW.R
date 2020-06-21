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

# Build correlation matrix
ggcorr(sub_data,label=TRUE)
# Interesting things to note - not many strong relationships except with walk

sub_data %>% gather() %>% ggplot(aes(value)) + facet_wrap(~ key, scales = "free") + geom_density()
# A few heavy tailed distributions - job distribution from census, average_number_jobs, total_jobs

hist(sub_data$average_ridership) # highly skewed - could do sweep with log instead but otherwise doesn't matter

sub_data %>%
  gather(-average_ridership, key = "var", value = "value") %>%
  ggplot(aes(x = value, y = average_ridership)) +
  geom_point() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
# Not much to note except maybe relationship with walk (and maybe mean_commute)

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

glm.1 = glm(treatment~.-rail_dist-subway_dist-ferry_dist-tram_dist-bus_dist-average_ridership, family="binomial",data=sub_data)
summary(glm.1)
glm.1.probs = predict(glm.1,type="response")
glm.1.pred = (glm.1.probs > 0.5)
table(glm.1.pred, sub_data$treatment)

# Use a subset of variables in logistic regression model
glm.2 = glm(treatment~income_per_cap+average_pop+average_num_jobs+total_pop+total_num_jobs+walk+percent_minority+percent_employed+percent_poverty, family="binomial",data=sub_data)
summary(glm.2)
glm.2.probs = predict(glm.2,type="response")
vif(glm.2)
glm.2.pred = (glm.2.probs > 0.5)
table(glm.2.pred, sub_data$treatment)
lrtest(glm.1,glm.2)

plot(glm.2) # Entry 27 sticks out
sub_data[27,] # Looks fairly poor but still high ridership - maybe look into where this is

roc.2 = roc(sub_data$treatment,glm.2.probs)
auc(roc.2)
ggroc(data=roc.2)

# Testing log - doesn't improve
glm.3 = glm(treatment~income_per_cap+average_pop+log(average_num_jobs)+total_pop+log(total_num_jobs)+walk+percent_minority+percent_employed+percent_poverty, family="binomial",data=sub_data)
summary(glm.3)
# Also tested interactions, polynomials of most significant variables 

range.rail = quantile(sub_data$rail_dist, probs = seq(0.25, 0.75, by= 0.01)) # decile
output.rail = sweep.boot(range.rail, sub_data, sub_data$rail_dist, "treatment~income_per_cap+average_pop+average_num_jobs+total_pop+total_num_jobs+walk+percent_minority+percent_employed+percent_poverty")
results.rail = sweep(range.rail, sub_data, sub_data$rail_dist, "treatment~income_per_cap+average_pop+average_num_jobs+total_pop+total_num_jobs+walk+percent_minority+percent_employed+percent_poverty")
output.rail.mean = results.rail$output
rail.max.weights = results.rail$max.weights
upper.rail.ci = output.rail.mean + 1.96 * apply(output.rail,1, sd)
lower.rail.ci = output.rail.mean - 1.96 * apply(output.rail,1, sd)

par(mfrow=c(1,2), cex=1.5, lwd=1, bg=NA)
plot(range.rail, output.rail.mean, type="l", xlab='Percentile (km)', ylab='Treatment effect (trips)', ylim = c(min(lower.rail.ci), max(upper.rail.ci)), lwd=2)
lines(range.rail, upper.rail.ci, lwd=2, lty=3)
lines(range.rail, lower.rail.ci, lwd=2, lty=3)
plot(range.rail,rail.max.weights, type="l", xlab='Percentile (km)', ylab = 'Max IPW Weight', lwd=2)

##################################### SUBWAY ###############################################
par(mfrow=c(1,1))
# Subway distance - start by building logistic regression for 50th percentile
quant_50 = quantile(sub_data$subway_dist, 0.50)
sub_data = sub_data %>% mutate(treatment=(pmin(subway_dist)<=quant_50))

sub_data %>%
  gather(-treatment, key = "var", value = "value") %>%
  ggplot(aes(x = treatment, y = value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

glm.1 = glm(treatment~.-rail_dist-subway_dist-ferry_dist-tram_dist-bus_dist-average_ridership, family="binomial",data=sub_data)
summary(glm.1)
glm.1.probs = predict(glm.1,type="response")
glm.1.pred = (glm.1.probs > 0.5)
table(glm.1.pred, sub_data$treatment)

glm.2 = glm(treatment~income_per_cap+average_pop+average_num_jobs+walk+percent_minority+percent_employed+percent_poverty+percent_transit+mean_commute+professional+office, family="binomial",data=sub_data)
summary(glm.2)
glm.2.probs = predict(glm.2,type="response")
vif(glm.2)
glm.2.pred = (glm.2.probs > 0.5)
table(glm.2.pred, sub_data$treatment)
lrtest(glm.1,glm.2)

plot(glm.2) # 27 again

roc.2 = roc(sub_data$treatment,glm.2.probs)
auc(roc.2) # pretty good AUC
ggroc(data=roc.2)

glm.3 = glm(treatment~income_per_cap+log(average_pop+1)+log(average_num_jobs+1)+walk+percent_minority+percent_employed+percent_poverty+percent_transit+mean_commute+professional+office, family="binomial",data=sub_data)
summary(glm.3)
# Also tested polynomials for most significant variables and interactions

range.subway = quantile(sub_data$subway_dist, probs = seq(0.25, 0.75, by= 0.01)) # decile
output.subway = sweep.boot(range.subway, sub_data, sub_data$subway_dist,"treatment~income_per_cap+average_pop+average_num_jobs+walk+percent_minority+percent_employed+percent_poverty+percent_transit+mean_commute+professional+office")
results.subway = sweep(range.subway, sub_data, sub_data$subway_dist,"treatment~income_per_cap+average_pop+average_num_jobs+walk+percent_minority+percent_employed+percent_poverty+percent_transit+mean_commute+professional+office")
output.subway.mean = results.subway$output
subway.max.weights = results.subway$max.weights
upper.subway.ci = output.subway.mean + 1.96 * apply(output.subway,1, sd)
lower.subway.ci = output.subway.mean - 1.96 * apply(output.subway,1, sd)

par(mfrow=c(1,2), cex=1.5, lwd=1, bg=NA)
plot(range.subway, output.subway.mean, type="l", xlab='Percentile (km)', ylab='Treatment effect (trips)', ylim = c(min(lower.subway.ci), max(upper.subway.ci)), lwd=2)
lines(range.subway, upper.subway.ci, lwd=2, lty=3)
lines(range.subway, lower.subway.ci, lwd=2, lty=3)
plot(range.subway,subway.max.weights, type="l", xlab='Percentile (km)', ylab = 'Max IPW Weight', lwd=2)

##################################### TRAM ###############################################
par(mfrow=c(1,1))
# Tram distance - start by building logistic regression for 50th percentile
quant_50 = quantile(sub_data$tram_dist, 0.50)
sub_data = sub_data %>% mutate(treatment=(pmin(tram_dist)<=quant_50))

sub_data %>%
  gather(-treatment, key = "var", value = "value") %>%
  ggplot(aes(x = treatment, y = value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

glm.1 = glm(treatment~.-rail_dist-subway_dist-ferry_dist-tram_dist-bus_dist-average_ridership, family="binomial",data=sub_data)
summary(glm.1)
glm.1.probs = predict(glm.1,type="response")
glm.1.pred = (glm.1.probs > 0.5)
table(glm.1.pred, sub_data$treatment)

roc.1 = roc(sub_data$treatment,glm.1.probs)
auc(roc.1) # really good AUC
ggroc(data=roc.1)

glm.2 = glm(treatment~income_per_cap+average_pop+average_num_jobs+total_pop+total_num_jobs+walk+percent_minority+percent_employed+percent_poverty+percent_transit+mean_commute+professional+office, family="binomial",data=sub_data)
summary(glm.2)
glm.2.probs = predict(glm.2,type="response")
glm.2.pred = (glm.2.probs > 0.5)
table(glm.2.pred, sub_data$treatment)
lrtest(glm.1,glm.2)

plot(glm.2) 

roc.2 = roc(sub_data$treatment,glm.2.probs)
auc(roc.2) # still very good AUC
ggroc(data=roc.2)

# Polynomials seem to improve
glm.3 = glm(treatment~income_per_cap+poly(average_pop,degree=2)+poly(average_num_jobs,degree=2)+poly(total_pop,degree=2)+poly(total_num_jobs,degree=2)+poly(walk,degree=2)+percent_minority+percent_employed+percent_poverty+percent_transit+poly(mean_commute,degree=2)+professional+office, family="binomial",data=sub_data)
summary(glm.3)
glm.3.probs = predict(glm.3,type="response")
glm.3.pred = (glm.3.probs > 0.5)
table(glm.3.pred, sub_data$treatment)
lrtest(glm.2,glm.3)

roc.3 = roc(sub_data$treatment,glm.3.probs)
auc(roc.3) 
ggroc(data=roc.3)

glm.4 = glm(treatment~income_per_cap+poly(average_pop,degree=2)+poly(total_num_jobs,degree=2)+poly(walk,degree=2)+percent_minority+percent_employed+percent_poverty+percent_transit+poly(mean_commute,degree=2)+professional+office, family="binomial",data=sub_data)
summary(glm.4)
lrtest(glm.3,glm.4)

# Interactions did not prove useful

#range.tram = quantile(sub_data$tram_dist, probs = seq(0.25, 0.75, by= 0.01)) # decile
#output.tram = sweep.boot(range.tram, sub_data, sub_data$tram_dist,"treatment~income_per_cap+poly(average_pop,degree=2)+poly(average_num_jobs,degree=2)+poly(total_pop,degree=2)+poly(total_num_jobs,degree=2)+poly(walk,degree=2)+percent_minority+percent_employed+percent_poverty+percent_transit+poly(mean_commute,degree=2)+professional+office")
#results.tram = sweep(range.tram, sub_data, sub_data$tram_dist,"treatment~income_per_cap+poly(average_pop,degree=2)+poly(average_num_jobs,degree=2)+poly(total_pop,degree=2)+poly(total_num_jobs,degree=2)+poly(walk,degree=2)+percent_minority+percent_employed+percent_poverty+percent_transit+poly(mean_commute,degree=2)+professional+office")
#output.tram.mean = results.tram$output
#tram.max.weights = results.tram$max.weights
#upper.tram.ci = output.tram.mean + 1.96 * apply(output.tram,1, sd)
#lower.tram.ci = output.tram.mean - 1.96 * apply(output.tram,1, sd)
#plot(range.tram, output.tram.mean, type="l", xlab='Percentile (km)', ylab='Treatment effect (trips)', main='Tram', ylim = c(min(lower.tram.ci), max(upper.tram.ci)))

par(mfrow=c(1,2), cex=1.5, lwd=1, bg=NA)
plot(range.tram, output.tram.mean, type="l", xlab='Percentile (km)', ylab='Treatment effect (trips)', ylim = c(min(lower.tram.ci), max(upper.tram.ci)), lwd=2)
lines(range.tram, upper.tram.ci, lwd=2, lty=3)
lines(range.tram, lower.tram.ci, lwd=2, lty=3)
plot(range.tram,tram.max.weights, type="l", xlab='Percentile (km)', ylab = 'Max IPW Weight', lwd=2)

# This shows that we cannot get useful information - logistic regression is too good and we do not 
# satisfy positivity

##################################### BUS AND FERRY ###############################################
par(mfrow=c(1,1))
# Bus distances
hist(sub_data$bus_dist)
quant_50 = quantile(sub_data$bus_dist, 0.50) ## Really small!

# Ferry distances
hist(sub_data$ferry_dist)
quant_50 = quantile(sub_data$ferry_dist, 0.50) ## Really huge!

# I don't think we can justify a causal effect for either

##################################### COMBO TRANSPORTATION ###############################################

sub_data = sub_data %>% mutate(treatment=(pmin(tram_dist,subway_dist,rail_dist)<=0.402))
summary(sub_data$treatment) # Not bad breakdown

sub_data %>%
  gather(-treatment, key = "var", value = "value") %>%
  ggplot(aes(x = treatment, y = value)) +
  geom_boxplot() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

glm.1 = glm(treatment~.-rail_dist-subway_dist-ferry_dist-tram_dist-bus_dist-average_ridership, family="binomial",data=sub_data)
summary(glm.1)
glm.1.probs = predict(glm.1,type="response")
glm.1.pred = (glm.1.probs > 0.5)
table(glm.1.pred, sub_data$treatment)

roc.1 = roc(sub_data$treatment,glm.1.probs)
auc(roc.1) # pretty good AUC
ggroc(data=roc.1)

glm.2 = glm(treatment~walk+percent_transit+mean_commute+income_per_cap+average_pop+average_num_jobs, family="binomial", data=sub_data)
summary(glm.2)
lrtest(glm.1,glm.2)
glm.2.probs = predict(glm.2,type="response")
glm.2.pred = (glm.2.probs > 0.5)
table(glm.2.pred, sub_data$treatment)

roc.2 = roc(sub_data$treatment,glm.2.probs)
auc(roc.2) # pretty good AUC
ggroc(data=roc.2)

glm.3 = glm(treatment~poly(walk,degree=2)+poly(percent_transit,degree=2)+mean_commute+income_per_cap+average_pop+average_num_jobs, family="binomial", data=sub_data)
summary(glm.3)
lrtest(glm.2,glm.3)
# Interactions not useful
glm.3.probs = predict(glm.3,type="response")
glm.3.pred = (glm.3.probs > 0.5)
table(glm.3.pred, sub_data$treatment)

roc.3 = roc(sub_data$treatment,glm.3.probs)
auc(roc.3) # pretty good AUC
ggroc(data=roc.3)

range.trans = seq(0.2, 0.804, by= 0.01) 
output.trans = sweep.boot(range.trans, sub_data, pmin(sub_data$tram_dist,sub_data$subway_dist,sub_data$rail_dist),"treatment~poly(walk,degree=2)+poly(percent_transit,degree=2)+mean_commute+income_per_cap+average_pop+average_num_jobs")
results.trans = sweep(range.trans, sub_data, pmin(sub_data$tram_dist,sub_data$subway_dist,sub_data$rail_dist),"treatment~poly(walk,degree=2)+poly(percent_transit,degree=2)+mean_commute+income_per_cap+average_pop+average_num_jobs")
output.trans.mean = results.trans$output
trans.max.weights = results.trans$max.weights
upper.trans.ci = output.trans.mean + 1.96 * apply(output.trans,1, sd)
lower.trans.ci = output.trans.mean - 1.96 * apply(output.trans,1, sd)

par(mfrow=c(1,2), cex=1.5, lwd=1, bg=NA)
plot(range.trans, output.trans.mean, type="l", xlab='Percentile (km)', ylab='Treatment effect (trips)', ylim = c(min(lower.trans.ci), max(upper.trans.ci)), lwd=2)
lines(range.trans, upper.trans.ci, lwd=2, lty=3)
lines(range.trans, lower.trans.ci, lwd=2, lty=3)
plot(range.trans,trans.max.weights, type="l", xlab='Percentile (km)', ylab = 'Max IPW Weight', lwd=2)

#write.csv(range.trans, 'Data/IPW/boston_combo_IPW_range.csv')
#write.csv(output.trans, 'Data/IPW/boston_combo_IPW_sweep.csv')
#write.csv(results.trans, 'Data/IPW/boston_combo_IPW_res.csv')
#write.csv(upper.trans.ci, 'Data/IPW/boston_combo_IPW_upper.csv')
#write.csv(lower.trans.ci, 'Data/IPW/boston_combo_IPW_lower.csv')
