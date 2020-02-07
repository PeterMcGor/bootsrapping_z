library(mgcv)
library(ggplot2)
library(dplyr)
empirical_confidence_interval <- function(statistic, alpha = 0.95){
  lower = quantile(statistic, max(0,((1.0-alpha)/2.0)))[[1]]
  upper = quantile(statistic, min((alpha+((1.0-alpha)/2.0))))[[1]]
  return (c(lower, upper))
}

#when se can be estimated as the standard deviation. To our pupose its the same to scale it 
get_z_se <- function(data_w_preds) {
  data_w_preds$z_se =  (data_w_preds$tcv - data_w_preds$fit) / data_w_preds$se.fit
  return(data_w_preds)
}


#the variance of the posterior distribution
get_z_frm_post_sd <- function(data_w_preds, gam_model) {
  data_w_preds$z_sd_post = (data_w_preds$tcv - data_w_preds$fit) / sqrt(gam_model$sig2)
  return(data_w_preds)
}

#Andre uses for the expected variance all the data (training and test), note really sure about how correct is that
get_z_pseudoandre <- function(data_w_preds, expected_sd) {
  data_w_preds$z_an = (data_w_preds$tcv - data_w_preds$fit) / sqrt(expected_sd**2*data_w_preds$se.fit) 
  return(data_w_preds)
}

get_all_z <- function(data_w_preds, gam_model, expected_sd) {
  return (data_w_preds %>% get_z_se() %>% get_z_frm_post_sd(gam_model) %>% get_z_pseudoandre(expected_sd) )
}

setwd("/home/pmacias/Projects/JanJo/bootsrapping_z/")

thedataF_pilar2 = read.csv("thedataF_pilar2.csv")
thedataM_pilar2 = read.csv("thedataM_pilar2.csv")


######### MODEL TCV #############
# tiene sentido meter el participant como un efecto más??
a = gam(tcv ~s(age, bs = "cs", k =4) +acode+s(participant, bs = "re"), method = "ML", data =thedataF_pilar2)
#Due to the number of mesures per patients the model does not really change much. It would be more correct the following (I think...)
#a' = gam(tcv ~s(age, participant, bs = "fs", k =4) +acode, method = "ML", data =thedataF_pilar2)
# or
#a' = gam(tcv ~ age + s(age, by = participant) +acode, method = "ML", data =thedataF_pilar2)
# I cannot applied a' due to the number of observations. Something in the middle woud be it similar to a 
b= gam(tcv ~ ti(age, participant, bs = c("cs", "re"), k = c(3,3)) + acode, method = "ML", data =thedataF_pilar2)

c = gam(tcv ~s(age, bs = "cs", k =4), method = "ML", data =thedataF_pilar2)

a_predictions_raw = predict(a, se.fit =T)  
a_predictions_corrected = predict(a, type = "response", se.fit =T, exclude = "s(participant)")  
b_predictions_raw = predict(b, se.fit =T)  
b_predictions_corrected = predict(b, type = "response", se.fit =T, exclude = "s(participant)")  

thedataF_pilar2 = within(thedataF_pilar2, {
  fit_a = a_predictions_raw$fit
  fit.se = a_predictions_raw$se.fit
  fit_a_corr = a_predictions_corrected$fit
  fit_a_corr.se = a_predictions_corrected$se.fit
  fit_b = b_predictions_raw$fit
  fit_b.se = b_predictions_raw$se.fit
  fit_b_corr = b_predictions_corrected$fit
  fit_b_corr.se = b_predictions_corrected$se.fit
})


#Model comparison
summary(a)
summary(b)
# sumaries reveals a slightly better r^2 for and slightly worse deviance explained
anova(a,b) #No differences
AIC(a,b) #The samllest the best fiting (overfitting could occurs)

#This visualization clear the overfitting problem when participant is included. 
# The number of unique subjects is not enough for modeling this effect 
#(we cannot even employ the needed formulas due to the unique number of subjects)
# Actually the c model (the simplest, in yellow) points to an asyntotically behaviour which probably makes more sense than the a model
# when acode is included, the curve shows some small overfitting problems
# WARNING: this does not mean that the model with acode and participant is incorrect. 
# Hopefully, a model with 10 or more repeated measures would work better than the simplest.
ggplot(thedataF_pilar2, aes(x = age, y = tcv, color = acode)) + 
  geom_point(show.legend = F) + 
  geom_line(aes(y = fit_a), color = female_color, alpha = 0.3) + 
  geom_line(aes(y = fit_a_corr), color = "green", alpha = 0.8) + 
  geom_line(aes(y = fit_a), color = "darkgrey", alpha = 0.3) +  
  geom_line(aes(y = fit_b_corr), color = "purple", alpha = 0.3)+ 
  geom_line(aes(y = predict(c)), color = male_color, alpha = 0.8)



thedataM_pilar2$participant <- rep(thedataF_pilar2$subID[1], nrow(thedataM_pilar2)) # esto para qué?
females_ss = thedataF_pilar2[c("age", "tcv", "acode", "participant")]
males_ss = thedataM_pilar2[c("age", "tcv", "acode", "participant")]
pred_vals_thedataM_pilar2 <- as.data.frame(predict(a, newdata =males_ss,type = "response", se.fit =T, exclude = "s(participant)"))

### SE vs sigma
## That standard error is the uncertainty in the estimate for the expected/fitted value. 
##### the variance of the posterior distribution, you want the estimated value of the variance (or scale parameter) for the model you fitted
#Creo que esto no tiene mucho sentido. Para tener la sigma, habría que plantear un modelo heterocedastico (una sigma diferente point/bin wise)
# Se podría calcular el intervalo de confianza con el SE como mu +- 1.96 fit.se
## Cuando Andre saca las Zs qué sigma utiliza??
#con el GAM la sigma es la misma para todos los puntos

#I use the expected variance of the training data
pred_vals_thedataM_pilar2 = get_all_z(cbind(males_ss, pred_vals_thedataM_pilar2), gam_model = a,
                                      expected_sd = sd(thedataF_pilar2$tcv))


# Fit models to 100 bootstrap replicates of the data
set.seed(69)
n_boots = 100
predictions = replicate(
  n_boots,{
    boot = females_ss [sample.int(nrow(females_ss), replace = T), ]
    boot = boot[c("age", "tcv", "acode", "participant")]
    males_ss$participant <- rep(boot$participant[1], nrow(males_ss))
    model = gam(tcv ~ s(age, bs="cs", k=4) +acode +s(participant, bs = "re"), method = "ML", data = boot)

    exp_sd = sd(boot$tcv) #Andre employs the whole sample not just training...
    type = "response"
    predictions_females = as.data.frame(predict(model, type = type, se.fit =T, exclude = "s(participant)"))
    predictions_females = get_all_z(cbind(boot, predictions_females), gam_model = model, expected_sd = exp_sd)
    
    predictions_males = as.data.frame(predict(model, newdata =males_ss, type = type, se.fit =T, exclude = "s(participant)"))
    predictions_males = get_all_z(cbind(males_ss, predictions_males), gam_model = model,expected_sd = exp_sd )
    return(list("predictions_females"=predictions_females, "predictions_males"=predictions_males))
  }
)


female_model = predictions[,1]$predictions_females
female_model$replicate = 1
male_model = predictions[,1]$predictions_males
male_model$replicate = 1
#female_predictions = as.data.frame(predictions[,1]$predictions_females)
#female_predictions$replicate = 1
for (i in 2:dim(predictions)[2]) {
  female_ss = predictions[,i]$predictions_females
  male_ss = predictions[,i]$predictions_males
  
  female_ss$replicate = i
  male_ss$replicate = i
  
  female_model = rbind(female_model, female_ss)
  male_model = rbind(male_model, male_ss)
}

female_color = "royalblue1"
male_color = "yellow"

#Just tcv fitting...must be the same
ggplot() +
  geom_line(data=female_model, aes(x = age, y = fit, group = replicate), color = female_color) +
  geom_line(data=female_model, aes(x = age, y = fit, group = replicate), color = male_color) + 
  geom_point(data = thedataF_pilar2, aes(x=age, y=tcv), color = female_color) +
  geom_point(data = thedataM_pilar2, aes(x=age, y=tcv), color = male_color)

#visualize Z differences betwwen sexs
ggplot(data = female_model, aes(x = age, y=z_an)) + 
  geom_point(color = female_color) + 
  geom_point(data = male_model, color = male_color) + 
  geom_smooth(color = female_color) + 
  geom_smooth(data = male_model,  color = male_color)

ggplot(thedataF_pilar2, aes(x=age, y =tcv)) + 
  geom_smooth(color = female_color) + 
  geom_point(color = female_color) + 
  geom_smooth(data=thedataM_pilar2, color = male_color) + 
  geom_point(data=thedataM_pilar2, color = male_color) + 
  geom_smooth(data = female_model, aes(y=fit), color = female_color, linetype = "dashed") +
  geom_smooth(data = male_model, aes(y=fit), color = male_color, linetype = "dashed")

male_model_ci = aggregate(z_an ~ age, FUN = empirical_confidence_interval, data=male_model) 
# aggragattes return a matrix
male_model_ci = cbind(male_model_ci$age, as.data.frame(male_model_ci$z_an) ) %>% setNames(c("age", "z_an_lower", "z_an_upper"))
pred_males_w_bootstrap = merge(pred_vals_thedataM_pilar2[!duplicated(pred_vals_thedataM_pilar2$age), ], male_model_ci, by = "age")
pred_males_w_bootstrap$z_an_in_ci <- 
  pred_males_w_bootstrap$z_an > pred_males_w_bootstrap$z_an_lower & pred_males_w_bootstrap$z_an < pred_males_w_bootstrap$z_an_upper

message("percentage of z inside the interval ", 100 *sum(pred_males_w_bootstrap$z_an_in_ci) / nrow(pred_males_w_bootstrap))

#Warning. The points are the z values when all samples are employed. the error bar, the ones obtained with the boostrapping 
ggplot(data = pred_males_w_bootstrap, aes(x = age)) + 
  geom_point(aes(y = z_an, color = factor(acode) )) + 
  geom_errorbar(aes(ymin = z_an_lower, ymax = z_an_upper, color = factor(acode)), width = 0.4) #+ 
  #geom_line(data = male_model,  color = male_color, aes(group=replicate, y=z_an))




 