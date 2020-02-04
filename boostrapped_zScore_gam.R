library(mgcv)
library(ggplot2)
library(itsadug)

empirical_confidence_interval <- function(statistic, alpha = 0.95){
  lower = quantile(statistic, max(0,((1.0-alpha)/2.0)))[[1]]
  upper = quantile(statistic, min((alpha+((1.0-alpha)/2.0))))[[1]]
  return (c(lower, upper))
}

setwd("/home/pmacias/Projects/JanJo/")

thedataF_pilar2 = read.csv("thedataF_pilar2.csv")
thedataM_pilar2 = read.csv("thedataM_pilar2.csv")


######### MODEL TCV #############
# tiene sentido meter el participant como un efecto más??
a = gam(tcv ~s(age, bs = "cs", k =4) +acode+s(participant, bs = "re"), method = "ML", data =thedataF_pilar2)
#Due to the number of mesures per patients the model does not really change much. It would be more correct the following (I think...)
#a' = gam(tcv ~s(age, participant, bs = "fs", k =4) +acode, method = "ML", data =thedataF_pilar2)
# or
#a' = gam(tcv ~ age + s(age, by = participant) +acode, method = "ML", data =thedataF_pilar2)
# I cannot applied a' due to the number of observations. Something in the middle woud be: 
b= gam(tcv ~ ti(age, participant, bs = c("cs", "re"), k = c(3,3)) + acode, method = "ML", data =thedataF_pilar2)

#Model comparison
summary(a)
summary(b)
# sumaries reveals a slightly better r^2 for and slightly worse deviance explained
anova(a,b) #No differences
AIC(a,b) #The samllest the best fiting (overfitting could occurs)

# visualize. ggplot and gam models are not really freindly together (itsadug library kind of ease the process)
pred_a = get_predictions(a, 
                         cond = list(age = seq(min(thedataF_pilar2$age, na.rm=T),
                                             max(thedataF_pilar2$age, na.rm=T), 
                                             length.out = nrow(thedataF_pilar2)), se=T))

pred_b = get_predictions(b, 
                         cond = list(age = seq(min(thedataF_pilar2$age, na.rm=T),
                                               max(thedataF_pilar2$age, na.rm=T), 
                                               length.out = nrow(thedataF_pilar2)), se=T))


pred_a$tcv  = pred_a$fit #create dummy tcv for consistance
pred_b$tcv = pred_b$fit

ggplot(thedataF_pilar2, aes(x = age, y = tcv, acode = acode, participant = participant)) + 
  geom_line(alpha=.3,aes(group=participant)) + 
  geom_point(alpha=.3) +
  geom_ribbon(data=pred_a, alpha=.6, aes(ymin=fit-CI, ymax=fit+CI), show.legend = F, fill='forestgreen') +  
  geom_line(data = pred_a, show.legend = F, color='forestgreen') + 
  geom_ribbon(data=pred_b, alpha=.2, aes(ymin=fit-CI, ymax=fit+CI), show.legend = F, fill='royalblue1') +  
  geom_line(data = pred_b, show.legend = F, color='royalblue1')


thedataM_pilar2$participant <- rep(thedataF_pilar2$subID[1], nrow(thedataM_pilar2)) # esto para qué?
pred_vals_thedataM_pilar2 <- as.data.frame(predict(a, newdata =thedataM_pilar2,type = "response", se.fit =T, exclude = "s(participant)"))

### SE vs sigma
## That standard error is the uncertainty in the estimate for the expected/fitted value. 
##### the variance of the posterior distribution, you want the estimated value of the variance (or scale parameter) for the model you fitted
#Creo que esto no tiene mucho sentido. Para tener la sigma, habría que plantear un modelo heterocedastico (una sigma diferente point/bin wise)
# Se podría calcular el intervalo de confianza con el SE como mu +- 1.96 fit.se
## Cuando Andre saca las Zs qué sigma utiliza??
pred_vals_thedataM_pilar2$sd.fit <- pred_vals_thedataM_pilar2$se.fit*sqrt(311) 
#con el GAM la sigma es la misma para todos los puntos
pred_vals_thedataF_pilar2$real_sd.fit = sqrt(a$sig2)

pred_vals_thedataM_pilar2$Z <- (thedataM_pilar2$tcv - pred_vals_thedataM_pilar2$fit) / pred_vals_thedataM_pilar2$sd.fit
pred_vals_thedataM_pilar2$Z_sig <- (thedataM_pilar2$tcv - pred_vals_thedataM_pilar2$fit) / sqrt(a$sig2)

z = as.data.frame(pred_vals_thedataM_pilar2$Z)
z_sig = as.data.frame(pred_vals_thedataM_pilar2$Z_sig)



# Fit models to 100 bootstrap replicates of the data
predictions = replicate(
  10,{
    boot = thedataF_pilar2[sample.int(nrow(thedataF_pilar2), replace = TRUE), ]
    thedataM_pilar2$participant <- rep(boot$subID[1], nrow(thedataM_pilar2))
    model = gam(tcv ~s(age, bs = "cs", k =4) +acode +s(participant, bs = "re"), method = "ML", data =boot)
    # Output predictions at each point that we'll want to plot later
    predictions = as.data.frame(predict(model,newdata =thedataM_pilar2,type = "response", se.fit =T, exclude = "s(participant)"))
    #predictions$real_sd.fit =  sqrt(model$sig2)
    return (thedataM_pilar2$tcv - predictions$fit) / sqrt(model$sig2)
  }
)
df <- data.frame(matrix(unlist(predictions), ncol=length(predictions), byrow=F))
#write.csv(df, "./Documentos/GM/sex differences docs/bootstrap_tcv.csv")
df = read.csv("./Documentos/GM/sex differences docs/bootstrap_tcv.csv")
df[1]<-NULL

# rename columns
even_indexes<-seq(2,2000,2)
odd_indexes<-seq(1,1999,2)
colnames(df)[odd_indexes] = "fit"
colnames(df)[even_indexes] = "se.fit"

# add suffix to differentiate between bootstraps
colnames(df) <- paste(colnames(df), rep(seq(1,1000), each =2), sep = "_")

# calculate sd
se = df[even_indexes]
sd = se*sqrt(389)

#calculate z
fit = df[odd_indexes]
n <- 1000
males<-as.data.frame(do.call("cbind", replicate(n, thedataM_pilar2$tcv, simplify = FALSE)))

#rename all cols to be the same
colnames(fit)<-seq(1:1000)
colnames(sd)<-seq(1:1000)
colnames(males)<-seq(1:1000)

# substract and divide to get z-scores
substraction=as.data.frame(Map("-", males, fit))
ZZ = as.data.frame(Map("/", substraction, sd))
#write.csv(ZZ, "./Documentos/GM/sex differences docs/boot_z_scores_tcv.csv")


######### MODEL cortex_sulcalwidth_brainvisa #############
a = gam(cortex_sulcalwidth_brainvisa ~s(age, bs = "cs", k =4) +acode+s(participant, bs = "re"), method = "ML", data =thedataF_pilar2)
thedataM_pilar2$participant <- rep(thedataF_pilar2$subID[1], nrow(thedataM_pilar2))
pred_vals_thedataM_pilar2 <- as.data.frame(predict(a, newdata =thedataM_pilar2,type = "response", se.fit =T, exclude = "s(participant)"))

pred_vals_thedataM_pilar2$sd.fit <- pred_vals_thedataM_pilar2$se.fit*sqrt(311)

pred_vals_thedataM_pilar2$Z <- (thedataM_pilar2$cortex_sulcalwidth_brainvisa - pred_vals_thedataM_pilar2$fit) / pred_vals_thedataM_pilar2$sd.fit
z = as.data.frame(pred_vals_thedataM_pilar2$Z)


# Fit models to 100 bootstrap replicates of the data
predictions = replicate(
  1000,{
    boot = thedataF_pilar2[sample.int(nrow(thedataF_pilar2), replace = TRUE), ]
    thedataM_pilar2$participant <- rep(boot$subID[1], nrow(thedataM_pilar2))
    model = gam(cortex_sulcalwidth_brainvisa ~s(age, bs = "cs", k =4) +acode +s(participant, bs = "re"), method = "ML", data =boot)
    # Output predictions at each point that we'll want to plot later
    predict(model,newdata =thedataM_pilar2,type = "response", se.fit =T, exclude = "s(participant)")
  }
)
df <- data.frame(matrix(unlist(predictions), ncol=length(predictions), byrow=F))
#write.csv(df, "./Documentos/GM/sex differences docs/bootstrap_cortex_sulcalwidth_brainvisa.csv")
df = read.csv("./Documentos/GM/sex differences docs/bootstrap_cortex_sulcalwidth_brainvisa.csv")
df[1]<-NULL

# rename columns
even_indexes<-seq(2,2000,2)
odd_indexes<-seq(1,1999,2)
colnames(df)[odd_indexes] = "fit"
colnames(df)[even_indexes] = "se.fit"

# add suffix to differentiate between bootstraps
colnames(df) <- paste(colnames(df), rep(seq(1,1000), each =2), sep = "_")

# calculate sd
se = df[even_indexes]
sd = se*sqrt(389)

#calculate z
fit = df[odd_indexes]
n <- 1000
males<-as.data.frame(do.call("cbind", replicate(n, thedataM_pilar2$cortex_sulcalwidth_brainvisa, simplify = FALSE)))

#rename all cols to be the same
colnames(fit)<-seq(1:1000)
colnames(sd)<-seq(1:1000)
colnames(males)<-seq(1:1000)

# substract and divide to get z-scores
substraction=as.data.frame(Map("-", males, fit))
ZZ = as.data.frame(Map("/", substraction, sd))
#write.csv(ZZ, "./Documentos/GM/sex differences docs/boot_z_scores_cortex_sulcalwidth_brainvisa.csv")




######### MODEL cortex_sulcallength_brainvisa #############
a = gam(cortex_sulcallength_brainvisa ~s(age, bs = "cs", k =4) +acode+s(participant, bs = "re"), method = "ML", data =thedataF_pilar2)
thedataM_pilar2$participant <- rep(thedataF_pilar2$subID[1], nrow(thedataM_pilar2))
pred_vals_thedataM_pilar2 <- as.data.frame(predict(a, newdata =thedataM_pilar2,type = "response", se.fit =T, exclude = "s(participant)"))

pred_vals_thedataM_pilar2$sd.fit <- pred_vals_thedataM_pilar2$se.fit*sqrt(311)

pred_vals_thedataM_pilar2$Z <- (thedataM_pilar2$cortex_sulcallength_brainvisa - pred_vals_thedataM_pilar2$fit) / pred_vals_thedataM_pilar2$sd.fit
z = as.data.frame(pred_vals_thedataM_pilar2$Z)


# Fit models to 100 bootstrap replicates of the data
predictions = replicate(
  1000,{
    boot = thedataF_pilar2[sample.int(nrow(thedataF_pilar2), replace = TRUE), ]
    thedataM_pilar2$participant <- rep(boot$subID[1], nrow(thedataM_pilar2))
    model = gam(cortex_sulcallength_brainvisa ~s(age, bs = "cs", k =4) +acode +s(participant, bs = "re"), method = "ML", data =boot)
    # Output predictions at each point that we'll want to plot later
    predict(model,newdata =thedataM_pilar2,type = "response", se.fit =T, exclude = "s(participant)")
  }
)

df <- data.frame(matrix(unlist(predictions), ncol=length(predictions), byrow=F))
write.csv(df, "./Documentos/GM/sex differences docs/bootstrap_cortex_cortex_sulcallength_brainvisa.csv")
df = read.csv("./Documentos/GM/sex differences docs/bootstrap_cortex_sulcallength_brainvisa.csv")
df[1]<-NULL

# rename columns
even_indexes<-seq(2,2000,2)
odd_indexes<-seq(1,1999,2)
colnames(df)[odd_indexes] = "fit"
colnames(df)[even_indexes] = "se.fit"

# add suffix to differentiate between bootstraps
colnames(df) <- paste(colnames(df), rep(seq(1,1000), each =2), sep = "_")

# calculate sd
se = df[even_indexes]
sd = se*sqrt(389)

#calculate z
fit = df[odd_indexes]
n <- 1000
males<-as.data.frame(do.call("cbind", replicate(n, thedataM_pilar2$cortex_sulcallength_brainvisa, simplify = FALSE)))

#rename all cols to be the same
colnames(fit)<-seq(1:1000)
colnames(sd)<-seq(1:1000)
colnames(males)<-seq(1:1000)

# substract and divide to get z-scores
substraction=as.data.frame(Map("-", males, fit))
ZZ = as.data.frame(Map("/", substraction, sd))
write.csv(ZZ, "./Documentos/GM/sex differences docs/boot_z_scores_cortex_sulcallength_brainvisa.csv")
