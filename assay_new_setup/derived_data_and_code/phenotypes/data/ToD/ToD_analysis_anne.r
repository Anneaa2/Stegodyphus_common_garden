setwd("/home/anneaa/spider2/faststorage/assay_study_anne/phenotypes/data/ToD/")
library(agricolae)

#Plot is further below

# temp of death
ToD_raw <- read.table("Complete_data_file_ToD_270218.txt", header = TRUE, sep = "", dec = ',')
	#write.table(ToD_raw,"Complete_data_file_ToD_270218A.txt",sep="\t",row.names=F,col.names=T,quote=F) # correct nest names
ToD <- ToD_raw[!grepl("E001", ToD_raw$Nest),] # remove nest from data since... ####Nest E001 has very low values in all acclimation temperatures, and is removed here########

modelToD <- lm(ToD ~ Population + Treatment + Population *Treatment, data=ToD)
summary(modelToD)





library(lmtest)

modelToD_final_pop <- lm(ToD ~ Population, data=ToD)
modelToD_final_temp <- lm(ToD ~ Treatment, data=ToD)
modelToD_final_add <- lm(ToD ~ Population + Treatment, data=ToD)
modelToD_final <- lm(ToD ~ Population + Treatment + Population * Treatment, data=ToD)
modelToD_final_M <- lm(ToD ~ Population + Treatment + Mass + Population * Treatment, data=ToD)
modelToD_final_MF <- lm(ToD ~ Population + Treatment + Mass + Fed + Population * Treatment, data=ToD)
modelToD_final_MFD <- lm(ToD ~ Population + Treatment + Mass + Fed + D_21 + Population * Treatment, data=ToD)

pdf("modeltest_raw_dat_full_model.pdf",width=14,height=10)
par(mfrow=c(4,4))
plot(modelToD_final_pop,main="pop")
plot(modelToD_final_temp,main="temp")
plot(modelToD_final_add,main="additive")
plot(modelToD_final,main="add+interaction")
dev.off()

pdf("modeltest_raw_dat_full_model_MFD.pdf",width=14,height=10)
par(mfrow=c(4,4))
plot(modelToD_final,main="add+interaction")
plot(modelToD_final_M,main="full+mass")
plot(modelToD_final_MF,main="full+mass+fed")
plot(modelToD_final_MFD,main="full+mass+fed+d21")
dev.off()


summary(modelToD_final)
'Call:
lm(formula = ToD ~ Population + Treatment + Population * Treatment,
    data = ToD)

Residuals:
    Min      1Q  Median      3Q     Max
-6.8358 -0.5820  0.0195  0.5672  4.7499

Coefficients:
                              Estimate Std. Error t value Pr(>|t|)
(Intercept)                   48.28377    0.16383 294.720  < 2e-16 ***
PopulationKarasburg            1.64401    0.26219   6.270 4.10e-10 ***
PopulationOtavi                0.53286    0.23072   2.310   0.0210 *
PopulationStampriet            1.22447    0.23920   5.119 3.26e-07 ***
Treatment                      0.09875    0.00713  13.850  < 2e-16 ***
PopulationKarasburg:Treatment -0.08705    0.01147  -7.589 4.23e-14 ***
PopulationOtavi:Treatment     -0.02413    0.01009  -2.391   0.0169 *
PopulationStampriet:Treatment -0.05683    0.01060  -5.363 8.76e-08 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.022 on 3141 degrees of freedom
  (511 observations deleted due to missingness)
Multiple R-squared:  0.1053,    Adjusted R-squared:  0.1033
F-statistic: 52.81 on 7 and 3141 DF,  p-value: < 2.2e-16
'
summary(modelToD_final_MFD)
'
Call:
lm(formula = ToD ~ Population + Treatment + Mass + Fed + D_21 +
    Population * Treatment, data = ToD)

Residuals:
    Min      1Q  Median      3Q     Max
-6.6747 -0.5901  0.0002  0.5859  4.7007

Coefficients:
                                Estimate Std. Error t value Pr(>|t|)
(Intercept)                   48.0096348  0.1712739 280.309  < 2e-16 ***
PopulationKarasburg            1.5910929  0.2629724   6.050 1.63e-09 ***
PopulationOtavi                0.5330953  0.2357226   2.262 0.023800 *
PopulationStampriet            1.2894807  0.2572158   5.013 5.68e-07 ***
Treatment                      0.1035439  0.0073749  14.040  < 2e-16 ***
Mass                          -0.0071516  0.0029113  -2.457 0.014088 *
Fed                            0.0577308  0.0175210   3.295 0.000996 ***
D_21                           0.0032503  0.0009998   3.251 0.001164 **
PopulationKarasburg:Treatment -0.0844533  0.0114370  -7.384 2.00e-13 ***
PopulationOtavi:Treatment     -0.0205297  0.0102966  -1.994 0.046264 *
PopulationStampriet:Treatment -0.0645528  0.0113071  -5.709 1.25e-08 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.018 on 2880 degrees of freedom
  (769 observations deleted due to missingness)
Multiple R-squared:  0.1191,    Adjusted R-squared:  0.1161
F-statistic: 38.95 on 10 and 2880 DF,  p-value: < 2.2e-16
'
	###---------- No Effect of D21
modelToDD_21 <- lm(ToD ~ D_21, data=ToD)
summary(modelToDD_21)
'
Call:
lm(formula = ToD ~ D_21, data = ToD)

Residuals:
    Min      1Q  Median      3Q     Max
-7.3191 -0.7180  0.0807  0.6817  4.4817

Coefficients:
              Estimate Std. Error  t value Pr(>|t|)
(Intercept)  5.042e+01  3.477e-02 1450.074   <2e-16 ***
D_21        -1.722e-05  6.791e-04   -0.025     0.98
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.08 on 3147 degrees of freedom
  (511 observations deleted due to missingness)
Multiple R-squared:  2.042e-07, Adjusted R-squared:  -0.0003176
F-statistic: 0.0006427 on 1 and 3147 DF,  p-value: 0.9798
'
	###---------- Effect of mass
modelToD_mass <- lm(ToD ~ Mass, data=ToD)
summary(modelToD_mass)
'Call:
lm(formula = ToD ~ Mass, data = ToD)

Residuals:
    Min      1Q  Median      3Q     Max
-7.3037 -0.6852  0.0422  0.7048  4.5294

Coefficients:
             Estimate Std. Error t value Pr(>|t|)
(Intercept) 50.292681   0.040968 1227.61  < 2e-16 ***
Mass         0.008185   0.002496    3.28  0.00105 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.081 on 2889 degrees of freedom
  (769 observations deleted due to missingness)
Multiple R-squared:  0.003709,  Adjusted R-squared:  0.003364
F-statistic: 10.76 on 1 and 2889 DF,  p-value: 0.001052
'
	###---------- Effect of Feeding
modelToD_fed <- lm(ToD ~ Fed, data=ToD)
summary(modelToD_fed)
'Call:
lm(formula = ToD ~ Fed, data = ToD)

Residuals:
    Min      1Q  Median      3Q     Max
-7.1901 -0.6901  0.0099  0.7003  4.5051

Coefficients:
            Estimate Std. Error  t value Pr(>|t|)
(Intercept) 50.18530    0.04108 1221.698  < 2e-16 ***
Fed          0.10479    0.01634    6.415 1.62e-10 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.073 on 3147 degrees of freedom
  (511 observations deleted due to missingness)
Multiple R-squared:  0.01291,   Adjusted R-squared:  0.01259
F-statistic: 41.15 on 1 and 3147 DF,  p-value: 1.619e-10
'


lrtest(modelToD_final_temp)
lrtest(modelToD_final_pop)
lrtest(modelToD_final_temp, modelToD_final_add)
lrtest(modelToD_final_pop, modelToD_final_add)
lrtest(modelToD_final_add, modelToD_final)




modelToDN_D_21_ <- lm(ToD ~ D_21, data=ToD)
summary(modelToDN_D_21_)
'
Call:
lm(formula = ToD ~ D_21, data = ToD)

Residuals:
    Min      1Q  Median      3Q     Max
-7.3191 -0.7180  0.0807  0.6817  4.4817

Coefficients:
              Estimate Std. Error  t value Pr(>|t|)
(Intercept)  5.042e+01  3.477e-02 1450.074   <2e-16 ***
D_21        -1.722e-05  6.791e-04   -0.025     0.98
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.08 on 3147 degrees of freedom
  (511 observations deleted due to missingness)
Multiple R-squared:  2.042e-07, Adjusted R-squared:  -0.0003176
F-statistic: 0.0006427 on 1 and 3147 DF,  p-value: 0.9798
'





##-----------reduced data ("ToD_reduced") --> getting mean for each nest/temp group

index <- ToD[,4]==1
ToD_mean <- colMeans(matrix(ToD$ToD,10),na.rm=T)
Mass_mean <- colMeans(matrix(ToD$Mass,10),na.rm=T)
Fed_mean <- colMeans(matrix(ToD$Fed,10),na.rm=T)
D_21_mean <- colMeans(matrix(ToD$D_21,10),na.rm=T)
ToD_reduced <- cbind(ToD[index,],ToD_mean,Mass_mean,Fed_mean,D_21_mean)

####Nest E001 has very low values in all acclimation temperatures, and is removed here########

ToD_final <- ToD_reduced[!grepl("E001", ToD_reduced$Nest),]
write.csv(ToD_final, file="/faststorage/project/spider2/assay_study_anne/phenotypes/data/ToD/ToD_nest-temp_mean.csv") 

# Plots with temperature
pdf("effects_of_temp_pop.pdf")
par(mfcol=c(2,2))
plot(ToD_final$Treatment, ToD_final$ToD_mean,col="grey", ylab="CTmax",xlab="Acclimation temperature", main="Effect of temperature")#, ylim=c(44,54))
abline(lm(ToD_final$ToD_mean ~ ToD_final$Treatment),col="grey")
plot(ToD_final$Treatment, ToD_final$Mass_mean, col="red", ylab="Mass",xlab="Acclimation temperature") #, ylim=c(44,54)
abline(lm(ToD_final$Mass_mean ~ ToD_final$Treatment),col="red")
#plots with pop
plot(ToD_final$Population, ToD_final$ToD_mean,col="grey", ylab="CTmax",xlab="Population", main="Effect of Population")
#summary(lm(ToD_final$ToD_mean ~ ToD_final$Population)) only K different to B. Ie. B,O and S --> a. K --> b
plot(ToD_final$Population, ToD_final$Mass_mean, col="red",ylab="Mass",xlab="Population")
#summary(lm(ToD_final$Mass_mean ~ ToD_final$Population)) #test for sign groups
	# no need to plot Fed abd Daysat21
dev.off()
#library(agricolae)
HSD.test(aov(Mass_mean ~ Population,data=ToD_final), "tx", group=TRUE)


modelToD_final <- lm(ToD_mean ~ Population + Treatment + Population * Treatment, data=ToD_final)
summary(modelToD_final)
ToDaov <- aov(ToD_mean ~ Population + Treatment + Population*Treatment, data=ToD_final)
summary(ToDaov)


#---------- Log likelihood ratio tests
modelToD_final <- lm(ToD_mean ~ Population + Treatment + Population * Treatment, data=ToD_final)
modelToD_final_add <- lm(ToD_mean ~ Population + Treatment, data=ToD_final)
modelToD_final_pop <- lm(ToD_mean ~ Population, data=ToD_final)
modelToD_final_temp <- lm(ToD_mean ~ Treatment, data=ToD_final)

library(lmtest)

lrtest(modelToD_final_temp)
lrtest(modelToD_final_pop)
lrtest(modelToD_final_temp, modelToD_final_add)
lrtest(modelToD_final_pop, modelToD_final_add)
lrtest(modelToD_final_add, modelToD_final)




modelToD_finalMonl <- lm(ToD_mean ~ Mass_mean, data=ToD_final)
summary(modelToD_finalMonl)	
	#effect of mass
modelToD_finalM <- lm(ToD_mean ~ Population + Treatment + Mass_mean + Population * Treatment, data=ToD_final)
summary(modelToD_finalM)	
	# additional effect of mass	(mean of mass)
modelToD_finalFonl <- lm(ToD_mean ~ Fed_mean, data=ToD_final)
summary(modelToD_finalFonl)	
	# effect of feeding
modelToD_finalF <- lm(ToD_mean ~ Population + Treatment + Fed_mean + Population * Treatment, data=ToD_final)
summary(modelToD_finalF)
	# additional effect of feeding
modelToD_finalDonl <- lm(ToD_mean ~ D_21_mean, data=ToD_final)
summary(modelToD_finalDonl)	
	# No effect of days at 21
modelToD_finalD <- lm(ToD_mean ~ Population + Treatment + D_21_mean + Population * Treatment, data=ToD_final)
summary(modelToD_finalD)
	# additional effect of days at 21
	
	
modelToD_final_ALL <- lm(ToD_mean ~ Population + Treatment + Mass_mean + Fed_mean + D_21_mean + Population * Treatment, data=ToD_final) # reducing dataframe to only use samples included in this model
modelToD_final_ALLM <- lm(ToD_mean ~ Population + Treatment + Mass_mean + Population * Treatment, data=modelToD_final_ALL$model)
modelToD_final_ALLF <- lm(ToD_mean ~ Population + Treatment + Fed_mean + Population * Treatment, data=modelToD_final_ALL$model)
modelToD_final_ALLD <- lm(ToD_mean ~ Population + Treatment + D_21_mean + Population * Treatment, data=modelToD_final_ALL$model)
modelToD_final <- lm(ToD_mean ~ Population + Treatment + Population * Treatment, data=modelToD_final_ALL$model)

lrtest(modelToD_final,modelToD_final_ALLM) 
	# additional effect of Mass
lrtest(modelToD_final,modelToD_final_ALLF) 
	# additional effect of Feeding
lrtest(modelToD_final,modelToD_final_ALLD) 
	# additional effect of days at 21

