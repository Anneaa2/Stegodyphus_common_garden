setwd("/home/anneaa/spider2/faststorage/assay_study_anne/phenotypes/data/ToD/")
library(agricolae)

# temp of death
ToD_raw <- read.table("Complete_data_file_ToD_270218.txt", header = TRUE, sep = "", dec = ',')

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









#subset - ToD - population
ToDN <- subset(ToD, Population == "Otavi")
ToDS <- subset(ToD, Population == "Stampriet")
ToDK <- subset(ToD, Population == "Karasburg")
ToDB <- subset(ToD, Population == "Betta")


par(mfrow=c(3,4))
plot(ToDN$D_21, ToDN$ToD, xlim =c(0,100), ylim=c(44,54))
abline(lm(ToDN$ToD ~ ToDN$D_21))
plot(ToDS$D_21, ToDS$ToD, xlim =c(0,100), ylim=c(44,54))
abline(lm(ToDS$ToD ~ ToDS$D_21))
plot(ToDK$D_21, ToDK$ToD, xlim =c(0,100), ylim=c(44,54))
abline(lm(ToDK$ToD ~ ToDK$D_21))
plot(ToDB$D_21, ToDB$ToD, xlim =c(0,100), ylim=c(44,54))
abline(lm(ToDB$ToD ~ ToDB$D_21))

#models - D_21
modelToDN_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDN)
statToDN_D_21 <- summary(modelToDN_D_21)
modelToDS_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDS)
statToDS_D_21 <- summary(modelToDS_D_21)
modelToDK_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDK)
statToDK_D_21 <- summary(modelToDK_D_21)
modelToDB_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDB)
statToDB_D_21 <- summary(modelToDB_D_21)

plot(ToDN$Mass, ToDN$ToD, xlim =c(0,50), ylim=c(44,54))
abline(lm(ToDN$ToD ~ ToDN$Mass))
plot(ToDS$Mass, ToDS$ToD, xlim =c(0,50), ylim=c(44,54))
abline(lm(ToDS$ToD ~ ToDS$Mass))
plot(ToDK$Mass, ToDK$ToD, xlim =c(0,50), ylim=c(44,54))
abline(lm(ToDK$ToD ~ ToDK$Mass))
plot(ToDB$Mass, ToDB$ToD, xlim =c(0,50), ylim=c(44,54))
abline(lm(ToDB$ToD ~ ToDB$Mass))

#models - Mass
modelToDN_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDN)
statToDN_Mass <- summary(modelToDN_Mass)
modelToDS_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDS)
statToDS_Mass <- summary(modelToDS_Mass)
modelToDK_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDK)
statToDK_Mass <- summary(modelToDK_Mass)
modelToDB_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDB)
statToDB_Mass <- summary(modelToDB_Mass)

plot(ToDN$Fed, ToDN$ToD, xlim =c(0,5), ylim=c(44,54))
abline(lm(ToDN$ToD ~ ToDN$Fed))
plot(ToDS$Fed, ToDS$ToD, xlim =c(0,5), ylim=c(44,54))
abline(lm(ToDS$ToD ~ ToDS$Fed))
plot(ToDK$Fed, ToDK$ToD, xlim =c(0,5), ylim=c(44,54))
abline(lm(ToDK$ToD ~ ToDK$Fed))
plot(ToDB$Fed, ToDB$ToD, xlim =c(0,5), ylim=c(44,54))
abline(lm(ToDB$ToD ~ ToDB$Fed))

#models - Fed
modelToDN_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDN)
statToDN_Fed <- summary(modelToDN_Fed)
modelToDS_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDS)
statToDS_Fed <- summary(modelToDS_Fed)
modelToDK_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDK)
statToDK_Fed <- summary(modelToDK_Fed)
modelToDB_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDB)
statToDB_Fed <- summary(modelToDB_Fed)


#subset - ToD - population - treatment
ToDN15 <- subset(ToDN, Treatment == "15")
ToDN19 <- subset(ToDN, Treatment == "19")
ToDN23 <- subset(ToDN, Treatment == "23")
ToDN25 <- subset(ToDN, Treatment == "25")
ToDN29 <- subset(ToDN, Treatment == "29")

###

pdf("Otavi_temp.pdf", width=18, height=15)
par(mfrow=c(3,5))
plot(ToDN15$D_21, ToDN15$ToD, xlim =c(0,100),main="15")
abline(lm(ToDN15$ToD ~ ToDN15$D_21))
modelToDN15_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDN15)
statToDN15_D_21 <- summary(modelToDN15_D_21)
plot(ToDN19$D_21, ToDN19$ToD, xlim =c(0,100),main="19")
abline(lm(ToDN19$ToD ~ ToDN19$D_21))
modelToDN19_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDN19)
statToDN19_D_21 <- summary(modelToDN19_D_21)
plot(ToDN23$D_21, ToDN23$ToD, xlim =c(0,100),main="23")
abline(lm(ToDN23$ToD ~ ToDN23$D_21))
modelToDN23_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDN23)
statToDN23_D_21 <- summary(modelToDN23_D_21)
plot(ToDN25$D_21, ToDN25$ToD, xlim =c(0,100),main="25")
abline(lm(ToDN25$ToD ~ ToDN25$D_21))
modelToDN25_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDN25)
statToDN25_D_21 <- summary(modelToDN25_D_21)
plot(ToDN29$D_21, ToDN29$ToD, xlim =c(0,100),main="29")
abline(lm(ToDN29$ToD ~ ToDN29$D_21))
modelToDN29_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDN29)
statToDN29_D_21 <- summary(modelToDN29_D_21)
####
plot(ToDN15$Mass, ToDN15$ToD, xlim =c(0,50))
abline(lm(ToDN15$ToD ~ ToDN15$Mass))
modelToDN15_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDN15)
statToDN15_Mass <- summary(modelToDN15_Mass)
plot(ToDN19$Mass, ToDN19$ToD, xlim =c(0,50))
abline(lm(ToDN19$ToD ~ ToDN19$Mass))
modelToDN19_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDN19)
statToDN19_Mass  <- summary(modelToDN19_Mass )
plot(ToDN23$Mass , ToDN23$ToD, xlim =c(0,50))
abline(lm(ToDN23$ToD ~ ToDN23$Mass ))
modelToDN23_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDN23)
statToDN23_Mass <- summary(modelToDN23_Mass)
plot(ToDN25$Mass, ToDN25$ToD, xlim =c(0,50))
abline(lm(ToDN25$ToD ~ ToDN25$Mass))
modelToDN25_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDN25)
statToDN25_Mass <- summary(modelToDN25_Mass)
plot(ToDN29$Mass, ToDN29$ToD, xlim =c(0,50))
abline(lm(ToDN29$ToD ~ ToDN29$Mass))
modelToDN29_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDN29)
statToDN29_Mass <- summary(modelToDN29_Mass)
####
plot(ToDN15$Fed, ToDN15$ToD, xlim =c(0,5))
abline(lm(ToDN15$ToD ~ ToDN15$Fed))
modelToDN15_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDN15)
statToDN15_Fed <- summary(modelToDN15_Fed)
plot(ToDN19$Fed, ToDN19$ToD, xlim =c(0,5))
abline(lm(ToDN19$ToD ~ ToDN19$Fed))
modelToDN19_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDN19)
statToDN19_Fed <- summary(modelToDN19_Fed)
plot(ToDN23$Fed, ToDN23$ToD, xlim =c(0,5))
abline(lm(ToDN23$ToD ~ ToDN23$Fed))
modelToDN23_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDN23)
statToDN23_Fed <- summary(modelToDN23_Fed)
plot(ToDN25$Fed, ToDN25$ToD, xlim =c(0,5))
abline(lm(ToDN25$ToD ~ ToDN25$Fed))
modelToDN25_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDN25)
statToDN25_Fed <- summary(modelToDN25_Fed)
plot(ToDN29$Fed, ToDN29$ToD, xlim =c(0,5))
abline(lm(ToDN29$ToD ~ ToDN29$Fed))
modelToDN29_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDN29)
statToDN29_Fed <- summary(modelToDN29_Fed)
dev.off()
#####

ToDS15 <- subset(ToDS, Treatment == "15")
ToDS19 <- subset(ToDS, Treatment == "19")
ToDS23 <- subset(ToDS, Treatment == "23")
ToDS25 <- subset(ToDS, Treatment == "25")
ToDS29 <- subset(ToDS, Treatment == "29")

###
pdf("Stampriet.pdf", width=18, height=15)
par(mfrow=c(3,5))
plot(ToDS15$D_21, ToDS15$ToD, xlim =c(0,100))
abline(lm(ToDS15$ToD ~ ToDS15$D_21))
modelToDS15_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDS15)
statToDS15_D_21 <- summary(modelToDS15_D_21)
plot(ToDS19$D_21, ToDS19$ToD, xlim =c(0,100))
abline(lm(ToDS19$ToD ~ ToDS19$D_21))
modelToDS19_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDS19)
statToDS19_D_21 <- summary(modelToDS19_D_21)
plot(ToDS23$D_21, ToDS23$ToD, xlim =c(0,100))
abline(lm(ToDS23$ToD ~ ToDS23$D_21))
modelToDS23_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDS23)
statToDS23_D_21 <- summary(modelToDS23_D_21)
plot(ToDS25$D_21, ToDS25$ToD, xlim =c(0,100))
abline(lm(ToDS25$ToD ~ ToDS25$D_21))
modelToDS25_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDS25)
statToDS25_D_21 <- summary(modelToDS25_D_21)
plot(ToDS29$D_21, ToDS29$ToD, xlim =c(0,100))
abline(lm(ToDS29$ToD ~ ToDS29$D_21))
modelToDS29_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDS29)
statToDS29_D_21 <- summary(modelToDS29_D_21)
####
plot(ToDS15$Mass, ToDS15$ToD, xlim =c(0,50))
abline(lm(ToDS15$ToD ~ ToDS15$Mass))
modelToDS15_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDS15)
statToDS15_Mass <- summary(modelToDS15_Mass)
plot(ToDS19$Mass, ToDS19$ToD, xlim =c(0,50))
abline(lm(ToDS19$ToD ~ ToDS19$Mass))
modelToDS19_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDS19)
statToDS19_Mass  <- summary(modelToDS19_Mass )
plot(ToDS23$Mass , ToDS23$ToD, xlim =c(0,50))
abline(lm(ToDS23$ToD ~ ToDS23$Mass ))
modelToDS23_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDS23)
statToDS23_Mass <- summary(modelToDS23_Mass)
plot(ToDS25$Mass, ToDS25$ToD, xlim =c(0,50))
abline(lm(ToDS25$ToD ~ ToDS25$Mass))
modelToDS25_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDS25)
statToDS25_Mass <- summary(modelToDS25_Mass)
plot(ToDS29$Mass, ToDS29$ToD, xlim =c(0,50))
abline(lm(ToDS29$ToD ~ ToDS29$Mass))
modelToDS29_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDS29)
statToDS29_Mass <- summary(modelToDS29_Mass)
####
plot(ToDS15$Fed, ToDS15$ToD, xlim =c(0,5))
abline(lm(ToDS15$ToD ~ ToDS15$Fed))
modelToDS15_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDS15)
statToDS15_Fed <- summary(modelToDS15_Fed)
plot(ToDS19$Fed, ToDS19$ToD, xlim =c(0,5))
abline(lm(ToDS19$ToD ~ ToDS19$Fed))
modelToDS19_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDS19)
statToDS19_Fed <- summary(modelToDS19_D_21)
plot(ToDS23$Fed, ToDS23$ToD, xlim =c(0,5))
abline(lm(ToDS23$ToD ~ ToDS23$Fed))
modelToDS23_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDS23)
statToDS23_Fed <- summary(modelToDS23_Fed)
plot(ToDS25$Fed, ToDS25$ToD, xlim =c(0,5))
abline(lm(ToDS25$ToD ~ ToDS25$Fed))
modelToDS25_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDS25)
statToDS25_Fed <- summary(modelToDS25_Fed)
plot(ToDS29$Fed, ToDS29$ToD, xlim =c(0,5))
abline(lm(ToDS29$ToD ~ ToDS29$Fed))
modelToDS29_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDS29)
statToDS29_Fed <- summary(modelToDS29_Fed)
dev.off()
#####

ToDK15 <- subset(ToDK, Treatment == "15")
ToDK19 <- subset(ToDK, Treatment == "19")
ToDK23 <- subset(ToDK, Treatment == "23")
ToDK25 <- subset(ToDK, Treatment == "25")
ToDK29 <- subset(ToDK, Treatment == "29")

###
pdf("Karasburg.pdf", width=18, height=15)
par(mfrow=c(3,5))
plot(ToDK15$D_21, ToDK15$ToD, xlim =c(0,100))
abline(lm(ToDK15$ToD ~ ToDK15$D_21))
modelToDK15_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDK15)
statToDK15_D_21 <- summary(modelToDK15_D_21)
plot(ToDK19$D_21, ToDK19$ToD, xlim =c(0,100))
abline(lm(ToDK19$ToD ~ ToDK19$D_21))
modelToDK19_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDK19)
statToDK19_D_21 <- summary(modelToDK19_D_21)
plot(ToDK23$D_21, ToDK23$ToD, xlim =c(0,100))
abline(lm(ToDK23$ToD ~ ToDK23$D_21))
modelToDK23_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDK23)
statToDK23_D_21 <- summary(modelToDK23_D_21)
plot(ToDK25$D_21, ToDK25$ToD, xlim =c(0,100))
abline(lm(ToDK25$ToD ~ ToDK25$D_21))
modelToDK25_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDK25)
statToDK25_D_21 <- summary(modelToDK25_D_21)
plot(ToDK29$D_21, ToDK29$ToD, xlim =c(0,100))
abline(lm(ToDK29$ToD ~ ToDK29$D_21))
modelToDK29_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDK29)
statToDK29_D_21 <- summary(modelToDK29_D_21)
####
plot(ToDK15$Mass, ToDK15$ToD, xlim =c(0,50))
abline(lm(ToDK15$ToD ~ ToDK15$Mass))
modelToDK15_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDK15)
statToDK15_Mass <- summary(modelToDK15_Mass)
plot(ToDK19$Mass, ToDK19$ToD, xlim =c(0,50))
abline(lm(ToDK19$ToD ~ ToDK19$Mass))
modelToDK19_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDK19)
statToDK19_Mass <- summary(modelToDK19_Mass)
plot(ToDK23$Mass, ToDK23$ToD, xlim =c(0,50))
abline(lm(ToDK23$ToD ~ ToDK23$Mass))
modelToDK23_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDK23)
statToDK23_Mass <- summary(modelToDK23_Mass)
plot(ToDK25$Mass, ToDK25$ToD, xlim =c(0,50))
abline(lm(ToDK25$ToD ~ ToDK25$Mass))
modelToDK25_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDK25)
statToDK25_Mass <- summary(modelToDK25_Mass)
plot(ToDK29$Mass, ToDK29$ToD, xlim =c(0,50))
abline(lm(ToDK29$ToD ~ ToDK29$Mass))
modelToDK29_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDK29)
statToDK29_Mass <- summary(modelToDK29_Mass)
####
plot(ToDK15$Fed, ToDK15$ToD, xlim =c(0,5))
abline(lm(ToDK15$ToD ~ ToDK15$Fed))
modelToDK15_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDK15)
statToDK15_Fed <- summary(modelToDK15_Fed)
plot(ToDK19$Fed, ToDK19$ToD, xlim =c(0,5))
abline(lm(ToDK19$ToD ~ ToDK19$Fed))
modelToDK19_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDK19)
statToDK19_Fed <- summary(modelToDK19_Fed)
plot(ToDK23$Fed, ToDK23$ToD, xlim =c(0,5))
abline(lm(ToDK23$ToD ~ ToDK23$Fed))
modelToDK23_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDK23)
statToDK23_Fed <- summary(modelToDK23_Fed)
plot(ToDK25$Fed, ToDK25$ToD, xlim =c(0,5))
abline(lm(ToDK25$ToD ~ ToDK25$Fed))
modelToDK25_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDK25)
statToDK25_Fed <- summary(modelToDK25_Fed)
plot(ToDK29$Fed, ToDK29$ToD, xlim =c(0,5))
abline(lm(ToDK29$ToD ~ ToDK29$Fed))
modelToDK29_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDK29)
statToDK29_Fed <- summary(modelToDK29_Fed)
dev.off()
#####

ToDB15 <- subset(ToDB, Treatment == "15")
ToDB19 <- subset(ToDB, Treatment == "19")
ToDB23 <- subset(ToDB, Treatment == "23")
ToDB25 <- subset(ToDB, Treatment == "25")
ToDB29 <- subset(ToDB, Treatment == "29")

###
pdf("Betta.pdf", width=18, height=15)
par(mfrow=c(3,5))
plot(ToDB15$D_21, ToDB15$ToD, xlim =c(0,100))
abline(lm(ToDB15$ToD ~ ToDB15$D_21))
modelToDB15_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDB15)
statToDB15_D_21 <- summary(modelToDB15_D_21)
plot(ToDB19$D_21, ToDB19$ToD, xlim =c(0,100))
abline(lm(ToDB19$ToD ~ ToDB19$D_21))
modelToDB19_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDB19)
statToDB19_D_21 <- summary(modelToDB19_D_21)
plot(ToDB23$D_21, ToDB23$ToD, xlim =c(0,100))
abline(lm(ToDB23$ToD ~ ToDB23$D_21))
modelToDB23_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDB23)
statToDB23_D_21 <- summary(modelToDB23_D_21)
plot(ToDB25$D_21, ToDB25$ToD, xlim =c(0,100))
abline(lm(ToDB25$ToD ~ ToDB25$D_21))
modelToDB25_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDB25)
statToDB25_D_21 <- summary(modelToDB25_D_21)
plot(ToDB29$D_21, ToDB29$ToD, xlim =c(0,100))
abline(lm(ToDB29$ToD ~ ToDB29$D_21))
modelToDB29_D_21 <- lm(ToD ~ D_21 + Treatment + D_21*Treatment, data=ToDB29)
statToDB29_D_21 <- summary(modelToDB29_D_21)
####
plot(ToDB15$Mass, ToDB15$ToD, xlim =c(0,50))
abline(lm(ToDB15$ToD ~ ToDB15$Mass))
modelToDB15_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDB15)
statToDB15_Mass <- summary(modelToDB15_Mass)
plot(ToDB19$Mass, ToDB19$ToD, xlim =c(0,50))
abline(lm(ToDB19$ToD ~ ToDB19$Mass))
modelToDB19_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDB19)
statToDB19_Mass  <- summary(modelToDB19_Mass)
plot(ToDB23$Mass , ToDB23$ToD, xlim =c(0,50))
abline(lm(ToDB23$ToD ~ ToDB23$Mass ))
modelToDB23_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDB23)
statToDB23_Mass <- summary(modelToDB23_Mass)
plot(ToDB25$Mass, ToDB25$ToD, xlim =c(0,50))
abline(lm(ToDB25$ToD ~ ToDB25$Mass))
modelToDB25_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDB25)
statToDB25_Mass <- summary(modelToDB25_Mass)
plot(ToDB29$Mass, ToDB29$ToD, xlim =c(0,50))
abline(lm(ToDB29$ToD ~ ToDB29$Mass))
modelToDB29_Mass <- lm(ToD ~ Mass + Treatment + Mass*Treatment, data=ToDB29)
statToDB29_Mass <- summary(modelToDB29_Mass)
####
plot(ToDB15$Fed, ToDB15$ToD, xlim =c(0,5))
abline(lm(ToDB15$ToD ~ ToDB15$Fed))
modelToDB15_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDB15)
statToDB15_Fed <- summary(modelToDB15_Fed)
plot(ToDB19$Fed, ToDB19$ToD, xlim =c(0,5))
abline(lm(ToDB19$ToD ~ ToDB19$Fed))
modelToDB19_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDB19)
statToDB19_Fed  <- summary(modelToDB19_Fed )
plot(ToDB23$Fed , ToDB23$ToD, xlim =c(0,5))
abline(lm(ToDB23$ToD ~ ToDB23$Fed ))
modelToDB23_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDB23)
statToDB23_Fed <- summary(modelToDB23_Fed)
plot(ToDB25$Fed, ToDB25$ToD, xlim =c(0,5))
abline(lm(ToDB25$ToD ~ ToDB25$Fed))
modelToDB25_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDB25)
statToDB25_Fed <- summary(modelToDB25_Fed)
plot(ToDB29$Fed, ToDB29$ToD, xlim =c(0,5))
abline(lm(ToDB29$ToD ~ ToDB29$Fed))
modelToDB29_Fed <- lm(ToD ~ Fed + Treatment + Fed*Treatment, data=ToDB29)
statToDB29_Fed <- summary(modelToDB29_Fed)
dev.off()
#####



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






#########

#subset - ToD_final - population
ToD_finalN <- subset(ToD_final, Population == "Otavi")
ToD_finalS <- subset(ToD_final, Population == "Stampriet")
ToD_finalK <- subset(ToD_final, Population == "Karasburg")
ToD_finalB <- subset(ToD_final, Population == "Betta")

par(mfrow=c(1,4))
plot(ToD_finalN$Treatment, ToD_finalN$ToD_mean, ylim=c(48,52), xlab = "Treatment", ylab = "CTmax", main = "Critical thermal maximum - Otavi")
abline(lm(ToD_finalN$ToD_mean ~ ToD_finalN$Treatment), col = "red")
plot(ToD_finalS$Treatment, ToD_finalS$ToD_mean, ylim=c(48,52), xlab = "Treatment", ylab = "CTmax", main = "Critical thermal maximum - Stampriet")
abline(lm(ToD_finalS$ToD_mean ~ ToD_finalS$Treatment), col = "blue")
plot(ToD_finalK$Treatment, ToD_finalK$ToD_mean, ylim=c(48,52), xlab = "Treatment", ylab = "CTmax", main = "Critical thermal maximum - Karasburg")
abline(lm(ToD_finalK$ToD_mean ~ ToD_finalK$Treatment), col = "green")
plot(ToD_finalB$Treatment, ToD_finalB$ToD_mean, ylim=c(48,52), xlab = "Treatment", ylab = "CTmax", main = "Critical thermal maximum - Betta")
abline(lm(ToD_finalB$ToD_mean ~ ToD_finalB$Treatment), col = "black")

#models
modelToD_finalN <- lm(ToD ~ Treatment, data=ToD_finalN)
statToD_finalN <- summary(modelToD_finalN)
modelToD_finalS <- lm(ToD ~ Treatment, data=ToD_finalS)
statToD_finalS <- summary(modelToD_finalS)
modelToD_finalK <- lm(ToD ~ Treatment, data=ToD_finalK)
statToD_finalK <- summary(modelToD_finalK)
modelToD_finalB <- lm(ToD ~ Treatment, data=ToD_finalB)
statToD_finalB <- summary(modelToD_finalB)

par(mfrow=c(4,4))
plot(modelToD_finalN)
plot(modelToD_finalS)
plot(modelToD_finalK)
plot(modelToD_finalB)

#subset - ToD_final - treatment
ToD_final15 <- subset(ToD_final, Treatment == "15")
ToD_final19 <- subset(ToD_final, Treatment == "19")
ToD_final23 <- subset(ToD_final, Treatment == "23")
ToD_final25 <- subset(ToD_final, Treatment == "25")
ToD_final29 <- subset(ToD_final, Treatment == "29")

par(mfrow=c(1,5))
plot(factor(ToD_final15$Population, levels = c("Betta", "S", "Karasburg", "North")), ToD_final15$ToD_mean, ylim=c(48,52), xlab = "Population", ylab = "CTmax", main = "Acclimation - 15")
plot(factor(ToD_final19$Population, levels = c("Betta", "S", "Karasburg", "North")), ToD_final19$ToD_mean, ylim=c(48,52), xlab = "Population", ylab = "CTmax", main = "Acclimation - 19")
plot(factor(ToD_final23$Population, levels = c("Betta", "S", "Karasburg", "North")), ToD_final23$ToD_mean, ylim=c(48,52), xlab = "Population", ylab = "CTmax", main = "Acclimation - 23")
plot(factor(ToD_final25$Population, levels = c("Betta", "S", "Karasburg", "North")), ToD_final25$ToD_mean, ylim=c(48,52), xlab = "Population", ylab = "CTmax", main = "Acclimation - 25")
plot(factor(ToD_final29$Population, levels = c("Betta", "S", "Karasburg", "North")), ToD_final29$ToD_mean, ylim=c(48,52), xlab = "Population", ylab = "CTmax", main = "Acclimation - 29")

#models
modelToD_final15 <- lm(ToD ~ Population, data=ToD_final15)
statToD_final15 <- summary(modelToD_final15)
modelToD_final19 <- lm(ToD ~ Population, data=ToD_final19)
statToD_final19 <- summary(modelToD_final19)
modelToD_final23 <- lm(ToD ~ Population, data=ToD_final23)
statToD_final23 <- summary(modelToD_final23)
modelToD_final25 <- lm(ToD ~ Population, data=ToD_final25)
statToD_final25 <- summary(modelToD_final25)
modelToD_final29 <- lm(ToD ~ Population, data=ToD_final29)
statToD_final29 <- summary(modelToD_final29)

par(mfrow=c(5,4))
plot(modelToD_final15)
plot(modelToD_final19)
plot(modelToD_final23)
plot(modelToD_final25)
plot(modelToD_final29)

#####################
#plots

ToD_finalN15 <- subset(ToD_final, Population == "North" & Treatment == "15")
ToD_finalN19 <- subset(ToD_final, Population == "North" & Treatment == "19")
ToD_finalN23 <- subset(ToD_final, Population == "North" & Treatment == "23")
ToD_finalN25 <- subset(ToD_final, Population == "North" & Treatment == "25")
ToD_finalN29 <- subset(ToD_final, Population == "North" & Treatment == "29")
ToD_finalS15 <- subset(ToD_final, Population == "S" & Treatment == "15")
ToD_finalS19 <- subset(ToD_final, Population == "S" & Treatment == "19")
ToD_finalS23 <- subset(ToD_final, Population == "S" & Treatment == "23")
ToD_finalS25 <- subset(ToD_final, Population == "S" & Treatment == "25")
ToD_finalS29 <- subset(ToD_final, Population == "S" & Treatment == "29")
ToD_finalK15 <- subset(ToD_final, Population == "Karasburg" & Treatment == "15")
ToD_finalK19 <- subset(ToD_final, Population == "Karasburg" & Treatment == "19")
ToD_finalK23 <- subset(ToD_final, Population == "Karasburg" & Treatment == "23")
ToD_finalK25 <- subset(ToD_final, Population == "Karasburg" & Treatment == "25")
ToD_finalK29 <- subset(ToD_final, Population == "Karasburg" & Treatment == "29")
ToD_finalB15 <- subset(ToD_final, Population == "Betta" & Treatment == "15")
ToD_finalB19 <- subset(ToD_final, Population == "Betta" & Treatment == "19")
ToD_finalB23 <- subset(ToD_final, Population == "Betta" & Treatment == "23")
ToD_finalB25 <- subset(ToD_final, Population == "Betta" & Treatment == "25")
ToD_finalB29 <- subset(ToD_final, Population == "Betta" & Treatment == "29")


Treat <- c(15, 19, 23, 25, 29)
NmCTmax <- c(colMeans(matrix(ToD_finalN15$ToD_mean),na.rm=T), colMeans(matrix(ToD_finalN19$ToD_mean),na.rm=T), colMeans(matrix(ToD_finalN23$ToD_mean),na.rm=T), colMeans(matrix(ToD_finalN25$ToD_mean),na.rm=T), colMeans(matrix(ToD_finalN29$ToD_mean),na.rm=T))
SmCTmax <- c(colMeans(matrix(ToD_finalS15$ToD_mean),na.rm=T), colMeans(matrix(ToD_finalS19$ToD_mean),na.rm=T), colMeans(matrix(ToD_finalS23$ToD_mean),na.rm=T), colMeans(matrix(ToD_finalS25$ToD_mean),na.rm=T), colMeans(matrix(ToD_finalS29$ToD_mean),na.rm=T))
KmCTmax <- c(colMeans(matrix(ToD_finalK15$ToD_mean),na.rm=T), colMeans(matrix(ToD_finalK19$ToD_mean),na.rm=T), colMeans(matrix(ToD_finalK23$ToD_mean),na.rm=T), colMeans(matrix(ToD_finalK25$ToD_mean),na.rm=T), colMeans(matrix(ToD_finalK29$ToD_mean),na.rm=T))
BmCTmax <- c(colMeans(matrix(ToD_finalB15$ToD_mean),na.rm=T), colMeans(matrix(ToD_finalB19$ToD_mean),na.rm=T), colMeans(matrix(ToD_finalB23$ToD_mean),na.rm=T), colMeans(matrix(ToD_finalB25$ToD_mean),na.rm=T), colMeans(matrix(ToD_finalB29$ToD_mean),na.rm=T))


plot(Treat, NmCTmax, , col = "red", ylim = c(48,52), ylab = "CTmax", xlab ="Acclimation temperature", main ="Critical thermal maximum")
abline(lm(NmCTmax ~Treat), col = "red")
legend(9,7,"Betta")
par(new=TRUE)
plot(Treat, SmCTmax, , col = "blue", ylim = c(48,52), ylab = "CTmax", xlab ="Acclimation temperature")
abline(lm(SmCTmax ~Treat), col = "blue")
par(new=TRUE)
plot(Treat, KmCTmax, , col = "green", ylim = c(48,52), ylab = "CTmax", xlab ="Acclimation temperature")
abline(lm(KmCTmax ~Treat), col = "green")
par(new=TRUE)
plot(Treat, BmCTmax, , col = "black", ylim = c(48,52), ylab = "CTmax", xlab ="Acclimation temperature")
abline(lm(BmCTmax ~Treat), col = "black")
legend(15.5, 52, legend=c("North", "Stampriet", "Karasburg", "Betta"),
       col=c("red", "blue", "green", "black"), lty=1, cex=0.8)

