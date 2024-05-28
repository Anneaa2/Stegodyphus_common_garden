
########## Growth ####

setwd("/home/anneaa/spider2/faststorage/assay_study_anne/phenotypes/data/growth_surv/")
library(agricolae)
library(betareg)

#plot(0,type='n',axes=FALSE,ann=FALSE)
 
Grow <- read.table("growth_and_survival.txt", header = TRUE, sep = "", dec = '.',stringsAsFactors=F,as.is=c(F,F,F,T,T,	T,T,T,T,T,	T))
par(mfrow=c(1,2))
plot(Grow$Population,Grow$growth,ylab="Growth factor",xlab="population")
plot(Grow$Treatment,Grow$growth)

modelGrow <- lm(growth ~ Population + Treatment + Population*Treatment, data=Grow)
summary(modelGrow)

a <- aov(growth ~ Population + Treatment + Population*Treatment, data=Grow)
summary(a)

## No effect of population, effect of treatment.

#subset - Grow - population
GrowN <- subset(Grow, Population == "Otavi")
GrowS <- subset(Grow, Population == "Stampriet")
GrowK <- subset(Grow, Population == "Karasburg")
GrowB <- subset(Grow, Population == "Betta")

pdf("Effect_of_treat_perpop_GROW.pdf")
par(mfrow=c(1,4))
plot(GrowN$Treatment,GrowN$growth, main="Otavi", ylab="Growth ratio")
abline(lm(growth ~ Treatment,data=GrowN))
plot(GrowS$Treatment,GrowS$growth, main="Stampriet",ylab=NA)
abline(lm(growth ~ Treatment,data=GrowS))
plot(GrowK$Treatment,GrowK$growth, main="Karasburg",ylab=NA)
abline(lm(growth ~ Treatment,data=GrowK))
plot(GrowB$Treatment,GrowB$growth, main="Betta",ylab=NA)
abline(lm(growth ~ Treatment,data=GrowB))
dev.off()

#subset - Grow - population - treatment
Grow15 <- subset(Grow, Treatment == "15")
Grow19 <- subset(Grow, Treatment == "19")
Grow23 <- subset(Grow, Treatment == "23")
Grow25 <- subset(Grow, Treatment == "25")
Grow29 <- subset(Grow, Treatment == "29")

pdf("Effect_of_pop_pertreat_GROW.pdf")
par(mfrow=c(1,5))
plot(Grow15$Population,Grow15$growth,main="15",las=2,ylim=c(0,4),ylab="Growth ratio",xlab=NA)
plot(Grow19$Population,Grow19$growth,main="19",las=2,ylim=c(0,4),ylab=NA,xlab=NA)
plot(Grow23$Population,Grow23$growth,main="23",las=2,ylim=c(0,4),ylab=NA,xlab=NA)
plot(Grow25$Population,Grow25$growth,main="25",las=2,ylim=c(0,4),ylab=NA,xlab=NA)
plot(Grow29$Population,Grow29$growth,main="29",las=2,ylim=c(0,4),ylab=NA,xlab=NA)
dev.off()


##### Survival



par(mfrow=c(1,2))
plot(Grow$Population,Grow$survival,ylab="Survival",xlab="population")
plot(Grow$Treatment,Grow$survival)

modelSURV <- glm(survival ~ Population + Treatment + Population*Treatment, data=Grow,family=quasibinomial(link="logit"))
summary(modelSURV)
plot(modelSURV)
# quasibinomial could be adequate regression and look ok. Not good, but okay. 

# Trying out betaregression with inflated 1'es
	# https://stats.stackexchange.com/questions/496050/why-is-betareg-giving-invalid-dependent-variable-error
library(gamlss)
library(gamlss.dist)

modelSURVb <- gamlss(survival ~ Population + Treatment + Population*Treatment, data=Grow,family=BEOI)
summary(modelSURVb)
#	 Looks good in theory. plots not so pretty.
'******************************************************************
Family:  c("BEOI", "One Inflated Beta")

Call:  gamlss(formula = survival ~ Population + Treatment +
    Population * Treatment, family = BEOI, data = Grow)

Fitting method: RS()

------------------------------------------------------------------
Mu link function:  logit
Mu Coefficients:
                               Estimate Std. Error t value Pr(>|t|)
(Intercept)                    0.500595   0.406585   1.231  0.21923
PopulationKarasburg            0.276512   0.653318   0.423  0.67243
PopulationOtavi               -0.017669   0.577900  -0.031  0.97563
PopulationStampriet           -0.120374   0.723580  -0.166  0.86799
Treatment                      0.059121   0.019219   3.076  0.00229 **
PopulationKarasburg:Treatment -0.003250   0.030478  -0.107  0.91515
PopulationOtavi:Treatment      0.002641   0.027521   0.096  0.92361
PopulationStampriet:Treatment  0.012205   0.033933   0.360  0.71934
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

------------------------------------------------------------------
Sigma link function:  log
Sigma Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  2.35515    0.09326   25.25   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

------------------------------------------------------------------
Nu link function:  logit
Nu Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  -1.2826     0.1391   -9.22   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

------------------------------------------------------------------
No. of observations in the fit:  304
Degrees of Freedom for the fit:  10
      Residual Deg. of Freedom:  294
                      at cycle:  8

Global Deviance:     -190.1661
            AIC:     -170.1661
            SBC:     -132.9958
******************************************************************
'	
plot(modelSURVb)




a_surv <- aov(survival ~ Population + Treatment + Population*Treatment, data=Grow)
summary(a_surv)
	## No effect of population, effect of treatment
	
	
#subset - Grow - population

pdf("Effect_of_treat_perpop_SURV.pdf")
par(mfrow=c(1,4))
plot(GrowN$Treatment,GrowN$survival, main="Otavi", ylab="Survival")
abline(lm(survival ~ Treatment,data=GrowN))
plot(GrowS$Treatment,GrowS$survival, main="Stampriet",ylab=NA)
abline(lm(survival ~ Treatment,data=GrowS))
plot(GrowK$Treatment,GrowK$survival, main="Karasburg",ylab=NA)
abline(lm(survival ~ Treatment,data=GrowK))
plot(GrowB$Treatment,GrowB$survival, main="Betta",ylab=NA)
abline(lm(survival ~ Treatment,data=GrowB))
dev.off()

#subset - Grow - population - treatment

pdf("Effect_of_pop_pertreat_SURV.pdf")
par(mfrow=c(1,5))
plot(Grow15$Population,Grow15$survival,main="15",las=2,ylim=c(0.2,1),ylab="Survival",xlab=NA)
plot(Grow19$Population,Grow19$survival,main="19",las=2,ylim=c(0.2,1),ylab=NA,xlab=NA)
plot(Grow23$Population,Grow23$survival,main="23",las=2,ylim=c(0.2,1),ylab=NA,xlab=NA)
plot(Grow25$Population,Grow25$survival,main="25",las=2,ylim=c(0.2,1),ylab=NA,xlab=NA)
plot(Grow29$Population,Grow29$survival,main="29",las=2,ylim=c(0.2,1),ylab=NA,xlab=NA)
dev.off()




