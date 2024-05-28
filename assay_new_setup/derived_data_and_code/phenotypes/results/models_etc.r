#!/home/anneaa/miniconda3/envs/assay/bin/Rscript

#source("/faststorage/project/spider2/assay_study_anne/phenotypes/results/models_etc.r")
#srun --mem=10g --pty /bin/bash

### modify to contain models. But datasets are being entered here, so it can be modified accordingly.

#conda create -n nice_tables R r-agricolae r-performance r-multcomp r-emmeans r-sjPlot r-TMB

### How do I print all this to nice tables?
#conda activate nice_tables; R

#conda activate assay_copy; R
#setwd("/faststorage/project/spider2/assay_study_anne/phenotypes/results/")

library(performance) # PSEUdo r2 for binomial model.
library(agricolae)
library(sjPlot) #for printing nicer tables
library(multcomp)
library(emmeans)
library(multcompView)
library(scales)
library(dplyr)
#library(tibble)
#library(flextable)

popcols=c("steelblue1","blue","tan4","green3","snow4","red") # B, K, N, S, G, E/O
popcols=c(alpha("lightskyblue1",.8),"blue3","tan4","green3","red2") # B, K, N, S, E/O    #### OBS new modified popcols
names(popcols) <- c("Betta", "Karasburg", "Ndumo", "Stampriet", "Otavi")
#popcols <- popcols[c(1,2,6,4)]	# B, K, E/O,	S
popcols2 <- popcols[c(grep("Betta",names(popcols)), grep("Karasburg",names(popcols)), grep("Otavi",names(popcols)),grep("Stampriet",names(popcols)))]	# B, K, E/O,	S

Mean_Center=T

### Functions for script
print_nice_anova_tables <- function(anova_output, filename, title_is){
  rnames = rownames(anova_output)
  cnames = colnames(anova_output)
  df_anova <- as.data.frame(anova_output)
  # Assign astericks according to p values
  
  df_anova$sign[df_anova[, which( cnames == cnames[ grep("Pr.*",cnames) ] ) ] < 0.05] <- "*" 
  df_anova$sign[df_anova[, which( cnames == cnames[ grep("Pr.*",cnames) ] ) ] < 0.01] <- "**"
  df_anova$sign[df_anova[, which( cnames == cnames[ grep("Pr.*",cnames) ] ) ] < 0.001] <- "***"
  df_anova$sign[df_anova[, which( cnames == cnames[ grep("Pr.*",cnames) ] ) ] > 0.05] <- ""
  colnames(df_anova) <- c(cnames,"")
  rownames(df_anova) <- rnames
  #sign_values <- capture.output(anova_output)[NROW(capture.output(anova_output))]
  #tab_df(df_anova,file=filename, title=title_is, footnote=sign_values, show.footnote = T,  show.rownames = T, col.header = colnames(df_anova))
  print( df_anova )
  tab_df(df_anova,file=filename, title=title_is, show.rownames = T, col.header = colnames(df_anova), encoding = "Windows-1252")
}

####################################
#------------- Growth--------------#
####################################
Grow <- read.table("assay_new_setup/derived_data_and_code/phenotypes/data/growth_surv/growth_and_survival.txt", header = TRUE, sep = "", dec = '.',stringsAsFactors=F,as.is=c(F,F,F,T,T,	T,T,T,T,T,	T))
#Grow$Treatment <- as.factor(Grow$Treatment)
# Mean center temperature data:
par(mfrow = c(1,2)); hist(Grow$Treatment, main="Non-centered")
if(Mean_Center) { 
  Grow$Treatment <- Grow$Treatment-mean(Grow$Treatment) 
  mean_center_var = "_gMeanC" 
}else{ mean_center_var = "" }
hist(Grow$Treatment,main="Gmean_centered"); dev.off()
# full model:
select_mod <- step(lm(growth ~ Population * Treatment,data=Grow))
 'Start:  AIC=-425.05
growth ~ Population * Treatment

                       Df Sum of Sq    RSS     AIC
<none>                              71.252 -425.05
- Population:Treatment  3    1.5597 72.811 -424.46
'
select_mod
summary(select_mod)
'Call:
lm(formula = growth ~ Population * Treatment, data = Grow)

Coefficients:
                  (Intercept)            PopulationKarasburg
                    -0.784306                      -0.161972
              PopulationOtavi            PopulationStampriet
                     0.522662                      -0.164497
                    Treatment  PopulationKarasburg:Treatment
                     0.104460                       0.005018
    PopulationOtavi:Treatment  PopulationStampriet:Treatment
                    -0.028690                       0.004328
'
anova(select_mod,test="LRT")	# main effects
print_nice_anova_tables(anova_output=anova(select_mod,test="LRT"), title_is="Anova on Growth model", filename=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/Growth_seleced_model", mean_center_var, "_Anova.doc"))
'Analysis of Variance Table

Response: growth
                      Df Sum Sq Mean Sq  F value  Pr(>F)
Population             3  0.837   0.279   1.1587 0.32577
Treatment              1 66.112  66.112 274.6483 < 2e-16 ***
Population:Treatment   3  1.560   0.520   2.1598 0.09286 .
Residuals            296 71.252   0.241
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
'
tab_model(select_mod, file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/Growth_seleced_model", mean_center_var, ".doc"),p.style="numeric_stars", show.se=T, show.r2=T,show.stat=T,show.aic=T)
 

summary(select_mod)			## export to table!
#tab_model(select_mod, file="phenotypes/results/Growth_seleced_model.doc",p.style="numeric_stars", show.se=T, show.r2=T,show.stat=T,show.aic=T)
tab_model(select_mod, file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/Growth_seleced_model", mean_center_var, ".doc"),p.style="numeric_stars", show.se=T, show.r2=T,show.stat=T,show.aic=T)

'
Call:
lm(formula = growth ~ Population * Treatment, data = Grow)

Residuals:
     Min       1Q   Median       3Q      Max
-1.40860 -0.28546 -0.01554  0.22298  2.32367

Coefficients:
                               Estimate Std. Error t value Pr(>|t|)
(Intercept)                   -0.784306   0.230634  -3.401 0.000765 ***
PopulationKarasburg           -0.161972   0.376624  -0.430 0.667464
PopulationOtavi                0.522662   0.323513   1.616 0.107250
PopulationStampriet           -0.164497   0.419652  -0.392 0.695351
Treatment                      0.104460   0.010151  10.291  < 2e-16 ***
PopulationKarasburg:Treatment  0.005018   0.016577   0.303 0.762330
PopulationOtavi:Treatment     -0.028690   0.014275  -2.010 0.045364 *
PopulationStampriet:Treatment  0.004328   0.018635   0.232 0.816502
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4906 on 296 degrees of freedom
Multiple R-squared:  0.4902,    Adjusted R-squared:  0.4781
F-statistic: 40.66 on 7 and 296 DF,  p-value: < 2.2e-16
'
 
 # Post hoc: Treatment
HSD.test(select_mod,"Treatment", group=T)$group  # Tests the data??
'      growth groups
29 2.0895684      a
25 2.0595457      a
23 1.3633561      b
19 0.9705685      c
15 0.9221065      c
'
CLD(emmeans (select_mod,"Treatment",at=list(Treatment=c(15,19,23,25,29))),Letters=letters)
## Tests the MODEL!
'NOTE: Results may be misleading due to involvement in interactions
 Treatment emmean     SE  df lower.CL upper.CL .group
        15  0.759 0.0533 296    0.654    0.864  a
        19  1.158 0.0357 296    1.087    1.228   b
        23  1.556 0.0307 296    1.496    1.616    c
        25  1.755 0.0352 296    1.686    1.825     d
        29  2.154 0.0526 296    2.050    2.257      e

Results are averaged over the levels of: Population
Confidence level used: 0.95
P value adjustment: tukey method for comparing a family of 5 estimates
significance level used: alpha = 0.05
'

# population
HSD.test(select_mod,"Population", group=T)$group  # Tests the data??
'            growth groups
Betta     1.534714      a
Karasburg 1.484138      a
Stampriet 1.434181      a
Otavi     1.411240      a'

## Testing the difference in DATA: (not model)
	# BETTA Only
HSD.test(lm(growth ~ Treatment,data=Grow[which(Grow$Population=="Betta"),]), "Treatment",group=T)$group
'
      growth groups
25 2.1840170      a
29 2.1715966      a
23 1.4272403      b
19 0.9475527      b
15 0.9431614      b
'
# Karasburg Only
HSD.test(lm(growth ~ Treatment,data=Grow[which(Grow$Population=="Karasburg"),]), "Treatment",group=T)$group
' growth groups
29 2.2145974      a
25 2.0485161      a
23 1.3140789      b
19 1.0534604     bc
15 0.7900347      c
'
# Stampriet Only
HSD.test(lm(growth ~ Treatment,data=Grow[which(Grow$Population=="Stampriet"),]), "Treatment",group=T)$group
'   growth groups
29 2.3909400      a
25 1.8656822      a
23 1.2710234      b
19 0.9233278      b
15 0.9221288      b
'
# Otavi Only
HSD.test(lm(growth ~ Treatment,data=Grow[which(Grow$Population=="Otavi"),]), "Treatment",group=T)$group
' 
      growth groups
25 2.0069128      a
29 1.8119742      a
23 1.3702436      b
15 0.9775144      c
19 0.9653675      c
'


####################################
#------------- Survival------------#
####################################

#Function for mcFaddens pseudo R2:
R2logit<- function(model){
    R2<- 1-(model$deviance/model$null.deviance)
    return(R2)
    }
	
surv_deaths_dat <- data.frame(Grow$spiders_3, Grow$spiders_1-Grow$spiders_3, (Grow$Treatment), Grow$Population)
colnames(surv_deaths_dat) <- c("N_alive","N_dead","Treatment","Population")

model_1 <- glm(as.matrix(surv_deaths_dat[,1:2]) ~ Population * Treatment, data=surv_deaths_dat,  family=binomial(link="logit"))
	
select_mod <- step(model_1)
summary(select_mod)
tab_model(select_mod, file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/Survival_seleced_model",mean_center_var,".doc"),p.style="numeric_stars", show.se=T, show.r2=T,show.stat=T,show.aic=T,transform=NULL)
tab_model(select_mod, file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/Survival_seleced_model",mean_center_var,"_transformed.doc"),p.style="numeric_stars", show.se=T, show.r2=T,show.stat=T,show.aic=T)
#tab_model(select_mod, file="Survival_seleced_model_fstat.doc",p.style="numeric_stars", show.se=T, show.r2=T,show.stat=T,show.aic=T)
	'Start:  AIC=1758.92
as.matrix(surv_deaths_dat[, 1:2]) ~ Population * Treatment

                       Df Deviance    AIC
- Population:Treatment  3   1074.2 1754.8
<none>                      1072.2 1758.9

Step:  AIC=1754.83
as.matrix(surv_deaths_dat[, 1:2]) ~ Population + Treatment

             Df Deviance    AIC
<none>            1074.2 1754.8
- Population  3   1092.2 1766.8
- Treatment   1   1307.4 1986.1
'

# Additive is the best model
	
	# calculate p-value of model by comparing to intercept model:
model_0 <- glm(as.matrix(surv_deaths_dat[,1:2]) ~ 1, data=surv_deaths_dat,  family=binomial(link="logit"))
anova(select_mod,model_0, test="Chisq")
print_nice_anova_tables(anova_output=anova(select_mod,model_0, test="Chisq"), title_is="Anova on Survival model", filename=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/Survival_seleced_model", mean_center_var, "_Anova_modelSignChi.doc"))
'Analysis of Deviance Table

Model 1: as.matrix(surv_deaths_dat[, 1:2]) ~ Population + Treatment
Model 2: as.matrix(surv_deaths_dat[, 1:2]) ~ 1
  Resid. Df Resid. Dev Df Deviance  Pr(>Chi)
1       299     1074.2
2       303     1325.4 -4  -251.24 < 2.2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
'

	# get pseudo R:
r2_mcfadden(select_mod)		# for multinomial and logistic regressions:
'# R2 for Generalized Linear Regression
       R2: 0.126
  adj. R2: 0.125
'
anova(select_mod,test="LRT")
print_nice_anova_tables(anova_output=anova(select_mod,test="LRT"), title_is="Anova on Survival model", filename=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/Survival_seleced_model", mean_center_var, "_Anova.doc"))
'Analysis of Deviance Table

Model: binomial, link: logit

Response: as.matrix(surv_deaths_dat[, 1:2])

Terms added sequentially (first to last)


           Df Deviance Resid. Df Resid. Dev  Pr(>Chi)
NULL                         303     1325.4
Population  3   17.992       300     1307.4 0.0004415 ***
Treatment   1  233.247       299     1074.2 < 2.2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
'
summary(select_mod)  ## export table
'Call:
glm(formula = as.matrix(surv_deaths_dat[, 1:2]) ~ Population +
    Treatment, family = binomial(link = "logit"), data = surv_deaths_dat)

Deviance Residuals:
    Min       1Q   Median       3Q      Max
-9.4375  -0.5888   0.6336   1.7769   3.0757

Coefficients:
                     Estimate Std. Error z value Pr(>|z|)
(Intercept)         -0.305976   0.160119  -1.911    0.056 .
PopulationKarasburg  0.418623   0.104222   4.017  5.9e-05 ***
PopulationOtavi      0.059449   0.081479   0.730    0.466
PopulationStampriet  0.118218   0.107960   1.095    0.274
Treatment            0.109126   0.007362  14.823  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 1325.4  on 303  degrees of freedom
Residual deviance: 1074.1  on 299  degrees of freedom
AIC: 1754.8

Number of Fisher Scoring iterations: 5
'
## Calculate pseudo R2

	#Temperature effects
CLD(emmeans (select_mod,specs="Treatment",at=list(Treatment=c(15,19,23,25,29))),Letters=letters)

'Treatment emmean     SE  df asymp.LCL asymp.UCL .group
        15   1.48 0.0533 Inf      1.38      1.58  a
        19   1.92 0.0385 Inf      1.84      1.99   b
        23   2.35 0.0430 Inf      2.27      2.44    c
        25   2.57 0.0519 Inf      2.47      2.67     d
        29   3.01 0.0752 Inf      2.86      3.16      e

Results are averaged over the levels of: Population
Results are given on the logit (not the response) scale.
Confidence level used: 0.95
Results are given on the log odds ratio (not the response) scale.
P value adjustment: tukey method for comparing a family of 5 estimates
significance level used: alpha = 0.05
'

	# Population Effects
CLD(emmeans (select_mod,specs=pairwise ~ Population)$emmeans,Letters=letters)
' Population emmean     SE  df asymp.LCL asymp.UCL .group
 Betta        2.11 0.0590 Inf      1.99      2.22  a
 Otavi        2.17 0.0596 Inf      2.05      2.28  a
 Stampriet    2.23 0.0927 Inf      2.04      2.41  ab
 Karasburg    2.53 0.0882 Inf      2.35      2.70   b

Results are given on the logit (not the response) scale.
Confidence level used: 0.95
Results are given on the log odds ratio (not the response) scale.
P value adjustment: tukey method for comparing a family of 4 estimates
significance level used: alpha = 0.05
'
		# Model treatment effects:
# Betta Only
surv_deaths_dat_sub <- surv_deaths_dat[which(surv_deaths_dat$Population == "Betta"),]
model_1 <- glm(as.matrix(surv_deaths_dat_sub[,1:2]) ~ Treatment + 0, data=surv_deaths_dat_sub ,  family=binomial(link="logit"))
CLD(emmeans (model_1,specs="Treatment",at=list(Treatment=c(15,19,23,25,29))),Letters=letters)
	# all different
# Karasburg Only
surv_deaths_dat_sub <- surv_deaths_dat[which(surv_deaths_dat$Population == "Karasburg"),]
model_1 <- glm(as.matrix(surv_deaths_dat_sub[,1:2]) ~ Treatment + 0, data=surv_deaths_dat_sub ,  family=binomial(link="logit"))
CLD(emmeans (model_1,specs="Treatment",at=list(Treatment=c(15,19,23,25,29))),Letters=letters)
	# all different
# Stampriet Only
surv_deaths_dat_sub <- surv_deaths_dat[which(surv_deaths_dat$Population == "Stampriet"),]
model_1 <- glm(as.matrix(surv_deaths_dat_sub[,1:2]) ~ Treatment + 0, data=surv_deaths_dat_sub ,  family=binomial(link="logit"))
CLD(emmeans (model_1,specs="Treatment",at=list(Treatment=c(15,19,23,25,29))),Letters=letters)
	# all different
# Otavi Only
surv_deaths_dat_sub <- surv_deaths_dat[which(surv_deaths_dat$Population == "Otavi"),]
model_1 <- glm(as.matrix(surv_deaths_dat_sub[,1:2]) ~ Treatment + 0, data=surv_deaths_dat_sub ,  family=binomial(link="logit"))
CLD(emmeans (model_1,specs="Treatment",at=list(Treatment=c(15,19,23,25,29))),Letters=letters)
	# all different


####################################
#------------- CTmax---------------#
####################################
ToD_final <- read.csv(file="assay_new_setup/derived_data_and_code/phenotypes/data/ToD/ToD_nest-temp_mean.csv", header=T,stringsAsFactors=T) 
selec_mod_estimates_noMeanC <- lm(ToD_mean ~ Population:Treatment + Population-1 + Mass_mean + Fed_mean + D_21_mean , data=ToD_final)
coef(selec_mod_estimates_noMeanC)
write.csv(file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/ctmax_select_mod_NOTMEANC_","intc_slope.csv"),coef(selec_mod_estimates_noMeanC))
'coef(selec_mod_estimates_noMeanC)
              PopulationBetta           PopulationKarasburg 
                 48.100131012                  49.629602900 
              PopulationOtavi           PopulationStampriet 
                 48.570277073                  49.313092129 
                    Mass_mean                      Fed_mean 
                 -0.011820391                   0.051216707 
                    D_21_mean     PopulationBetta:Treatment 
                  0.003073138                   0.104047911 
PopulationKarasburg:Treatment     PopulationOtavi:Treatment 
                  0.022100135                   0.086150406 
PopulationStampriet:Treatment 
                  0.041651792   '
if(Mean_Center) { 
  ToD_final$Treatment_not_centered <- ToD_final$Treatment
  ToD_final$Treatment <- ToD_final$Treatment-mean(ToD_final$Treatment) }
#ToD_final$Treatment <- as.factor(ToD_final$Treatment)
#nrow(ToD_final)
#[1] 366
#ToD_final <- na.omit(ToD_final[,c(3,4,6,12,13,14,15)])
#nrow(ToD_final)
#[1] 309

#newdat$pred = predict(fit, newdata = newdat)

full_mod <- lm(ToD_mean ~ Population*Treatment + Mass_mean + Fed_mean + D_21_mean , data=ToD_final)
select_mod <- step(full_mod)
summary(select_mod) # simplest with equal explanatory power (Rsq of 0.4088)

'Start:  AIC=-471.45
ToD_mean ~ Population * Treatment + Mass_mean + Fed_mean + D_21_mean

                       Df Sum of Sq    RSS     AIC
<none>                              62.578 -471.45
- Mass_mean             1    0.8134 63.392 -469.46
- Fed_mean              1    0.9173 63.496 -468.95
- D_21_mean             1    1.0563 63.635 -468.28
- Population:Treatment  3    7.7394 70.318 -441.42
'
anova(select_mod, test="LRT")
print_nice_anova_tables(anova_output=anova(select_mod,test="LRT"), title_is="Anova on CTmax model", filename=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/CTmax_seleced_model", mean_center_var, "_Anova.doc"))
'AAnalysis of Variance Table

Response: ToD_mean
                      Df Sum Sq Mean Sq  F value    Pr(>F)
Population             3  3.821  1.2737   6.0652 0.0005093 ***
Treatment              1 27.673 27.6730 131.7797 < 2.2e-16 ***
Mass_mean              1  1.813  1.8130   8.6337 0.0035578 **
Fed_mean               1  1.134  1.1340   5.4002 0.0208070 *
D_21_mean              1  1.091  1.0906   5.1936 0.0233779 *
Population:Treatment   3  7.739  2.5798  12.2852 1.339e-07 ***
Residuals            298 62.578  0.2100
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
'
summary(select_mod) 
tab_model(select_mod, file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/CTmax_seleced_model",mean_center_var,".doc"),p.style="numeric_stars", show.se=T, show.r2=T,show.stat=T,show.aic=T)
#tab_model(select_mod, file="CTmax_seleced_model_fstat.doc",p.style="numeric_stars", show.se=T, show.r2=T,show.stat=T,show.aic=T)
'Call:
lm(formula = ToD_mean ~ Population * Treatment + Mass_mean +
    Fed_mean + D_21_mean, data = ToD_final)

Residuals:
   Min     1Q Median     3Q    Max
-1.587 -0.320 -0.007  0.266  1.419

Coefficients:
                               Estimate Std. Error t value Pr(>|t|)
(Intercept)                   48.100131   0.225759 213.059  < 2e-16 ***
PopulationKarasburg            1.529472   0.349355   4.378 1.66e-05 ***
PopulationOtavi                0.470146   0.317210   1.482 0.139362
PopulationStampriet            1.212961   0.348519   3.480 0.000576 ***
Treatment                      0.104048   0.010346  10.057  < 2e-16 ***
Mass_mean                     -0.011820   0.006006  -1.968 0.049988 *
Fed_mean                       0.051217   0.024505   2.090 0.037464 *
D_21_mean                      0.003073   0.001370   2.243 0.025643 *
PopulationKarasburg:Treatment -0.081948   0.015328  -5.346 1.79e-07 ***
PopulationOtavi:Treatment     -0.017898   0.013988  -1.279 0.201729
PopulationStampriet:Treatment -0.062396   0.015330  -4.070 6.02e-05 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4583 on 298 degrees of freedom
  (57 observations deleted due to missingness)
Multiple R-squared:  0.4088,    Adjusted R-squared:  0.389
F-statistic: 20.61 on 10 and 298 DF,  p-value: < 2.2e-16
'

# output_type="png"
if (	output_type == "pdf"	)	{ 	
	pdf(paste0("assay_new_setup/derived_data_and_code/phenotypes/results/CT_max_confounding_fact_PLOTS.pdf"),width=4,height=7) 	
}else if ( 	output_type == "png"	){
	png(paste0("assay_new_setup/derived_data_and_code/phenotypes/results/CT_max_confounding_fact_PLOTS.png"),width=4,height=7,units="in",res=400) 	
}
par( mfrow=c(3,1), cex=.8, mar=c(4.5,4.1,1,1), oma=c(0,0,1,0))		

plot(ToD_final$Mass_mean, ToD_final$ToD_mean, xlab="Body mass (mg)",ylab="CTmax")
abline(lm(ToD_final$ToD_mean~ ToD_final$Mass_mean),col="blue")
mtext("a",side=3,line=0.5, adj=-.15,cex=1.2)

plot(ToD_final$Fed_mean, ToD_final$ToD_mean, xlab="Days since feeding",ylab="CTmax")
abline(lm(ToD_final$ToD_mean~ ToD_final$Fed_mean),col="blue")
mtext("b",side=3,line=0.5, adj=-.15,cex=1.2)

plot(ToD_final$D_21_mean, ToD_final$ToD_mean, xlab="Preacclimation Duration (days)",ylab="CTmax")
abline(lm(ToD_final$ToD_mean~ ToD_final$D_21_mean),col="blue")
mtext("c",side=3,line=0.5, adj=-.15,cex=1.2)

dev.off()


# Effect of population on mass:
anova(lm(Mass_mean~ Population, data=ToD_final))
'Analysis of Variance Table

Response: Mass_mean
            Df  Sum Sq Mean Sq F value    Pr(>F)
Population   3  2901.8  967.26  29.494 < 2.2e-16 ***
Residuals  305 10002.5   32.80
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
'
summary(lm(Mass_mean~ Population, data=ToD_final))
'Call:
lm(formula = Mass_mean ~ Population, data = ToD_final)

Residuals:
     Min       1Q   Median       3Q      Max
-11.3628  -4.7494  -0.5704   3.8436  17.3006

Coefficients:
                    Estimate Std. Error t value Pr(>|t|)
(Intercept)          12.7204     0.5727  22.212  < 2e-16 ***
PopulationKarasburg   1.9205     0.9257   2.075  0.03886 *
PopulationOtavi       5.7192     0.8422   6.791 5.83e-11 ***
PopulationStampriet  -2.8170     0.9304  -3.028  0.00267 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 5.727 on 305 degrees of freedom
  (57 observations deleted due to missingness)
Multiple R-squared:  0.2249,    Adjusted R-squared:  0.2172
F-statistic: 29.49 on 3 and 305 DF,  p-value: < 2.2e-16
'
HSD.test(lm(Mass_mean~ Population, data=ToD_final), "Population",group=T)$group
'   Mass_mean groups
Otavi     18.439524      a
Karasburg 14.640812      b
Betta     12.720355      b
Stampriet  9.903345      c
'



 #######
# Population effect on CTnax
cld(emmeans (select_mod,~ Population),Letters=letters)
'NOTE: Results may be misleading due to involvement in interactions
 Population emmean      SE  df lower.CL upper.CL .group
 Karasburg   50.19 0.06420 298    50.06    50.31  a
 Stampriet   50.30 0.07665 298    50.15    50.45  ab
 Betta       50.46 0.04919 298    50.37    50.56   b
 Otavi       50.54 0.06355 298    50.41    50.66   b

Results are given on the : (not the response) scale.
Confidence level used: 0.95
P value adjustment: tukey method for comparing a family of 4 estimates
significance level used: alpha = 0.05
'
write.csv(file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/ctmax_pop_means_different",mean_center_var,".csv"),cld(emmeans (select_mod,~ Population),Letters=letters))
HSD.test(select_mod, "Population",group=T)$group
'           ToD_mean groups
Betta     50.49424      a
Stampriet 50.45086      a
Otavi     50.40761      a
Karasburg 50.18863      b
'

# Treatment effect on CTmax
CLD(emmeans (select_mod,~ Treatment,at=list(Treatment=c(15,19,23,25,29))),Letters=letters)
'NOTE: Results may be misleading due to involvement in interactions
 Treatment emmean      SE  df lower.CL upper.CL .group
        15  49.93 0.05681 298    49.81    50.04  a
        19  50.18 0.03487 298    50.11    50.25   b
        23  50.43 0.02795 298    50.38    50.49    c
        25  50.56 0.03398 298    50.49    50.63     d
        29  50.82 0.05553 298    50.71    50.92      e

Results are averaged over the levels of: Population
Results are given on the : (not the response) scale.
Confidence level used: 0.95
P value adjustment: tukey method for comparing a family of 5 estimates
significance level used: alpha = 0.05
'
HSD.test(select_mod, "Treatment",group=T)$group
' ToD_mean groups
29 50.83891      a
25 50.50632      b
23 50.49014      b
19 50.31301      b
15 49.90298      c
'

# INT effects on CTmax
CLD(emtrends(select_mod, ~Population, var="Treatment"),Letters=letters)
cld(emtrends(select_mod, ~Population, var="Treatment"),Letters=letters)
tab_df(cld(emtrends(select_mod, ~Population,var="Treatment"),Letters=letters), digits=3,
        file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/ctmax_pop_trends_different",mean_center_var,".doc"),
        title="CTmax - Trendlines significance groups", footnote="Confidence level used: 0.95")
'  Population Treatment.trend     SE  df lower.CL upper.CL .group
 Karasburg           0.0221 0.0131 298 -0.00364   0.0478  a
 Stampriet           0.0417 0.0125 298  0.01696   0.0663  a
 Otavi               0.0862 0.0112 298  0.06406   0.1082   b
 Betta               0.1040 0.0103 298  0.08369   0.1244   b

Trends are based on the : (transformed) scale
Confidence level used: 0.95
P value adjustment: tukey method for comparing a family of 4 estimates
significance level used: alpha = 0.05
'
write.csv(file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/ctmax_pop_trends_different",mean_center_var,".csv"),CLD(emtrends(select_mod, ~Population,var="Treatment"),Letters=letters) )
# are trends for populations different from zero?
summary(emtrends(select_mod, ~Population,var="Treatment"),infer=c(TRUE,TRUE),null=0)
tab_df(summary(emtrends(select_mod, ~Population,var="Treatment"),infer=c(TRUE,TRUE),null=0), digits=3,
        file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/ctmax_pop_trends_diff0",mean_center_var,".doc"),
        title="CTmax - Trendlines diverging from zero", footnote="Confidence level used: 0.95")
write.csv(file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/ctmax_pop_trends_different_from_zero",mean_center_var,".csv"),summary(emtrends(select_mod, ~Population,var="Treatment"),infer=c(TRUE,TRUE),null=0))
'Population Treatment.trend     SE  df lower.CL upper.CL t.ratio p.value
 Betta               0.1040 0.0103 298  0.08369   0.1244 10.057  <.0001 ***
 Karasburg           0.0221 0.0131 298 -0.00364   0.0478  1.690  0.0921 
 Otavi               0.0862 0.0112 298  0.06406   0.1082  7.674  <.0001 ***
 Stampriet           0.0417 0.0125 298  0.01696   0.0663  3.319  0.0010 **

Trends are based on the : (transformed) scale 
Confidence level used: 0.95 '
pairs(emtrends(select_mod, ~Population,var="Treatment"),Letters=letters)
'contrast              estimate     SE  df t.ratio p.value
 Betta - Karasburg       0.0819 0.0153 298  5.346  <.0001 ***
 Betta - Otavi           0.0179 0.0140 298  1.279  0.5767 
 Betta - Stampriet       0.0624 0.0153 298  4.070  0.0003 **
 Karasburg - Otavi      -0.0641 0.0158 298 -4.053  0.0004 **
 Karasburg - Stampriet  -0.0196 0.0171 298 -1.145  0.6620 
 Otavi - Stampriet       0.0445 0.0159 298  2.806  0.0272 *

P value adjustment: tukey method for comparing a family of 4 estimates'


# pop effects of treatment
cld(emmeans (select_mod,~ Treatment,at=list(Treatment=c(15,19,23,25,29)),by="Population"),Letters=letters)
'Population = Betta:
 Treatment emmean      SE  df lower.CL upper.CL .group
        15  49.73 0.09159 298    49.55    49.91  a
        19  50.15 0.06079 298    50.03    50.27   b
        23  50.56 0.04926 298    50.47    50.66    c
        25  50.77 0.05547 298    50.66    50.88     d
        29  51.19 0.08335 298    51.02    51.35      e

Population = Karasburg:
 Treatment emmean      SE  df lower.CL upper.CL .group
        15  50.03 0.11216 298    49.81    50.25  a
        19  50.12 0.07546 298    49.97    50.27  a
        23  50.21 0.06544 298    50.08    50.34  a
        25  50.25 0.07502 298    50.11    50.40  a
        29  50.34 0.11147 298    50.12    50.56  a

Population = Otavi:
 Treatment emmean      SE  df lower.CL upper.CL .group
        15  49.93 0.09329 298    49.75    50.12  a
        19  50.28 0.06723 298    50.15    50.41   b
        23  50.62 0.06612 298    50.49    50.75    c
        25  50.80 0.07622 298    50.65    50.95     d
        29  51.14 0.10821 298    50.93    51.35      e

Population = Stampriet:
 Treatment emmean      SE  df lower.CL upper.CL .group
        15  50.01 0.11847 298    49.78    50.24  a
        19  50.18 0.08651 298    50.01    50.35   b
        23  50.34 0.07728 298    50.19    50.49    c
        25  50.43 0.08429 298    50.26    50.59     d
        29  50.59 0.11468 298    50.37    50.82      e

Confidence level used: 0.95
P value adjustment: tukey method for comparing a family of 5 estimates
significance level used: alpha = 0.05
NOTE: Compact letter displays can be misleading
      because they show NON-findings rather than findings.
      Consider using 'pairs()', 'pwpp()', or 'pwpm()' instead.
'

#Betta only ( on CTnax)
modde <- lm(ToD_mean ~ Treatment + Mass_mean + Fed_mean + D_21_mean,data=ToD_final[which(ToD_final$Population=="Betta"),])
HSD.test(modde, "Treatment",group=T)$group
'   ToD_mean groups
29 51.19486      a
25 50.67747      b
23 50.53587      b
19 50.32594      b
15 49.73707      c
'

#Stampriet only  ( on CTnax)
HSD.test(lm(ToD_mean ~ Treatment + Mass_mean + Fed_mean + D_21_mean,data=ToD_final[which(ToD_final$Population=="Stampriet"),]), "Treatment",group=T)$group
'ToD_mean groups
23 50.81370      a
29 50.55229      a
25 50.46736     ab
19 50.44929     ab
15 50.04738      b
'

#Karasburg only  ( on CTnax)
HSD.test(lm(ToD_mean ~ Treatment + Mass_mean + Fed_mean + D_21_mean,data=ToD_final[which(ToD_final$Population=="Karasburg"),]), "Treatment",group=T)$group
'ToD_mean groups
19 50.44506      a
23 50.25111     ab
25 50.20389     ab
29 50.15861     ab
15 49.88814      b
'

#Otavi only  ( on CTnax)
HSD.test(lm(ToD_mean ~ Treatment + Mass_mean + Fed_mean + D_21_mean,data=ToD_final[which(ToD_final$Population=="Otavi"),]), "Treatment",group=T)$group
'ToD_mean groups
29 51.10125      a
25 50.54842      b
23 50.38296     bc
19 50.11241     bc
15 49.98573      c
'


####################################
#------------- CCR-----------------#
####################################
CCR_final <- read.csv(file="assay_new_setup/derived_data_and_code/phenotypes/data/CCR/CCR_nest-temp_mean.csv", header=T,stringsAsFactors=T) 

selec_mod_estimates_noMeanC <- lm(CCR_mean ~ Population:Treatment + Population-1 + Mass_mean + Fed, data=CCR_final)
coef(selec_mod_estimates_noMeanC)
write.csv(file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/ccr_select_mod_NOTMEANC_","intc_slope.csv"),coef(selec_mod_estimates_noMeanC))
'    coef(selec_mod_estimates_noMeanC)
              PopulationBetta           PopulationKarasburg 
                  12.30775464                    7.63258203 
              PopulationOtavi           PopulationStampriet 
                  12.78086642                    9.34830560 
                    Mass_mean                           Fed 
                  -0.04147644                   -0.30967295 
    PopulationBetta:Treatment PopulationKarasburg:Treatment 
                   0.02775537                    0.20059241 
    PopulationOtavi:Treatment PopulationStampriet:Treatment 
                   0.03271787                    0.13833830 '


if(Mean_Center) { 
  CCR_final$Treatment_not_centered <- CCR_final$Treatment
  CCR_final$Treatment <- CCR_final$Treatment-mean(CCR_final$Treatment) }


full_mod <- lm(CCR_mean ~ Population*Treatment + Mass_mean + Fed + D_21 , data=CCR_final)
select_mod <- step(full_mod)
'Start:  AIC=490.97
CCR_mean ~ Population * Treatment + Mass_mean + Fed + D_21

                       Df Sum of Sq    RSS    AIC
- D_21                  1     0.754 1367.5 489.15
<none>                              1366.8 490.97
- Mass_mean             1     9.268 1376.1 491.20
- Population:Treatment  3    38.024 1404.8 494.03
- Fed                   1    37.179 1404.0 497.83

Step:  AIC=489.15
CCR_mean ~ Population + Treatment + Mass_mean + Fed + Population:Treatment

                       Df Sum of Sq    RSS    AIC
<none>                              1367.5 489.15
- Mass_mean             1     9.632 1377.2 489.47
- Population:Treatment  3    37.906 1405.5 492.18
- Fed                   1    36.648 1404.2 495.88
'


anova(select_mod, test="LRT")
print_nice_anova_tables(anova_output=anova(select_mod,test="LRT"), title_is="Anova on CCRTemp model", filename=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/CCR_seleced_model", mean_center_var, "_Anova.doc"))
'Analysis of Variance Table

Response: CCR_mean
                      Df  Sum Sq Mean Sq F value   Pr(>F)
Population             3   48.43  16.144  3.7777 0.010920 *
Treatment              1   22.64  22.638  5.2973 0.022001 *
Mass_mean              1    2.77   2.769  0.6480 0.421433
Fed                    1   36.79  36.789  8.6086 0.003588 **
Population:Treatment   3   37.91  12.635  2.9567 0.032615 *
Residuals            320 1367.54   4.274
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
'
summary(select_mod)
tab_model(select_mod, file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/CCR_seleced_model",mean_center_var,".doc"),p.style="numeric_stars", show.se=T, show.r2=T,show.stat=T,show.aic=T,show.fstat=T)
#tab_model(select_mod, file="CCR_seleced_model_fstat.doc",p.style="numeric_stars", show.se=T, show.r2=T,show.stat=T,show.aic=T,show.fstat=T)
'Call:
lm(formula = CCR_mean ~ Population + Treatment + Mass_mean +
    Fed + Population:Treatment, data = CCR_final)

Residuals:
    Min      1Q  Median      3Q     Max
-6.4318 -1.3124 -0.0169  1.1538  7.2381

Coefficients:
                               Estimate Std. Error t value Pr(>|t|)
(Intercept)                   12.307755   1.003979  12.259  < 2e-16 ***
PopulationKarasburg           -4.675173   1.563557  -2.990  0.00301 **
PopulationOtavi                0.473112   1.375415   0.344  0.73109
PopulationStampriet           -2.959449   1.543165  -1.918  0.05603 .
Treatment                      0.027755   0.047151   0.589  0.55651
Mass_mean                     -0.041476   0.027628  -1.501  0.13427
Fed                           -0.309673   0.105749  -2.928  0.00365 **
PopulationKarasburg:Treatment  0.172837   0.069032   2.504  0.01279 *
PopulationOtavi:Treatment      0.004962   0.060463   0.082  0.93464
PopulationStampriet:Treatment  0.110583   0.068070   1.625  0.10524
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 2.067 on 320 degrees of freedom
  (41 observations deleted due to missingness)
Multiple R-squared:  0.09797,   Adjusted R-squared:  0.0726
F-statistic: 3.862 on 9 and 320 DF,  p-value: 0.0001159
'

#output_type="pdf"
if (	output_type == "pdf"	)	{ 	
	pdf(paste0("assay_new_setup/derived_data_and_code/phenotypes/results/CCR_confounding_fact_PLOTS.pdf"),width=4,height=7) 	
}else if ( 	output_type == "png"	){
	png(paste0("assay_new_setup/derived_data_and_code/phenotypes/results/CCR_confounding_fact_PLOTS.png"),width=4,height=7,units="in",res=400) 	
}
par( mfrow=c(3,1), cex=.8, mar=c(4.5,4.1,1,1), oma=c(0,0,1,0))		

plot(CCR_final$Mass_mean, CCR_final$CCR_mean, xlab="Body mass (mg)",ylab="CCRTemp")
abline(lm(CCR_final$CCR_mean~ CCR_final$Mass_mean),col="blue")
mtext("a",side=3,line=0.5, adj=-.15,cex=1.2)

plot(CCR_final$Fed, CCR_final$CCR_mean, xlab="Days since feeding",ylab="CCRTemp")
abline(lm(CCR_final$CCR_mean~ CCR_final$Fed),col="blue")
mtext("b",side=3,line=0.5, adj=-.15,cex=1.2)

plot(CCR_final$D_21, CCR_final$CCR_mean, xlab="Preacclimation Duration (days)",ylab="CCRTemp")
abline(lm(CCR_final$CCR_mean~ CCR_final$D_21),col="blue")
mtext("c",side=3,line=0.5, adj=-.15,cex=1.2)

dev.off()

anova(lm(Mass_mean~ Population, data=CCR_final))
'Analysis of Variance Table

Response: Mass_mean
            Df Sum Sq Mean Sq F value    Pr(>F)
Population   3 3024.9 1008.31  33.409 < 2.2e-16 ***
Residuals  328 9899.4   30.18
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
'


 #######
# Population effect on CCR
cld(emmeans (select_mod,~ Population),Letters=letters)
'NOTE: Results may be misleading due to involvement in interactions
 Population emmean    SE  df lower.CL upper.CL .group
 Karasburg    10.8 0.275 320     10.2     11.3  a
 Stampriet    11.1 0.278 320     10.6     11.7  a
 Betta        11.6 0.212 320     11.2     12.0  ab
 Otavi        12.2 0.235 320     11.7     12.7   b

Confidence level used: 0.95
P value adjustment: tukey method for comparing a family of 4 estimates
significance level used: alpha = 0.05
'
write.csv(file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/ccr_pop_means_different",mean_center_var,".csv"),cld(emmeans (select_mod,~ Population),Letters=letters))

HSD.test(select_mod, "Population",group=T)$group
'              CCR_mean groups
Otavi     12.02960      a
Betta     11.57844     ab
Stampriet 11.24236     ab
Karasburg 10.99911      b

'

# Treatment effect on CTmax
cld(emmeans (select_mod,~ Treatment,at=list(Treatment=c(15,19,23,25,29))),Letters=letters)
'NOTE: Results may be misleading due to involvement in interactions
 Treatment emmean    SE  df lower.CL upper.CL .group
        15   10.7 0.258 320     10.2     11.2  a
        19   11.1 0.157 320     10.8     11.4   b
        23   11.5 0.119 320     11.3     11.8    c
        25   11.7 0.144 320     11.4     12.0     d
        29   12.1 0.240 320     11.6     12.6      e

Results are averaged over the levels of: Population
Confidence level used: 0.95
P value adjustment: tukey method for comparing a family of 5 estimates
significance level used: alpha = 0.05
'
HSD.test(select_mod, "Treatment",group=T)$group
' CCR_mean groups
23 12.53415      a
25 12.08058     ab
29 11.17136     bc
19 11.03899      c
15 10.88735      c
'

# INT effects on CTmax
cld(emtrends(select_mod, ~Population, var="Treatment"),Letters=letters)
tab_df(cld(emtrends(select_mod, ~Population, var="Treatment"),Letters=letters), digits=3,
        file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/ccr_pop_trends_different",mean_center_var,".doc"),
        title="CCRTemp - Trendlines significance groups", footnote="Confidence level used: 0.95")
'  Population Treatment.trend     SE  df lower.CL upper.CL .group
 Betta               0.0278 0.0472 320  -0.0650    0.121  a
 Otavi               0.0327 0.0471 320  -0.0599    0.125  a
 Stampriet           0.1383 0.0555 320   0.0292    0.247  a
 Karasburg           0.2006 0.0589 320   0.0848    0.316  a

Confidence level used: 0.95
P value adjustment: tukey method for comparing a family of 4 estimates
significance level used: alpha = 0.05
'
# are trends for populations different from zero?
summary(emtrends(select_mod, ~Population,var="Treatment"),infer=c(TRUE,TRUE),null=0)
write.csv(file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/ccr_pop_trends_different_from_zero",mean_center_var,".csv"),summary(emtrends(select_mod, ~Population,var="Treatment"),infer=c(TRUE,TRUE),null=0))
tab_df(summary(emtrends(select_mod, ~Population,var="Treatment"),infer=c(TRUE,TRUE),null=0), digits=3,
        file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/ccr_pop_trends_diff0",mean_center_var,".doc"),
        title="CCRTemp - Trendlines diverging from zero", footnote="Confidence level used: 0.95")
'Population Treatment.trend     SE  df lower.CL upper.CL t.ratio p.value
 Betta               0.0278 0.0472 320  -0.0650    0.121 0.589   0.5565 
 Karasburg           0.2006 0.0589 320   0.0848    0.316 3.408   0.0007 **
 Otavi               0.0327 0.0471 320  -0.0599    0.125 0.695   0.4877 
 Stampriet           0.1383 0.0555 320   0.0292    0.247 2.494   0.0131 *

Confidence level used: 0.95 '
pairs(emtrends(select_mod, ~Population,var="Treatment"),Letters=letters)
' contrast              estimate     SE  df t.ratio p.value
 Betta - Karasburg     -0.17284 0.0690 320 -2.504  0.0612
 Betta - Otavi         -0.00496 0.0605 320 -0.082  0.9998
 Betta - Stampriet     -0.11058 0.0681 320 -1.625  0.3662
 Karasburg - Otavi      0.16787 0.0690 320  2.433  0.0730
 Karasburg - Stampriet  0.06225 0.0759 320  0.820  0.8449
 Otavi - Stampriet     -0.10562 0.0680 320 -1.552  0.4076
'
write.csv(file=paste0("assay_new_setup/derived_data_and_code/phenotypes/results/ccr_pop_trends_different",mean_center_var,".csv"),cld(emtrends(select_mod, ~Population,var="Treatment"),Letters=letters))


#populationwise acclimation effects
cld(emmeans (select_mod,~ Treatment,at=list(Treatment=c(15,19,23,25,29)),by="Population"),Letters=letters)
'Population = Betta:
 Treatment emmean    SE  df lower.CL upper.CL .group
        15  11.43 0.409 320    10.62     12.2  a
        19  11.54 0.266 320    11.02     12.1  a
        23  11.65 0.213 320    11.23     12.1  a
        25  11.71 0.243 320    11.23     12.2  a
        29  11.82 0.375 320    11.08     12.6  a

Population = Karasburg:
 Treatment emmean    SE  df lower.CL upper.CL .group
        15   9.35 0.502 320     8.36     10.3  a
        19  10.15 0.331 320     9.50     10.8   b
        23  10.95 0.280 320    10.40     11.5    c
        25  11.35 0.323 320    10.72     12.0     d
        29  12.16 0.490 320    11.19     13.1      e

Population = Otavi:
 Treatment emmean    SE  df lower.CL upper.CL .group
        15  11.98 0.364 320    11.26     12.7  a
        19  12.11 0.248 320    11.62     12.6  a
        23  12.24 0.248 320    11.75     12.7  a
        25  12.30 0.296 320    11.72     12.9  a
        29  12.44 0.440 320    11.57     13.3  a

Population = Stampriet:
 Treatment emmean    SE  df lower.CL upper.CL .group
        15  10.13 0.512 320     9.12     11.1  a
        19  10.68 0.346 320    10.00     11.4  a
        23  11.24 0.275 320    10.69     11.8  a
        25  11.51 0.301 320    10.92     12.1  a
        29  12.07 0.441 320    11.20     12.9  a

Confidence level used: 0.95
P value adjustment: tukey method for comparing a family of 5 estimates
significance level used: alpha = 0.05
NOTE: Compact letter displays can be misleading
      because they show NON-findings rather than findings.
      Consider using "pairs(), pwpp(), or pwpm()" instead.
'	  
#Betta only
HSD.test(lm(CCR_mean ~ Treatment + Mass_mean + Fed,data=CCR_final[which(CCR_final$Population=="Betta"),]), "Treatment",group=T)$group
' CCR_mean groups
23 12.18571      a
19 11.69224      a
29 11.36637      a
15 11.34620      a
25 11.30170      a
'
#Stampriet only
HSD.test(lm(CCR_mean ~ Treatment + Mass_mean + Fed,data=CCR_final[which(CCR_final$Population=="Stampriet"),]), "Treatment",group=T)$group
'    CCR_mean groups
23 12.681539      a
25 12.235220      a
19 10.809379     ab
29 10.727436     ab
15  9.680809      b
'
#Karasburg only
HSD.test(lm(CCR_mean ~ Treatment + Mass_mean + Fed,data=CCR_final[which(CCR_final$Population=="Karasburg"),]), "Treatment",group=T)$group
'    CCR_mean groups
25 12.774699      a
23 12.194087     ab
29 10.766852    abc
19 10.135079     bc
15  9.335449      c
'

#Otavi only
HSD.test(lm(CCR_mean ~ Treatment + Mass_mean + Fed,data=CCR_final[which(CCR_final$Population=="Otavi"),]), "Treatment",group=T)$group
'   CCR_mean groups
23 12.96206      a
25 12.35626     ab
15 12.15796     ab
29 11.52531     ab
19 11.12949      b
'

# plot for mass, pop and mass, treatment
output_type = "png"
ymax=round_any(max(c(ToD_final$Mass_mean,CCR_final$Mass_mean), na.rm=T),6,f=ceiling)
if (	output_type == "pdf"	)	{ 	
	pdf(paste0("assay_new_setup/derived_data_and_code/phenotypes/results/Mass_pop_effect1.pdf"),width=7,height=5) 	
}else if ( 	output_type == "png"	){
	png(paste0("assay_new_setup/derived_data_and_code/phenotypes/results/Mass_pop_effect1.png"),width=7,height=5,units="in",res=400) 	
}
par( mfrow=c(2,2),cex=.7, mar=c(2.5,2.5,1,1),oma=c(2,4,1,0), xpd=NA)		
plot(ToD_final$Population, ToD_final$Mass_mean, ylab="Body mass",xlab=NA,col=popcols2, ylim=c(0,ymax), yaxt='n')
axis(2, las=2)
mtext("CTmax", las=3, side=2, line=4.5, xpd=NA)
mtext("a",side=3,line=0.5, adj=-.15,cex=1.2)
post_hoc <- HSD.test(lm(Mass_mean ~ Population ,data=ToD_final), "Population",group=T)$group
post_hoc_letters <- post_hoc[order(rownames(post_hoc)),]$groups
text(x=seq(1,4,1), y=ymax, labels=post_hoc_letters)

plot(as.factor(ToD_final$Treatment_not_centered), ToD_final$Mass_mean, ylab=NA, xlab=NA, col=heat.colors(5,rev=T), ylim=c(0,ymax), yaxt='n')
axis(2, las=2)
mtext("b",side=3,line=0.5, adj=-.15,cex=1.2)
post_hoc <- HSD.test(lm(Mass_mean ~ Treatment_not_centered ,data=ToD_final), "Treatment",group=T)$group
post_hoc_letters <- post_hoc[order(rownames(post_hoc)),]$groups
text(x=seq(1,5,1), y=ymax, labels=post_hoc_letters)

plot(CCR_final$Population, CCR_final$Mass_mean, ylab="Body mass",xlab="Population",col=popcols2, ylim=c(0,ymax), yaxt='n')
axis(2, las=2)
mtext("CCRTemp", las=3, side=2, line=4.5, xpd=NA)
mtext("c",side=3,line=0.5, adj=-.15,cex=1.2)
post_hoc <- HSD.test(lm(Mass_mean ~ Population ,data=CCR_final), "Population",group=T)$group
post_hoc_letters <- post_hoc[order(rownames(post_hoc)),]$groups
text(x=seq(1,4,1), y=ymax, labels=post_hoc_letters)

plot(as.factor(CCR_final$Treatment_not_centered), CCR_final$Mass_mean, ylab=NA, xlab=expression(paste("Acclimation Temperature (", degree,"C)")),col=heat.colors(5,rev=T), ylim=c(0,ymax), yaxt='n')
axis(2, las=2)
post_hoc <- HSD.test(lm(Mass_mean ~ Treatment_not_centered ,data=CCR_final), "Treatment",group=T)$group
post_hoc_letters <- post_hoc[order(rownames(post_hoc)),]$groups
mtext("d",side=3,line=0.5, adj=-.15,cex=1.2)
text(x=seq(1,5,1), y=ymax, labels=post_hoc_letters)

dev.off()

