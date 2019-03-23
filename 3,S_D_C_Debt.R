getwd()
setwd("C:/temp/Debt")
community<-read.csv("Ca_for_R_total.csv", header =TRUE)

library(nlme)#求pvalue，但是语法不一样，注意
library(lme4)#lme4包，不提供p-value

#计算R2的方程
r2.corr.mer <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}

#ABG
ABG_Ra_1=lmer(ABG_Ra ~ tb + (tb | local), community) #截距和斜率都有随机项
ABG_Ra_2=lmer(ABG_Ra ~ tb + (1 | local), community)#仅截距有随机项
ABG_Ra_3=lmer(ABG_Ra ~ tb + (0 + tb | local), community)#仅斜率有随机项
summary(ABG_Ra_1)
summary(ABG_Ra_2)
summary(ABG_Ra_3)
anova(ABG_Ra_1,ABG_Ra_2,ABG_Ra_3)#结果表明仅考虑截距随机项最佳 model:ABG_Ra_2
coef(ABG_Ra_2) #intercept for each local

ABG_Ra_2_P= lme (ABG_Ra ~ tb, random= ~ 1|local, data=community,method="REML")#求最优方程的p-value
summary(ABG_Ra_2_P)
anova(ABG_Ra_2_P)

r2.corr.mer(ABG_Ra_2)

#PR+PB
PRB_Ra_1=lmer(PRB_Ra ~ tb + (tb | local), community) #截距和斜率都有随机项
PRB_Ra_2=lmer(PRB_Ra ~ tb + (1 | local), community)#仅截距有随机项
PRB_Ra_3=lmer(PRB_Ra ~ tb + (0 + tb | local), community)#仅斜率有随机项
summary(PRB_Ra_1)
summary(PRB_Ra_2)
summary(PRB_Ra_3)
anova(PRB_Ra_1,PRB_Ra_2,PRB_Ra_3)#结果表明仅考虑截距随机项最佳 model:ABG_Ra_2
coef(PRB_Ra_1) #intercept and slope for each local

PRB_Ra_1_P= lme (PRB_Ra ~ tb, random= ~ 1+tb|local, data=community,method="REML")#求最优方程的p-value
summary(PRB_Ra_1_P)
anova(PRB_Ra_1_P)

r2.corr.mer(PRB_Ra_1)

#Leymus+Stipa
LS_Ra_1=lmer(LS_Ra ~ tb + (tb | local), community) #截距和斜率都有随机项
LS_Ra_2=lmer(LS_Ra ~ tb + (1 | local), community)#仅截距有随机项
LS_Ra_3=lmer(LS_Ra ~ tb + (0 + tb | local), community)#仅斜率有随机项
summary(LS_Ra_1)
summary(LS_Ra_2)
summary(LS_Ra_3)
anova(LS_Ra_1,LS_Ra_2,LS_Ra_3)#结果表明仅考虑截距随机项最佳 model:ABG_Ra_2
coef(LS_Ra_1) #intercept for each local

LS_Ra_2_P= lme (LS_Ra ~ tb, random= ~ 1|local, data=community,method="REML")#求最优方程的p-value
summary(LS_Ra_2_P)
anova(LS_Ra_2_P)

r2.corr.mer(LS_Ra_2)

#community height
CH_Ra_1=lmer(CH_Ra ~ tb + (tb | local), community) #截距和斜率都有随机项
CH_Ra_2=lmer(CH_Ra ~ tb + (1 | local), community)#仅截距有随机项
CH_Ra_3=lmer(CH_Ra ~ tb + (0 + tb | local), community)#仅斜率有随机项
summary(CH_Ra_1)
summary(CH_Ra_2)
summary(CH_Ra_3)
anova(CH_Ra_1,CH_Ra_2,CH_Ra_3)#结果表明仅考虑截距随机项最佳 model:ABG_Ra_2
coef(CH_Ra_1) #intercept for each local

CH_Ra_1_P= lme (CH_Ra ~ tb, random= ~ 1+tb|local, data=community,method="REML")#求最优方程的p-value
summary(CH_Ra_1_P)
anova(CH_Ra_1_P)

r2.corr.mer(CH_Ra_1)

#species number  结果不显著,暂定生物量的上述结果
SN_Ra_1=lmer(SN_Ra ~ tb + (tb | local), community) #截距和斜率都有随机项
SN_Ra_2=lmer(SN_Ra ~ tb + (1 | local), community)#仅截距有随机项
SN_Ra_3=lmer(SN_Ra ~ tb + (0 + tb | local), community)#仅斜率有随机项
summary(SN_Ra_1)
summary(SN_Ra_2)
summary(SN_Ra_3)
anova(SN_Ra_1,SN_Ra_2,SN_Ra_3)#结果表明仅考虑截距随机项最佳 model:ABG_Ra_2
coef(SN_Ra_3) #intercept for each local

SN_Ra_1_P= lme (SN_Ra ~ tb, random= ~ 1+tb|local, data=community,method="REML")#求最优方程的p-value
summary(SN_Ra_1_P)
anova(SN_Ra_1_P)

r2.corr.mer(SN_Ra_1)

#root 0-10cm
Ra_Ra_1=lmer(Ra_ra ~ tb + (tb | local), community) #截距和斜率都有随机项
Ra_Ra_2=lmer(Ra_ra ~ tb + (1 | local), community)#仅截距有随机项
Ra_Ra_3=lmer(Ra_ra ~ tb + (0 + tb | local), community)#仅斜率有随机项
summary(Ra_Ra_1)
summary(Ra_Ra_2)
summary(Ra_Ra_3)
anova(Ra_Ra_1,Ra_Ra_2,Ra_Ra_3)#结果表明仅考虑截距随机项最佳 model:ABG_Ra_2
coef(Ra_Ra_2) #intercept for each local

Ra_Ra_2_P= lme (Ra_ra ~ tb, random= ~ tb|local, data=community,method="REML")#求最优方程的p-value
summary(Ra_Ra_2_P)
anova(Ra_Ra_2_P)

r2.corr.mer(Ra_Ra_2)

#root 10-30cm
Rb_Ra_1=lmer(Rb_ra ~ tb + (tb | local), community) #截距和斜率都有随机项
Rb_Ra_2=lmer(Rb_ra ~ tb + (1 | local), community)#仅截距有随机项
Rb_Ra_3=lmer(Rb_ra ~ tb + (0 + tb | local), community)#仅斜率有随机项
summary(Rb_Ra_1)
summary(Rb_Ra_2)
summary(Rb_Ra_3)
anova(Rb_Ra_1,Rb_Ra_2,Rb_Ra_3)#结果表明仅考虑截距随机项最佳 model:ABG_Ra_2
coef(Rb_Ra_2) #intercept for each local

Rb_Ra_2_P= lme (Rb_ra ~ tb, random= ~ 1|local, data=community,method="REML")#求最优方程的p-value
summary(Rb_Ra_2_P)
anova(Rb_Ra_2_P)

r2.corr.mer(Rb_Ra_2)

#SOC 0-10cm
Soc_a_Ra_1=lmer(Soc_a_ra ~ tb + (tb | local), community) #截距和斜率都有随机项
Soc_a_Ra_2=lmer(Soc_a_ra ~ tb + (1 | local), community)#仅截距有随机项
Soc_a_Ra_3=lmer(Soc_a_ra ~ tb + (0 + tb | local), community)#仅斜率有随机项
summary(Soc_a_Ra_1)
summary(Soc_a_Ra_2)
summary(Soc_a_Ra_3)
anova(Soc_a_Ra_1,Soc_a_Ra_2,Soc_a_Ra_3)#结果表明仅考虑截距随机项最佳 model:ABG_Ra_2
coef(Soc_a_Ra_2) #intercept for each local

Soc_a_Ra_2_P= lme (Soc_a_ra ~ tb, random= ~ 1|local, data=community,method="REML")#求最优方程的p-value
summary(Soc_a_Ra_2_P)
anova(Soc_a_Ra_2_P)

r2.corr.mer(Soc_a_Ra_2)

#SOC 10-30cm
Soc_b_Ra_1=lmer(Soc_b_ra ~ tb + (tb | local), community) #截距和斜率都有随机项
Soc_b_Ra_2=lmer(Soc_b_ra ~ tb + (1 | local), community)#仅截距有随机项
Soc_b_Ra_3=lmer(Soc_b_ra ~ tb + (0 + tb | local), community)#仅斜率有随机项
summary(Soc_b_Ra_1)
summary(Soc_b_Ra_2)
summary(Soc_b_Ra_3)
anova(Soc_b_Ra_1,Soc_b_Ra_2,Soc_b_Ra_3)#结果表明仅考虑截距随机项最佳 model:ABG_Ra_2
coef(Soc_b_Ra_2) #intercept for each local

Soc_b_Ra_2_P= lme (Soc_b_ra ~ tb, random= ~ 1|local, data=community,method="REML")#求最优方程的p-value
summary(Soc_b_Ra_2_P)
anova(Soc_b_Ra_2_P)

r2.corr.mer(Soc_b_Ra_2)

#TN 0-10cm
TN_a_Ra_1=lmer(TN_a_Ra ~ tb + (tb | local), community) #截距和斜率都有随机项
TN_a_Ra_2=lmer(TN_a_Ra ~ tb + (1 | local), community)#仅截距有随机项
TN_a_Ra_3=lmer(TN_a_Ra ~ tb + (0 + tb | local), community)#仅斜率有随机项
summary(TN_a_Ra_1)
summary(TN_a_Ra_2)
summary(TN_a_Ra_3)
anova(TN_a_Ra_1,TN_a_Ra_2,TN_a_Ra_3)#结果表明仅考虑截距随机项最佳 model:ABG_Ra_2
coef(TN_a_Ra_1) #intercept for each local

TN_a_Ra_2_P= lme (TN_a_Ra ~ tb, random= ~ 1|local, data=community,method="REML")#求最优方程的p-value
summary(TN_a_Ra_2_P)
anova(TN_a_Ra_2_P)

r2.corr.mer(TN_a_Ra_2)

#TN 10-30cm
TN_b_Ra_1=lmer(TN_b_Ra ~ tb + (tb | local), community) #截距和斜率都有随机项
TN_b_Ra_2=lmer(TN_b_Ra ~ tb + (1 | local), community)#仅截距有随机项
TN_b_Ra_3=lmer(TN_b_Ra ~ tb + (0 + tb | local), community)#仅斜率有随机项
summary(TN_b_Ra_1)
summary(TN_b_Ra_2)
summary(TN_b_Ra_3)
anova(TN_b_Ra_1,TN_b_Ra_2,TN_b_Ra_3)#结果表明仅考虑截距随机项最佳 model:ABG_Ra_2
coef(TN_b_Ra_2) #intercept for each local

TN_b_Ra_2_P= lme (TN_b_Ra ~ tb, random= ~ 1|local, data=community,method="REML")#求最优方程的p-value
summary(TN_b_Ra_2_P)
anova(TN_b_Ra_2_P)

r2.corr.mer(TN_b_Ra_2)