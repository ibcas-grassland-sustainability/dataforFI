getwd()
setwd("C:/temp/Debt")
community<-read.csv("Ca_for_R_total_1.csv", header =TRUE)

library(nlme)#��pvalue�������﷨��һ����ע��
library(lme4)#lme4�������ṩp-value

#����R2�ķ���
r2.corr.mer <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}

#AGB
AGB_1=lmer(AGB ~ tb + (tb | local), community) #�ؾ��б�ʶ��������
AGB_2=lmer(AGB ~ tb + (1 | local), community)#���ؾ��������
AGB_3=lmer(AGB ~ tb + (0 + tb | local), community)#��б���������
summary(AGB_1)
summary(AGB_2)
summary(AGB_3)
anova(AGB_1,AGB_2,AGB_3)#������������ǽؾ��������� model:ABG_Ra_2
coef(AGB_2) #intercept for each local

AGB_2_P= lme (AGB ~ tb, random= ~ 1|local, data=community,method="REML")#�����ŷ��̵�p-value
summary(AGB_2_P)
anova(AGB_2_P)

r2.corr.mer(AGB_2)

#PR+PB
PRB_1=lmer(PRB ~ tb + (tb | local), community) #�ؾ��б�ʶ��������
PRB_2=lmer(PRB ~ tb + (1 | local), community)#���ؾ��������
PRB_3=lmer(PRB ~ tb + (0 + tb | local), community)#��б���������
summary(PRB_1)
summary(PRB_2)
summary(PRB_3)
anova(PRB_1,PRB_2,PRB_3)#������������ǽؾ��������� model:ABG_Ra_2
coef(PRB_1) #intercept and slope for each local

PRB_1_P= lme (PRB ~ tb, random= ~ 1+tb|local, data=community,method="REML")#�����ŷ��̵�p-value
summary(PRB_1_P)
anova(PRB_1_P)

r2.corr.mer(PRB_1)

#Leymus+Stipa
LS_1=lmer(LS ~ tb + (tb | local), community) #�ؾ��б�ʶ��������
LS_2=lmer(LS ~ tb + (1 | local), community)#���ؾ��������
LS_3=lmer(LS ~ tb + (0 + tb | local), community)#��б���������
summary(LS_1)
summary(LS_2)
summary(LS_3)
anova(LS_1,LS_2,LS_3)#������������ǽؾ��������� model:ABG_Ra_2
coef(LS_1) #intercept for each local

LS_1_P= lme (LS ~ tb, random= ~ 1 |local, data=community,method="REML")#�����ŷ��̵�p-value
summary(LS_1_P)
anova(LS_1_P)

r2.corr.mer(LS_1)

#community height
CH_1=lmer(CH ~ tb + (tb | local), community) #�ؾ��б�ʶ��������
CH_2=lmer(CH ~ tb + (1 | local), community)#���ؾ��������
CH_3=lmer(CH ~ tb + (0 + tb | local), community)#��б���������
summary(CH_1)
summary(CH_2)
summary(CH_3)
anova(CH_1,CH_2,CH_3)#������������ǽؾ��������� model:ABG_Ra_2
coef(CH_1) #intercept for each local

CH_1_P= lme (CH ~ 1+tb, random= ~ 1|local, data=community,method="REML")#�����ŷ��̵�p-value
summary(CH_1_P)
anova(CH_1_P)

r2.corr.mer(CH_1)

#species number  ���������,�ݶ����������������
SN_1=lmer(SN ~ tb + (tb | local), community) #�ؾ��б�ʶ��������
SN_2=lmer(SN ~ tb + (1 | local), community)#���ؾ��������
SN_3=lmer(SN ~ tb + (0 + tb | local), community)#��б���������
summary(SN_1)
summary(SN_2)
summary(SN_3)
anova(SN_1,SN_2,SN_3)#������������ǽؾ��������� model:ABG_Ra_2
coef(SN_3) #intercept for each local

SN_1_P= lme (SN ~ tb, random= ~ 1+tb|local, data=community,method="REML")#�����ŷ��̵�p-value
summary(SN_1_P)
anova(SN_1_P)

r2.corr.mer(SN_2)

#root 0-10cm
Ra_1=lmer(roota ~ tb + (tb | local), community) #�ؾ��б�ʶ��������
Ra_2=lmer(roota ~ tb + (1 | local), community)#���ؾ��������
Ra_3=lmer(roota ~ tb + (0 + tb | local), community)#��б���������
summary(Ra_1)
summary(Ra_2)
summary(Ra_3)
anova(Ra_1,Ra_2,Ra_3)#������������ǽؾ��������� model:ABG_Ra_2
coef(Ra_2) #intercept for each local

Ra_1_P= lme (roota ~ tb, random= ~ 1|local, data=community,method="REML")#�����ŷ��̵�p-value
summary(Ra_1_P)
anova(Ra_1_P)

r2.corr.mer(Ra_1)

#root 10-30cm
Rb_1=lmer(rootb ~ tb + (tb | local), community) #�ؾ��б�ʶ��������
Rb_2=lmer(rootb ~ tb + (1 | local), community)#���ؾ��������
Rb_3=lmer(rootb ~ tb + (0 + tb | local), community)#��б���������
summary(Rb_1)
summary(Rb_2)
summary(Rb_3)
anova(Rb_1,Rb_2,Rb_3)#������������ǽؾ��������� model:ABG_Ra_2
coef(Rb_2) #intercept for each local

Rb_2_P= lme (rootb ~ tb, random= ~ 1|local, data=community,method="REML")#�����ŷ��̵�p-value
summary(Rb_2_P)
anova(Rb_2_P)

r2.corr.mer(Rb_2)

#SOC 0-10cm
Soc_a_1=lmer(soca ~ tb + (tb | local), community) #�ؾ��б�ʶ��������
Soc_a_2=lmer(soca ~ tb + (1 | local), community)#���ؾ��������
Soc_a_3=lmer(soca ~ tb + (0 + tb | local), community)#��б���������
summary(Soc_a_1)
summary(Soc_a_2)
summary(Soc_a_3)
anova(Soc_a_1,Soc_a_2,Soc_a_3)#������������ǽؾ��������� model:ABG_Ra_2
coef(Soc_a_2) #intercept for each local

Soc_a_2_P= lme (soca ~ tb, random= ~ 1|local, data=community,method="REML")#�����ŷ��̵�p-value
summary(Soc_a_2_P)
anova(Soc_a_2_P)

r2.corr.mer(Soc_a_2)

#SOC 10-30cm
Soc_b_1=lmer(socb ~ tb + (tb | local), community) #�ؾ��б�ʶ��������
Soc_b_2=lmer(socb ~ tb + (1 | local), community)#���ؾ��������
Soc_b_3=lmer(socb ~ tb + (0 + tb | local), community)#��б���������
summary(Soc_b_1)
summary(Soc_b_2)
summary(Soc_b_3)
anova(Soc_b_1,Soc_b_2,Soc_b_3)#������������ǽؾ��������� model:ABG_Ra_2
coef(Soc_b_2) #intercept for each local

Soc_b_2_P= lme (socb~ tb, random= ~ 1|local, data=community,method="REML")#�����ŷ��̵�p-value
summary(Soc_b_2_P)
anova(Soc_b_2_P)

r2.corr.mer(Soc_b_2)

#TN 0-10cm
TN_a_1=lmer(Tna ~ tb + (tb | local), community) #�ؾ��б�ʶ��������
TN_a_2=lmer(Tna ~ tb + (1 | local), community)#���ؾ��������
TN_a_3=lmer(Tna ~ tb + (0 + tb | local), community)#��б���������
summary(TN_a_1)
summary(TN_a_2)
summary(TN_a_3)
anova(TN_a_1,TN_a_2,TN_a_3)#������������ǽؾ��������� model:ABG_Ra_2
coef(TN_a_1) #intercept for each local

TN_a_2_P= lme (Tna ~ tb, random= ~ 1|local, data=community,method="REML")#�����ŷ��̵�p-value
summary(TN_a_2_P)
anova(TN_a_2_P)

r2.corr.mer(TN_a_2)

#TN 10-30cm
TN_b_1=lmer(Tnb ~ tb + (tb | local), community) #�ؾ��б�ʶ��������
TN_b_2=lmer(Tnb ~ tb + (1 | local), community)#���ؾ��������
TN_b_3=lmer(Tnb ~ tb + (0 + tb | local), community)#��б���������
summary(TN_b_1)
summary(TN_b_2)
summary(TN_b_3)
anova(TN_b_1,TN_b_2,TN_b_3)#������������ǽؾ��������� model:ABG_Ra_2
coef(TN_b_2) #intercept for each local

TN_b_1_P= lme (Tnb ~ tb, random= ~ 1|local, data=community,method="REML")#�����ŷ��̵�p-value
summary(TN_b_1_P)
anova(TN_b_1_P)

r2.corr.mer(TN_b_1)