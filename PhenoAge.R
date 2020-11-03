
#Adaptive responses to SARS-CoV-2 infection linked to accelerated aging measures predict adverse outcomes in patients with severe COVID-19 
#Data Analysis: Alejandro Márquez-Salinas (carlosfermin38@gmail.com); Carlos A. Fermín-Martínez (carlosfermin38@gmail.com);Omar Yaxmehen Bello-Chavolla (oyaxbell@yahoo.com.mx)
#Latest version of Analysis 31-Oct-2020
#Any question regarding analysis contact Omar Yaxmehen Bello-Chavolla


#### Database ####
library(haven); library(tidyverse); library(pROC); library(factoextra); library(rgl); library(gridExtra)
library(FactoMineR);library(survival); library(fpc); library(NbClust); library(ggimage);library(glmnet)
library(ggsci);library(ggpubr); library(survminer); library(cluster); library(ggplotify); library(UpSetR)
library(nortest); library(bestNormalize); library(flextable); library(viridis); library(bioage)

setwd("~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/COVID-Nutricion/PhenoAge")
setwd("C:/Users/Usuario Lab. Datos/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/COVID-Nutricion/PhenoAge")
setwd("C:/Users/facmed/UNIVERSIDAD NACIONAL AUT?NOMA DE M?XICO/OMAR YAXMEHEN BELLO CHAVOLLA - PhenoAge")
nutri<-read.csv("nutri.csv")
nutri$edad_cat[nutri$edad<50]<-0;nutri$edad_cat[nutri$edad>=50 & nutri$edad<=70]<-1;nutri$edad_cat[nutri$edad>70]<-2
nutri$sexo<-ifelse(nutri$sexo=="H",1,0)
nutri$sexo <- factor(nutri$sexo, levels= c(0,1),labels= c("Women", "Male"))

#Comorbidities
nutri$comorb[(nutri$obesidad+nutri$dm2+nutri$epoc+nutri$irc+nutri$vih+
                nutri$hipertension+nutri$cardiovascular+nutri$insuficiencia_hepatica+
                nutri$tabaquismo)>0]<-1;nutri$comorb[(nutri$obesidad+nutri$dm2+nutri$epoc+nutri$irc+nutri$vih+
                                                        nutri$hipertension+nutri$cardiovascular+nutri$insuficiencia_hepatica+
                                                        nutri$tabaquismo)==0]<-0
nutri$comorb1 <- as.numeric(nutri$obesidad+nutri$dm2+nutri$hipertension+
                              nutri$cardiovascular+nutri$insuficiencia_hepatica+nutri$tabaquismo)
nutri$comorb2<-as.numeric(nutri$obesidad+nutri$dm2+nutri$epoc+nutri$irc+nutri$vih+nutri$hipertension+
  nutri$cardiovascular+nutri$insuficiencia_hepatica+nutri$tabaquismo)
nutri$comorb2[nutri$comorb2==0]<-0;nutri$comorb2[nutri$comorb2==1]<-1
nutri$comorb2[nutri$comorb2==2]<-2;nutri$comorb2[nutri$comorb2>=3]<-3

nutri$comorb<-factor(nutri$comorb, labels = c("No-comorb", ">=1 comorb"))
nutri$comorb2<-factor(nutri$comorb2, labels = c("No-comorb", "1 comorb", "2 comorb", ">=3 comorb"))
table(nutri$comorb);table(nutri$comorb2)
sum(is.na(nutri$PhenoAge))
nrow(nutri)

#PhenoAge and PhenoAge Accel
nutri$xb<- (-19.907-0.0336*nutri$alb0+0.0095*(nutri$cr0*88.42)+0.1953*(nutri$GLU/18)+
  0.0954*log(nutri$pcr0)-0.0120*nutri$lin0+0.0268*nutri$VCM+
  0.3306*nutri$ADE+0.00188*nutri$fa0+0.0554*(nutri$leu0/1000)+
  0.0804*nutri$edad)
nutri$M<-1-exp((-1.51714*exp(nutri$xb))/0.0076927)
nutri$PhenoAge <- 141.5 + ((log(-0.00553*log(1-nutri$M)))/(0.09165))
nutri<-nutri[!is.infinite(nutri$PhenoAge),]
m0<-lm(PhenoAge~edad, data=nutri)
nutri$pred<-predict(m0, nutri)
nutri$PhenoAccelAge<-nutri$PhenoAge-nutri$pred
nutri$accel1[nutri$PhenoAccelAge>0]<-1;nutri$accel1[nutri$PhenoAccelAge<=0]<-0

#Outcomes
nutri$critical<-NULL
nutri$critical<-as.numeric((nutri$intubado+nutri$uci)>=1)
nutri$variable2<-nutri$defuncion; nutri$variable2[nutri$variable2==1]<-2
nutri$groups<-nutri$critical+nutri$variable2; nutri$groups[nutri$groups==3]<-2
nutri$def1<-factor(nutri$defuncion, labels = c("Non-lethal", "Lethal"))
nutri$groups1<-factor(nutri$groups, labels = c("Severe", "Critical", "Lethal"))

#Additional variables
nutri$dm2_age40<-as.numeric((as.numeric(nutri$edad<=40)+nutri$dm2)==2)
nutri$pulsepr<-nutri$tas-nutri$tad
nutri$TyG<-log((nutri$GLU*nutri$TG)/2)


#### Study population - Characteristics ####
nutri1<-nutri[!is.na(nutri$PhenoAccelAge),]
nutri2 <- nutri1 %>% filter(nutri$accel1==0)
nutri3 <- nutri1 %>% filter(nutri$accel1==1)
nutri4 <- nutri1 %>% filter(nutri$accel1>=0)
table(nutri4$groups>=1)%>%prop.table()
table(nutri4$defuncion)%>%prop.table()
table(nutri4$groups)%>%prop.table()

#Fisher's exact test should be carried out for the following variables
table(nutri$accel1,nutri$dm2_age40)
table(nutri$accel1,nutri$irc)
table(nutri$accel1,nutri$asma)
table(nutri$accel1,nutri$vih)

# Sexo
Male0 <- table(nutri2$sexo)[2]
pMale0 <- round(((table(nutri2$sexo)%>%prop.table())[2])*100,1)
Male1 <- table(nutri3$sexo)[2]
pMale1 <- round(((table(nutri3$sexo)%>%prop.table())[2])*100,1)
Male2 <- table(nutri4$sexo)[2]
pMale2 <- round(((table(nutri4$sexo)%>%prop.table())[2])*100,1)
p1 <- format.pval(prop.test(x=c(Male0,Male1), n=c(nrow(nutri2),nrow(nutri3)))$p.value, eps = .001, digits = 3) 

# Edad
hist(nutri4$edad);nortest::ad.test(nutri4$edad)
Age0 <- round(summary(nutri2$edad)[3],2)
Age0.1 <- round(summary(nutri2$edad)[2],2)
Age0.3 <- round(summary(nutri2$edad)[5],2)
Age1 <- round(summary(nutri3$edad)[3],2)
Age1.1 <- round(summary(nutri3$edad)[2],2)
Age1.3 <- round(summary(nutri3$edad)[5],2)
Age2 <- round(summary(nutri4$edad)[3],2)
Age2.1 <- round(summary(nutri4$edad)[2],2)
Age2.3 <- round(summary(nutri4$edad)[5],2)
p2 <- format.pval(wilcox.test(nutri4$edad~nutri4$accel1)$p.value, eps = .001, digits = 3) 

tapply(nutri4$edad,nutri4$accel1,summary)


# PhenoAge
hist(nutri4$PhenoAge);nortest::ad.test(nutri4$PhenoAge)
PhenoAge0 <- round(summary(nutri2$PhenoAge)[3],2)
PhenoAge0.1 <- round(summary(nutri2$PhenoAge)[2],2)
PhenoAge0.3 <- round(summary(nutri2$PhenoAge)[5],2)
PhenoAge1 <- round(summary(nutri3$PhenoAge)[3],2)
PhenoAge1.1 <- round(summary(nutri3$PhenoAge)[2],2)
PhenoAge1.3 <- round(summary(nutri3$PhenoAge)[5],2)
PhenoAge2 <- round(summary(nutri4$PhenoAge)[3],2)
PhenoAge2.1 <- round(summary(nutri4$PhenoAge)[2],2)
PhenoAge2.3 <- round(summary(nutri4$PhenoAge)[5],2)
p3<- format.pval(wilcox.test(nutri2$PhenoAge, nutri3$PhenoAge)$p.value, eps = .001, digits = 3)

# Patient status
Severe0 <- table(nutri2$groups1)[1]
pSevere0 <- round(((table(nutri2$groups1)%>%prop.table())[1])*100,1)
Critial0 <- table(nutri2$groups1)[2]
pCritial0 <- round(((table(nutri2$groups1)%>%prop.table())[2])*100,1)
Lethal0 <- table(nutri2$groups1)[3]
pLethal0 <- round(((table(nutri2$groups1)%>%prop.table())[3])*100,1)
Severe1 <- table(nutri3$groups1)[1]
pSevere1 <- round(((table(nutri3$groups1)%>%prop.table())[1])*100,1)
Critial1 <- table(nutri3$groups1)[2]
pCritial1 <- round(((table(nutri3$groups1)%>%prop.table())[2])*100,1)
Lethal1 <- table(nutri3$groups1)[3]
pLethal1 <- round(((table(nutri3$groups1)%>%prop.table())[3])*100,1)
Severe2 <- table(nutri4$groups1)[1]
pSevere2 <- round(((table(nutri4$groups1)%>%prop.table())[1])*100,1)
Critial2 <- table(nutri4$groups1)[2]
pCritial2 <- round(((table(nutri4$groups1)%>%prop.table())[2])*100,1)
Lethal2 <- table(nutri4$groups1)[3]
pLethal2 <- round(((table(nutri4$groups1)%>%prop.table())[3])*100,1)
p4 <- format.pval(prop.test(x=c(Severe0,Severe1), n=c(nrow(nutri2),nrow(nutri3)))$p.value, eps = .001, digits = 3) 
p5 <- format.pval(prop.test(x=c(Critial0,Critial1), n=c(nrow(nutri2),nrow(nutri3)))$p.value, eps = .001, digits = 3) 
p6 <- format.pval(prop.test(x=c(Lethal0,Lethal1), n=c(nrow(nutri2),nrow(nutri3)))$p.value, eps = .001, digits = 3) 

# >= 1 comorb
Comorb0 <- table(nutri2$comorb)[2]
pComorb0 <- round(((table(nutri2$comorb)%>%prop.table())[2])*100,1)
Comorb1 <- table(nutri3$comorb)[2]
pComorb1 <- round(((table(nutri3$comorb)%>%prop.table())[2])*100,1)
Comorb2 <- table(nutri4$comorb)[2]
pComorb2 <- round(((table(nutri4$comorb)%>%prop.table())[2])*100,1)
p7 <- format.pval(prop.test(x=c(Comorb0,Comorb1), n=c(nrow(nutri2),nrow(nutri3)))$p.value, eps = .001, digits = 3) 

# Diabetes
Diabetes0 <- table(nutri2$dm2)[2]
pDiabetes0 <- round(((table(nutri2$dm2)%>%prop.table())[2])*100,1)
Diabetes1 <- table(nutri3$dm2)[2]
pDiabetes1 <- round(((table(nutri3$dm2)%>%prop.table())[2])*100,1)
Diabetes2 <- table(nutri4$dm2)[2]
pDiabetes2 <- round(((table(nutri4$dm2)%>%prop.table())[2])*100,1)
p8 <- format.pval(prop.test(x=c(Diabetes0,Diabetes1), n=c(nrow(nutri2),nrow(nutri3)))$p.value, eps = .001, digits = 3) 

# Early-onset diabetes
EDiabetes0 <- table(nutri2$dm2_age40)[2]
pEDiabetes0 <- round(((table(nutri2$dm2_age40)%>%prop.table())[2])*100,1)
EDiabetes1 <- table(nutri3$dm2_age40)[2]
pEDiabetes1 <- round(((table(nutri3$dm2_age40)%>%prop.table())[2])*100,1)
EDiabetes2 <- table(nutri4$dm2_age40)[2]
pEDiabetes2 <- round(((table(nutri4$dm2_age40)%>%prop.table())[2])*100,1)
p9 <- format.pval(fisher.test(table(nutri$accel1,nutri$dm2_age40))$p.value, eps = .001, digits = 3) 

# Obesity
Obesidad0 <- table(nutri2$obesidad)[2]
pObesidad0 <- round(((table(nutri2$obesidad)%>%prop.table())[2])*100,1)
Obesidad1 <- table(nutri3$obesidad)[2]
pObesidad1 <- round(((table(nutri3$obesidad)%>%prop.table())[2])*100,1)
Obesidad2 <- table(nutri4$obesidad)[2]
pObesidad2 <- round(((table(nutri4$obesidad)%>%prop.table())[2])*100,1)
p10 <- format.pval(prop.test(x=c(Obesidad0,Obesidad1), n=c(nrow(nutri2),nrow(nutri3)))$p.value, eps = .001, digits = 3) 

# Cardiovascular disease
Cardiovascular0 <- table(nutri2$cardiovascular)[2]
pCardiovascular0 <- round(((table(nutri2$cardiovascular)%>%prop.table())[2])*100,1)
Cardiovascular1 <- table(nutri3$cardiovascular)[2]
pCardiovascular1 <- round(((table(nutri3$cardiovascular)%>%prop.table())[2])*100,1)
Cardiovascular2 <- table(nutri4$cardiovascular)[2]
pCardiovascular2 <- round(((table(nutri4$cardiovascular)%>%prop.table())[2])*100,1)
p11 <- format.pval(prop.test(x=c(Cardiovascular0,Cardiovascular1), n=c(nrow(nutri2),nrow(nutri3)))$p.value, eps = .001, digits = 3) 

# Hipertension
Hipertension0 <- table(nutri2$hipertension)[2]
pHipertension0 <- round(((table(nutri2$hipertension)%>%prop.table())[2])*100,1)
Hipertension1 <- table(nutri3$hipertension)[2]
pHipertension1 <- round(((table(nutri3$hipertension)%>%prop.table())[2])*100,1)
Hipertension2 <- table(nutri4$hipertension)[2]
pHipertension2 <- round(((table(nutri4$hipertension)%>%prop.table())[2])*100,1)
p12 <- format.pval(prop.test(x=c(Hipertension0,Hipertension1), n=c(nrow(nutri2),nrow(nutri3)))$p.value, eps = .001, digits = 3) 

# Chronic kidney disease
IRC0 <- table(nutri2$irc)[2]
pIRC0 <- round(((table(nutri2$irc)%>%prop.table())[2])*100,1)
IRC1 <- table(nutri3$irc)[2]
pIRC1 <- round(((table(nutri3$irc)%>%prop.table())[2])*100,1)
IRC2 <- table(nutri4$irc)[2]
pIRC2 <- round(((table(nutri4$irc)%>%prop.table())[2])*100,1)
p13 <- format.pval(fisher.test(table(nutri$accel1,nutri$irc))$p.value, eps = .001, digits = 3) 

# COPD
EPOC0 <- table(nutri2$epoc)[2]
pEPOC0 <- round(((table(nutri2$epoc)%>%prop.table())[2])*100,1)
EPOC1 <- table(nutri3$epoc)[2]
pEPOC1 <- round(((table(nutri3$epoc)%>%prop.table())[2])*100,1)
EPOC2 <- table(nutri4$epoc)[2]
pEPOC2 <- round(((table(nutri4$epoc)%>%prop.table())[2])*100,1)
p14 <- format.pval(fisher.test(table(nutri$accel1,nutri$epoc))$p.value, eps = .001, digits = 3) 

# Asthma
Asma0 <- table(nutri2$asma)[2]
pAsma0 <- round(((table(nutri2$asma)%>%prop.table())[2])*100,1)
Asma1 <- table(nutri3$asma)[2]
pAsma1 <- round(((table(nutri3$asma)%>%prop.table())[2])*100,1)
Asma2 <- table(nutri4$asma)[2]
pAsma2 <- round(((table(nutri4$asma)%>%prop.table())[2])*100,1)
p15 <- format.pval(fisher.test(table(nutri$accel1,nutri$asma))$p.value, eps = .001, digits = 3) 

# Inmmunosupresion
Inmunosupresion0 <- table(nutri2$inmunoupresion)[2]
pInmunosupresion0 <- round(((table(nutri2$inmunoupresion)%>%prop.table())[2])*100,1)
Inmunosupresion1 <- table(nutri3$inmunoupresion)[2]
pInmunosupresion1 <- round(((table(nutri3$inmunoupresion)%>%prop.table())[2])*100,1)
Inmunosupresion2 <- table(nutri4$inmunoupresion)[2]
pInmunosupresion2 <- round(((table(nutri4$inmunoupresion)%>%prop.table())[2])*100,1)
p16 <- format.pval(prop.test(x=c(Inmunosupresion0,Inmunosupresion1), n=c(nrow(nutri2),nrow(nutri3)))$p.value, eps = .001, digits = 3) 

# VIH
VIH0 <- table(nutri2$vih)[2]
pVIH0 <- round(((table(nutri2$vih)%>%prop.table())[2])*100,1)
VIH1 <- table(nutri3$vih)[2]
pVIH1 <- round(((table(nutri3$vih)%>%prop.table())[2])*100,1)
VIH2 <- table(nutri4$vih)[2]
pVIH2 <- round(((table(nutri4$vih)%>%prop.table())[2])*100,1)
p17 <- format.pval(fisher.test(table(nutri$accel1,nutri$vih))$p.value, eps = .001, digits = 3) 

# Smoking
Tabaquismo0 <- table(nutri2$tabaquismo)[2]
pTabaquismo0 <- round(((table(nutri2$tabaquismo)%>%prop.table())[2])*100,1)
Tabaquismo1 <- table(nutri3$tabaquismo)[2]
pTabaquismo1 <- round(((table(nutri3$tabaquismo)%>%prop.table())[2])*100,1)
Tabaquismo2 <- table(nutri4$tabaquismo)[2]
pTabaquismo2 <- round(((table(nutri4$tabaquismo)%>%prop.table())[2])*100,1)
p18 <- format.pval(prop.test(x=c(Tabaquismo0,Tabaquismo1), n=c(nrow(nutri2),nrow(nutri3)))$p.value, eps = .001, digits = 3) 

tab0<-data.frame("Characteristics"=c("Male (%)","Age (years)", "PhenoAge (years)", "Severe (%)", "Critical (%)", "Lethal (%)", ">=1 comorbidity (%)",
                                     "Diabetes Mellitus (%)", "Early-onset Diabetes (%)", "Obesity (%)", "Cardiovascular disease (%)",
                                     "Hipertension (%)", "Chronic kidney disease (%)", "COPD (%)", "Asthma (%)", "Immunosuppression (%)","HIV (%)","Smoking (%)"),
                 "Overall (n=703)"=c(paste0(Male2," (",pMale2,")"),paste0(Age2," (",Age2.1,"-",Age2.3,")"),paste0(PhenoAge2," (",PhenoAge2.1,"-",PhenoAge2.3,")"),
                             paste0(Severe2," (",pSevere2,")"),paste0(Critial2," (",pCritial2,")"),paste0(Lethal2," (",pLethal2,")"),
                             paste0(Comorb2," (",pComorb2,")"),paste0(Diabetes2," (",pDiabetes2,")"),paste0(EDiabetes2," (",pEDiabetes2,")") ,paste0(Obesidad2," (",pObesidad2,")"),
                             paste0(Cardiovascular2," (",pCardiovascular2,")"),paste0(Hipertension2," (",pHipertension2,")"),paste0(IRC2," (",pIRC2,")"),
                             paste0(EPOC2," (",pEPOC2,")"),paste0(Asma2," (",pAsma2,")"),paste0(Inmunosupresion2," (",pInmunosupresion2,")"),
                             paste0(VIH2," (",pVIH2,")"),paste0(Tabaquismo2," (",pTabaquismo2,")")),
                 "PhenoAccelAge>0 (n=420)"=c(paste0(Male0," (",pMale0,")"),paste0(Age0," (",Age0.1,"-",Age0.3,")"),paste0(PhenoAge0," (",PhenoAge0.1,"-",PhenoAge0.3,")"),
                                             paste0(Severe0," (",pSevere0,")"),paste0(Critial0," (",pCritial0,")"),paste0(Lethal0," (",pLethal0,")"),
                                             paste0(Comorb0," (",pComorb0,")"),paste0(Diabetes0," (",pDiabetes0,")"),paste0(EDiabetes0," (",pEDiabetes0,")") ,paste0(Obesidad0," (",pObesidad0,")"),
                                             paste0(Cardiovascular0," (",pCardiovascular0,")"),paste0(Hipertension0," (",pHipertension0,")"),paste0(IRC0," (",pIRC0,")"),
                                             paste0(EPOC0," (",pEPOC0,")"),paste0(Asma0," (",pAsma0,")"),paste0(Inmunosupresion0," (",pInmunosupresion0,")"),
                                             paste0(VIH0," (",pVIH0,")"),paste0(Tabaquismo0," (",pTabaquismo0,")")),
                 "PhenoAccelAge<=0    (n=607)"=c(paste0(Male1," (",pMale1,")"),paste0(Age1," (",Age1.1,"-",Age1.3,")"),paste0(PhenoAge1," (",PhenoAge1.1,"-",PhenoAge1.3,")"),
                                             paste0(Severe1," (",pSevere1,")"),paste0(Critial1," (",pCritial1,")"),paste0(Lethal1," (",pLethal1,")"),
                                             paste0(Comorb1," (",pComorb1,")"),paste0(Diabetes1," (",pDiabetes1,")"),paste0(EDiabetes1," (",pEDiabetes1,")"),paste0(Obesidad1," (",pObesidad1,")"),
                                             paste0(Cardiovascular1," (",pCardiovascular1,")"),paste0(Hipertension1," (",pHipertension1,")"),paste0(IRC1," (",pIRC1,")"),
                                             paste0(EPOC1," (",pEPOC1,")"),paste0(Asma1," (",pAsma1,")"),paste0(Inmunosupresion1," (",pInmunosupresion1,")"),
                                             paste0(VIH1," (",pVIH1,")"),paste0(Tabaquismo1," (",pTabaquismo1,")")),
                 "P-value"=c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18))
tab0<-`names<-`(tab0,c("Characteristics","Overall (n=1068)","Physiological aging (n=631)","Accelerated aging   (n=438)","P-value"))
tab0<-align(flextable(tab0,cwidth = c(2,1.5,1.5,2,1)),align = "center",part = "all")
save_as_docx(tab0,path="tabla0.docx")


#### Cox regression analyses - Lethality ####
nutri1<-nutri[!is.na(nutri$PhenoAccelAge),]
nutri1<-nutri1%>%dplyr::select(critical, FU_time, defuncion, edad,
                              alb0, cr0, fa0, leu0, VCM, ADE, lin0, pcr0, GLU, PhenoAge, PhenoAccelAge, accel1, sexo, comorb,comorb1,comorb2,groups1)
set.seed(123); imp<-mice::mice(nutri1, m=1, maxit=1); nutri1<-complete(imp, "long")
nrow(nutri1)
m1<-coxph(Surv(FU_time, defuncion)~edad, data=nutri1)
cox.zph(m1)
summary(m1)
BIC1<-round(BIC(m1)-BIC(m1),2)
C1<-round(summary(m1)$concordance[1],3)
B1 <- round(summary(m1)$coefficients[1],3)
p1 <- format.pval(summary(m1)$coefficients[5], eps = .001, digits = 3)
HR1 <- round(summary(m1)$conf.int[1],3)
HR1l <- round(summary(m1)$conf.int[3],3)
HR1u <- round(summary(m1)$conf.int[4],3)
HR1 <- paste0(HR1," ","(",HR1l,"-",HR1u,")")

m2<-coxph(Surv(FU_time, defuncion)~GLU, data=nutri1)
cox.zph(m2)
summary(m2)
BIC2<-round(BIC(m1)-BIC(m2),2)
C2<-round(summary(m2)$concordance[1],3)
B2 <- round(summary(m2)$coefficients[1],3)
p2 <- format.pval(summary(m2)$coefficients[5], eps = .001, digits = 3)
HR2 <- round(summary(m2)$conf.int[1],3)
HR2l <- round(summary(m2)$conf.int[3],3)
HR2u <- round(summary(m2)$conf.int[4],3)
HR2 <- paste0(HR2," ","(",HR2l,"-",HR2u,")")

m3<-coxph(Surv(FU_time, defuncion)~pcr0, data=nutri1)
cox.zph(m3)
summary(m3)
BIC3<-round(BIC(m1)-BIC(m3),2)
C3<-round(summary(m3)$concordance[1],3)
B3 <- round(summary(m3)$coefficients[1],3)
p3 <- format.pval(summary(m3)$coefficients[5], eps = .001, digits = 3)
HR3 <- round(summary(m3)$conf.int[1],3)
HR3l <- round(summary(m3)$conf.int[3],3)
HR3u <- round(summary(m3)$conf.int[4],3)
HR3 <- paste0(HR3," ","(",HR3l,"-",HR3u,")")

m4<-coxph(Surv(FU_time, defuncion)~fa0, data=nutri1)
cox.zph(m4)
summary(m4)
BIC4<-round(BIC(m1)-BIC(m4),2)
C4<-round(summary(m4)$concordance[1],3)
B4 <- round(summary(m4)$coefficients[1],3)
p4 <- format.pval(summary(m4)$coefficients[5], eps = .001, digits = 3)
HR4 <- round(summary(m4)$conf.int[1],3)
HR4l <- round(summary(m4)$conf.int[3],3)
HR4u <- round(summary(m4)$conf.int[4],3)
HR4 <- paste0(HR4," ","(",HR4l,"-",HR4u,")")

m5<-coxph(Surv(FU_time, defuncion)~ADE, data=nutri1)
cox.zph(m5)
summary(m5)
BIC5<-round(BIC(m1)-BIC(m5),2)
C5<-round(summary(m5)$concordance[1],3)
B5 <- round(summary(m5)$coefficients[1],3)
p5 <- format.pval(summary(m5)$coefficients[5], eps = .001, digits = 3)
HR5 <- round(summary(m5)$conf.int[1],3)
HR5l <- round(summary(m5)$conf.int[3],3)
HR5u <- round(summary(m5)$conf.int[4],3)
HR5 <- paste0(HR5," ","(",HR5l,"-",HR5u,")")

m6<-coxph(Surv(FU_time, defuncion)~VCM, data=nutri1)
cox.zph(m6)
summary(m6)
BIC6<-round(BIC(m1)-BIC(m6),2)
C6<-round(summary(m6)$concordance[1],3)
B6 <- round(summary(m6)$coefficients[1],3)
p6 <- format.pval(summary(m6)$coefficients[5], eps = .001, digits = 3)
HR6 <- round(summary(m6)$conf.int[1],3)
HR6l <- round(summary(m6)$conf.int[3],3)
HR6u <- round(summary(m6)$conf.int[4],3)
HR6 <- paste0(HR6," ","(",HR6l,"-",HR6u,")")

m7<-coxph(Surv(FU_time, defuncion)~lin0, data=nutri1)
cox.zph(m7)
summary(m7)
BIC7<-round(BIC(m1)-BIC(m7),2)
C7<-round(summary(m7)$concordance[1],3)
B7 <- round(summary(m7)$coefficients[1],3)
p7 <- format.pval(summary(m7)$coefficients[5], eps = .001, digits = 3)
HR7 <- round(summary(m7)$conf.int[1],3)
HR7l <- round(summary(m7)$conf.int[3],3)
HR7u <- round(summary(m7)$conf.int[4],3)
HR7 <- paste0(HR7," ","(",HR7l,"-",HR7u,")")

nutri1$leu_2<-nutri1$leu0/1000
m8<-coxph(Surv(FU_time, defuncion)~leu_2, data=nutri1)
cox.zph(m8)
summary(m8)
BIC8<-round(BIC(m1)-BIC(m8),2)
C8<-round(summary(m8)$concordance[1],3)
B8 <- round(summary(m8)$coefficients[1],3)
p8 <- format.pval(summary(m8)$coefficients[5], eps = .001, digits = 3)
HR8 <- round(summary(m8)$conf.int[1],3)
HR8l <- round(summary(m8)$conf.int[3],3)
HR8u <- round(summary(m8)$conf.int[4],3)
HR8 <- paste0(HR8," ","(",HR8l,"-",HR8u,")")

m9<-coxph(Surv(FU_time, defuncion)~cr0, data=nutri1)
cox.zph(m9)
summary(m9)
BIC9<-round(BIC(m1)-BIC(m9),2)
C9<-round(summary(m9)$concordance[1],3)
B9 <- round(summary(m9)$coefficients[1],3)
p9 <- format.pval(summary(m9)$coefficients[5], eps = .001, digits = 3)
HR9 <- round(summary(m9)$conf.int[1],3)
HR9l <- round(summary(m9)$conf.int[3],3)
HR9u <- round(summary(m9)$conf.int[4],3)
HR9 <- paste0(HR9," ","(",HR9l,"-",HR9u,")")

m10<-coxph(Surv(FU_time, defuncion)~alb0, data=nutri1)
cox.zph(m10)
summary(m10)
BIC10<-round(BIC(m1)-BIC(m10),2)
C10<-round(summary(m10)$concordance[1],3)
B10 <- round(summary(m10)$coefficients[1],3)
p10 <- format.pval(summary(m10)$coefficients[5], eps = .001, digits = 3)
HR10 <- round(summary(m10)$conf.int[1],3)
HR10l <- round(summary(m10)$conf.int[3],3)
HR10u <- round(summary(m10)$conf.int[4],3)
HR10 <- paste0(HR10," ","(",HR10l,"-",HR10u,")")

mSa02<-coxph(Surv(FU_time, defuncion)~sao2aa, data=nutri1)
cox.zph(mSa02)
summary(mSa02)
BICSaO2<-round(BIC(m1)-BIC(mSa02),2)
CSaO2<-round(summary(mSa02)$concordance[1],3)
BSaO2 <- round(summary(mSa02)$coefficients[1],3)
pSaO2 <- format.pval(summary(mSa02)$coefficients[5], eps = .001, digits = 3)
HRSaO2 <- round(summary(mSa02)$conf.int[1],3)
HRSaO2l <- round(summary(mSa02)$conf.int[3],3)
HRSaO2u <- round(summary(mSa02)$conf.int[4],3)
HRSaO2 <- paste0(HRSaO2," ","(",HRSaO2l,"-",HRSaO2u,")")

m11<-coxph(Surv(FU_time, defuncion)~PhenoAge, data=nutri1)
cox.zph(m11)
summary(m11)
BIC11<-round(BIC(m1)-BIC(m11),2)
C11<-round(summary(m11)$concordance[1],3)
B11 <- round(summary(m11)$coefficients[1],3)
p11 <- format.pval(summary(m11)$coefficients[5], eps = .001, digits = 3)
HR11 <- round(summary(m11)$conf.int[1],3)
HR11l <- round(summary(m11)$conf.int[3],3)
HR11u <- round(summary(m11)$conf.int[4],3)
HR11 <- paste0(HR11," ","(",HR11l,"-",HR11u,")")

m12<-coxph(Surv(FU_time, defuncion)~PhenoAccelAge, data=nutri1)
cox.zph(m12)
summary(m12)
BIC12<-round(BIC(m1)-BIC(m12),2)
C12<-round(summary(m12)$concordance[1],3)
B12 <- round(summary(m12)$coefficients[1],3)
p12 <- format.pval(summary(m12)$coefficients[5], eps = .001, digits = 3)
HR12 <- round(summary(m12)$conf.int[1],3)
HR12l <- round(summary(m12)$conf.int[3],3)
HR12u <- round(summary(m12)$conf.int[4],3)
HR12 <- paste0(HR12," ","(",HR12l,"-",HR12u,")")

tab1<-data.frame("Parameter"=c("Age (years)","Glucose (mg/dL)","CRP","Alkaline Phosphatase","RDW","MCV","Lymphocytes (%)",
                         "Leucocytes (x1000)","Creatinine","Albumin","SaO2","PhenoAge (years)", "PhenoAgeAccel"),
           "DelBIC"=c(BIC1,BIC2,BIC3,BIC4,BIC5,BIC6,BIC7,BIC8,BIC9,BIC10,BICSaO2,BIC11,BIC12),
           "C-statistic"=c(C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,CSaO2,C11,C12),
           "Beta-coefficient"=c(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,BSaO2,B11,B12),
           "HR (95%CI)"=c(HR1,HR2,HR3,HR4,HR5,HR6,HR7,HR8,HR9,HR10,HRSaO2,HR11,HR12),
           "P-value"=c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,pSaO2,p11,p12))
tab1<-`names<-`(tab1,c("Parameter","Delta BIC","C-statistic","B-coefficient","HR (95%CI)","P-value"))
tab1<-align(flextable(tab1,cwidth = c(3,1,1,1,3,1)),align = "center",part = "all")
save_as_docx(tab1,path="tabla1.docx")#YA EST? GUARDADA


### Lethality models ###
m1<-coxph(Surv(FU_time, defuncion)~lin0+VCM+GLU+pcr0+ADE+edad, data=nutri1)
cox.zph(m2)
summary(m2)
  C1<-round(summary(m1)$concordance[1],3)
  B1 <- round(summary(m1)$coefficients[1],3)
  B2 <- round(summary(m1)$coefficients[2],3)
  B3 <- round(summary(m1)$coefficients[3],3)
  B4 <- round(summary(m1)$coefficients[4],3)
  B5 <- round(summary(m1)$coefficients[5],3)
  B6 <- round(summary(m1)$coefficients[6],3)
  HR1 <- round(summary(m1)$conf.int[1],3)
  HR1l <- round(summary(m1)$conf.int[13],3)
  HR1u <- round(summary(m1)$conf.int[19],3)
  HR1 <- paste0(HR1," ","(",HR1l,"-",HR1u,")")
  HR2 <- round(summary(m1)$conf.int[2],3)
  HR2l <- round(summary(m1)$conf.int[14],3)
  HR2u <- round(summary(m1)$conf.int[20],3)
  HR2 <- paste0(HR2," ","(",HR2l,"-",HR2u,")")
  HR3 <- round(summary(m1)$conf.int[3],3)
  HR3l <- round(summary(m1)$conf.int[15],3)
  HR3u <- round(summary(m1)$conf.int[21],3)
  HR3 <- paste0(HR3," ","(",HR3l,"-",HR3u,")")
  HR4 <- round(summary(m1)$conf.int[4],3)
  HR4l <- round(summary(m1)$conf.int[16],3)
  HR4u <- round(summary(m1)$conf.int[22],3)
  HR4 <- paste0(HR4," ","(",HR4l,"-",HR4u,")")
  HR5 <- round(summary(m1)$conf.int[5],3)
  HR5l <- round(summary(m1)$conf.int[17],3)
  HR5u <- round(summary(m1)$conf.int[23],3)
  HR5 <- paste0(HR5," ","(",HR5l,"-",HR5u,")")
  HR6 <- round(summary(m1)$conf.int[6],3)
  HR6l <- round(summary(m1)$conf.int[18],3)
  HR6u <- round(summary(m1)$conf.int[24],3)
  HR6 <- paste0(HR6," ","(",HR6l,"-",HR6u,")")
  p1 <- format.pval(summary(m1)$coefficients[25], eps = .001, digits = 3) 
  p2 <- format.pval(summary(m1)$coefficients[26], eps = .001, digits = 3)
  p3 <- format.pval(summary(m1)$coefficients[27], eps = .001, digits = 3)
  p4 <- format.pval(summary(m1)$coefficients[28], eps = .001, digits = 3)
  p5 <- format.pval(summary(m1)$coefficients[29], eps = .001, digits = 3)
  p6 <- format.pval(summary(m1)$coefficients[30], eps = .001, digits = 3)
  
m2<-coxph(Surv(FU_time, defuncion)~PhenoAccelAge+edad, data=nutri1)
cox.zph(m2)
summary(m2)
  C2<-round(summary(m2)$concordance[1],3)
  B7 <- round(summary(m2)$coefficients[1],3)
  B8 <- round(summary(m2)$coefficients[2],3)
  HR7 <- round(summary(m2)$conf.int[1],3)
  HR7l <- round(summary(m2)$conf.int[5],3)
  HR7u <- round(summary(m2)$conf.int[7],3)
  HR7 <- paste0(HR7," ","(",HR7l,"-",HR7u,")")
  HR8 <- round(summary(m2)$conf.int[2],3)
  HR8l <- round(summary(m2)$conf.int[6],3)
  HR8u <- round(summary(m2)$conf.int[8],3)
  HR8 <- paste0(HR8," ","(",HR8l,"-",HR8u,")")
  p7 <- format.pval(summary(m2)$coefficients[9], eps = .001, digits = 3)
  p8 <- format.pval(summary(m2)$coefficients[10], eps = .001, digits = 3)

tab2 <- data.frame("Model"=c(" "," ","PhenoAge Components",paste0("C-Statistic"," ",C1)," "," ",
                             "PhenoAgeAccel + Edad", paste0("C-Statistic"," ",C2)),
                   "Parameter"=c("Lymphocytes (%)","MCV","Glucose (mg/dL)","CRP","RDW","Age","PhenoAgeAccel","Age"),
                   "B-coefficient"=c(B1,B2,B3,B4,B5,B6,B7,B8),
                   "HR (95%CI)"=c(HR1,HR2,HR3,HR4,HR5,HR6,HR7,HR8),
                   "p-value"=c(p1,p2,p3,p4,p5,p6,p7,p8))

tab2 <-`names<-`(tab2,c("Model","Parameter","B-coefficient","HR (95%CI)","p-value"))
tab2 <-align(flextable(tab2,cwidth = c(2,1.5,1,2,1)),align = "center",part = "all")
save_as_docx(tab2,path="tabla2.docx")#YA EST? GUARDADA


#### Cox regression analyses - Lethality adjusted by sex and comorbidities ####
nutri1<-nutri[!is.na(nutri$PhenoAccelAge),]
nutri1<-nutri1%>%dplyr::select(critical, FU_time, defuncion, edad,
                              alb0, cr0, fa0, leu0, VCM, ADE, lin0, pcr0, GLU, PhenoAge, PhenoAccelAge, accel1, sao2aa, sexo, comorb,comorb1,comorb2,groups1)
set.seed(123); imp<-mice::mice(nutri1, m=1, maxit=1); nutri1<-complete(imp, "long")
nrow(nutri1)
m1<-coxph(Surv(FU_time, defuncion)~edad+comorb1+sexo, data=nutri1)
cox.zph(m1)
summary(m1)
BIC1<-round(BIC(m1)-BIC(m1),2)
C1<-round(summary(m1)$concordance[1],3)
B1 <- round(summary(m1)$coefficients[1],3)
p1 <- format.pval(summary(m1)$coefficients[13], eps = .001, digits = 3)
HR1 <- round(summary(m1)$conf.int[1],3)
HR1l <- round(summary(m1)$conf.int[7],3)
HR1u <- round(summary(m1)$conf.int[10],3)
HR1 <- paste0(HR1," ","(",HR1l,"-",HR1u,")")

m2<-coxph(Surv(FU_time, defuncion)~GLU+comorb1+sexo, data=nutri1)
cox.zph(m2)
summary(m2)
BIC2<-round(BIC(m1)-BIC(m2),2)
C2<-round(summary(m2)$concordance[1],3)
B2 <- round(summary(m2)$coefficients[1],3)
p2 <- format.pval(summary(m2)$coefficients[13], eps = .001, digits = 3)
HR2 <- round(summary(m2)$conf.int[1],3)
HR2l <- round(summary(m2)$conf.int[7],3)
HR2u <- round(summary(m2)$conf.int[10],3)
HR2 <- paste0(HR2," ","(",HR2l,"-",HR2u,")")

m3<-coxph(Surv(FU_time, defuncion)~pcr0+comorb1+sexo, data=nutri1)
cox.zph(m3)
summary(m3)
BIC3<-round(BIC(m1)-BIC(m3),2)
C3<-round(summary(m3)$concordance[1],3)
B3 <- round(summary(m3)$coefficients[1],3)
p3 <- format.pval(summary(m3)$coefficients[13], eps = .001, digits = 3)
HR3 <- round(summary(m3)$conf.int[1],3)
HR3l <- round(summary(m3)$conf.int[7],3)
HR3u <- round(summary(m3)$conf.int[10],3)
HR3 <- paste0(HR3," ","(",HR3l,"-",HR3u,")")

m4<-coxph(Surv(FU_time, defuncion)~fa0+comorb1+sexo, data=nutri1)
cox.zph(m4)
summary(m4)
BIC4<-round(BIC(m1)-BIC(m4),2)
C4<-round(summary(m4)$concordance[1],3)
B4 <- round(summary(m4)$coefficients[1],3)
p4 <- format.pval(summary(m4)$coefficients[13], eps = .001, digits = 3)
HR4 <- round(summary(m4)$conf.int[1],3)
HR4l <- round(summary(m4)$conf.int[7],3)
HR4u <- round(summary(m4)$conf.int[10],3)
HR4 <- paste0(HR4," ","(",HR4l,"-",HR4u,")")

m5<-coxph(Surv(FU_time, defuncion)~ADE+comorb1+sexo, data=nutri1)
cox.zph(m5)
summary(m5)
BIC5<-round(BIC(m1)-BIC(m5),2)
C5<-round(summary(m5)$concordance[1],3)
B5 <- round(summary(m5)$coefficients[1],3)
p5 <- format.pval(summary(m5)$coefficients[13], eps = .001, digits = 3)
HR5 <- round(summary(m5)$conf.int[1],3)
HR5l <- round(summary(m5)$conf.int[7],3)
HR5u <- round(summary(m5)$conf.int[10],3)
HR5 <- paste0(HR5," ","(",HR5l,"-",HR5u,")")

m6<-coxph(Surv(FU_time, defuncion)~VCM+comorb1+sexo, data=nutri1)
cox.zph(m6)
summary(m6)
BIC6<-round(BIC(m1)-BIC(m6),2)
C6<-round(summary(m6)$concordance[1],3)
B6 <- round(summary(m6)$coefficients[1],3)
p6 <- format.pval(summary(m6)$coefficients[13], eps = .001, digits = 3)
HR6 <- round(summary(m6)$conf.int[1],3)
HR6l <- round(summary(m6)$conf.int[7],3)
HR6u <- round(summary(m6)$conf.int[10],3)
HR6 <- paste0(HR6," ","(",HR6l,"-",HR6u,")")

m7<-coxph(Surv(FU_time, defuncion)~lin0+comorb1+sexo, data=nutri1)
cox.zph(m7)
summary(m7)
BIC7<-round(BIC(m1)-BIC(m7),2)
C7<-round(summary(m7)$concordance[1],3)
B7 <- round(summary(m7)$coefficients[1],3)
p7 <- format.pval(summary(m7)$coefficients[13], eps = .001, digits = 3)
HR7 <- round(summary(m7)$conf.int[1],3)
HR7l <- round(summary(m7)$conf.int[7],3)
HR7u <- round(summary(m7)$conf.int[10],3)
HR7 <- paste0(HR7," ","(",HR7l,"-",HR7u,")")

nutri1$leu_2<-nutri1$leu0/1000
m8<-coxph(Surv(FU_time, defuncion)~leu_2+comorb1+sexo, data=nutri1)
cox.zph(m8)
summary(m8)
BIC8<-round(BIC(m1)-BIC(m8),2)
C8<-round(summary(m8)$concordance[1],3)
B8 <- round(summary(m8)$coefficients[1],3)
p8 <- format.pval(summary(m8)$coefficients[13], eps = .001, digits = 3)
HR8 <- round(summary(m8)$conf.int[1],3)
HR8l <- round(summary(m8)$conf.int[7],3)
HR8u <- round(summary(m8)$conf.int[10],3)
HR8 <- paste0(HR8," ","(",HR8l,"-",HR8u,")")

m9<-coxph(Surv(FU_time, defuncion)~cr0+comorb1+sexo, data=nutri1)
cox.zph(m9)
summary(m9)
BIC9<-round(BIC(m1)-BIC(m9),2)
C9<-round(summary(m9)$concordance[1],3)
B9 <- round(summary(m9)$coefficients[1],3)
p9 <- format.pval(summary(m9)$coefficients[13], eps = .001, digits = 3)
HR9 <- round(summary(m9)$conf.int[1],3)
HR9l <- round(summary(m9)$conf.int[7],3)
HR9u <- round(summary(m9)$conf.int[10],3)
HR9 <- paste0(HR9," ","(",HR9l,"-",HR9u,")")

m10<-coxph(Surv(FU_time, defuncion)~alb0+comorb1+sexo, data=nutri1)
cox.zph(m10)
summary(m10)
BIC10<-round(BIC(m1)-BIC(m10),2)
C10<-round(summary(m10)$concordance[1],3)
B10 <- round(summary(m10)$coefficients[1],3)
p10 <- format.pval(summary(m10)$coefficients[13], eps = .001, digits = 3)
HR10 <- round(summary(m10)$conf.int[1],3)
HR10l <- round(summary(m10)$conf.int[7],3)
HR10u <- round(summary(m10)$conf.int[10],3)
HR10 <- paste0(HR10," ","(",HR10l,"-",HR10u,")")

mSa02<-coxph(Surv(FU_time, defuncion)~sao2aa+comorb1+sexo, data=nutri1)
cox.zph(mSa02)
summary(mSa02)
BICSaO2<-round(BIC(m1)-BIC(mSa02),2)
CSaO2<-round(summary(mSa02)$concordance[1],3)
BSaO2 <- round(summary(mSa02)$coefficients[1],3)
pSaO2 <- format.pval(summary(mSa02)$coefficients[13], eps = .001, digits = 3)
HRSaO2 <- round(summary(mSa02)$conf.int[1],3)
HRSaO2l <- round(summary(mSa02)$conf.int[7],3)
HRSaO2u <- round(summary(mSa02)$conf.int[10],3)
HRSaO2 <- paste0(HRSaO2," ","(",HRSaO2l,"-",HRSaO2u,")")

m11<-coxph(Surv(FU_time, defuncion)~PhenoAge+comorb1+sexo, data=nutri1)
cox.zph(m11)
summary(m11)
BIC11<-round(BIC(m1)-BIC(m11),2)
C11<-round(summary(m11)$concordance[1],3)
B11 <- round(summary(m11)$coefficients[1],3)
p11 <- format.pval(summary(m11)$coefficients[13], eps = .001, digits = 3)
HR11 <- round(summary(m11)$conf.int[1],3)
HR11l <- round(summary(m11)$conf.int[7],3)
HR11u <- round(summary(m11)$conf.int[10],3)
HR11 <- paste0(HR11," ","(",HR11l,"-",HR11u,")")

m12<-coxph(Surv(FU_time, defuncion)~PhenoAccelAge+comorb1+sexo, data=nutri1)
cox.zph(m12)
summary(m12)
BIC12<-round(BIC(m1)-BIC(m12),2)
C12<-round(summary(m12)$concordance[1],3)
B12 <- round(summary(m12)$coefficients[1],3)
p12 <- format.pval(summary(m12)$coefficients[13], eps = .001, digits = 3)
HR12 <- round(summary(m12)$conf.int[1],3)
HR12l <- round(summary(m12)$conf.int[7],3)
HR12u <- round(summary(m12)$conf.int[10],3)
HR12 <- paste0(HR12," ","(",HR12l,"-",HR12u,")")

tab1_comorb<-data.frame("Parameter"=c("Age (years)","Glucose (mg/dL)","CRP","Alkaline Phosphatase","RDW","MCV","Lymphocytes (%)",
                               "Leucocytes (x1000)","Creatinine","Albumin","SaO2","PhenoAge (years)", "PhenoAgeAccel"),
                 "DelBIC"=c(BIC1,BIC2,BIC3,BIC4,BIC5,BIC6,BIC7,BIC8,BIC9,BIC10,BICSaO2,BIC11,BIC12),
                 "C-statistic"=c(C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,CSaO2,C11,C12),
                 "Beta-coefficient"=c(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,BSaO2,B11,B12),
                 "HR (95%CI)"=c(HR1,HR2,HR3,HR4,HR5,HR6,HR7,HR8,HR9,HR10,HRSaO2,HR11,HR12),
                 "P-value"=c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,pSaO2,p11,p12))
tab1_comorb<-`names<-`(tab1_comorb,c("Parameter","Delta BIC","C-statistic","B-coefficient","HR (95%CI)","P-value"))
tab1_comorb<-align(flextable(tab1_comorb,cwidth = c(3,1,1,1,3,1)),align = "center",part = "all")
save_as_docx(tab1_comorb,path="tabla1_comorb.docx")#YA EST? GUARDADA


### Lethality models adjusted by sex and comorbidities ###
m1<-coxph(Surv(FU_time, defuncion)~lin0+VCM+GLU+pcr0+ADE+edad+sexo+comorb1, data=nutri1)
cox.zph(m1)
summary(m1)
C1<-round(summary(m1)$concordance[1],3)
B1 <- round(summary(m1)$coefficients[1],3)
B2 <- round(summary(m1)$coefficients[2],3)
B3 <- round(summary(m1)$coefficients[3],3)
B4 <- round(summary(m1)$coefficients[4],3)
B5 <- round(summary(m1)$coefficients[5],3)
B6 <- round(summary(m1)$coefficients[6],3)
B7 <- round(summary(m1)$coefficients[7],3)
B8 <- round(summary(m1)$coefficients[8],3)
HR1 <- round(summary(m1)$conf.int[1],3)
HR1l <- round(summary(m1)$conf.int[17],3)
HR1u <- round(summary(m1)$conf.int[25],3)
HR1 <- paste0(HR1," ","(",HR1l,"-",HR1u,")")
HR2 <- round(summary(m1)$conf.int[2],3)
HR2l <- round(summary(m1)$conf.int[18],3)
HR2u <- round(summary(m1)$conf.int[26],3)
HR2 <- paste0(HR2," ","(",HR2l,"-",HR2u,")")
HR3 <- round(summary(m1)$conf.int[3],3)
HR3l <- round(summary(m1)$conf.int[19],3)
HR3u <- round(summary(m1)$conf.int[27],3)
HR3 <- paste0(HR3," ","(",HR3l,"-",HR3u,")")
HR4 <- round(summary(m1)$conf.int[4],3)
HR4l <- round(summary(m1)$conf.int[20],3)
HR4u <- round(summary(m1)$conf.int[28],3)
HR4 <- paste0(HR4," ","(",HR4l,"-",HR4u,")")
HR5 <- round(summary(m1)$conf.int[5],3)
HR5l <- round(summary(m1)$conf.int[21],3)
HR5u <- round(summary(m1)$conf.int[29],3)
HR5 <- paste0(HR5," ","(",HR5l,"-",HR5u,")")
HR6 <- round(summary(m1)$conf.int[6],3)
HR6l <- round(summary(m1)$conf.int[22],3)
HR6u <- round(summary(m1)$conf.int[30],3)
HR6 <- paste0(HR6," ","(",HR6l,"-",HR6u,")")
HR7 <- round(summary(m1)$conf.int[7],3)
HR7l <- round(summary(m1)$conf.int[23],3)
HR7u <- round(summary(m1)$conf.int[31],3)
HR7 <- paste0(HR7," ","(",HR7l,"-",HR7u,")")
HR8 <- round(summary(m1)$conf.int[8],3)
HR8l <- round(summary(m1)$conf.int[24],3)
HR8u <- round(summary(m1)$conf.int[32],3)
HR8 <- paste0(HR8," ","(",HR8l,"-",HR8u,")")
p1 <- format.pval(summary(m1)$coefficients[33], eps = .001, digits = 3) 
p2 <- format.pval(summary(m1)$coefficients[34], eps = .001, digits = 3)
p3 <- format.pval(summary(m1)$coefficients[35], eps = .001, digits = 3)
p4 <- format.pval(summary(m1)$coefficients[36], eps = .001, digits = 3)
p5 <- format.pval(summary(m1)$coefficients[37], eps = .001, digits = 3)
p6 <- format.pval(summary(m1)$coefficients[38], eps = .001, digits = 3)
p7 <- format.pval(summary(m1)$coefficients[39], eps = .001, digits = 3)
p8 <- format.pval(summary(m1)$coefficients[40], eps = .001, digits = 3)

m2<-coxph(Surv(FU_time, defuncion)~PhenoAccelAge+edad+sexo+comorb1, data=nutri1)
BIC(m2)
cox.zph(m2)
summary(m2)
C2<-round(summary(m2)$concordance[1],3)
B9 <- round(summary(m2)$coefficients[1],3)
B10 <- round(summary(m2)$coefficients[2],3)
B11 <- round(summary(m2)$coefficients[3],3)
B12 <- round(summary(m2)$coefficients[4],3)
HR9 <- round(summary(m2)$conf.int[1],3)
HR9l <- round(summary(m2)$conf.int[9],3)
HR9u <- round(summary(m2)$conf.int[13],3)
HR9 <- paste0(HR9," ","(",HR9l,"-",HR9u,")")
HR10 <- round(summary(m2)$conf.int[2],3)
HR10l <- round(summary(m2)$conf.int[10],3)
HR10u <- round(summary(m2)$conf.int[14],3)
HR10 <- paste0(HR10," ","(",HR10l,"-",HR10u,")")
HR11 <- round(summary(m2)$conf.int[3],3)
HR11l <- round(summary(m2)$conf.int[11],3)
HR11u <- round(summary(m2)$conf.int[15],3)
HR11 <- paste0(HR11," ","(",HR11l,"-",HR11u,")")
HR12 <- round(summary(m2)$conf.int[4],3)
HR12l <- round(summary(m2)$conf.int[12],3)
HR12u <- round(summary(m2)$conf.int[16],3)
HR12 <- paste0(HR12," ","(",HR12l,"-",HR12u,")")
p9 <- format.pval(summary(m2)$coefficients[17], eps = .001, digits = 3)
p10 <- format.pval(summary(m2)$coefficients[18], eps = .001, digits = 3)
p11 <- format.pval(summary(m2)$coefficients[19], eps = .001, digits = 3)
p12 <- format.pval(summary(m2)$coefficients[20], eps = .001, digits = 3)

tab2_comorb <- data.frame("Model"=c(" "," "," ", "PhenoAge","Components",paste0("C-Statistic"," ",C1)," "," ",
                             " ","PhenoAgeAccel + Age", paste0("C-Statistic"," ",C2)," "), 
                   "Parameter"=c("Lymphocytes (%)","MCV","Glucose (mg/dL)","C-reactive protein","RDW","Age","Male","Comorbidities","PhenoAgeAccel","Age","Male","Comorbidities"),
                   "B-coefficient"=c(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12),
                   "HR (95%CI)"=c(HR1,HR2,HR3,HR4,HR5,HR6,HR7,HR8,HR9,HR10,HR11,HR12),
                   "p-value"=c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12))

tab2_comorb <-`names<-`(tab2_comorb,c("Model","Parameter","B-coefficient","HR (95%CI)","p-value"))
tab2_comorb <-align(flextable(tab2_comorb,cwidth = c(2,1.5,1,2,1)),align = "center",part = "all")
save_as_docx(tab2_comorb,path="tabla2_comorb.docx")#YA EST? GUARDADA


#### Cox regression analyses - Adverse Outcomes ####
nutri1<-nutri[!is.na(nutri$PhenoAccelAge),]
nutri1<-nutri1%>%dplyr::select(critical, FU_time, defuncion, edad,
                              alb0, cr0, fa0, leu0, VCM, ADE, lin0, pcr0, GLU, PhenoAge, PhenoAccelAge, accel1, sao2aa, sexo, comorb,comorb2,groups1)
set.seed(123); imp<-mice::mice(nutri1, m=1, maxit=1); nutri1<-complete(imp, "long")
nrow(nutri1)

m1<-coxph(Surv(FU_time, critical)~edad, data=nutri1)
cox.zph(m1)
summary(m1)
BIC1<-round(BIC(m1)-BIC(m1),2)
C1<-round(summary(m1)$concordance[1],3)
B1 <- round(summary(m1)$coefficients[1],3)
p1 <- format.pval(summary(m1)$coefficients[5], eps = .001, digits = 3)
HR1 <- round(summary(m1)$conf.int[1],3)
HR1l <- round(summary(m1)$conf.int[3],3)
HR1u <- round(summary(m1)$conf.int[4],3)
HR1 <- paste0(HR1," ","(",HR1l,"-",HR1u,")")

m2<-coxph(Surv(FU_time, critical)~GLU, data=nutri1)
cox.zph(m2)
summary(m2)
BIC2<-round(BIC(m1)-BIC(m2),2)
C2<-round(summary(m2)$concordance[1],3)
B2 <- round(summary(m2)$coefficients[1],3)
p2 <- format.pval(summary(m2)$coefficients[5], eps = .001, digits = 3)
HR2 <- round(summary(m2)$conf.int[1],3)
HR2l <- round(summary(m2)$conf.int[3],3)
HR2u <- round(summary(m2)$conf.int[4],3)
HR2 <- paste0(HR2," ","(",HR2l,"-",HR2u,")")

m3<-coxph(Surv(FU_time, critical)~pcr0, data=nutri1)
cox.zph(m3)
summary(m3)
BIC3<-round(BIC(m1)-BIC(m3),2)
C3<-round(summary(m3)$concordance[1],3)
B3 <- round(summary(m3)$coefficients[1],3)
p3 <- format.pval(summary(m3)$coefficients[5], eps = .001, digits = 3)
HR3 <- round(summary(m3)$conf.int[1],3)
HR3l <- round(summary(m3)$conf.int[3],3)
HR3u <- round(summary(m3)$conf.int[4],3)
HR3 <- paste0(HR3," ","(",HR3l,"-",HR3u,")")

m4<-coxph(Surv(FU_time, critical)~fa0, data=nutri1)
cox.zph(m4)
summary(m4)
BIC4<-round(BIC(m1)-BIC(m4),2)
C4<-round(summary(m4)$concordance[1],3)
B4 <- round(summary(m4)$coefficients[1],3)
p4 <- format.pval(summary(m4)$coefficients[5], eps = .001, digits = 3)
HR4 <- round(summary(m4)$conf.int[1],3)
HR4l <- round(summary(m4)$conf.int[3],3)
HR4u <- round(summary(m4)$conf.int[4],3)
HR4 <- paste0(HR4," ","(",HR4l,"-",HR4u,")")

m5<-coxph(Surv(FU_time, critical)~ADE, data=nutri1)
cox.zph(m5)
summary(m5)
BIC5<-round(BIC(m1)-BIC(m5),2)
C5<-round(summary(m5)$concordance[1],3)
B5 <- round(summary(m5)$coefficients[1],3)
p5 <- format.pval(summary(m5)$coefficients[5], eps = .001, digits = 3)
HR5 <- round(summary(m5)$conf.int[1],3)
HR5l <- round(summary(m5)$conf.int[3],3)
HR5u <- round(summary(m5)$conf.int[4],3)
HR5 <- paste0(HR5," ","(",HR5l,"-",HR5u,")")

m6<-coxph(Surv(FU_time, critical)~VCM, data=nutri1)
cox.zph(m6)
summary(m6)
BIC6<-round(BIC(m1)-BIC(m6),2)
C6<-round(summary(m6)$concordance[1],3)
B6 <- round(summary(m6)$coefficients[1],3)
p6 <- format.pval(summary(m6)$coefficients[5], eps = .001, digits = 3)
HR6 <- round(summary(m6)$conf.int[1],3)
HR6l <- round(summary(m6)$conf.int[3],3)
HR6u <- round(summary(m6)$conf.int[4],3)
HR6 <- paste0(HR6," ","(",HR6l,"-",HR6u,")")

m7<-coxph(Surv(FU_time, critical)~lin0, data=nutri1)
cox.zph(m7)
summary(m7)
BIC7<-round(BIC(m1)-BIC(m7),2)
C7<-round(summary(m7)$concordance[1],3)
B7 <- round(summary(m7)$coefficients[1],3)
p7 <- format.pval(summary(m7)$coefficients[5], eps = .001, digits = 3)
HR7 <- round(summary(m7)$conf.int[1],3)
HR7l <- round(summary(m7)$conf.int[3],3)
HR7u <- round(summary(m7)$conf.int[4],3)
HR7 <- paste0(HR7," ","(",HR7l,"-",HR7u,")")

nutri1$leu_2<-nutri1$leu0/1000
m8<-coxph(Surv(FU_time, critical)~leu_2, data=nutri1)
cox.zph(m8)
summary(m8)
BIC8<-round(BIC(m1)-BIC(m8),2)
C8<-round(summary(m8)$concordance[1],3)
B8 <- round(summary(m8)$coefficients[1],3)
p8 <- format.pval(summary(m8)$coefficients[5], eps = .001, digits = 3)
HR8 <- round(summary(m8)$conf.int[1],3)
HR8l <- round(summary(m8)$conf.int[3],3)
HR8u <- round(summary(m8)$conf.int[4],3)
HR8 <- paste0(HR8," ","(",HR8l,"-",HR8u,")")

m9<-coxph(Surv(FU_time, critical)~cr0, data=nutri1)
cox.zph(m9)
summary(m9)
BIC9<-round(BIC(m1)-BIC(m9),2)
C9<-round(summary(m9)$concordance[1],3)
B9 <- round(summary(m9)$coefficients[1],3)
p9 <- format.pval(summary(m9)$coefficients[5], eps = .001, digits = 3)
HR9 <- round(summary(m9)$conf.int[1],3)
HR9l <- round(summary(m9)$conf.int[3],3)
HR9u <- round(summary(m9)$conf.int[4],3)
HR9 <- paste0(HR9," ","(",HR9l,"-",HR9u,")")

m10<-coxph(Surv(FU_time, critical)~alb0, data=nutri1)
cox.zph(m10)
summary(m10)
BIC10<-round(BIC(m1)-BIC(m10),2)
C10<-round(summary(m10)$concordance[1],3)
B10 <- round(summary(m10)$coefficients[1],3)
p10 <- format.pval(summary(m10)$coefficients[5], eps = .001, digits = 3)
HR10 <- round(summary(m10)$conf.int[1],3)
HR10l <- round(summary(m10)$conf.int[3],3)
HR10u <- round(summary(m10)$conf.int[4],3)
HR10 <- paste0(HR10," ","(",HR10l,"-",HR10u,")")

mSa02<-coxph(Surv(FU_time, critical)~sao2aa, data=nutri1)
cox.zph(mSa02)
summary(mSa02)
BICSaO2<-round(BIC(m1)-BIC(mSa02),2)
CSaO2<-round(summary(mSa02)$concordance[1],3)
BSaO2 <- round(summary(mSa02)$coefficients[1],3)
pSaO2 <- format.pval(summary(mSa02)$coefficients[5], eps = .001, digits = 3)
HRSaO2 <- round(summary(mSa02)$conf.int[1],3)
HRSaO2l <- round(summary(mSa02)$conf.int[3],3)
HRSaO2u <- round(summary(mSa02)$conf.int[4],3)
HRSaO2 <- paste0(HRSaO2," ","(",HRSaO2l,"-",HRSaO2u,")")

m11<-coxph(Surv(FU_time, critical)~PhenoAge, data=nutri1)
cox.zph(m11)
summary(m11)
BIC11<-round(BIC(m1)-BIC(m11),2)
C11<-round(summary(m11)$concordance[1],3)
B11 <- round(summary(m11)$coefficients[1],3)
p11 <- format.pval(summary(m11)$coefficients[5], eps = .001, digits = 3)
HR11 <- round(summary(m11)$conf.int[1],3)
HR11l <- round(summary(m11)$conf.int[3],3)
HR11u <- round(summary(m11)$conf.int[4],3)
HR11 <- paste0(HR11," ","(",HR11l,"-",HR11u,")")

m12<-coxph(Surv(FU_time, critical)~PhenoAccelAge, data=nutri1)
cox.zph(m12)
summary(m12)
BIC12<-round(BIC(m1)-BIC(m12),2)
C12<-round(summary(m12)$concordance[1],3)
B12 <- round(summary(m12)$coefficients[1],3)
p12 <- format.pval(summary(m12)$coefficients[5], eps = .001, digits = 3)
HR12 <- round(summary(m12)$conf.int[1],3)
HR12l <- round(summary(m12)$conf.int[3],3)
HR12u <- round(summary(m12)$conf.int[4],3)
HR12 <- paste0(HR12," ","(",HR12l,"-",HR12u,")")

tab3<-data.frame("Parameter"=c("Age (years)","Glucose (mg/dL)","CRP","Alkaline Phosphatase","RDW","MCV","Lymphocytes (%)",
                               "Leucocytes (x1000)","Creatinine","Albumin","SaO2","PhenoAge (years)", "PhenoAgeAccel"),
                 "DelBIC"=c(BIC1,BIC2,BIC3,BIC4,BIC5,BIC6,BIC7,BIC8,BIC9,BIC10,BICSaO2,BIC11,BIC12),
                 "C-statistic"=c(C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,CSaO2,C11,C12),
                 "Beta-coefficient"=c(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,BSaO2,B11,B12),
                 "HR (95%CI)"=c(HR1,HR2,HR3,HR4,HR5,HR6,HR7,HR8,HR9,HR10,HRSaO2,HR11,HR12),
                 "P-value"=c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,pSaO2,p11,p12))
tab3<-`names<-`(tab3,c("Parameter","Delta BIC","C-statistic","B-coefficient","HR (95%CI)","P-value"))

tab3<-align(flextable(tab3,cwidth = c(3,1,1,1,3,1)),align = "center",part = "all")
save_as_docx(tab3,path="tabla3.docx")#YA EST? GUARDADA


### Adverse outcomes models ###
m1<-coxph(Surv(FU_time, critical)~lin0+GLU+pcr0+edad, data=nutri1)
cox.zph(m1)
summary(m1)
C1<-round(summary(m1)$concordance[1],3)
B1 <- round(summary(m1)$coefficients[1],3)
B2 <- round(summary(m1)$coefficients[2],3)
B3 <- round(summary(m1)$coefficients[3],3)
B4 <- round(summary(m1)$coefficients[4],3)
HR1 <- round(summary(m1)$conf.int[1],3)
HR1l <- round(summary(m1)$conf.int[9],3)
HR1u <- round(summary(m1)$conf.int[13],3)
HR1 <- paste0(HR1," ","(",HR1l,"-",HR1u,")")
HR2 <- round(summary(m1)$conf.int[2],3)
HR2l <- round(summary(m1)$conf.int[10],3)
HR2u <- round(summary(m1)$conf.int[14],3)
HR2 <- paste0(HR2," ","(",HR2l,"-",HR2u,")")
HR3 <- round(summary(m1)$conf.int[3],3)
HR3l <- round(summary(m1)$conf.int[11],3)
HR3u <- round(summary(m1)$conf.int[15],3)
HR3 <- paste0(HR3," ","(",HR3l,"-",HR3u,")")
HR4 <- round(summary(m1)$conf.int[4],3)
HR4l <- round(summary(m1)$conf.int[12],3)
HR4u <- round(summary(m1)$conf.int[16],3)
HR4 <- paste0(HR4," ","(",HR4l,"-",HR4u,")")
p1 <- format.pval(summary(m1)$coefficients[17], eps = .001, digits = 3)
p2 <- format.pval(summary(m1)$coefficients[18], eps = .001, digits = 3)
p3 <- format.pval(summary(m1)$coefficients[19], eps = .001, digits = 3)
p4 <- format.pval(summary(m1)$coefficients[20], eps = .001, digits = 3)

m2<-coxph(Surv(FU_time, critical)~PhenoAccelAge+edad, data=nutri1)
cox.zph(m2)
summary(m2)
  C2<-round(summary(m2)$concordance[1],3)
  B5 <- round(summary(m2)$coefficients[1],3)
  B6 <- round(summary(m2)$coefficients[2],3)
  HR5 <- round(summary(m2)$conf.int[1],3)
  HR5l <- round(summary(m2)$conf.int[5],3)
  HR5u <- round(summary(m2)$conf.int[7],3)
  HR5 <- paste0(HR5," ","(",HR5l,"-",HR5u,")")
  HR6 <- round(summary(m2)$conf.int[2],3)
  HR6l <- round(summary(m2)$conf.int[6],3)
  HR6u <- round(summary(m2)$conf.int[8],3)
  HR6 <- paste0(HR6," ","(",HR6l,"-",HR6u,")")
  p5 <- format.pval(summary(m2)$coefficients[9], eps = .001, digits = 3)
  p6 <- format.pval(summary(m2)$coefficients[10], eps = .001, digits = 3)
  
tab4 <- data.frame("Model"=c(" ","PhenoAge Components",paste0("C-Statistic"," ",C1)," ",
                               "PhenoAgeAccel + Edad", paste0("C-Statistic"," ",C2)),
                   "Parameter"=c("Lymphocytes (%)","Glucose (mg/dL)","CRP","Age (years)","PhenoAgeAccel","Age (years)"),
                   "B-coefficient"=c(B1,B2,B3,B4,B5,B6),
                   "HR (95%CI)"=c(HR1,HR2,HR3,HR4,HR5,HR6),
                   "p-value"=c(p1,p2,p3,p4,p5,p6))
  
tab4 <-`names<-`(tab4,c("Model","Parameter","B-coefficient","HR (95%CI)","p-value"))
tab4 <-align(flextable(tab4,cwidth = c(2,1.5,1,2,1)),align = "center",part = "all")
save_as_docx(tab4,path="tabla4.docx")#YA EST? GUARDADA

#### Cox regression analyses - Adverse outcomes adjusted by sex and comorbidities ####
nutri1<-nutri[!is.na(nutri$PhenoAccelAge),]
nutri1<-nutri1%>%dplyr::select(critical, FU_time, defuncion, edad,
                              alb0, cr0, fa0, leu0, VCM, ADE, lin0, pcr0, GLU, PhenoAge, PhenoAccelAge, accel1, sao2aa, sexo, comorb,comorb1,comorb2,groups1, edad_cat)
set.seed(123); imp<-mice::mice(nutri1, m=1, maxit=1); nutri1<-complete(imp, "long")
nrow(nutri1)
m1<-coxph(Surv(FU_time, critical)~edad+comorb1+sexo, data=nutri1)
cox.zph(m1)
summary(m1)
BIC1<-round(BIC(m1)-BIC(m1),2)
C1<-round(summary(m1)$concordance[1],3)
B1 <- round(summary(m1)$coefficients[1],3)
p1 <- format.pval(summary(m1)$coefficients[13], eps = .001, digits = 3)
HR1 <- round(summary(m1)$conf.int[1],3)
HR1l <- round(summary(m1)$conf.int[7],3)
HR1u <- round(summary(m1)$conf.int[10],3)
HR1 <- paste0(HR1," ","(",HR1l,"-",HR1u,")")

m2<-coxph(Surv(FU_time, critical)~GLU+comorb1+sexo, data=nutri1)
cox.zph(m2)
summary(m2)
BIC2<-round(BIC(m1)-BIC(m2),2)
C2<-round(summary(m2)$concordance[1],3)
B2 <- round(summary(m2)$coefficients[1],3)
p2 <- format.pval(summary(m2)$coefficients[13], eps = .001, digits = 3)
HR2 <- round(summary(m2)$conf.int[1],3)
HR2l <- round(summary(m2)$conf.int[7],3)
HR2u <- round(summary(m2)$conf.int[10],3)
HR2 <- paste0(HR2," ","(",HR2l,"-",HR2u,")")

m3<-coxph(Surv(FU_time, critical)~pcr0+comorb1+sexo, data=nutri1)
cox.zph(m3)
summary(m3)
BIC3<-round(BIC(m1)-BIC(m3),2)
C3<-round(summary(m3)$concordance[1],3)
B3 <- round(summary(m3)$coefficients[1],3)
p3 <- format.pval(summary(m3)$coefficients[13], eps = .001, digits = 3)
HR3 <- round(summary(m3)$conf.int[1],3)
HR3l <- round(summary(m3)$conf.int[7],3)
HR3u <- round(summary(m3)$conf.int[10],3)
HR3 <- paste0(HR3," ","(",HR3l,"-",HR3u,")")

m4<-coxph(Surv(FU_time, critical)~fa0+comorb1+sexo, data=nutri1)
cox.zph(m4)
summary(m4)
BIC4<-round(BIC(m1)-BIC(m4),2)
C4<-round(summary(m4)$concordance[1],3)
B4 <- round(summary(m4)$coefficients[1],3)
p4 <- format.pval(summary(m4)$coefficients[13], eps = .001, digits = 3)
HR4 <- round(summary(m4)$conf.int[1],3)
HR4l <- round(summary(m4)$conf.int[7],3)
HR4u <- round(summary(m4)$conf.int[10],3)
HR4 <- paste0(HR4," ","(",HR4l,"-",HR4u,")")

m5<-coxph(Surv(FU_time, critical)~ADE+comorb1+sexo, data=nutri1)
cox.zph(m5)
summary(m5)
BIC5<-round(BIC(m1)-BIC(m5),2)
C5<-round(summary(m5)$concordance[1],3)
B5 <- round(summary(m5)$coefficients[1],3)
p5 <- format.pval(summary(m5)$coefficients[13], eps = .001, digits = 3)
HR5 <- round(summary(m5)$conf.int[1],3)
HR5l <- round(summary(m5)$conf.int[7],3)
HR5u <- round(summary(m5)$conf.int[10],3)
HR5 <- paste0(HR5," ","(",HR5l,"-",HR5u,")")

m6<-coxph(Surv(FU_time, critical)~VCM+comorb1+sexo, data=nutri1)
cox.zph(m6)
summary(m6)
BIC6<-round(BIC(m1)-BIC(m6),2)
C6<-round(summary(m6)$concordance[1],3)
B6 <- round(summary(m6)$coefficients[1],3)
p6 <- format.pval(summary(m6)$coefficients[13], eps = .001, digits = 3)
HR6 <- round(summary(m6)$conf.int[1],3)
HR6l <- round(summary(m6)$conf.int[7],3)
HR6u <- round(summary(m6)$conf.int[10],3)
HR6 <- paste0(HR6," ","(",HR6l,"-",HR6u,")")

m7<-coxph(Surv(FU_time, critical)~lin0+comorb1+sexo, data=nutri1)
cox.zph(m7)
summary(m7)
BIC7<-round(BIC(m1)-BIC(m7),2)
C7<-round(summary(m7)$concordance[1],3)
B7 <- round(summary(m7)$coefficients[1],3)
p7 <- format.pval(summary(m7)$coefficients[13], eps = .001, digits = 3)
HR7 <- round(summary(m7)$conf.int[1],3)
HR7l <- round(summary(m7)$conf.int[7],3)
HR7u <- round(summary(m7)$conf.int[10],3)
HR7 <- paste0(HR7," ","(",HR7l,"-",HR7u,")")

nutri1$leu_2<-nutri1$leu0/1000
m8<-coxph(Surv(FU_time, critical)~leu_2+comorb1+sexo, data=nutri1)
cox.zph(m8)
summary(m8)
BIC8<-round(BIC(m1)-BIC(m8),2)
C8<-round(summary(m8)$concordance[1],3)
B8 <- round(summary(m8)$coefficients[1],3)
p8 <- format.pval(summary(m8)$coefficients[13], eps = .001, digits = 3)
HR8 <- round(summary(m8)$conf.int[1],3)
HR8l <- round(summary(m8)$conf.int[7],3)
HR8u <- round(summary(m8)$conf.int[10],3)
HR8 <- paste0(HR8," ","(",HR8l,"-",HR8u,")")

m9<-coxph(Surv(FU_time, critical)~cr0+comorb1+sexo, data=nutri1)
cox.zph(m9)
summary(m9)
BIC9<-round(BIC(m1)-BIC(m9),2)
C9<-round(summary(m9)$concordance[1],3)
B9 <- round(summary(m9)$coefficients[1],3)
p9 <- format.pval(summary(m9)$coefficients[13], eps = .001, digits = 3)
HR9 <- round(summary(m9)$conf.int[1],3)
HR9l <- round(summary(m9)$conf.int[7],3)
HR9u <- round(summary(m9)$conf.int[10],3)
HR9 <- paste0(HR9," ","(",HR9l,"-",HR9u,")")

m10<-coxph(Surv(FU_time, critical)~alb0+comorb1+sexo, data=nutri1)
cox.zph(m10)
summary(m10)
BIC10<-round(BIC(m1)-BIC(m10),2)
C10<-round(summary(m10)$concordance[1],3)
B10 <- round(summary(m10)$coefficients[1],3)
p10 <- format.pval(summary(m10)$coefficients[13], eps = .001, digits = 3)
HR10 <- round(summary(m10)$conf.int[1],3)
HR10l <- round(summary(m10)$conf.int[7],3)
HR10u <- round(summary(m10)$conf.int[10],3)
HR10 <- paste0(HR10," ","(",HR10l,"-",HR10u,")")

mSa02<-coxph(Surv(FU_time, critical)~sao2aa+comorb1+sexo, data=nutri1)
cox.zph(mSa02)
summary(mSa02)
BICSaO2<-round(BIC(m1)-BIC(mSa02),2)
CSaO2<-round(summary(mSa02)$concordance[1],3)
BSaO2 <- round(summary(mSa02)$coefficients[1],3)
pSaO2 <- format.pval(summary(mSa02)$coefficients[13], eps = .001, digits = 3)
HRSaO2 <- round(summary(mSa02)$conf.int[1],3)
HRSaO2l <- round(summary(mSa02)$conf.int[7],3)
HRSaO2u <- round(summary(mSa02)$conf.int[10],3)
HRSaO2 <- paste0(HRSaO2," ","(",HRSaO2l,"-",HRSaO2u,")")

m11<-coxph(Surv(FU_time, critical)~PhenoAge+comorb1+sexo, data=nutri1)
cox.zph(m11)
summary(m11)
BIC11<-round(BIC(m1)-BIC(m11),2)
C11<-round(summary(m11)$concordance[1],3)
B11 <- round(summary(m11)$coefficients[1],3)
p11 <- format.pval(summary(m11)$coefficients[13], eps = .001, digits = 3)
HR11 <- round(summary(m11)$conf.int[1],3)
HR11l <- round(summary(m11)$conf.int[7],3)
HR11u <- round(summary(m11)$conf.int[10],3)
HR11 <- paste0(HR11," ","(",HR11l,"-",HR11u,")")

m12<-coxph(Surv(FU_time, critical)~PhenoAccelAge+comorb1+sexo, data=nutri1)
cox.zph(m12)
summary(m12)
BIC12<-round(BIC(m1)-BIC(m12),2)
C12<-round(summary(m12)$concordance[1],3)
B12 <- round(summary(m12)$coefficients[1],3)
p12 <- format.pval(summary(m12)$coefficients[13], eps = .001, digits = 3)
HR12 <- round(summary(m12)$conf.int[1],3)
HR12l <- round(summary(m12)$conf.int[7],3)
HR12u <- round(summary(m12)$conf.int[10],3)
HR12 <- paste0(HR12," ","(",HR12l,"-",HR12u,")")

tab3_comorb<-data.frame("Parameter"=c("Age (years)","Glucose (mg/dL)","CRP","Alkaline Phosphatase","RDW","MCV","Lymphocytes (%)",
                                      "Leucocytes (x1000)","Creatinine","Albumin","SaO2","PhenoAge (years)", "PhenoAgeAccel"),
                        "DelBIC"=c(BIC1,BIC2,BIC3,BIC4,BIC5,BIC6,BIC7,BIC8,BIC9,BIC10,BICSaO2,BIC11,BIC12),
                        "C-statistic"=c(C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,CSaO2,C11,C12),
                        "Beta-coefficient"=c(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,BSaO2,B11,B12),
                        "HR (95%CI)"=c(HR1,HR2,HR3,HR4,HR5,HR6,HR7,HR8,HR9,HR10,HRSaO2,HR11,HR12),
                        "P-value"=c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,pSaO2,p11,p12))
tab3_comorb<-`names<-`(tab3_comorb,c("Parameter","Delta BIC","C-statistic","B-coefficient","HR (95%CI)","P-value"))
tab3_comorb<-align(flextable(tab3_comorb,cwidth = c(3,1,1,1,3,1)),align = "center",part = "all")
save_as_docx(tab3_comorb,path="tabla3_comorb.docx")#YA EST? GUARDADA


#Adverse outcomes models adjusted by sex and comorbidities

m1<-coxph(Surv(FU_time, critical)~lin0+GLU+pcr0+edad+sexo+comorb1, data=nutri1)
cox.zph(m1)
summary(m1)
C1<-round(summary(m1)$concordance[1],3)
B1 <- round(summary(m1)$coefficients[1],3)
B2 <- round(summary(m1)$coefficients[2],3)
B3 <- round(summary(m1)$coefficients[3],3)
B4 <- round(summary(m1)$coefficients[4],3)
B5 <- round(summary(m1)$coefficients[5],3)
B6 <- round(summary(m1)$coefficients[6],3)
HR1 <- round(summary(m1)$conf.int[1],3)
HR1l <- round(summary(m1)$conf.int[13],3)
HR1u <- round(summary(m1)$conf.int[19],3)
HR1 <- paste0(HR1," ","(",HR1l,"-",HR1u,")")
HR2 <- round(summary(m1)$conf.int[2],3)
HR2l <- round(summary(m1)$conf.int[14],3)
HR2u <- round(summary(m1)$conf.int[20],3)
HR2 <- paste0(HR2," ","(",HR2l,"-",HR2u,")")
HR3 <- round(summary(m1)$conf.int[3],3)
HR3l <- round(summary(m1)$conf.int[15],3)
HR3u <- round(summary(m1)$conf.int[21],3)
HR3 <- paste0(HR3," ","(",HR3l,"-",HR3u,")")
HR4 <- round(summary(m1)$conf.int[4],3)
HR4l <- round(summary(m1)$conf.int[16],3)
HR4u <- round(summary(m1)$conf.int[22],3)
HR4 <- paste0(HR4," ","(",HR4l,"-",HR4u,")")
HR5 <- round(summary(m1)$conf.int[5],3)
HR5l <- round(summary(m1)$conf.int[17],3)
HR5u <- round(summary(m1)$conf.int[23],3)
HR5 <- paste0(HR5," ","(",HR5l,"-",HR5u,")")
HR6 <- round(summary(m1)$conf.int[6],3)
HR6l <- round(summary(m1)$conf.int[18],3)
HR6u <- round(summary(m1)$conf.int[24],3)
HR6 <- paste0(HR6," ","(",HR6l,"-",HR6u,")")
p1 <- format.pval(summary(m1)$coefficients[25], eps = .001, digits = 3)
p2 <- format.pval(summary(m1)$coefficients[26], eps = .001, digits = 3)
p3 <- format.pval(summary(m1)$coefficients[27], eps = .001, digits = 3)
p4 <- format.pval(summary(m1)$coefficients[28], eps = .001, digits = 3)
p5 <- format.pval(summary(m1)$coefficients[29], eps = .001, digits = 3)
p6 <- format.pval(summary(m1)$coefficients[30], eps = .001, digits = 3)

m2<-coxph(Surv(FU_time, critical)~PhenoAccelAge+edad+sexo+comorb1, data=nutri1)
cox.zph(m2)
summary(m2)
C2<-round(summary(m2)$concordance[1],3)
B7 <- round(summary(m2)$coefficients[1],3)
B8 <- round(summary(m2)$coefficients[2],3)
B9 <- round(summary(m2)$coefficients[3],3)
B10 <- round(summary(m2)$coefficients[4],3)
HR7 <- round(summary(m2)$conf.int[1],3)
HR7l <- round(summary(m2)$conf.int[9],3)
HR7u <- round(summary(m2)$conf.int[13],3)
HR7 <- paste0(HR7," ","(",HR7l,"-",HR7u,")")
HR8 <- round(summary(m2)$conf.int[2],3)
HR8l <- round(summary(m2)$conf.int[10],3)
HR8u <- round(summary(m2)$conf.int[14],3)
HR8 <- paste0(HR8," ","(",HR8l,"-",HR8u,")")
HR9 <- round(summary(m2)$conf.int[3],3)
HR9l <- round(summary(m2)$conf.int[11],3)
HR9u <- round(summary(m2)$conf.int[15],3)
HR9 <- paste0(HR9," ","(",HR9l,"-",HR9u,")")
HR10 <- round(summary(m2)$conf.int[4],3)
HR10l <- round(summary(m2)$conf.int[12],3)
HR10u <- round(summary(m2)$conf.int[16],3)
HR10 <- paste0(HR10," ","(",HR10l,"-",HR10u,")")
p7 <- format.pval(summary(m2)$coefficients[17], eps = .001, digits = 3)
p8 <- format.pval(summary(m2)$coefficients[18], eps = .001, digits = 3)
p9 <- format.pval(summary(m2)$coefficients[19], eps = .001, digits = 3)
p10 <- format.pval(summary(m2)$coefficients[20], eps = .001, digits = 3)

tab4_comorb <- data.frame("Model"=c(" "," ","PhenoAge Components",paste0("C-Statistic"," ",C1)," "," "," ",
                                    "PhenoAgeAccel + Edad", paste0("C-Statistic"," ",C2), " "),
                          "Parameter"=c("Lymphocytes (%)","Glucose (mg/dL)","CRP","Age (years)","Male Sex", "Comorbidities",
                                        "PhenoAgeAccel","Age (years)","Male Sex", "Comorbidities"),
                          "B-coefficient"=c(B1,B2,B3,B4,B5,B6,B7,B8,B9,B10),
                          "HR (95%CI)"=c(HR1,HR2,HR3,HR4,HR5,HR6,HR7,HR8,HR9,HR10),
                          "p-value"=c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10))

tab4_comorb <-`names<-`(tab4_comorb,c("Model","Parameter","B-coefficient","HR (95%CI)","p-value"))
tab4_comorb <-align(flextable(tab4_comorb,cwidth = c(2,1.5,1,2,1)),align = "center",part = "all")
save_as_docx(tab4_comorb,path="tabla4_adv_aj.docx")#YA EST? GUARDADA


#### Cox PhenoAge ####
nutri1<-nutri[!is.na(nutri$PhenoAccelAge),]
nutri1<-nutri1%>%dplyr::select(critical, FU_time, defuncion, edad,
                               alb0, cr0, fa0, leu0, VCM, ADE, lin0, pcr0, GLU, PhenoAge, PhenoAccelAge, accel1, sexo, comorb,comorb1,comorb2,groups1, edad_cat)
set.seed(123); imp<-mice::mice(nutri1, m=1, maxit=1); nutri1<-complete(imp, "long")
nrow(nutri1)

nutri1$edad_cat<-factor(nutri1$edad_cat, labels = c("<50 years","50-70 years",">70 years"))

m12<-coxph(Surv(FU_time, defuncion)~PhenoAge*edad_cat+sexo+comorb1, data=nutri1)
cox.zph(m12)
summary(m12)


#### PhenoAge across severity spectrum ####
nutri1<-nutri[!is.na(nutri$PhenoAccelAge),]
nutri1<-nutri1%>%dplyr::select(critical, FU_time, defuncion, edad,
                              alb0, cr0, fa0, leu0, VCM, ADE, lin0, pcr0, GLU, PhenoAge, PhenoAccelAge, accel1, sexo, comorb,comorb2,groups1, edad_cat)
set.seed(123); imp<-mice::mice(nutri1, m=1, maxit=1); nutri1<-complete(imp, "long")
nrow(nutri1)

summary(lm(nutri$PhenoAge~nutri$edad))

nutri1$edad_cat<-factor(nutri1$edad_cat, labels = c("<50 years","50-70 years",">70 years"))

f1a<-ggplot(nutri1, aes(x=groups1, y=PhenoAge, fill=groups1))+geom_boxplot()+
  theme_classic()+ylab("PhenoAge")+labs(fill="COVID-19 severity")+
  ggpubr::stat_compare_means(comparisons = list(c("Severe", "Critical"), c("Lethal", "Critical"), c("Severe", "Lethal")))+
  facet_wrap(~comorb2,ncol=4,nrow=1)+xlab("")+scale_fill_manual(values=c("gold","darkorange2","red4"))

f1b<-ggplot(nutri1, aes(x=groups1, y=PhenoAccelAge, fill=groups1))+geom_boxplot()+
  theme_classic()+ylab("PhenoAgeAccel")+labs(fill="COVID-19 severity")+
  ggpubr::stat_compare_means(comparisons = list(c("Severe", "Critical"), c("Lethal", "Critical"), c("Severe", "Lethal")))+
  facet_wrap(~comorb2,ncol=4,nrow=1)+xlab("")+scale_fill_manual(values=c("gold","darkorange2","red4"))

f1c<-ggplot(nutri1, aes(x=groups1, y=PhenoAccelAge, fill=groups1))+geom_boxplot()+
  theme_classic()+ylab("PhenoAgeAccel")+labs(fill="COVID-19 severity")+
  ggpubr::stat_compare_means(comparisons = list(c("Severe", "Critical"), c("Lethal", "Critical"), c("Severe", "Lethal")))+
  facet_wrap(~edad_cat,ncol=4,nrow=1)+xlab("")+scale_fill_manual(values=c("gold","darkorange2","red4"))

f1d<-ggplot(nutri1, aes(x=edad_cat, y=PhenoAccelAge, fill=edad_cat))+geom_boxplot()+
  theme_classic()+ylab("PhenoAgeAccel")+labs(fill="Age group")+
  ggpubr::stat_compare_means()+
  facet_wrap(~groups1,ncol=4,nrow=1)+xlab("")+scale_fill_manual(values=c("gold","darkorange2","red4"))


f1e<-ggplot(nutri, aes(x=edad, y=PhenoAge, col=groups))+geom_point()+
  geom_smooth(method="lm",color="black")+theme_classic()+labs(col="Severity")+ylab("PhenoAge (years)")+xlab("Chronological age (years)")+
  scale_color_gradient2(midpoint=1, low="gold", mid="darkorange2",high="red4")

f1f<-ggplot(nutri, aes(x=edad, y=PhenoAccelAge, col=groups))+geom_point()+
  geom_smooth(method="lm", color="black")+theme_classic()+labs(col="Severity")+ylab("PhenoAgeAccel (years)")+xlab("Chronological age (years)")+
  scale_color_gradient2(midpoint=1, low="gold", mid="darkorange2",high="red4")


fig1<-ggarrange(f1a,f1b,f1c,f1d,f1e, f1f,ncol=2,nrow=3,labels = LETTERS[1:6])

ggsave(filename = "Figure1.jpg", 
       fig1,
       width = 45, 
       height = 30,
       units=c("cm"),
       dpi = 600,
       limitsize = FALSE)


f1.1a<-ggplot(nutri1, aes(x=comorb2, y=PhenoAge, fill=comorb2))+geom_boxplot()+
  theme_classic()+ylab("PhenoAge")+labs(fill="Comorbidities")+
  ggpubr::stat_compare_means()+
  facet_wrap(~groups1,ncol=3,nrow=1)+xlab("")+
  scale_fill_manual(values=c("paleturquoise2","deepskyblue3","dodgerblue4","midnightblue"))

f1.1b<-ggplot(nutri1, aes(x=comorb2, y=PhenoAccelAge, fill=comorb2))+geom_boxplot()+
  theme_classic()+ylab("PhenoAgeAccel")+labs(fill="Comorbidities")+
  ggpubr::stat_compare_means()+
  facet_wrap(~groups1,ncol=3,nrow=1)+xlab("")+
  scale_fill_manual(values=c("paleturquoise2","deepskyblue3","dodgerblue4","midnightblue"))

fig1.1<-ggarrange(f1.1a,f1.1b,ncol=1,nrow=2,labels = c("A", "B"))

ggsave(filename = "Figure1.1.jpg", 
       fig1.1,
       width = 40, 
       height = 25,
       units=c("cm"),
       dpi = 600,
       limitsize = FALSE)

f1a1<-ggplot(nutri1, aes(x=groups1, y=PhenoAge, fill=groups1))+geom_boxplot()+
  theme_classic()+ylab("PhenoAge")+labs(fill="COVID-19 severity")+
  ggpubr::stat_compare_means(comparisons = list(c("Severe", "Critical"), c("Lethal", "Critical"), c("Severe", "Lethal")))+
  facet_wrap(~edad_65,ncol=4,nrow=1)+xlab("")+scale_fill_manual(values=c("gold","darkorange2","red4"))

f1b1<-ggplot(nutri1, aes(x=groups1, y=PhenoAccelAge, fill=groups1))+geom_boxplot()+
  theme_classic()+ylab("PhenoAgeAccel")+labs(fill="COVID-19 severity")+
  ggpubr::stat_compare_means(comparisons = list(c("Severe", "Critical"), c("Lethal", "Critical"), c("Severe", "Lethal")))+
  facet_wrap(~edad_65,ncol=4,nrow=1)+xlab("")+scale_fill_manual(values=c("gold","darkorange2","red4"))

f1c1<-ggplot(nutri1, aes(x=edad_65, y=PhenoAccelAge, fill=edad_65))+geom_boxplot()+
  theme_classic()+ylab("PhenoAgeAccel")+labs(fill="Age group")+
  ggpubr::stat_compare_means()+
  facet_wrap(~groups1,ncol=4,nrow=1)+xlab("")+scale_fill_manual(values=c("gold","red4"))


fig1.2<-ggarrange(f1a1,f1b1,f1c1,labels = c("A", "B", "C"), nrow=1)

ggsave(filename = "Figure1.2.jpg", 
       fig1.2,
       width = 55, 
       height = 25,
       units=c("cm"),
       dpi = 600,
       limitsize = FALSE)


#### ROC curves ####
nutri1<-nutri[!is.na(nutri$PhenoAccelAge),]
nutri1<-nutri1%>%dplyr::select(critical, FU_time, defuncion, edad,
                              alb0, cr0, fa0, leu0, VCM, ADE, lin0, pcr0, GLU, PhenoAge, PhenoAccelAge, accel1, sexo, comorb,comorb2,groups1)
set.seed(123); imp<-mice::mice(nutri1, m=1, maxit=1); nutri1<-complete(imp, "long")
nrow(nutri1)

r1<-pROC::roc(nutri1$critical, nutri1$PhenoAge, ci=T)
r2<-pROC::roc(nutri1$critical, nutri1$PhenoAccelAge, ci=T)
r3<-pROC::roc(nutri1$critical, nutri1$edad, ci=T)
r4<-pROC::roc(nutri1$critical, nutri1$pcr0, ci=T)
r5<-pROC::roc(nutri1$critical, nutri1$alb0, ci=T)
r6<-pROC::roc(nutri1$critical, nutri1$GLU, ci=T)
r7<-pROC::roc(nutri1$critical, nutri1$lin0, ci=T)

roc.test(r1, r2, method = "boot")
roc.test(r3, r1, method = "boot")
roc.test(r4, r1, method = "boot")
roc.test(r5, r1, method = "boot")
roc.test(r6, r1, method = "boot")
roc.test(r7, r1, method = "boot")

roc.list<-list(r1, r2, r3, r4, r5, r6, r7)
print(roc.list)
table1<-matrix(c(round(r1$ci[1],3), round(r1$auc,3), round(r1$ci[2],3),
                 round(r2$ci[1],3), round(r2$auc,3), round(r2$ci[2],3),
                 round(r3$ci[1],3), round(r3$auc,3), round(r3$ci[2],3),
                 round(r4$ci[1],3), round(r4$auc,3), round(r4$ci[2],3),
                 round(r5$ci[1],3), round(r5$auc,3), round(r5$ci[2],3),
                 round(r6$ci[1],3), round(r6$auc,3), round(r6$ci[2],3),
                 round(r7$ci[1],3), round(r7$auc,3), round(r7$ci[2],3)),ncol=3,byrow=T)
table1<-table1[,c(2,1,3)]

colnames(table1)<-c("AUC", "Low-CI","Up-CI")
rownames(table1)<-c("PhenoAge", "PhenoAgeAccel", "CA", "CRP", "Albumin", "Glucose", "Lymph (%)")
names(roc.list)<-c("PhenoAge", "PhenoAgeAccel", "CA", "CRP", "Albumin", "Glucose", "Lymph (%)")
names(roc.list)<-c("PhenoAge", "PhenoAgeAccel", "CA", "CRP", "Albumin", "Glucose", "Lymph (%)")

fig4a<-ggroc(roc.list, aes("linetype", "color"), size=1.5)+
  theme(legend.position="none")+
  theme_minimal() + 
  ggtitle("") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="solid")+
  xlab("False Positive Rate (1-Specificity)") + 
  ylab("True Positive Rate (Sensitivity)") +
  labs(col="Scores", linetype="Scores")+
  theme(
    panel.background = element_blank(), 
    axis.title.x = element_text(size =14, face = 'plain'),
    axis.title.y = element_text(size =14, face = 'plain'),
    axis.text.x = element_text(size = 12, face ='plain'),
    axis.text.y = element_text(size = 12, face ='plain'))+theme_classic()+
  annotation_custom(tableGrob(table1), xmin=-0.5, xmax=0, ymin=0.2, ymax=0.4)

r1<-pROC::roc(nutri1$defuncion, nutri1$PhenoAge, ci=T)
r2<-pROC::roc(nutri1$defuncion, nutri1$PhenoAccelAge, ci=T)
r3<-pROC::roc(nutri1$defuncion, nutri1$edad, ci=T)
r4<-pROC::roc(nutri1$defuncion, nutri1$pcr0, ci=T)
r5<-pROC::roc(nutri1$defuncion, nutri1$ADE, ci=T)
r6<-pROC::roc(nutri1$defuncion, nutri1$VCM, ci=T)
r7<-pROC::roc(nutri1$defuncion, nutri1$GLU, ci=T)
r8<-pROC::roc(nutri1$defuncion, nutri1$lin0, ci=T)

roc.test(r1, r2, method = "boot")
roc.test(r3, r2, method = "boot")
roc.test(r4, r2, method = "boot")
roc.test(r5, r2, method = "boot")
roc.test(r6, r2, method = "boot")
roc.test(r7, r2, method = "boot")
roc.test(r8, r2, method = "boot")


roc.list<-list(r1, r2, r3, r4, r5, r6, r7, r8)
print(roc.list)
table1<-matrix(c(round(r1$ci[1],3), round(r1$auc,3), round(r1$ci[2],3),
                 round(r2$ci[1],3), round(r2$auc,3), round(r2$ci[2],3),
                 round(r3$ci[1],3), round(r3$auc,3), round(r3$ci[2],3),
                 round(r4$ci[1],3), round(r4$auc,3), round(r4$ci[2],3),
                 round(r5$ci[1],3), round(r5$auc,3), round(r5$ci[2],3),
                 round(r6$ci[1],3), round(r6$auc,3), round(r6$ci[2],3),
                 round(r7$ci[1],3), round(r7$auc,3), round(r7$ci[2],3),
                 round(r8$ci[1],3), round(r8$auc,3), round(r8$ci[2],3)),ncol=3,byrow=T)
table1<-table1[,c(2,1,3)]

colnames(table1)<-c("AUC", "Low-CI","Up-CI")
rownames(table1)<-c("PhenoAge", "PhenoAgeAccel", "CA", "CRP", "RDW","MCV", "Glucose", "Lymph (%)")
names(roc.list)<-c("PhenoAge", "PhenoAgeAccel", "CA", "CRP", "RDW","MCV", "Glucose", "Lymph (%)")
names(roc.list)<-c("PhenoAge", "PhenoAgeAccel", "CA", "CRP", "RDW","MCV", "Glucose", "Lymph (%)")

fig4b<-ggroc(roc.list, aes("linetype", "color"), size=1.5)+
  theme(legend.position="none")+
  theme_minimal() + 
  ggtitle("") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="black", linetype="solid")+
  xlab("False Positive Rate (1-Specificity)") + 
  ylab("True Positive Rate (Sensitivity)") +
  labs(col="Scores", linetype="Scores")+
  theme(
    panel.background = element_blank(), 
    axis.title.x = element_text(size =14, face = 'plain'),
    axis.title.y = element_text(size =14, face = 'plain'),
    axis.text.x = element_text(size = 12, face ='plain'),
    axis.text.y = element_text(size = 12, face ='plain'))+theme_classic()+
  annotation_custom(tableGrob(table1), xmin=-0.5, xmax=0, ymin=0.2, ymax=0.4)

#### Decision analyses ####
library(rmda)
nutri1<-nutri[!is.na(nutri$PhenoAccelAge),]
nutri1<-nutri1%>%dplyr::select(critical, FU_time, defuncion, edad,sao2aa,
                               alb0, cr0, fa0, leu0, VCM, ADE, lin0, pcr0, GLU, PhenoAge, PhenoAccelAge)%>% drop_na()
set.seed(123); imp<-mice::mice(nutri1, m=1, maxit=1); nutri1<-complete(imp, "long")


full.model_apparent1 <- decision_curve(critical~PhenoAge,
                                       data = nutri1,
                                       thresholds = seq(0, 1, by = .1),
                                       confidence.intervals = 'none')
full.model_apparent2 <- decision_curve(critical~PhenoAccelAge,
                                       data = nutri1,
                                       thresholds = seq(0, 1, by = .1),
                                       confidence.intervals = 'none')
full.model_apparent3 <- decision_curve(critical~edad,
                                       data = nutri1,
                                       thresholds = seq(0, 1, by = .1),
                                       confidence.intervals = 'none')
full.model_apparent4 <- decision_curve(critical~pcr0,
                                       data = nutri1,
                                       thresholds = seq(0,1, by = .1),
                                       confidence.intervals = 'none')
full.model_apparent5 <- decision_curve(critical~alb0,
                                       data = nutri1,
                                       thresholds = seq(0, 1, by = .1),
                                       confidence.intervals = 'none')


fig4c<-as.ggplot(~plot_decision_curve(list(full.model_apparent1, full.model_apparent2,full.model_apparent3, full.model_apparent4, 
                         full.model_apparent5),
                    curve.names = c('PhenoAge', 'PhenoAccelAge', 'Chronological Age', 'CRP', 'Albumin'),
                    col = c('red', 'blue', 'green', 'orange', 'purple', "yellow", "pink"),
                    lty = c(1,2,1,2,1,2,1),
                    lwd = c(3,2,3,2,3,2,3),
                    legend.position = 'topright'))+theme(plot.margin=unit(c(-1,-1,-1,-2.5),"cm"))


full.model_apparent1 <- decision_curve(defuncion~PhenoAge,
                                       data = nutri1,
                                       thresholds = seq(0, 0.7, by = .1),
                                       confidence.intervals = 'none')
full.model_apparent2 <- decision_curve(defuncion~PhenoAccelAge,
                                       data = nutri1,
                                       thresholds = seq(0, 0.7, by = .1),
                                       confidence.intervals = 'none')
full.model_apparent3 <- decision_curve(defuncion~edad,
                                       data = nutri1,
                                       thresholds = seq(0, 0.7, by = .1),
                                       confidence.intervals = 'none')
full.model_apparent4 <- decision_curve(defuncion~pcr0,
                                       data = nutri1,
                                       thresholds = seq(0, 0.7, by = .1),
                                       confidence.intervals = 'none')
full.model_apparent5 <- decision_curve(defuncion~alb0,
                                       data = nutri1,
                                       thresholds = seq(0, 0.7, by = .1),
                                       confidence.intervals = 'none')

fig4d<-as.ggplot(~plot_decision_curve(list(full.model_apparent1, full.model_apparent2,full.model_apparent3, full.model_apparent4, 
                         full.model_apparent5),
                    curve.names = c('PhenoAge', 'PhenoAccelAge', 'Chronological Age', 'SpO2', 'CRP', 'Albumin', 'PhenoAge+SpO2'),
                    col = c('red', 'blue', 'green', 'orange', 'purple', "yellow", "pink"),
                    lty = c(1,2,1,2,1,2,1),
                    lwd = c(3,2,3,2,3,2,3),
                    legend.position = 'topright'))+theme(plot.margin=unit(c(-1,-1,-1,-2.5),"cm"))

f3<-ggarrange(fig4a, fig4b, fig4c, fig4d,labels = LETTERS[1:4])

ggsave(file = "Figure4.jpg", 
       f3,
       bg = "transparent",
       width = 60, 
       height = 30,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)

#### Accelerated aging ####
library(fmsb)
colors_border<-c("#06C2D7","#9E3932")
colors_in<- rgb(col2rgb(colors_border)[1,],col2rgb(colors_border)[2,],
                col2rgb(colors_border)[3,],max=255,alpha=(50)*(255/100))

# Kaplan-meier curves
nutri1<-nutri[!is.na(nutri$PhenoAccelAge),]
nutri1<-nutri1%>%dplyr::select(critical, FU_time, defuncion, edad,alb0, cr0, fa0, leu0, VCM, ADE, lin0, pcr0, GLU, accel1)
set.seed(123); imp<-mice::mice(nutri1, m=1, maxit=1); nutri1<-complete(imp, "long")
nrow(nutri1)

mod2_km<-survfit(Surv(FU_time, critical) ~ accel1, data = nutri1);mod2_km
summary(mod2_km)
f3d<-ggsurvplot(mod2_km, data = nutri1, size = 1,palette = colors_border,conf.int = T,
                risk.table = T,pval = TRUE,ggtheme = theme_classic(),xlab="Time (Days)",
                ylab="Survival probability",
                legend.labs = c("Physiological aging","Accelerated aging"),
                ylim= c(0.3,1.0), xlim=c(0, 30),
                break.y.by= c(0.1), break.x.by= c(3),
                pval.coord = c(0, 0.45))+theme_survminer(base_size = 9, base_family = "Arial")
fig3d<-ggarrange(f3d$plot, f3d$table, heights = c(2, 0.7), ncol = 1, nrow = 2);fig3d


mod1_km<-survfit(Surv(FU_time, defuncion) ~ accel1, data = nutri1);mod1_km
summary(mod1_km)
f3e<-ggsurvplot(mod1_km, data = nutri1, size = 1,palette = colors_border,conf.int = T,
                risk.table = T,pval = TRUE,ggtheme = theme_classic(),xlab="Time (Days)",
                ylab="Survival probability",
                legend.labs = c("Physiological aging","Accelerated aging"),
                ylim= c(0.3,1.0),xlim=c(0, 30),
                break.y.by= c(0.1),break.x.by= c(3),
                pval.coord = c(0, 0.45))+theme_survminer(base_size = 9,base_family = "Arial")
fig3e<-ggarrange(f3e$plot, f3e$table, heights = c(2, 0.7),ncol = 1, nrow = 2);fig3e
fig3_bottom<-ggarrange(fig3d, fig3e,labels = c("D","E"))


#Spiderplot: Phenoage components
accel1<-nutri%>%dplyr::select(edad,alb0,cr0,fa0,leu0,VCM,ADE,lin0,pcr0,GLU,accel1, groups1)
#set.seed(123);bestNormalize(accel1$edad)$chosen_transform #BoxCox
#set.seed(123);bestNormalize(accel1$alb0)$chosen_transform #BoxCox
#set.seed(123);bestNormalize(accel1$cr0)$chosen_transform #Yeo-Johnson
#set.seed(123);bestNormalize(accel1$fa0)$chosen_transform #BoxCox
#set.seed(123);bestNormalize(accel1$leu0)$chosen_transform #asinh(x)
#set.seed(123);bestNormalize(accel1$VCM)$chosen_transform #orderNorm
#set.seed(123);bestNormalize(accel1$ADE)$chosen_transform #Yeo-Johnson
#set.seed(123);bestNormalize(accel1$lin0)$chosen_transform #Yeo-Johnson
#set.seed(123);bestNormalize(accel1$pcr0)$chosen_transform #orderNorm
#set.seed(123);bestNormalize(accel1$GLU)$chosen_transform #orderNorm

accel1<-accel1%>%dplyr::transmute(edad=boxcox(edad)$x.t,alb0=boxcox(alb0)$x.t,cr0=yeojohnson(cr0)$x.t,
                              fa0=boxcox(fa0)$x.t,leu0=arcsinh_x(leu0)$x.t,VCM=orderNorm(VCM)$x.t,
                              ADE=yeojohnson(ADE)$x.t,lin0=yeojohnson(lin0)$x.t,pcr0=orderNorm(pcr0)$x.t,
                              GLU=orderNorm(GLU)$x.t,accel1, groups1)%>%drop_na()

accel2<-as.data.frame(scale(accel1[,1:10]))
accel2$accel1<-accel1$accel1
accel2$groups<-accel1$groups1

t.test(accel2$edad~accel2$accel1)$p.value;wilcox.test(accel2$alb0~accel2$accel1)$p.value #age is not significant
wilcox.test(accel2$cr0~accel2$accel1)$p.value;wilcox.test(accel2$fa0~accel2$accel1)$p.value
t.test(accel2$leu0~accel2$accel1)$p.value;t.test(accel2$VCM~accel2$accel1)$p.value #MCV is not significant
wilcox.test(accel2$ADE~accel2$accel1)$p.value;t.test(accel2$lin0~accel2$accel1)$p.value
t.test(accel2$pcr0~accel2$accel1)$p.value;t.test(accel2$GLU~accel2$accel1)$p.value

data2<-aggregate(scale(accel2[,1:10]), list(accel2$accel1), median)[,-1]
colnames(data2) <- c("CA" , "Albumin***" , "Creatinin***" , "Alk-Ph***", "Leukocytes***" , "MCV", "RDW***", "Lymph***", "CRP***", "Glucose***" )
rownames(data2) <- c("Phys-Aging", "Accel-Aging")
data2 <- rbind(rep(2,10) , rep(-2,10) , data2)

fig3a<-ggplotify::as.ggplot(~radarchart(data2 , axistype=1 , seg=2,pcol=colors_border,pfcol=colors_in,
                      plwd=4 , plty=1, title="", cex.main=1.5,cglcol="grey20", cglty=1, axislabcol="gray35",
                      caxislabels=seq(c(-2),2, 2), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))


#Spiderplot: respiratory and metabolic function
accel3<-nutri%>%dplyr::select(fr,sao2aa,pafi0,ph0,lact0,imc,TG,TyG,pulsepr,fc,accel1,groups1)%>%drop_na()
#set.seed(123);bestNormalize(accel3$fr)$chosen_transform #Box Cox
#set.seed(123);bestNormalize(accel3$sao2aa)$chosen_transform #orderNorm
#set.seed(123);bestNormalize(accel3$pafi0)$chosen_transform #orderNorm
#set.seed(123);bestNormalize(accel3$ph0)$chosen_transform #orderNorm
#set.seed(123);bestNormalize(accel3$lact0)$chosen_transform #Box Cox
#set.seed(123);bestNormalize(accel3$imc)$chosen_transform #Yeo-Johnson
#set.seed(123);bestNormalize(accel3$TG)$chosen_transform #Yeo-Johnson
#set.seed(123);bestNormalize(accel3$TyG)$chosen_transform #orderNorm
#set.seed(123);bestNormalize(accel3$pulsepr)$chosen_transform #orderNorm
#set.seed(123);bestNormalize(accel3$fc)$chosen_transform #orderNorm

accel3<-accel3%>%transmute(fr=boxcox(fr)$x.t, sao2aa=orderNorm(sao2aa)$x.t, pafi0=orderNorm(pafi0)$x.t,
                           ph0=orderNorm(ph0)$x.t,lact0=boxcox(lact0)$x.t,imc=yeojohnson(imc)$x.t,
                           TG=yeojohnson(TG)$x.t,pulsepr=orderNorm(pulsepr)$x.t,
                           fc=orderNorm(fc)$x.t,accel1, groups1)

accel4<-as.data.frame(scale(accel3[,1:10]))
accel4$accel1<-accel3$accel1
accel4$groups<-accel3$groups1

wilcox.test(accel4$fr~accel4$accel1)$p.value;t.test(accel4$sao2aa~accel4$accel1)$p.value
t.test(accel4$pafi0~accel4$accel1)$p.value;wilcox.test(accel4$ph0~accel4$accel1)$p.value
wilcox.test(accel4$lact0~accel4$accel1)$p.value;t.test(accel4$imc~accel4$accel1)$p.value #BMI is not significant
wilcox.test(accel4$TG~accel4$accel1)$p.value
wilcox.test(accel4$pulsepr~accel4$accel1)$p.value;t.test(accel4$fc~accel4$accel1)$p.value #Pulsepr is not significant, HeartRate p=0.036

data4<-aggregate(accel4[,1:9], list(accel4$accel1), median)[,-c(1)]
colnames(data4) <- c("Respiratory rate***" , "SpO2***" , "PaFi***" , "Blood pH***", "Lactate***" , "BMI",
                    "Triglyicerides***", "Pulse\npressure", "Heart rate*")
rownames(data4) <- c("Phys-Aging", "Accel-Aging")
data4 <- rbind(rep(2,9) , rep(-2,9) , data4)

fig3b<-ggplotify::as.ggplot(~radarchart(data4 , axistype=1 , seg=2,pcol=colors_border,pfcol=colors_in,
                             plwd=4 , plty=1, title="", cex.main=1.5,cglcol="grey20", cglty=1, axislabcol="gray35",
                             caxislabels=seq(c(-2),2, 2), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))


#Spiderplot: Inflammatory biomarkers
accel5<-nutri%>%dplyr::select(fibri0ge00,pla0,dd0,ferritina0,dhl0,tpni0,bun0,alt0,ast0,bt0,accel1,groups1)%>%drop_na()
#set.seed(123);bestNormalize(accel5$fibrinogeno0)$chosen_transform #orderNorm
#set.seed(123);bestNormalize(accel5$pla0)$chosen_transform #orderNorm
#set.seed(123);bestNormalize(accel5$dd0)$chosen_transform #orderNorm
#set.seed(123);bestNormalize(accel5$ferritina0)$chosen_transform #orderNorm
#set.seed(123);bestNormalize(accel5$dhl0)$chosen_transform #asinh(x)
#set.seed(123);bestNormalize(accel5$tpni0)$chosen_transform #BoxCox
#set.seed(123);bestNormalize(accel5$bun0)$chosen_transform #Yeo-Johnson
#set.seed(123);bestNormalize(accel5$alt0)$chosen_transform #Yeo-Johnson
#set.seed(123);bestNormalize(accel5$ast0)$chosen_transform #BoxCox
#set.seed(123);bestNormalize(accel5$bt0)$chosen_transform #BoxCox

accel5<-accel5%>%transmute(fibrinogeno0=orderNorm(fibri0ge00)$x.t, pla0=orderNorm(pla0)$x.t, dd0=orderNorm(dd0)$x.t,
                           ferritina0=orderNorm(ferritina0)$x.t,dhl0=arcsinh_x(dhl0)$x.t,tpni0=boxcox(tpni0)$x.t,
                           bun0=yeojohnson(bun0)$x.t,alt0=yeojohnson(alt0)$x.t,ast0=boxcox(ast0)$x.t,
                           bt0=boxcox(bt0)$x.t,accel1, groups1)

accel6<-as.data.frame(scale(accel5[,1:10]))
accel6$accel1<-accel5$accel1
accel6$groups<-accel5$groups1

t.test(accel6$fibrinogeno0~accel6$accel1)$p.value;t.test(accel6$pla0~accel6$accel1)$p.value #Platelets p=0.02 (*)
t.test(accel6$dd0~accel6$accel1)$p.value;t.test(accel6$ferritina0~accel6$accel1)$p.value #D-dimer p=0.0015 (**), ferritin p=0.008 (**)
t.test(accel6$dhl0~accel6$accel1)$p.value;wilcox.test(accel6$tpni0~accel6$accel1)$p.value
t.test(accel6$bun0~accel6$accel1)$p.value;t.test(accel6$alt0~accel6$accel1)$p.value #ALT not significant
wilcox.test(accel6$ast0~accel6$accel1)$p.value;wilcox.test(accel6$bt0~accel6$accel1)$p.value #AST not significant, TBr p=0.00105 (**)

data6<-aggregate(accel6[,1:10], list(accel6$accel1), median)[,-c(1)]
colnames(data6) <- c("Fibrinogen***" , "Platelets*" , "D-Dimer**" , "Ferritin**", "LDH***" , "TPNI***", "BUN***", "ALT", "AST","Total\nbilirubin**")
rownames(data6) <- c("Phys-Aging", "Accel-Aging")
data6 <- rbind(rep(2,10) , rep(-2,10) , data6)

fig3c<-ggplotify::as.ggplot(~radarchart(data6 , axistype=1 , seg=2,pcol=colors_border,pfcol=colors_in,
                             plwd=4 , plty=1, title="", cex.main=1.5,cglcol="grey20", cglty=1, axislabcol="gray35",
                             caxislabels=seq(c(-2),2, 2), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(-1.5,-0.75,-2.5,-2.25),"cm"))

fig3_top<-ggarrange(fig3a, fig3b, fig3c, labels = c("A", "B", "C"),ncol=3)
fig2<-ggarrange(fig3_top,fig3_bottom, ncol=1, nrow=2, common.legend = T)
ggsave(filename = "Figure2.jpg", fig2, width = 40, height = 25, units=c("cm"), dpi = 300, limitsize = FALSE)



##### Clustering ####
#LASSO penalization
#Lethal
nutri3<-nutri%>%dplyr::transmute(edad, alb0, cr0, fa0, leu0, VCM, ADE, lin0, pcr0, GLU,FU_time,defuncion,critical, PhenoAge)%>%drop_na()
set.seed(123); imp<-mice::mice(nutri3, m=1, maxit=1); nutri3<-complete(imp, "long")
nrow(nutri1)

set.seed(123)
x <- as.matrix(nutri3%>%dplyr::select(edad, alb0, cr0, fa0, leu0, VCM, ADE, lin0, pcr0, GLU))
x<-scale(x)
y <- Surv(nutri3$FU_time+1, nutri3$defuncion)
n2 <- cv.glmnet(x, y, nfolds = 5,family="cox", alpha=1)
coef <- coef(n2, s = "lambda.min")
n3<-data.frame(name = coef@Dimnames[[1]][coef@i + 1], coefficient = coef@x)
coef
#Critical
set.seed(123)
x <- as.matrix(nutri3%>%dplyr::select(edad, alb0, cr0, fa0, leu0, VCM, ADE, lin0, pcr0, GLU))
x<-scale(x)
y <- Surv(nutri3$FU_time+1, nutri3$critical)
n2 <- cv.glmnet(x, y, nfolds = 5,family="cox", alpha=1)
coef <- coef(n2, s = "lambda.min")
n3<-data.frame(name = coef@Dimnames[[1]][coef@i + 1], coefficient = coef@x)
coef

#Correlation with age
data.frame(Variable=colnames(x)[-1],
           cor=round(as.numeric((sapply(as.data.frame(x[,2:10]),cor.test,x[,1]))[c((9*(1:9))-5)]),3),
           p.value=round(as.numeric((sapply(as.data.frame(x[,2:10]),cor.test,x[,1]))[c((9*(1:9))-6)]),3))


#Clusters
#We removed leucocytes (Cox regression for mortality with a LASSO regularization parameter) and
#alkaline phosphatase (correlation with age), we kept ADE to preserve the clusters' stability
set.seed(123)
nutri1<-nutri[!is.na(nutri$PhenoAccelAge),]
pheno<-nutri1 %>% dplyr::select(pin, GLU, lin0, pcr0, edad, alb0, cr0, ADE, VCM)
set.seed(123); imp<-mice::mice(pheno, m=1, maxit=1); pheno<-complete(imp, "long")
pheno1<-sapply(pheno, as.numeric)
pheno1<-scale(pheno1)[,-c(1)]

#Optimal number of clusters
fviz_nbclust(pheno1, kmeans, method = "silhouette")

res.nbclust <- NbClust(pheno1, distance = "euclidean",min.nc = 3, max.nc = 9, method = "complete", index ="all")
factoextra::fviz_nbclust(res.nbclust)+theme_minimal()+ggtitle("NbClust's optimal number of clusters")+theme_classic()


#Cluster stability
#The number of clusters with the most stability is 4
cf1 <- clusterboot(pheno1,B=1000,bootmethod=c("jitter"),clustermethod=kmeansCBI,krange=4,seed=123)
print(cf1)

#Cluster visualization
p1<-princomp(pheno1, cor=T)
set.seed(123);k1<-kmeansruns(pheno1, runs = 100, k=4)
fviz_cluster(k1,pheno1,palette=c("#9E0064","#AE480D","#2387CA","green3"))+theme_classic()
plot3d(p1$scores[,1:3], col=factor(k1$cluster,labels=c("#9E0064","#AE480D","#2387CA","green3")), size=1,  type='s')
table(k1$cluster)

#Addition of clusters to the database
pheno$cluster<-NULL
pheno$cluster<-k1$cluster
nutriCL<-merge(nutri, pheno[,c(1,10)], by="pin")
nutriCL$cluster<-factor(nutriCL$cluster)
nutriCL$cluster_cat<-relevel(nutriCL$cluster, ref=1)
nutriCL$cluster2<-factor(nutriCL$cluster,labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"))



##### Cluster characterization 1: age, Phenoage and PhenoAgeAccel ####
tapply(nutriCL$PhenoAge, nutriCL$cluster, median, na.rm=T)
tapply(nutriCL$edad, nutriCL$cluster, median, na.rm=T)
tapply(nutriCL$PhenoAccelAge, nutriCL$cluster, median, na.rm=T)

f10<-nutriCL%>%dplyr::select(cluster2, edad, cluster2) %>%
  ggplot(aes(x=cluster2, y=edad, fill=cluster2))+
  geom_boxplot()+ylab("Chronological Age (years)")+xlab(" ")+
  labs(fill="Cluster")+theme_classic()+stat_compare_means(size=7)+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  scale_fill_manual(values=c("#9E0064","#AE480D","#2387CA","green3"))+
  theme(text=element_text(size=25),axis.text=element_text(size=15))

f11<-nutriCL%>%dplyr::select(cluster2, PhenoAge, cluster2) %>%
  ggplot(aes(x=cluster2, y=PhenoAge, fill=cluster2))+
  geom_boxplot()+ylab("PhenoAge (years)")+xlab(" ")+
  labs(fill="Cluster")+theme_classic()+stat_compare_means(size=7)+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  scale_fill_manual(values=c("#9E0064","#AE480D","#2387CA","green3"))+
  theme(text=element_text(size=25),axis.text=element_text(size=15))

f12<-nutriCL%>%dplyr::select(cluster2, PhenoAccelAge, cluster2) %>%
  ggplot(aes(x=cluster2, y=PhenoAccelAge, fill=cluster2))+
  geom_boxplot()+ylab("PhenoAgeAccel (years)")+xlab(" ")+
  labs(fill="Cluster")+theme_classic()+stat_compare_means(size=7)+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  scale_fill_manual(values=c("#9E0064","#AE480D","#2387CA","green3"))+
  theme(text=element_text(size=25),axis.text=element_text(size=15),legend.position = "bottom",
        legend.key.size=unit(12,"mm"))

fig5b<-ggarrange(f10,f11,f12,nrow=1,common.legend = T,labels = LETTERS[1:3],legend="bottom",
                 legend.grob=get_legend(f12),font.label = list(size=25))

##### Cluster characterization 2: adverse outcomes ####

#Barplots
cov_ao<-nutriCL %>% group_by(cluster2)%>%
  summarise(critical=sum(critical, na.rm=T),uci=sum(uci, na.rm=T), intub=sum(intubado, na.rm=T), mort=sum(defuncion),n=n()) %>% 
  mutate(freq1 = critical/n, freq2 = uci/n, freq3 = intub/n, freq4=mort/n)%>%drop_na()
cov_ao

p.ao1<-format.pval(((table(nutriCL$cluster,nutriCL$critical)%>%prop.test())$p.value),eps=0.001,digits=3)
p.ao2<-format.pval(((table(nutriCL$cluster,nutriCL$uci)%>%prop.test())$p.value),eps=0.001,digits=3)
p.ao3<-format.pval(((table(nutriCL$cluster,nutriCL$intubado)%>%prop.test())$p.value),eps=0.001,digits=3)
p.ao4<-format.pval(((table(nutriCL$cluster,nutriCL$defuncion)%>%prop.test())$p.value),eps=0.001,digits=3)

clust<-c(rep(cov_ao$cluster2, 4))
freq<-c(cov_ao$freq1, cov_ao$freq2, cov_ao$freq3, cov_ao$freq4)
des<-c(rep(paste0("Adverse outcomes,"," ","p",p.ao1), 4), rep(paste0("ICU admission,"," ","p",p.ao2), 4),
       rep(paste0("IMV requirement,"," ","p",p.ao3), 4), rep(paste0("Lethality,"," ","p",p.ao4), 4))
cov1<-data.frame(clust, freq, des)

fig5a<-ggplot(cov1%>% arrange(), aes(y=freq, x=factor(clust),fill=factor(clust))) +
  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme_pubr()+ylab("Frequency (%)")+facet_wrap(~des, scales = "free")+
  scale_y_continuous(labels = scales::percent_format())+xlab("")+
  scale_fill_manual(values=c("#9E0064","#AE480D","#2387CA","green3"))+labs(fill="Cluster")+
  theme(text=element_text(size=25),axis.text=element_text(size=15),legend.position="bottom",legend.key.size = unit(12,"mm"))+
  geom_text(position = position_dodge(width= 0.8),aes(label=round(freq*100,2)), vjust=1.6,color="white", size=6.5)

Fig5<-ggarrange(fig5b,NULL,fig5a,nrow=3,heights=c(1,0.075,1.5),labels=c("","","D"),font.label = list(size=25))
ggsave(Fig5, filename = "Figure5.jpeg",
       width = 25, 
       height = 20,
       units=c("in"),
       dpi = 300,
       limitsize = FALSE)

#Cox regressions
nutriCL$cluster1<-relevel(nutriCL$cluster_cat, ref=4)
mk1<-coxph(Surv(FU_time, defuncion)~cluster1+edad+sexo+comorb, data=nutriCL)
summary(mk1)
mk2<-coxph(Surv(FU_time, critical)~cluster1+edad+sexo+comorb, data=nutriCL)
summary(mk2)

#Kaplan-Meier curves
km_clust1<-survfit(Surv(FU_time, defuncion) ~ cluster, data = nutriCL)
km_clust1
km1<-ggsurvplot(km_clust1, data = nutriCL, size = 1,palette = c("#9E0064","#AE480D","#2387CA","green3"),
                conf.int = F,risk.table = T,pval = TRUE,ggtheme = theme_classic(),xlab="Time (Days)",
                ylab="Survival probability",legend.labs = c("Cluster1","Cluster2","Cluster 3", "Cluster 4"),
                ylim= c(0,1.0),xlim=c(0, 30),break.y.by= c(0.1),break.x.by= c(3),
                pval.coord = c(0, 0.45))+theme_survminer(base_size = 9,base_family = "Arial")


km_clust2<-survfit(Surv(FU_time, critical) ~ cluster, data = nutriCL)
km_clust2
km2<-ggsurvplot(km_clust2, data = nutriCL, size = 1,palette = c("#9E0064","#AE480D","#2387CA","green3"),
                conf.int = F,risk.table = T,pval = TRUE,ggtheme = theme_classic(),xlab="Time (Days)",
                ylab="Survival probability",legend.labs = c("Cluster1","Cluster2","Cluster 3", "Cluster 4"),
                ylim= c(0,1.0),xlim=c(0, 30),break.y.by= c(0.1),break.x.by= c(3),
                pval.coord = c(0, 0.45))+theme_survminer(base_size = 9,base_family = "Arial")

km1b<-ggarrange(km1$plot, km1$table, heights = c(2, 0.7),ncol = 1, nrow = 2)
km2b<-ggarrange(km2$plot, km2$table, heights = c(2, 0.7),ncol = 1, nrow = 2)
ggarrange(km2b,km1b,nrow = 1,ncol=2)



#### Cluster characterization 3: Comorbidities ####

nutriCL$comorb3<-as.numeric(nutriCL$comorb)-1
nutriCL$sexo2<-as.numeric(nutriCL$sexo)-1
nutriCL$dm2_age40<-as.numeric((as.numeric(nutriCL$edad<=40)+nutriCL$dm2)>=2)
nutriCL$alta<-as.numeric(nutriCL$defuncion==0)

cov_comorb<-nutriCL %>% group_by(cluster2)%>%
  summarise(sexo=sum(sexo2, na.rm = T), comorb=sum(comorb3, na.rm=T), ob=sum(obesidad,na.rm=T), dm2=sum(dm2,na.rm=T), eod=sum(dm2_age40, na.rm = T),
            hipert=sum(hipertension,na.rm=T), cardio= sum(cardiovascular, na.rm=T), epoc=sum(epoc,na.rm=T), asma=sum(asma,na.rm=T),
            irc=sum(irc,na.rm=T), hep=sum(insuficiencia_hepatica, na.rm=T), inmuno=sum(inmunoupresion,na.rm=T), vih=sum(vih,na.rm=T),
            tbq=sum(tabaquismo, na.rm=T), alta=sum(alta, na.rm=T), n=n())%>%
  mutate(p1=sexo/n, p2=comorb/n, p3=ob/n, p4=dm2/n, p5=eod/n, p6=hipert/n, p7=cardio/n, p8=epoc/n, p9=asma/n, p10=irc/n, p11=hep/n, p12=inmuno/n,
         p13=vih/n, p14=tbq/n, p15=alta/n)%>%drop_na()

pc1<-format.pval(((table(nutriCL$cluster,nutriCL$sexo2)%>%prop.test())$p.value),eps=0.001,digits=3)
pc2<-format.pval(((table(nutriCL$cluster,nutriCL$comorb3)%>%prop.test())$p.value),eps=0.001,digits=3)
pc3<-format.pval(((table(nutriCL$cluster,nutriCL$obesidad)%>%prop.test())$p.value),eps=0.001,digits=3)
pc4<-format.pval(((table(nutriCL$cluster,nutriCL$dm2)%>%prop.test())$p.value),eps=0.001,digits=3)
pc5<-format.pval(((table(nutriCL$cluster,nutriCL$dm2_age40)%>%fisher.test())$p.value),eps=0.001,digits=3)
pc6<-format.pval(((table(nutriCL$cluster,nutriCL$hipertension)%>%prop.test())$p.value),eps=0.001,digits=3)
pc7<-format.pval(((table(nutriCL$cluster,nutriCL$cardiovascular)%>%fisher.test())$p.value),eps=0.001,digits=3)
pc8<-format.pval(((table(nutriCL$cluster,nutriCL$epoc)%>%fisher.test())$p.value),eps=0.001,digits=3)
pc9<-format.pval(((table(nutriCL$cluster,nutriCL$asma)%>%fisher.test())$p.value),eps=0.001,digits=3)
pc10<-format.pval(((table(nutriCL$cluster,nutriCL$irc)%>%fisher.test())$p.value),eps=0.001,digits=3)
pc11<-format.pval(((table(nutriCL$cluster,nutriCL$insuficiencia_hepatica)%>%fisher.test())$p.value),eps=0.001,digits=3)
pc12<-format.pval(((table(nutriCL$cluster,nutriCL$inmunoupresion)%>%fisher.test())$p.value),eps=0.001,digits=3)
pc13<-format.pval(((table(nutriCL$cluster,nutriCL$vih)%>%fisher.test())$p.value),eps=0.001,digits=3)
pc14<-format.pval(((table(nutriCL$cluster,nutriCL$tabaquismo)%>%fisher.test())$p.value),eps=0.001,digits=3)
pc15<-format.pval(((table(nutriCL$cluster,nutriCL$alta)%>%prop.test())$p.value),eps=0.001,digits=3)

data.frame(X=names(cov_comorb[,2:16]),Y=c(pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10,pc11,pc12,pc13,pc14,pc15))

clust2<-c(rep(cov_comorb$cluster2, 15))
freq2<-c(cov_comorb$p1, cov_comorb$p2, cov_comorb$p3, cov_comorb$p4,cov_comorb$p5,cov_comorb$p6,cov_comorb$p7,cov_comorb$p8,cov_comorb$p9,
        cov_comorb$p10, cov_comorb$p11, cov_comorb$p12, cov_comorb$p13, cov_comorb$p14, cov_comorb$p15)
des2<-c(rep("Male sex***", 4), rep(">=1 comorb**", 4), rep("Obesity", 4), rep("Type 2 diabetes***", 4), rep("Early-onset diabetes***", 4),
       rep("Hipertension***",4),rep("Cardiovascular disease",4),rep("COPD",4),rep("Asthma",4), rep("Chronic kidney disease",4),
       rep("Chronic liver disease",4), rep("Immunosuppression",4), rep("HIV infection",4), rep("Smoking",4), rep("Clinical improvement***",4))
facet2<-factor(des2, levels=des2[4*(1:15)])
cov2<-data.frame(clust2, freq2, des2,facet2)

fig6a<-ggplot(cov2%>% arrange(), aes(y=freq2, x=factor(clust2), fill=factor(clust2) )) +
  geom_bar(stat="identity", width=0.7, position=position_dodge(width=0.8))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  theme_pubr()+ylab("Frequency (%)")+facet_wrap(~facet2, nrow=3, ncol=5, as.table=T, scales = "free")+
  scale_y_continuous(labels = scales::percent_format())+xlab("")+
  scale_fill_manual(values=c("#9E0064","#AE480D","#2387CA","green3"))+labs(fill="Cluster")+
  theme(text=element_text(size=25),axis.text=element_text(size=15),legend.position="bottom",legend.key.size = unit(12,"mm"))+
  geom_text(position = position_dodge(width= 0.8),aes(label=round(freq2*100,2)), vjust=1.6,color="white", size=6.5)
ggsave(fig6a, filename = "Figure6.jpeg",width = 25,height = 15,units=c("in"),dpi = 300,limitsize = FALSE)

#Increasing number of comorbidities
nutriCL$comorb_o2<-as.ordered(nutriCL$comorb2)
cov_cb <- nutriCL %>% count(cluster, comorb_o2) %>% drop_na() %>%
  mutate(prop = as.numeric(table(nutriCL$comorb_o2, nutriCL$cluster)%>%prop.table(2)))
ggplot(data = cov_cb, aes(x = cluster, y = prop, fill = comorb_o2))+geom_bar(stat = 'identity', position = 'dodge', alpha = 2/3)+
  labs(x = 'Clusters', y = NULL, fill = 'Number of \ncomorbidities',title = 'Proportions by cluster')+
  theme(legend.position = "bottom",legend.title.align = 0.5)+scale_fill_manual(values=c("deepskyblue2","deepskyblue3","dodgerblue4","midnightblue"))



#### Cluster characterization 4: Clinical parameters ####
clust_col<-c("#9E0064","#AE480D","#2387CA","green3")
clust_fill<- rgb(col2rgb(clust_col)[1,],col2rgb(clust_col)[2,],
                 col2rgb(clust_col)[3,],max=255,alpha=(50)*(255/100))


#PhenoAge Components
CL1<-nutriCL%>%dplyr::transmute(edad=boxcox(edad)$x.t,alb0=boxcox(alb0)$x.t,cr0=yeojohnson(cr0)$x.t,
                                fa0=boxcox(fa0)$x.t,leu0=arcsinh_x(leu0)$x.t,VCM=orderNorm(VCM)$x.t,
                                ADE=yeojohnson(ADE)$x.t,lin0=yeojohnson(lin0)$x.t,pcr0=orderNorm(pcr0)$x.t,
                                GLU=orderNorm(GLU)$x.t,cluster)

CL2<-as.data.frame(scale(CL1[,1:10]))
CL2$cluster<-CL1$cluster
data.CL2<-aggregate(CL2[,1:10], list(CL2$cluster), median, na.rm=T)[,-c(1)]
colnames(data.CL2) <- c("CA" , "Albumin" , "Creatinin" , "Alk-Ph", "Leukocytes" , "MCV", "RDW", "Lymph", "CRP", "Glucose")
rownames(data.CL2) <- c("Cluster 1", "Cluster 2","Cluster 3","Cluster 4")
data.CL2 <- rbind(rep(2,10) , rep(-2,10) , data.CL2)

fig7.a1<-as.ggplot(~radarchart(data.CL2[c(1:2,3),] , axistype=1 , seg=2,pcol=clust_col[1],pfcol=clust_fill[1],
                               plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                               caxislabels=seq(c(-2),2, 2), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))
fig7.a2<-as.ggplot(~radarchart(data.CL2[c(1:2,4),] , axistype=1 , seg=2,pcol=clust_col[2],pfcol=clust_fill[2],
                               plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                               caxislabels=seq(c(-2),2, 2), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))
fig7.a3<-as.ggplot(~radarchart(data.CL2[c(1:2,5),] , axistype=1 , seg=2,pcol=clust_col[3],pfcol=clust_fill[3],
                               plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                               caxislabels=seq(c(-2),2, 2), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))
fig7.a4<-as.ggplot(~radarchart(data.CL2[c(1:2,6),] , axistype=1 , seg=2,pcol=clust_col[4],pfcol=clust_fill[4],
                               plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                               caxislabels=seq(c(-2),2, 2), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))


#Metabolic and respiratory function
CL3<-nutriCL%>%transmute(fr=boxcox(fr)$x.t, sao2aa=orderNorm(sao2aa)$x.t, pafi0=orderNorm(pafi0)$x.t,
                         ph0=orderNorm(ph0)$x.t,lact0=boxcox(lact0)$x.t,imc=yeojohnson(imc)$x.t,
                         TG=yeojohnson(TG)$x.t,pulsepr=orderNorm(pulsepr)$x.t,
                         fc=orderNorm(fc)$x.t,cluster)

CL4<-as.data.frame(scale(CL3[,1:9]))
CL4$cluster<-CL3$cluster
data.CL4<-aggregate(CL4[,1:9], list(CL4$cluster), median, na.rm=T)[,-c(1)]
colnames(data.CL4) <- c("Respiratory rate" , "SpO2" , "PaFi" , "Blood pH", "Lactate" , "BMI", "Triglyicerides", "Pulse\npressure", "Heart rate")
rownames(data.CL4) <- c("Cluster 1", "Cluster 2","Cluster 3","Cluster 4")
data.CL4 <- rbind(rep(2,9) , rep(-2,9) , data.CL4)

fig7.a5<-as.ggplot(~radarchart(data.CL4[c(1:2,3),] , axistype=1 , seg=2,pcol=clust_col[1],pfcol=clust_fill[1],
                               plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                               caxislabels=seq(c(-2),2, 2), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))
fig7.a6<-as.ggplot(~radarchart(data.CL4[c(1:2,4),] , axistype=1 , seg=2,pcol=clust_col[2],pfcol=clust_fill[2],
                               plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                               caxislabels=seq(c(-2),2, 2), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))
fig7.a7<-as.ggplot(~radarchart(data.CL4[c(1:2,5),] , axistype=1 , seg=2,pcol=clust_col[3],pfcol=clust_fill[3],
                               plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                               caxislabels=seq(c(-2),2, 2), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))
fig7.a8<-as.ggplot(~radarchart(data.CL4[c(1:2,6),] , axistype=1 , seg=2,pcol=clust_col[4],pfcol=clust_fill[4],
                               plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                               caxislabels=seq(c(-2),2, 2), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))


#Inflamatory biomarkers
CL5<-nutriCL%>%transmute(fibrinogeno0=orderNorm(fibri0ge00)$x.t, pla0=orderNorm(pla0)$x.t, dd0=orderNorm(dd0)$x.t,
                         ferritina0=orderNorm(ferritina0)$x.t,dhl0=arcsinh_x(dhl0)$x.t,tpni0=boxcox(tpni0)$x.t,
                         bun0=yeojohnson(bun0)$x.t,alt0=yeojohnson(alt0)$x.t,ast0=boxcox(ast0)$x.t,
                         bt0=boxcox(bt0)$x.t, CV=sqrt(cargaviralct),cluster)
data.CL6$`Viral load`
CL6<-as.data.frame(scale(CL5[,1:11]))
CL6$cluster<-CL5$cluster
data.CL6<-aggregate(CL6[,1:11], list(CL6$cluster), median, na.rm=T)[,-c(1)]
colnames(data.CL6) <- c("Fibrinogen" , "Platelets" , "D-Dimer" , "Ferritin", "LDH" , "TPNI", "BUN", "ALT", "AST","Total\nbilirubin", "Viral load CT")
rownames(data.CL6) <- c("Cluster 1", "Cluster 2","Cluster 3","Cluster 4")
data.CL6 <- rbind(rep(2,10) , rep(-2,10) , data.CL6)

fig7.a9<-as.ggplot(~radarchart(data.CL6[c(1:2,3),] , axistype=1 , seg=2,pcol=clust_col[1],pfcol=clust_fill[1],
                               plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                               caxislabels=seq(c(-2),2, 2), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))
fig7.a10<-as.ggplot(~radarchart(data.CL6[c(1:2,4),] , axistype=1 , seg=2,pcol=clust_col[2],pfcol=clust_fill[2],
                                plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                                caxislabels=seq(c(-2),2, 2), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))
fig7.a11<-as.ggplot(~radarchart(data.CL6[c(1:2,5),] , axistype=1 , seg=2,pcol=clust_col[3],pfcol=clust_fill[3],
                                plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                                caxislabels=seq(c(-2),2, 2), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))
fig7.a12<-as.ggplot(~radarchart(data.CL6[c(1:2,6),] , axistype=1 , seg=2,pcol=clust_col[4],pfcol=clust_fill[4],
                                plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                                caxislabels=seq(c(-2),2, 2), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))


## Symptoms

x1<-prop.table(table(nutriCL$fiebre,nutriCL$cluster),2)*100
x2<-prop.table(table(nutriCL$tos,nutriCL$cluster),2)*100
x3<-prop.table(table(nutriCL$disnea,nutriCL$cluster),2)*100
x4<-prop.table(table(nutriCL$cefalea,nutriCL$cluster),2)*100
x5<-prop.table(table(nutriCL$mialgias,nutriCL$cluster),2)*100
x6<-prop.table(table(nutriCL$artralgias,nutriCL$cluster),2)*100
x7<-prop.table(table(nutriCL$malestargeneral,nutriCL$cluster),2)*100

p1<-x1[c(2),]
p2<-x2[c(2),]
p3<-x3[c(2),]
p4<-x4[c(2),]
p5<-x5[c(2),]
p6<-x6[c(2),]
p7<-x7[c(2),]

table_sintomas<-matrix(c(p1,p2,p3,p4,p5,p6,p7),ncol=7,byrow=F)

rownames(table_sintomas) <- c("Cluster 1","Cluster 2","Cluster 3","Cluster 4")
colnames(table_sintomas) <- c("Fever","Cough", "Dyspnea",
                              "Headache","Myalgias","Arthralgias","General malaïse")
table_sintomas_2<-as.data.frame(table_sintomas)
table_sintomas_3<- rbind(rep(75,5) , rep(0,5) , table_sintomas_2)
table_sintomas_3<-as.data.frame(table_sintomas_3)[3:6,]
rownames(table_sintomas_3) <- c("Cluster 1", "Cluster 2","Cluster 3","Cluster 4")
table_sintomas_3 <- rbind(rep(100,10) , rep(0,10) , table_sintomas_3)


fig7.a13<-as.ggplot(~radarchart(table_sintomas_3[c(1:2,3),] , axistype=1 , seg=2,pcol=clust_col[1],pfcol=clust_fill[1],
                                plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                                caxislabels=seq(0,100, 50), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))
fig7.a14<-as.ggplot(~radarchart(table_sintomas_3[c(1:2,4),] , axistype=1 , seg=2,pcol=clust_col[2],pfcol=clust_fill[2],
                                plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                                caxislabels=seq(0,100, 50), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))
fig7.a15<-as.ggplot(~radarchart(table_sintomas_3[c(1:2,5),] , axistype=1 , seg=2,pcol=clust_col[3],pfcol=clust_fill[3],
                                plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                                caxislabels=seq(0,100, 50), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))
fig7.a16<-as.ggplot(~radarchart(table_sintomas_3[c(1:2,6),] , axistype=1 , seg=2,pcol=clust_col[4],pfcol=clust_fill[4],
                                plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                                caxislabels=seq(0,100, 50), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))


x1<-prop.table(table(nutriCL$diarrea,nutriCL$cluster),2)*100
x2<-prop.table(table(nutriCL$dolorcardiotoracico,nutriCL$cluster),2)*100
x3<-prop.table(table(nutriCL$rinorrea,nutriCL$cluster),2)*100
x4<-prop.table(table(nutriCL$anosmia,nutriCL$cluster),2)*100
x5<-prop.table(table(nutriCL$vomito,nutriCL$cluster),2)*100
x6<-prop.table(table(nutriCL$conjuntivitis,nutriCL$cluster),2)*100
x7<-prop.table(table(nutriCL$cianosis,nutriCL$cluster),2)*100
x8<-prop.table(table(nutriCL$disgueusia,nutriCL$cluster),2)*100
x9<-prop.table(table(nutriCL$dolorabdominal,nutriCL$cluster),2)*100


#Grafica Agrupada
p1<-x1[c(2),]
p2<-x2[c(2),]
p3<-x3[c(2),]
p4<-x4[c(2),]
p5<-x5[c(2),]
p6<-x6[c(2),]
p7<-x7[c(2),]
p8<-x8[c(2),]
p9<-x9[c(2),]

table_sintomas<-matrix(c(p1,p2,p3,p4,p5,p6,p7,p8,p9),ncol=9,byrow=F)

rownames(table_sintomas) <- c("Cluster 1","Cluster 2","Cluster 3","Cluster 4")
colnames(table_sintomas) <- c("Diarrhea","Chest Pain","Rhinorrhea","Anosmia","Vomiting",
                              "Conjunctivitis","Cyanosis","Dysgeusia", "Abdominal pain")
table_sintomas_2<-as.data.frame(table_sintomas)
table_sintomas_3<- rbind(rep(75,5) , rep(0,5) , table_sintomas_2)
table_sintomas_3<-as.data.frame(table_sintomas_3)[3:6,]
rownames(table_sintomas_3) <- c("Cluster 1", "Cluster 2","Cluster 3","Cluster 4")
table_sintomas_3 <- rbind(rep(50,10) , rep(0,10) , table_sintomas_3)


fig7.a17<-as.ggplot(~radarchart(table_sintomas_3[c(1:2,3),] , axistype=1 , seg=2,pcol=clust_col[1],pfcol=clust_fill[1],
                               plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                               caxislabels=seq(0,50, 25), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))
fig7.a18<-as.ggplot(~radarchart(table_sintomas_3[c(1:2,4),] , axistype=1 , seg=2,pcol=clust_col[2],pfcol=clust_fill[2],
                                plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                                caxislabels=seq(0,50, 25), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))
fig7.a19<-as.ggplot(~radarchart(table_sintomas_3[c(1:2,5),] , axistype=1 , seg=2,pcol=clust_col[3],pfcol=clust_fill[3],
                                plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                                caxislabels=seq(0,50, 25), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))
fig7.a20<-as.ggplot(~radarchart(table_sintomas_3[c(1:2,6),] , axistype=1 , seg=2,pcol=clust_col[4],pfcol=clust_fill[4],
                                plwd=4 , plty=1, title="", cex.main=2,cglcol="grey2", cglty=1, axislabcol="gray35",
                                caxislabels=seq(0,50, 25), cglwd=0.8,vlcex=1))+theme(plot.margin=unit(c(0,0,-4,0),"cm"))
nutriCL$Cluster<-nutriCL$cluster2

g1<-ggplot(nutriCL, aes(x=Cluster, y=edad, fill=Cluster))+scale_fill_manual(values=clust_fill)+geom_boxplot()+
  theme(legend.position = "top",legend.title=element_text(size=20), 
        legend.text=element_text(size=19))
legend<-get_legend(g1)

fig7<-ggarrange(fig7.a1,fig7.a5,fig7.a9,fig7.a13,fig7.a17,
                fig7.a2,fig7.a6,fig7.a10,fig7.a14,fig7.a18,
                fig7.a3,fig7.a7,fig7.a11,fig7.a15,fig7.a19,
                fig7.a4,fig7.a8,fig7.a12,fig7.a16,fig7.a20,ncol = 5,nrow = 4)+geom_subview(x=0.5, y=0.985, subview=legend)
                
ggsave(fig7, filename = "Figure7.jpeg",width = 35,height = 24,units=c("in"),dpi = 400,limitsize = FALSE)


#### Cluster characterization 5: Symptoms ####

#Respiratory
nutriCL$sym<-(nutriCL$fr_cat+nutriCL$fiebre+nutriCL$disnea+nutriCL$tos+nutriCL$cianosis)
nutriCL$sym[nutriCL$sym==0]<-0;nutriCL$sym[nutriCL$sym==1]<-1
nutriCL$sym[nutriCL$sym==2]<-2;nutriCL$sym[nutriCL$sym>=3]<-3
nutriCL$sym<-factor(nutriCL$sym,labels=c("0","1","2",">=3"))
nutriCL$sym<-as.ordered(nutriCL$sym);table(nutriCL$sym)

#### PhenoAge ####

table(nutriCL$sym,nutriCL$cluster)
cov_sym <- nutriCL %>% count(cluster, sym) %>% drop_na() %>%
  mutate(prop = as.numeric(table(nutriCL$sym,nutriCL$cluster)%>%prop.table(2)))
ggplot(data = cov_sym, aes(x = cluster, y = prop, fill = sym))+geom_bar(stat = 'identity', position = 'dodge', alpha = 2/3) +
  labs(x = 'Clusters', y = NULL, fill = 'Number of \nrespiratory symptoms',title = 'Proportions by clusters')+
  theme(legend.position = "bottom",legend.title.align = 0.5)+scale_fill_brewer(palette = "Oranges")

#Non-respiratory
nutriCL$sym2<-apply((nutriCL[,30:46]),1,sum,na.rm=T)
nutriCL$sym2[nutriCL$sym2>=1&nutriCL$sym2<=5]<-1;nutriCL$sym2[nutriCL$sym2>=6&nutriCL$sym2<=9]<-2
nutriCL$sym2[nutriCL$sym2>=10]<-3;nutriCL$sym2<-factor(nutriCL$sym2,labels=c("1-5","6-10",">=10"))
nutriCL$sym2<-as.ordered(nutriCL$sym2);table(nutriCL$sym2)

table(nutriCL$sym2,nutriCL$cluster)
cov_sym2 <- nutriCL %>% count(cluster, sym2) %>% drop_na() %>%
  mutate(prop = as.numeric(table(nutriCL$sym2,nutriCL$cluster)%>%prop.table(2)))
ggplot(data = cov_sym2, aes(x = cluster, y = prop, fill = sym2))+geom_bar(stat = 'identity', position = 'dodge', alpha = 2/3) +
  labs(x = 'Clusters', y = NULL, fill = 'Number of \nnon-respiratory symptoms',title = 'Proportions by clusters')+
  theme(legend.position = "bottom",legend.title.align = 0.5)+scale_fill_manual(values=c("tan1","tan3","tan4"))


