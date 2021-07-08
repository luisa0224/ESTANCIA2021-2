####PAQUETES####
library(haven); library(tidyverse); library(dplyr); library(mice); library(nhanesA); 
library(VIM); library(bestNormalize); library(factoextra); library(MAMI);
library(glmnetUtils); library(corrplot); library(miceadds); library(factoextra); library(NbClust);
library(FactoMineR); library(ggpubr); library(cluster); library(fpc); library(MASS);
library(caret); library(rms); library(jtools); library(ggiraphExtra); library(mediation);
library(ggeffects); library(finalfit); library(Gmisc); library(Hmisc); library(ggdag);

#setwd("C:/Users/B4-SUR-4/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/ESTANCIA")
setwd("~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/ESTANCIA")

####FUNCIONES####
categorize<-function(x){
  #Categorizes discrete variable into 1 and 0, preserving 1. Missing data substitutes as 0.
  c1<-c()
  name<-deparse(substitute(x))
  for (value in x){
    if (is.na(value)){
      cat=0
    }else if (value==1){
      cat=1
    }else{
      cat=0
    }
    c1<-c(c1, cat)
  }
  return(name<-c1)
}

categorizeinv<-function(x){
  #Categorizes discrete variable into 1 and 0, preserving 0. Missing data substitutes as 0.
  c1<-c()
  name<-deparse(substitute(x))
  for (value in x){
    if (is.na(value)){
      cat=0
    }else if (value==0){
      cat=0
    }else{
      cat=1
    }
    c1<-c(c1, cat)
  }
  return(name<-c1)
}

categ<-function(x, lim){
  #With a given limit, dichotomizes continuous variables, with greater-than criteria. 
  #Missing data substitutes as 0.
  category<-c()
  for(value in x){
    if(is.na(value)){
      cat=0
    }else if (value<lim){
      cat=0
    }else{
      cat=1
    }
    category<-c(category,cat)
  }
  return(category)
}

inv.categ<-function(x, lim){
  #With a given limit, dichotomizes continuous variables, with less-than criteria.
  #missing data substitutes as 0.
  category<-c()
  for(value in x){
    if(is.na(value)){
      cat=0
    }else if (value>lim){
      cat=0
    }else{
      cat=1
    }
    category<-c(category,cat)
  }
  return(category)
}

cap<-function(x){
  #Substitutes outliers in continuous variable with value in third or first quartile, respectively.
  qnt <- quantile(x, probs=c(.25, .75), na.rm = T)
  caps <- quantile(x, probs=c(.05, .95), na.rm = T)
  H <- 1.5 * IQR(x, na.rm = T)
  x[x < (qnt[1] - H)] <- caps[1]
  x[x > (qnt[2] + H)] <- caps[2]
  return(x)
}

pool.BIC<-function(data, x1, x2){
  #Averages BIC obtained from multiple models given by mice-imputed datasets and determines best model.
  BIC1<-c()
  m<-max(data$.imp)
  for(i in 1:m){
    model.imputed1 <- glm((x1), data = data[which(data$.imp == m),])
    BIC1[i] <- BIC(model.imputed1)
  }
  B1<-c(mean(BIC1))
  BIC2<-c()
  m<-max(data$.imp)
  for(i in 1:m){
    model.imputed2 <- glm((x2), data = data[which(data$.imp == m),])
    BIC2[i] <- BIC(model.imputed2)
  }
  B2<-c(mean(BIC2)) 
  if(abs(B2-B1)>2){
    if(B2>B1){
      return(c("Select first model. BIC: ",B1))
    }else{
      return(c("Select second model. BIC: ",B2))
    }
  }else{
    return(c("Select simpler model. BIC1: ",B1,". BIC2: ",B2))
  }
}

pool.BIC.r<-function(data, x1, x2){
  #Averages BIC obtained from multiple models with random effects 
  #given by mice-imputed datasets and determines best model.
  BIC1<-c()
  m<-max(data$.imp)
  for(i in 1:m){
    model.imputed1 <- lmer((x1), data = data[which(data$.imp == m),])
    BIC1[i] <- BIC(model.imputed1)
  }
  B1<-c(mean(BIC1))
  BIC2<-c()
  m<-max(data$.imp)
  for(i in 1:m){
    model.imputed2 <- lmer((x2), data = data[which(data$.imp == m),])
    BIC2[i] <- BIC(model.imputed2)
  }
  B2<-c(mean(BIC2)) 
  if(abs(B2-B1)>2){
    if(B2>B1){
      return(c("Select first model. BIC: ",B1))
    }else{
      return(c("Select second model. BIC: ",B2))
    }
  }else{
    return(c("Select simpler model. BIC1: ",B1,". BIC2: ",B2))
  }
}

BICS<-function(x){
  bics=c()
  mod=c()
  name=deparse(substitute(x))
  for(i in 1:length(x)){
    b=mean(pool(x[[i]])$glanced$BIC) 
    bics=c(bics,b)
  }
  name=bics 
  return(data.frame("Model"=1:length(x),"BIC"=name))}


####ARMADO DE BASE####
####1. LABS####
INS<-nhanes("INS_J")%>%dplyr::select(SEQN, LBXIN)
GLU<-nhanes("GLU_J")%>%dplyr::select(SEQN, LBXGLU)
CRP<-nhanes("HSCRP_J")%>%dplyr::select(SEQN,LBXHSCRP)
HDL<-nhanes("HDL_J")%>%dplyr::select(SEQN, LBDHDD)
COL<-nhanes("TCHOL_J")%>%dplyr::select(SEQN, LBXTC)
GHG<-nhanes("GHB_J")%>%dplyr::select(SEQN, LBXGH)
BQC<-nhanes("BIOPRO_J")%>%dplyr::select(SEQN, LBXSGTSI, LBXSATSI, LBXSAL, LBXSASSI, LBXSCK, LBXSCR, LBXSTB,
                                        LBXSTP, LBXSTR, LBXSUA)

labs.2018<-merge(CRP,
                 merge(COL,
                       merge(HDL, 
                             merge(GHG,merge(BQC, 
                                             merge(INS, GLU, by="SEQN", all.x=T)),
                                             by="SEQN", all.x=T),
                                   by="SEQN", all.x=T), 
                             by="SEQN", all.x=T),
                       by="SEQN", all.x=T)

####2. ANTROPOMETRÍA####

BMX.2018<-nhanes("BMX_J")%>%dplyr::select(SEQN, BMXWT, BMXHT, BMXWAIST, BMXHIP)

DXA.2018<-nhanes("DXX_J")%>%dplyr::select(SEQN, DXXTRFAT, DXDTRPF, DXDTOFAT, DXDTOPF)

LUX.2018<-nhanes("LUX_J")%>%dplyr::select(SEQN,  LUXSMED, LUXSIQR, LUXCAPM, LUXCPIQR)

antro.2018<-merge(BMX.2018, merge(LUX.2018, DXA.2018, by="SEQN", all.x=T), by="SEQN", all.x=T)

####3. CUESTIONARIO####
diet.2018<-nhanes("DR2TOT_J")%>%dplyr::select(SEQN, DR2TKCAL, DR2TPROT, DR2TCARB, DR2TSUGR, DR2TFIBE, DR2TTFAT, 
                                              DR2TSFAT, DR2TMFAT, DR2TPFAT, DR2TCHOL)

BPQ.2018<-nhanes("BPQ_J")%>%dplyr::select(SEQN, BPQ020, BPQ040A)

MCQ.2018<-nhanes("MCQ_J")%>%dplyr::select(SEQN, MCQ160B, MCQ160C, MCQ160D, MCQ160E, MCQ160F)

CDQ.2018<-nhanes("CDQ_J")%>%dplyr::select(SEQN, CDQ010)

DIQ.2018<-nhanes("DIQ_J")%>%dplyr::select(SEQN, DIQ010, DID040, DIQ050, DIQ070)

quest.2018<-merge(MCQ.2018, 
                  merge(DIQ.2018, 
                        merge(diet.2018, 
                              merge(BPQ.2018, CDQ.2018, by="SEQN", all=T), 
                              by="SEQN", all.x=T), 
                        by="SEQN", all.x=T), 
                  by="SEQN", all.x=T)

####4. DEMOGRÁFICOS####

DEMO.2018<-nhanes("DEMO_J")%>%dplyr::select(SEQN, RIAGENDR, RIDAGEYR, RIDRETH1)

####5. MANEJO DE BASE####

d.nhanes.2018<-merge(antro.2018, 
                           merge(labs.2018, 
                                 merge(quest.2018, DEMO.2018, by="SEQN", all.x=T),
                                 by="SEQN", all.x=T),
                           by="SEQN", all.x=T)

d.nhanes.2018.0<-d.nhanes.2018[!duplicated(d.nhanes.2018$SEQN),]
nhanes.2018.1<-d.nhanes.2018.0[!is.na(d.nhanes.2018.0$LUXSMED),]
nhanes.2018.2<-nhanes.2018.1[!is.na(nhanes.2018.1$LBXSUA),]
nhanes.2018.3<-nhanes.2018.2[!is.na(nhanes.2018.2$DXDTOPF)&nhanes.2018.2$RIDAGEYR>=18,]

data0<-nhanes.2018.3

data0[, c(1:54)]<-sapply(data0[, c(1:54)], as.numeric)

data0$GENDER<-categorize(data0$RIAGENDR)
data0$BPQ020<-categorize(data0$BPQ020)
data0$BPQ040A<-categorize(data0$BPQ040A)
data0$BPQ020<-categorize(data0$BPQ020)
data0$MCQ160B<-categorize(data0$MCQ160B)
data0$MCQ160C<-categorize(data0$MCQ160C)
data0$MCQ160D<-categorize(data0$MCQ160D)
data0$MCQ160E<-categorize(data0$MCQ160E)
data0$MCQ160F<-categorize(data0$MCQ160F)
data0$MCQ220<-categorize(data0$MCQ220)
data0$MCQ010<-categorize(data0$MCQ010)
data0$MCQ030<-categorize(data0$MCQ030)
data0$DIQ010<-categorize(data0$DIQ010)
data0$DID040<-categorize(data0$DID040)
data0$DIQ050<-categorize(data0$DIQ050)
data0$DIQ070<-categorize(data0$DIQ070)
data0$CDQ010<-categorize(data0$CDQ010)

data0$BMI<-(data0$BMXWT/((data0$BMXHT)/100)^2)
data0$METSIR<-(log(2*data0$LBXGLU+data0$LBXSTR)*data0$BMI)/log(data0$LBDHDD)
data0$WHR<-(data0$BMXWAIST)/(data0$BMXHT)
data0$METSVF<-(4.466+0.011*(log(data0$METSIR)^3)+3.239*(log(data0$WHR)^3)+0.319*(data0$GENDER)
                +0.594*(log(data0$RIDAGEYR)))
data0$METSVF.C<-categ(data0$METSVF, 7.18)
data0$FLI<-(exp((0.953*log(data0$LBXSTR))+(0.139*data0$BMI)+(0.718*log(data0$LBXSGTSI))+(0.053*data0$BMXWAIST)-15.745))/
  (1+exp((0.953*log(data0$LBXSTR))+(0.139*data0$BMI)+(0.718*log(data0$LBXSGTSI))+(0.053*data0$BMXWAIST-15.745)))*100
data0$ALTAST<-data0$LBXSASSI/data0$LBXSATSI

data0$t_SUA<-orderNorm(data0$LBXSUA)$x.t
data0$t_LUXSMED<-orderNorm(data0$LUXSMED)$x.t
data0$t_LUXCAPM<-orderNorm(data0$LUXCAPM)$x.t
data0$t_FLI<-orderNorm(data0$FLI)$x.t
data0$t_ALTAST<-orderNorm(data0$ALTAST)$x.t
data0$t_METSVF<-orderNorm(data0$METSVF)$x.t
data0$t_WHR<-orderNorm(data0$WHR)$x.t
data0$t_DR2TKCAL<-orderNorm(data0$DR2TKCAL)$x.t
data0$t_DR2TPROT<-orderNorm(data0$DR2TPROT)$x.t
data0$t_DR2TCARB<-orderNorm(data0$DR2TCARB)$x.t
data0$t_DR2TSUGR<-orderNorm(data0$DR2TSUGR)$x.t
data0$t_DR2TFIBE<-orderNorm(data0$DR2TFIBE)$x.t
data0$t_DR2TTFAT<-orderNorm(data0$DR2TTFAT)$x.t
data0$t_DR2TSFAT<-orderNorm(data0$DR2TSFAT)$x.t
data0$t_DR2TMFAT<-orderNorm(data0$DR2TMFAT)$x.t
data0$t_DR2TPFAT<-orderNorm(data0$DR2TPFAT)$x.t
data0$t_DR2TCHOL<-orderNorm(data0$DR2TCHOL)$x.t
data0$t_LBXSCR<-orderNorm(data0$LBXSCR)$x.t

data0$t_med<-data0$t_SUA*data0$t_ALTAST
data0$c_med<-categ(data0$t_med,quantile(data0$t_med, 0.8, na.rm=T))

data0$c_SUA<-NULL
data0$c_SUA[data0$GENDER==0&data0$LBXSUA<4.85]<-0
data0$c_SUA[data0$GENDER==0&data0$LBXSUA>=4.85]<-1
data0$c_SUA[data0$GENDER==1&data0$LBXSUA<6]<-0
data0$c_SUA[data0$GENDER==1&data0$LBXSUA>=6]<-1

data0$c_LUXCAPM<-NULL
data0$c_LUXCAPM[data0$LUXCAPM<=238]<-0
data0$c_LUXCAPM[data0$LUXCAPM>238&data0$LUXCAPM<=260]<-1
data0$c_LUXCAPM[data0$LUXCAPM>260&data0$LUXCAPM<=290]<-2
data0$c_LUXCAPM[data0$LUXCAPM>290]<-3

data0$p_LUXCAPM<-NULL
data0$p_LUXCAPM<-ifelse(data0$c_LUXCAPM!=0, 1, 0)

data0$c_LUXSMED<-NULL
data0$c_LUXSMED[data0$LUXSMED<=2]<-0
data0$c_LUXSMED[data0$LUXSMED>2&data0$LUXSMED<=7.5]<-0
data0$c_LUXSMED[data0$LUXSMED>7.5&data0$LUXSMED<=10]<-1
data0$c_LUXSMED[data0$LUXSMED>10&data0$LUXSMED<=14]<-2
data0$c_LUXSMED[data0$LUXSMED>14]<-2

data0$p_LUXSMED<-NULL
data0$p_LUXSMED<-ifelse(data0$c_LUXSMED==0, 0, 1)

data0$c_FLI<-NULL
data0$c_FLI[data0$FLI<=20]<-0
data0$c_FLI[data0$FLI>20&data0$FLI<=59.5]<-1
data0$c_FLI[data0$FLI>60]<-2

data0$p_FLI<-NULL
data0$p_FLI<-ifelse(data0$c_FLI==2, 1, 0)

set.seed(123)
md.pattern(data0, rotate.names=TRUE)

aggr_plot<-aggr(data0, col=c('navyblue','red'), 
                numbers=TRUE, sortVars=TRUE, labels=names(data0), 
                cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
print(aggr_plot)

imp1<-mice(data=data0, m=5, maxit=5, printFlag = FALSE, set.seed=123)
data1<-complete(imp1, "long")

####CORRELACIÓN####
df1<-cbind(data0$t_SUA, data0$t_METSVF, data0$t_DR2TCARB, data0$t_DR2TKCAL, 
           data0$t_DR2TCHOL, data0$t_DR2TFIBE,data0$t_DR2TSUGR, data0$t_DR2TTFAT, 
           data0$t_DR2TSFAT, data0$t_DR2TMFAT,data0$t_DR2TPFAT, data0$t_DR2TPROT, 
           data0$t_LUXSMED, data0$t_LUXCAPM, data0$t_FLI)
res1<-cor.mtest(df1, conf.level = .95)
colnames(df1)<-c("Uric\nacid", "Visceral\nfat", "Carbs", "Calories", "Cholesterol", 
                 "Fiber", "Sugar", "Total\nfat","Sat\nfat", "Mono\nfat", "Poly\nfat", 
                 "Protein", "Fibrosis", "Steatosis", "FLI")
cor1<-cor(df1, method="pearson", use = "complete.obs")
corrplot.mixed(cor1,lower.col = "black", number.cex = 1, upper = "circle",
               p.mat=res1$p, sig.level=.05)

df2<-cbind(data0$t_SUA, data0$t_METSVF, data0$WHR, data0$DXDTRPF, 
           data0$t_LUXSMED, data0$t_LUXCAPM, data0$t_FLI, data0$t_ALTAST)
res2<-cor.mtest(df2, conf.level = .95)
colnames(df2)<-c("Uric\nacid", "METSVF", "WHR", "Trunk\nfat", "Fibrosis", "Steatosis", "FLI", "AST/ALT")
cor2<-cor(df2, method="pearson", use = "complete.obs")
corrplot.mixed(cor2,lower.col = "black", number.cex = 1, upper = "circle",
               p.mat=res2$p, sig.level=.05)

res.mi.cor<-micombine.cor(mi.res=data1, variables = c(63:65, 10)); res.mi.cor
attr(res.mi.cor,"r_matrix")

####SELECCIÓN DE VARIABLES####

####Regresión penalizada
set.seed(123)

d1<-complete(imp1, 1)
d2<-complete(imp1, 2)
d3<-complete(imp1, 3)
d4<-complete(imp1, 4)
d5<-complete(imp1, 5)

data3<-list(d1, d2, d3, d4, d5)

m1 <- mami(imp1, model.family="binomial", outcome=c("p_LUXSMED"), method="LASSO", alpha=0.5, kfold=10,
           report.exp=TRUE, missing.data = "imputed", var.remove=c("SEQN","BMXWT","BMXHT","BMXWAIST","BMXHIP",
                                                                    "LUXSIQR","LUXCPIQR","DXXTRFAT","DXDTRPF",
                                                                    "DXDTOFAT","DXDTOPF","LBXHSCRP","LBXTC",
                                                                    "LBDHDD","LBXGH","LBXSGTSI","LBXSATSI",
                                                                    "LBXSAL","LBXSASSI","LBXSCK","LBXSCR",
                                                                    "LBXSTB","LBXSTP","LBXSTR","LBXSUA","LBXIN",
                                                                    "LBXGLU","MCQ160B","MCQ160C","MCQ160D",
                                                                    "MCQ160E","MCQ160F","DIQ010","DID040","DIQ050",
                                                                    "DIQ070","BPQ020","BPQ040A","CDQ010",
                                                                    "RIAGENDR","RIDAGEYR","RIDRETH1", "GENDER",
                                                                    "BMI","METSIR","WHR","METSVF.C",
                                                                    "t_SUA","t_METSVF", "t_LUXSMED", "LUXCAPM",
                                                                    "c_LUXCAPM", "LUXSMED", "DR2TKCAL", "DR2TPROT", "DR2TCARB", "DR2TSUGR", 
                                                                    "DR2TFIBE", "DR2TTFAT", "DR2TSFAT", "DR2TMFAT", 
                                                                    "DR2TPFAT", "DR2TCHOL", "METSVF", "p_LUXCAPM", "c_LUXSMED",
                                                                    "t_FLI", "FLI", "ALTAST", "t_ALTAST", "t_LUXCAPM",
                                                                    "t_WHR", "t_med", "c_med", "c_FLI", "p_FLI", "c_SUA"))
summary(m1)

m2 <- mami(imp1, model.family="binomial", outcome=c("p_LUXCAPM"), method="LASSO", alpha=0.5, kfold=10, 
           report.exp=TRUE, missing.data = "imputed", var.remove=c("SEQN","BMXWT","BMXHT","BMXWAIST","BMXHIP",
                                                                    "LUXSIQR","LUXCPIQR","DXXTRFAT","DXDTRPF",
                                                                    "DXDTOFAT","DXDTOPF","LBXHSCRP","LBXTC",
                                                                    "LBDHDD","LBXGH","LBXSGTSI","LBXSATSI",
                                                                    "LBXSAL","LBXSASSI","LBXSCK","LBXSCR",
                                                                    "LBXSTB","LBXSTP","LBXSTR","LBXSUA","LBXIN",
                                                                    "LBXGLU","MCQ160B","MCQ160C","MCQ160D",
                                                                    "MCQ160E","MCQ160F","DIQ010","DID040","DIQ050",
                                                                    "DIQ070","BPQ020","BPQ040A","CDQ010",
                                                                    "RIAGENDR","RIDAGEYR","RIDRETH1","GENDER",
                                                                    "BMI","METSIR","WHR","METSVF.C",
                                                                    "t_SUA","t_METSVF", "t_LUXSMED", "LUXCAPM",
                                                                    "c_LUXSMED", "DR2TKCAL", "DR2TPROT", "DR2TCARB", "DR2TSUGR", 
                                                                    "DR2TFIBE", "DR2TTFAT", "DR2TSFAT", "DR2TMFAT", 
                                                                    "DR2TPFAT", "DR2TCHOL", "METSVF", "p_LUXSMED", "LUXSMED",
                                                                    "c_LUXCAPM", "t_FLI", "FLI", "ALTAST", "t_ALTAST", "t_LUXCAPM",
                                                                   "t_WHR", "t_med", "c_med", "c_FLI", "p_FLI", "c_SUA"))
summary(m2)


mod <- glmnet(x=matrix(c(c(d1$t_DR2TPROT, d1$t_DR2TCARB, d1$t_DR2TSUGR, d1$t_DR2TFIBE, d1$t_DR2TTFAT, 
                                   d1$t_DR2TSFAT, d1$t_DR2TMFAT, d1$t_DR2TPFAT, d1$t_DR2TCHOL)), ncol=9), 
                      y=matrix(c(d1$p_LUXSMED), ncol=1), 
                      family="binomial",alpha = 1, k=5)

mod_cva <- cv.glmnet(x=matrix(c(c(d1$t_DR2TPROT, d1$t_DR2TCARB, d1$t_DR2TSUGR, d1$t_DR2TFIBE, d1$t_DR2TTFAT, 
                                 d1$t_DR2TSFAT, d1$t_DR2TMFAT, d1$t_DR2TPFAT, d1$t_DR2TCHOL)), ncol=9), 
                      y=matrix(c(d1$p_LUXSMED), ncol=1), 
                      family="binomial",alpha = 1, k=5)

coef <- coef(mod_cva, s = "lambda.min")
n3<-data.frame(name = coef@Dimnames[[1]][coef@i + 1], coefficient = coef@x)
n3

mod2 <- glmnet(x=matrix(c(c(d1$t_DR2TPROT, d1$t_DR2TCARB, d1$t_DR2TSUGR, d1$t_DR2TFIBE, d1$t_DR2TTFAT, 
                           d1$t_DR2TSFAT, d1$t_DR2TMFAT, d1$t_DR2TPFAT, d1$t_DR2TCHOL)), ncol=9), 
              y=matrix(c(d1$p_LUXCAPM), ncol=1), 
              family="binomial",alpha = 1, k=5)

mod_cva2 <- cv.glmnet(x=matrix(c(c(d1$t_DR2TPROT, d1$t_DR2TCARB, d1$t_DR2TSUGR, d1$t_DR2TFIBE, d1$t_DR2TTFAT, 
                                  d1$t_DR2TSFAT, d1$t_DR2TMFAT, d1$t_DR2TPFAT, d1$t_DR2TCHOL)), ncol=9), 
                     y=matrix(c(d1$p_LUXCAPM), ncol=1), 
                     family="binomial",alpha = 1, k=5)

coef <- coef(mod_cva2, s = "lambda.min")
n4<-data.frame(name = coef@Dimnames[[1]][coef@i + 1], coefficient = coef@x)
n4

####PCA/CLUSTERING####
data2<-complete(imp1, 1)

df_diet<-scale(data2%>%dplyr::select(DR2TPROT, DR2TCARB, DR2TSUGR, DR2TFIBE, DR2TTFAT, 
                                     DR2TSFAT, DR2TMFAT, DR2TPFAT, DR2TCHOL))

clusters1<-NbClust(data = df_diet, method = "kmeans", distance = "euclidean", min.nc = 2, max.nc = 15, 
                   index = "all", alphaBeale = 0.1)
summary(clusters1)

data2$clusters1<-clusters1$Best.partition

res.km <- eclust(df_diet, "kmeans", nstart = 25)
print(res.km)
fviz_gap_stat(res.km$gap_stat)

set.seed(123)
k2 <- kmeansruns(df_diet, k = 2, runs=100)
summary(k2)
data2$cluster<-k2$cluster

fviz_cluster(k2, data = df_diet)

cf1<-clusterboot(df_diet, B=1000, bootmethod=c("boot"), clustermethod=kmeansCBI,
                 krange=2, seed=123)
print(cf1)

pca<-PCA(df_diet, scale = TRUE)

#Scree plot
fviz_eig(pca)

#Score plot
fviz_pca_ind(pca,
             col.var = "cos2",gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE)

fviz_pca_ind(pca, repel = FALSE, addEllipses=TRUE,ellipse.level=0.95)

dim1<-fviz_contrib(pca, choice = "var", axes = 1, top = 10)
dim2<-fviz_contrib(pca, choice = "var", axes = 2, top = 10)

plot1<-ggarrange(dim1, dim2, ncol=2, nrow=1)
print(plot1)

##Loading plot
fviz_pca_var(pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = FALSE)

##Biplot
fviz_pca_biplot(pca, repel = FALSE, col.var = "#2E9FDF", col.ind = "#696969")

# Eigenvalues
eig.val <- get_eigenvalue(pca); eig.val

# Results for Variables
res.var <- get_pca_var(pca)
res.var$coord          # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

# Results for individuals
res.ind <- get_pca_ind(pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation 

####CARACTERIZACIÓN DE LOS CLUSTERS####

sex<-c("Female", "Male")
names(sex)<-c("0","1")

cluster_labs<-c("High", "Low")

plot2.0<-ggplot(data=data2,aes(x=factor(clusters1), y=DR2TTFAT, fill=factor(clusters1)))+
  geom_boxplot()+
  stat_compare_means(label.y=350)+
  xlab("Diet cluster")+
  ylab("Total fat")+
  facet_wrap(~GENDER, labeller=labeller(GENDER=sex))+
  scale_x_discrete(labels=cluster_labs)+
  scale_fill_discrete(name=" Consumption\n cluster" ,labels=cluster_labs)+
  theme_minimal()

print(plot2.0)

plot2.1<-ggplot(data=data2,aes(x=factor(clusters1), y=DR2TSUGR, fill=factor(clusters1)))+
  geom_boxplot()+
  stat_compare_means(label.y=800)+
  xlab("Diet cluster")+
  ylab("Total sugar")+
  facet_wrap(~GENDER, labeller=labeller(GENDER=sex))+
  scale_x_discrete(labels=cluster_labs)+
  scale_fill_discrete(name=" Consumption\n cluster" ,labels=cluster_labs)+
  theme_minimal()

print(plot2.1)

plot2.2<-ggplot(data=data2,aes(x=factor(clusters1), y=DR2TCARB, fill=factor(clusters1)))+
  geom_boxplot()+
  stat_compare_means(label.y=900)+
  xlab("Diet cluster")+
  ylab("Total carbohydrates")+
  facet_wrap(~GENDER, labeller=labeller(GENDER=sex))+
  scale_x_discrete(labels=cluster_labs)+
  scale_fill_discrete(name=" Consumption\n cluster" ,labels=cluster_labs)+
  theme_minimal()

print(plot2.2)


plot2.3<-ggplot(data=data2,aes(x=factor(clusters1), y=DR2TMFAT, fill=factor(clusters1)))+
  geom_boxplot()+
  stat_compare_means(label.y=125)+
  xlab("Diet cluster")+
  ylab("Monoinsaturated fat")+
  facet_wrap(~GENDER, labeller=labeller(GENDER=sex))+
  scale_x_discrete(labels=cluster_labs)+
  scale_fill_discrete(name=" Consumption\n cluster" ,labels=cluster_labs)+
  theme_minimal()

print(plot2.3)

plot2<-ggarrange(plot2.0, plot2.3, plot2.1, plot2.2, ncol=2, nrow=2, common.legend = T)
print(plot2)

ggsave("caractclusters.png", plot=plot2, width=30, height=20, units="cm", dpi=300)

#Clusters-> alto y bajo 

obese<-c("Non-obese","Obese")
names(obese)<-c("0","1")

data2$c_METSVF<-factor(data2$METSVF.C, levels=c("0","1"), labels=c("Non-obese","Obese"))

plot3<-ggplot(data=data2,aes(x=factor(c_LUXCAPM), y=LBXSUA, fill=factor(c_LUXCAPM)))+
  geom_boxplot()+
  stat_compare_means(label.y=15)+
  xlab("Steatosis degree")+
  ylab("Serum uric acid [mg/dL]")+
  scale_x_discrete(labels=c("S0","S1", "S2", "S3"))+
  facet_wrap(~c_METSVF, labeller=labeller(METSVF.C=obese))+
  scale_fill_discrete(name="Degree of\nsteatosis" ,labels=c("S0","S1", "S2", "S3"))+
  theme_minimal()

print(plot3)

plot4<-ggplot(data=data2,aes(x=factor(c_LUXSMED), y=LBXSUA, fill=factor(c_LUXSMED)))+
  geom_boxplot()+
  stat_compare_means(label.y=15)+
  xlab("Fibrosis degree")+
  ylab("Serum uric acid [mg/dL]")+
  scale_x_discrete(labels=c("F0", "F1", "F2", "F3", "F4"))+
  scale_fill_discrete(name="Degree of\nfibrosis" ,labels=c("F0", "F1", "F2", "F3", "F4"))+
  facet_wrap(~c_METSVF, labeller=labeller(METSVF.C=obese))+
  theme_minimal()

print(plot4)

plot5<-ggarrange(plot3, plot4, ncol=1, common.legend = F)
ggsave("SUAynafld.png", plot=plot5, width=30, height=20, units="cm", dpi=300)

####MODELOS DE REGRESIÓN####

#Esteatosis
data2$sex<-factor(data2$GENDER, levels=c(0,1),labels=c("Female", "Male"))

m<-polr(as.factor(c_LUXCAPM)~METSVF+t_SUA+t_LBXSCR+RIDAGEYR, data=d1, Hess=T)
summary(m)

(ctable <- coef(summary(m)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

(ctable <- cbind(ctable, "p value" = p))
data2$p_LUXCAPM<-as.numeric(data2$p_LUXCAPM)

r1_1<-with(data=imp1, exp=polr(as.factor(c_LUXCAPM)~t_METSVF+t_SUA+t_LBXSCR, Hess=T))
summary(pool(r1_1), conf.int = TRUE, exponentiate = TRUE)
summary(polr(as.factor(c_LUXCAPM)~t_METSVF+t_SUA+t_LBXSCR, data=data1, Hess=T))
BIC(polr(as.factor(c_LUXCAPM)~t_METSVF+t_SUA+t_LBXSCR, Hess=T, data=data2))

r1<-glm(p_LUXCAPM~METSVF+t_SUA, family="binomial", data=data2)
summary(r1, exponentiate=TRUE)

r1_2<-with(data=imp1, exp=polr(as.factor(c_LUXCAPM)~WHR*DXDTOPF+t_SUA+GENDER+t_LBXSCR+RIDAGEYR, Hess=T))
summary(pool(r1_2), conf.int = TRUE, exponentiate = TRUE)
BIC(polr(as.factor(c_LUXCAPM)~t_METSVF+t_SUA+t_LBXSCR, Hess=T, data=data2))

r2<-glm(p_LUXCAPM~WHR*DXDTOPF+t_SUA+sex, family="binomial", data=data2)
summary(r2, exponentiate=TRUE)

r1_3<-with(data=imp1, exp=polr(as.factor(c_LUXCAPM)~WHR+t_SUA+GENDER+t_LBXSCR+RIDAGEYR, Hess=T))
summary(pool(r1_3), conf.int = TRUE, exponentiate = TRUE)
BIC(polr(as.factor(c_LUXCAPM)~t_METSVF+t_SUA+t_LBXSCR, Hess=T, data=data2))

r3<-glm(p_LUXCAPM~WHR+t_SUA+sex, family="binomial", data=data2)
summary(r3, exponentiate=TRUE)

r1_4<-with(data=imp1, exp=polr(as.factor(c_LUXCAPM)~BMI+t_SUA+GENDER+t_LBXSCR+RIDAGEYR, Hess=T))
summary(pool(r1_4), conf.int = TRUE, exponentiate = TRUE)

r4<-glm(p_LUXCAPM~BMI+t_SUA+sex, family="binomial", data=data2)
summary(r4, exponentiate=TRUE)

r1_5<-with(data=imp1, exp=polr(as.factor(c_LUXCAPM)~DXDTRPF+t_SUA+GENDER+t_LBXSCR+RIDAGEYR, Hess=T))
summary(pool(r1_5), conf.int = TRUE, exponentiate = TRUE)

r5<-glm(p_LUXCAPM~DXDTRPF+t_SUA+sex, family="binomial", data=data2)
summary(r5, exponentiate=TRUE)

p_st<-ggplot(data=data2, aes(y=t_SUA, x=t_LUXCAPM, color=factor(cluster)))+
  geom_smooth(method="loess")+
  ylab("Uric acid quartile")+
  xlab("Steatosis quartile")+
  scale_color_discrete(name=" Consumption\n cluster" ,labels=cluster_labs)+
  theme_minimal()
  theme_minimal()

p1<-ggPredict(r1, se=TRUE,interactive=TRUE)
#p2<-ggPredict(r2, se=TRUE,interactive=TRUE)
p3<-ggPredict(r3, se=TRUE,interactive=TRUE)
p4<-ggPredict(r4, se=TRUE,interactive=TRUE)
p5<-ggPredict(r5, se=TRUE,interactive=TRUE)

pp1<-ggPredict(r1, se=TRUE,interactive=F)+
  ylab("Esteatosis probability")+
  xlab("METS-VF")+
  theme_minimal()
#pp2<-ggPredict(r2, se=TRUE,interactive=F)
pp3<-ggPredict(r3, se=TRUE,interactive=F)+
  ylab("Esteatosis probability")+
  xlab("Waist-to-height ratio")+
  theme_minimal()
pp4<-ggPredict(r4, se=TRUE,interactive=F)+
  ylab("Esteatosis probability")+
  xlab("BMI")+
  theme_minimal()
pp5<-ggPredict(r5, se=TRUE,interactive=F)+
  ylab("Esteatosis probability")+
  xlab("Trunk fat [%]")+
  theme_minimal()

f2<-ggarrange(pp1,pp3,pp4,pp5, ncol=2, nrow=2, common.legend = T)

ggsave("steatosis.png", plot=f2, width=30, height=30, units="cm", dpi=300)

r1_6<-with(data=imp1, exp=polr(as.factor(c_LUXCAPM)~rcs(t_METSVF)+t_SUA+t_LBXSCR+RIDAGEYR, Hess=T))
summary(pool(r1_6), conf.int = TRUE, exponentiate = TRUE)

data2$s_METSVF<-rcs(data2$t_METSVF)

r6<-glm(p_LUXCAPM~s_METSVF+t_SUA+t_LBXSCR+RIDAGEYR, family="binomial", data=data2)
summary(r6, conf.int = TRUE, exponentiate = TRUE)

#Fatty liver index

m2<-polr(as.factor(c_FLI)~METSVF+t_SUA, data=d1, Hess=T)
summary(m2)

(ctable <- coef(summary(m2)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

(ctable <- cbind(ctable, "p value" = p))

r2_1<-with(data=imp1, exp=polr(as.factor(c_FLI)~t_METSVF+t_SUA+t_LBXSCR, Hess=T))
summary(pool(r2_1), conf.int = TRUE, exponentiate = TRUE)

r7<-glm(p_FLI~METSVF+t_SUA, family="binomial", data=data2)
summary(r7, exponentiate=TRUE)

r2_2<-with(data=imp1, exp=polr(as.factor(c_FLI)~WHR*DXDTOPF+t_SUA+GENDER+t_LBXSCR+RIDAGEYR, Hess=T))
summary(pool(r2_2), conf.int = TRUE, exponentiate = TRUE)

r8<-glm(p_FLI~WHR*DXDTOPF+t_SUA+sex, family="binomial", data=data2)
summary(r8, exponentiate=TRUE)

r2_3<-with(data=imp1, exp=polr(as.factor(c_FLI)~WHR+t_SUA+GENDER+t_LBXSCR+RIDAGEYR, Hess=T))
summary(pool(r2_3), conf.int = TRUE, exponentiate = TRUE)

r9<-glm(p_FLI~WHR+t_SUA+sex, family="binomial", data=data2)
summary(r9, exponentiate=TRUE)

r2_4<-with(data=imp1, exp=polr(as.factor(c_FLI)~BMI+t_SUA+GENDER+t_LBXSCR+RIDAGEYR, Hess=T))
summary(pool(r2_4), conf.int = TRUE, exponentiate = TRUE)

r10<-glm(p_FLI~BMI+t_SUA+sex, family="binomial", data=data2)
summary(r10, exponentiate=TRUE)

r2_5<-with(data=imp1, exp=polr(as.factor(c_FLI)~DXDTRPF+t_SUA+GENDER+t_LBXSCR+RIDAGEYR, Hess=T))
summary(pool(r2_5), conf.int = TRUE, exponentiate = TRUE)

r11<-glm(p_FLI~DXDTRPF+t_SUA+sex, family="binomial", data=data2)
summary(r11, exponentiate=TRUE)

p_fli<-ggplot(data=data2, aes(y=t_SUA, x=t_FLI, color=factor(cluster)))+
  geom_smooth(method="loess")+
  ylab("Uric acid quartile")+
  xlab("FLI quartile")+
  scale_color_discrete(name=" Consumption\n cluster" ,labels=cluster_labs)+
  theme_minimal()
  theme_minimal()

p7<-ggPredict(r7, se=TRUE,interactive=TRUE)
#p8<-ggPredict(r8, se=TRUE,interactive=TRUE)
p9<-ggPredict(r9, se=TRUE,interactive=TRUE)
p10<-ggPredict(r10, se=TRUE,interactive=TRUE)
p11<-ggPredict(r11, se=TRUE,interactive=TRUE)

pp7<-ggPredict(r7, se=TRUE,interactive=F)+
  ylab("FLI probability")+
  xlab("METS-VF")+
  theme_minimal()
#pp8<-ggPredict(r8, se=TRUE,interactive=F)
pp9<-ggPredict(r9, se=TRUE,interactive=F)+
  ylab("FLI probability")+
  xlab("Waist-to-height ratio")+
  theme_minimal()
pp10<-ggPredict(r10, se=TRUE,interactive=F)+
  ylab("FLI probability")+
  xlab("BMI")+
  theme_minimal()
pp11<-ggPredict(r11, se=TRUE,interactive=F)+
  ylab("FLI probability")+
  xlab("Trunk fat [%]")+
  theme_minimal()

 
f3<-ggarrange(pp7,pp9,pp10,pp11, ncol=2, nrow=2, common.legend = T)

ggsave("fli.png", plot=f3, width=30, height=30, units="cm", dpi=300)

r6<-glm(p_FLI~s_METSVF+t_SUA, family="binomial", data=data2)
summary(r6, conf.int = TRUE, exponentiate = TRUE)

#Fibrosis

m3<-polr(as.factor(c_LUXSMED)~METSVF+t_SUA+t_LBXSCR+RIDAGEYR, data=d1, Hess=T)
summary(m3)

(ctable <- coef(summary(m2)))
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

(ctable <- cbind(ctable, "p value" = p))

r3_1<-with(data=imp1, exp=polr(as.factor(c_LUXSMED)~t_METSVF+t_SUA+t_LBXSCR+RIDAGEYR, Hess=T))
summary(pool(r3_1), conf.int = TRUE, exponentiate = TRUE)

r13<-glm(p_LUXSMED~METSVF+t_SUA, family="binomial", data=data2)
summary(r13, exponentiate=TRUE)

r3_2<-with(data=imp1, exp=polr(as.factor(c_LUXSMED)~WHR*DXDTOPF+t_SUA+GENDER+t_LBXSCR+RIDAGEYR, Hess=T))
summary(pool(r3_2), conf.int = TRUE, exponentiate = TRUE)

r14<-glm(p_LUXSMED~WHR*DXDTOPF+t_SUA+sex, family="binomial", data=data2)
summary(r14, exponentiate=TRUE)

r3_3<-with(data=imp1, exp=polr(as.factor(c_LUXSMED)~WHR+t_SUA+GENDER+t_LBXSCR+RIDAGEYR, Hess=T))
summary(pool(r3_3), conf.int = TRUE, exponentiate = TRUE)

r15<-glm(p_LUXSMED~WHR+t_SUA+sex, family="binomial", data=data2)
summary(r15, exponentiate=TRUE)

r3_4<-with(data=imp1, exp=polr(as.factor(c_LUXSMED)~BMI+t_SUA+GENDER+t_LBXSCR+RIDAGEYR, Hess=T))
summary(pool(r3_4), conf.int = TRUE, exponentiate = TRUE)

r16<-glm(p_LUXSMED~BMI+t_SUA+sex, family="binomial", data=data2)
summary(r16, exponentiate=TRUE)

r3_5<-with(data=imp1, exp=polr(as.factor(c_LUXSMED)~DXDTRPF+t_SUA+GENDER+t_LBXSCR+RIDAGEYR, Hess=T))
summary(pool(r3_5), conf.int = TRUE, exponentiate = TRUE)

r17<-glm(p_LUXSMED~DXDTRPF+t_SUA+sex, family="binomial", data=data2)
summary(r17, exponentiate=TRUE)

p_fib<-ggplot(data=data2, aes(y=t_SUA, x=t_LUXSMED, color=factor(cluster)))+
  geom_smooth(method="loess")+
  ylab("Uric acid quartile")+
  xlab("Fibrosis quartile")+
  scale_color_discrete(name=" Consumption\n cluster" ,labels=cluster_labs)+
  theme_minimal()

l_p<-ggarrange(p_st, p_fli, p_fib, ncol=3, common.legend = T)

ggsave("lines.png", plot=l_p, width=30, height=10, units="cm", dpi=300)

p13<-ggPredict(r13, se=TRUE,interactive=TRUE)
#p14<-ggPredict(r14, se=TRUE,interactive=TRUE)
p15<-ggPredict(r15, se=TRUE,interactive=TRUE)
p16<-ggPredict(r16, se=TRUE,interactive=TRUE)
p17<-ggPredict(r17, se=TRUE,interactive=TRUE)

pp13<-ggPredict(r13, se=TRUE,interactive=F)+
  ylab("Fibrosis probability")+
  xlab("METS-VF")+
  theme_minimal()
#pp14<-ggPredict(r14, se=TRUE,interactive=F)
pp15<-ggPredict(r15, se=TRUE,interactive=F)+
  ylab("Fibrosis probability")+
  xlab("Waist-to-height ratio")+
  theme_minimal()
pp16<-ggPredict(r16, se=TRUE,interactive=F)+
  ylab("Fibrosis probability")+
  xlab("BMI")+
  theme_minimal()
pp17<-ggPredict(r17, se=TRUE,interactive=F)+
  ylab("Fibrosis probability")+
  xlab("Trunk fat [%]")+
  theme_minimal()

f4<-ggarrange(pp13,pp15, pp16, pp17, ncol=2, nrow=2, common.legend = T)

ggsave("fibrosis.png", plot=f4, width=30, height=30, units="cm", dpi=300)

mediator<-with(data=imp1, exp=lm(as.factor(c_LUXCAPM)~as.factor(c_LUXSMED)+t_SUA+RIDAGEYR+t_LBXSCR+GENDER))

model1<-with(data=imp1, exp=polr(as.factor(c_LUXCAPM)~as.factor(c_LUXSMED)+t_SUA+RIDAGEYR+t_LBXSCR+GENDER, Hess=T))
summary(pool(model1), conf.int=TRUE, exponenciate=TRUE)

model2<-with(data=imp1, exp=polr(as.factor(c_LUXCAPM)~as.factor(c_LUXSMED)+t_ALTAST+RIDAGEYR+t_LBXSCR+GENDER, Hess=T))
summary(pool(model2), conf.int=TRUE, exponenciate=TRUE)

## Regression tables ##
e_names<-c("Parameter","Beta","SE","p-value","95% CI")

se_r.1.1<-summary(pool(r1_1), conf.int = TRUE, exponenciate=TRUE)
e_model1<-`names<-`(data.frame(c("Visceral fat", "SUA", "E0-E1", "E1-E2", "E2-E3"),
                               round(se_r.1.1[c(1,2,4,5,6),c(2,3)],3),
                               format.pval(se_r.1.1[c(1,2,4,5,6),c(6)],digits=2,eps=0.001),
                               paste0(round(se_r.1.1[c(1,2,4,5,6),c(7)],3)," - ",round(se_r.1.1[c(1,2,4,5,6),c(8)],3))),e_names)
e_model1[nrow(e_model1)+1,]<-NA

se_r.1.2<-summary(pool(r1_2), conf.int = TRUE, exponenciate=TRUE)
e_model2<-`names<-`(data.frame(c("Visceral fat", "SUA", "E0-E1", "E1-E2", "E2-E3"),
                               round(se_r.1.2[c(7,3,8,9,10),c(2,3)],3),
                               format.pval(se_r.1.2[c(7,3,8,9,10),c(6)],digits=2,eps=0.001),
                               paste0(round(se_r.1.2[c(7,3,8,9,10),c(7)],3)," - ",round(se_r.1.2[c(7,3,8,9,10),c(8)],3))),e_names)
e_model2[nrow(e_model2)+1,]<-NA

se_r.1.3<-summary(pool(r1_3), conf.int = TRUE, exponenciate=TRUE)
e_model3<-`names<-`(data.frame(c("Visceral fat", "SUA", "E0-E1", "E1-E2", "E2-E3"),
                               round(se_r.1.3[c(1,2,6,7,8),c(2,3)],3),
                               format.pval(se_r.1.3[c(1,2,6,7,8),c(6)],digits=2,eps=0.001),
                               paste0(round(se_r.1.3[c(1,2,6,7,8),c(7)],3)," - ",round(se_r.1.3[c(1,2,6,7,8),c(8)],3))),e_names)
e_model3[nrow(e_model3)+1,]<-NA

se_r.1.4<-summary(pool(r1_4), conf.int = TRUE, exponenciate=TRUE)
e_model4<-`names<-`(data.frame(c("Visceral fat", "SUA", "E0-E1", "E1-E2", "E2-E3"),
                               round(se_r.1.4[c(1,2,6,7,8),c(2,3)],3),
                               format.pval(se_r.1.4[c(1,2,6,7,8),c(6)],digits=2,eps=0.001),
                               paste0(round(se_r.1.4[c(1,2,6,7,8),c(7)],3)," - ",round(se_r.1.4[c(1,2,6,7,8),c(8)],3))),e_names)
e_model4[nrow(e_model4)+1,]<-NA

se_r.1.5<-summary(pool(r1_5), conf.int = TRUE, exponenciate=TRUE)
e_model5<-`names<-`(data.frame(c("Visceral fat", "SUA", "E0-E1", "E1-E2", "E2-E3"),
                               round(se_r.1.5[c(1,2,6,7,8),c(2,3)],3),
                               format.pval(se_r.1.5[c(1,2,6,7,8),c(6)],digits=2,eps=0.001),
                               paste0(round(se_r.1.5[c(1,2,6,7,8),c(7)],3)," - ",round(se_r.1.5[c(1,2,6,7,8),c(8)],3))),e_names)
e_model5[nrow(e_model5)+1,]<-NA

e_names2<-c("Parameter2","Beta2","SE2","p-value2","95% CI2")

se_r.2.1<-summary(pool(r2_1), conf.int = TRUE, exponenciate=TRUE)
e_model6<-`names<-`(data.frame(c("Visceral fat", "SUA", "FLI0-FLI1", "FLI1-FLI2"),
                               round(se_r.2.1[c(1,2,4,5),c(2,3)],3),
                               format.pval(se_r.2.1[c(1,2,4,5),c(6)],digits=2,eps=0.001),
                               paste0(round(se_r.2.1[c(1,2,4,5),c(7)],3)," - ",round(se_r.2.1[c(1,2,4,5),c(8)],3))),e_names2)
e_model6[nrow(e_model6)+1,]<-NA
e_model6[nrow(e_model6)+1,]<-NA

se_r.2.2<-summary(pool(r2_2), conf.int = TRUE, exponenciate=TRUE)
e_model7<-`names<-`(data.frame(c("Visceral fat", "SUA", "FLI0-FLI1", "FLI1-FLI2"),
                               round(se_r.2.2[c(7,3,8,9),c(2,3)],3),
                               format.pval(se_r.2.2[c(7,3,8,9),c(6)],digits=2,eps=0.001),
                               paste0(round(se_r.2.2[c(7,3,8,9),c(7)],3)," - ",round(se_r.2.2[c(7,3,8,9),c(8)],3))),e_names2)
e_model7[nrow(e_model7)+1,]<-NA
e_model7[nrow(e_model7)+1,]<-NA

se_r.2.3<-summary(pool(r2_3), conf.int = TRUE, exponenciate=TRUE)
e_model8<-`names<-`(data.frame(c("Visceral fat", "SUA", "FLI0-FLI1", "FLI1-FLI2"),
                               round(se_r.2.3[c(1,2,6,7),c(2,3)],3),
                               format.pval(se_r.2.3[c(1,2,6,7),c(6)],digits=2,eps=0.001),
                               paste0(round(se_r.2.3[c(1,2,6,7),c(7)],3)," - ",round(se_r.2.3[c(1,2,6,7),c(8)],3))),e_names2)
e_model8[nrow(e_model8)+1,]<-NA
e_model8[nrow(e_model8)+1,]<-NA

se_r.2.4<-summary(pool(r2_4), conf.int = TRUE, exponenciate=TRUE)
e_model9<-`names<-`(data.frame(c("Visceral fat", "SUA", "FLI0-FLI1", "FLI1-FLI2"),
                               round(se_r.2.4[c(1,2,6,7),c(2,3)],3),
                               format.pval(se_r.2.4[c(1,2,6,7),c(6)],digits=2,eps=0.001),
                               paste0(round(se_r.2.4[c(1,2,6,7),c(7)],3)," - ",round(se_r.2.4[c(1,2,6,7),c(8)],3))),e_names2)
e_model9[nrow(e_model9)+1,]<-NA
e_model9[nrow(e_model9)+1,]<-NA

se_r.2.5<-summary(pool(r2_5), conf.int = TRUE, exponenciate=TRUE)
e_model10<-`names<-`(data.frame(c("Visceral fat", "SUA", "FLI0-FLI1", "FLI1-FLI2"),
                               round(se_r.2.5[c(1,2,6,7),c(2,3)],3),
                               format.pval(se_r.2.5[c(1,2,6,7),c(6)],digits=2,eps=0.001),
                               paste0(round(se_r.2.5[c(1,2,6,7),c(7)],3)," - ",round(se_r.2.5[c(1,2,6,7),c(8)],3))),e_names2)
e_model10[nrow(e_model10)+1,]<-NA
e_model10[nrow(e_model10)+1,]<-NA

e_names3<-c("Parameter3","Beta3","SE3","p-value3","95% CI3")

se_r.3.1<-summary(pool(r3_1), conf.int = TRUE, exponenciate=TRUE)
e_model11<-`names<-`(data.frame(c("Visceral fat", "SUA", "F0-F1", "F1-F2", "F2-F3", "F3-F4"),
                               round(se_r.3.1[c(1,2,5,6,7,8),c(2,3)],3),
                               format.pval(se_r.3.1[c(1,2,5,6,7,8),c(6)],digits=2,eps=0.001),
                               paste0(round(se_r.3.1[c(1,2,5,6,7,8),c(7)],3)," - ",round(se_r.3.1[c(1,2,5,6,7,8),c(8)],3))),e_names3)

se_r.3.2<-summary(pool(r3_2), conf.int = TRUE, exponenciate=TRUE)
e_model12<-`names<-`(data.frame(c("Visceral fat", "SUA", "F0-F1", "F1-F2", "F2-F3", "F3-F4"),
                               round(se_r.3.2[c(7,3,8,9,10,11),c(2,3)],3),
                               format.pval(se_r.3.2[c(7,3,8,9,10,11),c(6)],digits=2,eps=0.001),
                               paste0(round(se_r.3.2[c(7,3,8,9,10,11),c(7)],3)," - ",round(se_r.3.2[c(7,3,8,9,10,11),c(8)],3))),e_names3)

se_r.3.3<-summary(pool(r3_3), conf.int = TRUE, exponenciate=TRUE)
e_model13<-`names<-`(data.frame(c("Visceral fat", "SUA", "F0-F1", "F1-F2", "F2-F3", "F3-F4"),
                               round(se_r.3.3[c(1,2,6,7,8,9),c(2,3)],3),
                               format.pval(se_r.3.3[c(1,2,6,7,8,9),c(6)],digits=2,eps=0.001),
                               paste0(round(se_r.3.3[c(1,2,6,7,8,9),c(7)],3)," - ",round(se_r.3.3[c(1,2,6,7,8,9),c(8)],3))),e_names3)

se_r.3.4<-summary(pool(r3_4), conf.int = TRUE, exponenciate=TRUE)
e_model14<-`names<-`(data.frame(c("Visceral fat", "SUA", "F0-F1", "F1-F2", "F2-F3", "F3-F4"),
                               round(se_r.3.4[c(1,2,6,7,8,9),c(2,3)],3),
                               format.pval(se_r.3.4[c(1,2,6,7,8,9),c(6)],digits=2,eps=0.001),
                               paste0(round(se_r.3.4[c(1,2,6,7,8,9),c(7)],3)," - ",round(se_r.3.4[c(1,2,6,7,8,9),c(8)],3))),e_names3)

se_r.3.5<-summary(pool(r3_5), conf.int = TRUE, exponenciate=TRUE)
e_model15<-`names<-`(data.frame(c("Visceral fat", "SUA", "F0-F1", "F1-F2", "F2-F3", "F3-F4"),
                                round(se_r.3.5[c(1,2,6,7,8,9),c(2,3)],3),
                                format.pval(se_r.3.5[c(1,2,6,7,8,9),c(6)],digits=2,eps=0.001),
                                paste0(round(se_r.3.5[c(1,2,6,7,8,9),c(7)],3)," - ",round(se_r.3.5[c(1,2,6,7,8,9),c(8)],3))),e_names3)

e_nmod3<-c("METS-VF","", "", "", "", "", "Interaction", "", "", "", "", "", "WHR", "", "", "", "","", "BMI", "", "", "", "","", "Trunk fat","", "", "", "", "")

e_tabla3 <- cbind("Model"=e_nmod3, rbind(e_model1,e_model2,e_model3,e_model4,e_model5),
                  rbind(e_model11,e_model12,e_model13,e_model14,e_model15),
                  rbind(e_model6,e_model7,e_model8,e_model9,e_model10))

#TAB1<-flextable::align(flextable::flextable(e_tabla3,cwidth=2),align="center",part="all")%>%flextable::autofit()
#flextable::save_as_docx(TAB1,path="~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/ESTANCIA/tab1.docx")

####MEDIACIONES####

datasets<-list(imp1=d1, imp2=d2, imp3=d3, imp4=d4, imp5=d5)
mediators<-c("c_SUA")
covariates<-c("t_LBXSCR+GENDER+RIDAGEYR+BMI")

treatment1<-c("t_LUXSMED","t_LUXSMED","t_LUXSMED","t_LUXSMED","t_LUXSMED")
outcome1<-c("p_LUXCAPM","p_LUXCAPM","p_LUXCAPM","p_LUXCAPM","p_LUXCAPM")

mediators2<-c("c_med")

outcome2<-c("p_LUXSMED","p_LUXSMED","p_LUXSMED","p_LUXSMED","p_LUXSMED")
treatment2<-c("t_LUXCAPM","t_LUXCAPM","t_LUXCAPM","t_LUXCAPM","t_LUXCAPM")


outcome3<-c("p_LUXCAPM","p_LUXCAPM","p_LUXCAPM","p_LUXCAPM","p_LUXCAPM")
treatment3<-c("t_WHR", "t_WHR", "t_WHR", "t_WHR", "t_WHR")
mediators3<-"WHR"

med3_1<-mediations(datasets, treatment3, mediators, outcome3, covariates,
                 families=c("binomial","gaussian"), interaction=FALSE, 
                 conf.level=.95, sims=100) 

summary(amelidiate(med3_1))
summary(med3_1)

#SUA debe ir después de grasa visceral

med4<-mediations(datasets, treatment2, mediators2, outcome2, covariates,
                 families=c("binomial","binomial"), interaction=FALSE, 
                 conf.level=.95, sims=100) 

summary(amelidiate(med4))

#ac
mod1<-lm(t_SUA~t_WHR+t_LBXSCR+GENDER+RIDAGEYR+BMI, data=data1)
summary(mod1)

#ab+c
mod2<-glm(p_LUXCAPM~t_WHR+t_SUA+t_LBXSCR+GENDER+RIDAGEYR+BMI, data=data1, family="binomial")
summary(mod2)

med5<-mediate(mod1, mod2, sims=1000, treat="t_WHR", outcome="p_LUXCAPM", mediator="t_SUA",
              conf.level=.95, boot=T)

s_med5<-summary(med5)


#ac
mod3<-lm(t_SUA~p_LUXCAPM+t_LBXSCR+GENDER+RIDAGEYR+BMI, data=data1)
summary(mod3)

#ab+c
mod4<-lm(t_WHR~p_LUXCAPM+t_SUA+t_LBXSCR+GENDER+RIDAGEYR+BMI, data=data1)
summary(mod4)

med6<-mediate(mod3, mod4, sims=1000, treat="p_LUXCAPM", outcome="t_WHR", mediator="t_SUA",
              conf.level=.95, boot=F)


s_med6<-summary(med6)

#ac
mod5<-lm(t_WHR~t_SUA+t_LBXSCR+GENDER+RIDAGEYR+BMI, data=data1)
summary(mod5)

#ab+c
mod6<-glm(p_LUXCAPM~t_SUA+t_WHR+t_LBXSCR+GENDER+RIDAGEYR+BMI, data=data1, family="binomial")
summary(mod6)

med7<-mediate(mod5, mod6, sims=1000, treat="t_SUA", outcome="p_LUXCAPM", mediator="t_WHR",
              conf.level=.95, boot=T)
s_med7<-summary(med7)


ggplot(aes(x=t_SUA, y=p_LUXCAPM), data=data1)+
  geom_smooth(method="lm")+
  geom_point(alpha=0.5)


m_names<-c("Model", "Treatment", "Mediator", "Outcome", "ACME", "95%IC1", "ADE", "95%IC2","Total Effect", "95%IC3", "% Mediated","95%IC4")

m_model1<-`names<-`(data.frame("1", "WHR", "Uric acid", "Steatosis", round(s_med5$d0, 3), paste0("(",round(s_med5$d0.ci[1], 3),"-",round(s_med5$d0.ci[2], 3),")"), 
                               round(s_med5$z0, 3), paste0("(", round(s_med5$z0.ci[1], 3),"-", round(s_med5$z0.ci[2], 3),")"),
                               round(s_med5$tau.coef, 3), paste0("(",round(s_med5$tau.ci[1], 3),"-", round(s_med5$tau.ci[2], 3),")"),
                               paste0((round(s_med5$n0,4)*100),"%"), paste0("(", (round(s_med5$n0.ci[1], 4)*100),"%-", (round(s_med5$n0.ci[2], 4)*100),"%)")), m_names)


m_model2<-`names<-`(data.frame("2", "Steatosis", "", "WHR", round(s_med6$d0, 3), paste0("(",round(s_med6$d0.ci[1], 3),"-",round(s_med6$d0.ci[2], 3),")"), 
                               round(s_med6$z0, 3), paste0("(", round(s_med6$z0.ci[1], 3),"-", round(s_med6$z0.ci[2], 3),")"),
                               round(s_med6$tau.coef, 3), paste0("(",round(s_med6$tau.ci[1], 3),"-", round(s_med6$tau.ci[2], 3),")"),
                               paste0((round(s_med6$n0,4)*100),"%"), paste0("(", (round(s_med6$n0.ci[1], 4)*100),"%-", (round(s_med6$n0.ci[2], 4)*100),"%)")), m_names)

m_model3<-`names<-`(data.frame("3", "Uric Acid", "WHR", "Steatosis", round(s_med7$d0, 3), paste0("(",round(s_med7$d0.ci[1], 3),"-",round(s_med7$d0.ci[2], 3),")"), 
                               round(s_med7$z0, 3), paste0("(", round(s_med7$z0.ci[1], 3),"-", round(s_med7$z0.ci[2], 3),")"),
                               round(s_med7$tau.coef, 3), paste0("(",round(s_med7$tau.ci[1], 3),"-", round(s_med7$tau.ci[2], 3),")"),
                               paste0((round(s_med7$n0,4)*100),"%"), paste0("(", (round(s_med7$n0.ci[1], 4)*100),"%-", (round(s_med7$n0.ci[2], 4)*100),"%)")), m_names)  

med_table<-rbind(m_model1, m_model2, m_model3)
  
med_tab<-flextable::align(flextable::flextable(med_table,cwidth=2),align="center",part="all")%>%flextable::autofit()
flextable::save_as_docx(med_tab,path="~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/ESTANCIA/med_tab.docx")

dag_sua_whr <- dagify(WHR ~ SUA,  SUA ~ steatosis, SUA ~ WHR, WHR~SUA, WHR~steatosis+SUA,
                      SUA~steatosis+WHR, fibrosis~steatosis, labels = c("WHR" = "Waist-to-\nheight\nratio", 
                                                     "steatosis" = "Steatosis",
                                                     "SUA" = "Uric\nacid",
                                                     "fibrosis" = "Fibrosis"),
                      exposure = c("WHR","SUA"),
                      outcome = "fibrosis")

dag1<-dag_sua_whr%>%tidy_dagitty()%>%
  dag_label(labels=c("WHR" = "Waist-to-\nheight\nratio", "steatosis" = "Steatosis", "SUA" = "Uric\nacid",
                     "fibrosis" = "Fibrosis"))%>%  
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_dag_edges() +
  geom_dag_point()+
  geom_dag_label_repel(aes(label = label, fill = label), col = "white", show.legend = FALSE) +
  theme_dag()

print(dag1)

ggsave(dag1, file="dag_nafld.png", height=15, width=20, units= "cm", dpi=300)


####TABLA DE POBLACIÓN####

data_pop<-data0%>%dplyr::select("LUXSMED","LUXCAPM","DXDTRPF", "DXDTOPF","LBXHSCRP","LBXTC","LBDHDD","LBXGH",
                                "LBXSGTSI","LBXSATSI","LBXSAL","LBXSASSI", "LBXSCR", "LBXSTB","LBXSTR",
                                "LBXSUA", "LBXIN","LBXGLU", "DIQ010","DR2TKCAL","DR2TPROT","DR2TCARB",
                                "DR2TSUGR", "DR2TFIBE","DR2TTFAT","DR2TSFAT", "DR2TMFAT" ,"DR2TPFAT","DR2TCHOL",
                                "RIAGENDR","RIDAGEYR", "BMI","METSIR","WHR","METSVF",
                                "METSVF.C","FLI","ALTAST", "p_LUXCAPM")

names<-c("Fibrosis [kPa]", "Steatosis [dB/m]", "Tunk fat %", "Total fat %", "CPR", "Total cholesterol [mg/dL]", "cHDL [mg/dL]",
         "Glycohemoglobin %", "GGT", "ALT", "Albumin", "AST", "Creatinine [mg/dL]", "Bilirrubin [mg/dL]", "Triglycerides [mg/dL]", 
         "Uric acid [mg/dL]", "Insulin [uU/mL]", "Glucose [mg/dL]", "Diabetes [%]", "Total calories", 
         "Protein [g]", "Carbohidrates [g]", "Sugar [g]", "Fiber [g]", "Total fat [g]", "Saturated fat [g]", 
         "Monounsaturated fat [g]", "Polyunsaturated fat [g]", "Cholesterol [g]", "Male [%]", 
         "Age [years]", "BMI []", "METS-IR","WHR", "METS-VF", "V_obese", "FLI", "AST/ALT", "steatosis")

colnames(data_pop)<-names

data_pop$`Male [%]`<-factor(data_pop$`Male [%]`, levels=c(1,2), labels=c("Male", "Female"))
data_pop$V_obese<-factor(data_pop$V_obese, labels=c("Non-obese", "Obese"))
data_pop$steatosis<-factor(data_pop$steatosis, labels=c("E0/E1", ">E1"))
data_pop$`Diabetes [%]`<-factor(data_pop$`Diabetes [%]`, levels=c(1,2), labels=c("Non-diabetic", "Diabetic"))

data_pop1<-data_pop%>%dplyr::filter(V_obese=="Non-obese")%>%dplyr::select(!("V_obese"))

htmlTable::setHtmlTableTheme(css.rgroup="")

getTable1Stats <- function(x, digits = 2, ...){
  getDescriptionStatsBy(x = x, 
                        by = data_pop1$steatosis,
                        digits = digits,
                        continuous_fn = describeMedian,
                        factor_fn = describeFactors,
                        header_count = TRUE,
                        statistics = list(continuous = getPvalKruskal,
                                          factor = getPvalChiSq,
                                          proportion = getPvalFisher))
}

getTable2Stats <- function(x, digits = 2, ...){
  getDescriptionStatsBy(x = x, 
                        by = data_pop2$steatosis,
                        digits = digits,
                        continuous_fn = describeMedian,
                        factor_fn = describeFactors,
                        header_count = TRUE,
                        statistics = list(continuous = getPvalKruskal,
                                          factor = getPvalChiSq,
                                          proportion = getPvalFisher))
}

tab1<-sapply(FUN=getTable1Stats, data_pop1, na.rm=T)

data_pop2<-data_pop%>%dplyr::filter(V_obese=="Obese")%>%dplyr::select(!("V_obese"))

tab2<-sapply(FUN=getTable2Stats, data_pop2, na.rm=T)

p_tab1<-data.frame(mergeDesc(tab1))

p_tab2<-data.frame(mergeDesc(tab2))

p_tab<-bind_rows(p_tab1, p_tab2)%>%rownames_to_column()

p_tab3<-flextable::align(flextable::flextable(p_tab,cwidth=2),align="center",part="all")%>%flextable::autofit()
flextable::save_as_docx(p_tab3,path="~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/ESTANCIA/ptab3.docx")

getTable3Stats <- function(x, digits = 2, ...){
  getDescriptionStatsBy(x = x, 
                        by = data_pop$V_obese,
                        digits = digits,
                        continuous_fn = describeMedian,
                        factor_fn = describeFactors,
                        header_count = TRUE,
                        statistics = list(continuous = getPvalKruskal,
                                          factor = getPvalChiSq,
                                          proportion = getPvalFisher))
}

tab4<-sapply(FUN=getTable3Stats, data_pop, na.rm=T)

mergeDesc(tab4)

#Tabla de dos columnas: con y sin grasa visceral 
#Dos subcolumnas con y sin esteatosis: alteraciones metabólicas y hepáticas deben ser más graves en con grasa visceral y esteatosis
#Sin grasa visceral y esteatosis fenotipo es más metabólico y no antropométrico


