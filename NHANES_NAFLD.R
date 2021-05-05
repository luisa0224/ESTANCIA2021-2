####PAQUETES####
library(haven); library(tidyverse); library(dplyr); library(mice); library(nhanesA); 
library(qpcR); library(VIM); library(bestNormalize); library(factoextra); library(MAMI);
library(glmnetUtils); library(corrplot); library(miceadds); library(factoextra); library(NbClust);
library(FactoMineR); library(ggpubr); library(cluster); library(fpc)

setwd("C:/Users/B4-SUR-4/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/ESTANCIA")

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
    model.imputed1 <- lm((x1), data = data[which(data$.imp == m),])
    BIC1[i] <- BIC(model.imputed1)
  }
  B1<-c(mean(BIC1))
  BIC2<-c()
  m<-max(data$.imp)
  for(i in 1:m){
    model.imputed2 <- lm((x2), data = data[which(data$.imp == m),])
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
data0$METSIR<-(log(2*data0$LBXGLU+data0$LBXTC)*data0$BMI)/log(data0$LBDHDD)
data0$WHR<-(data0$BMXWAIST)/(data0$BMXHT)
data0$METSVF<-(4.466+0.011*(log(data0$METSIR)^3)+3.239*(log(data0$WHR)^3)+0.319*(data0$GENDER)
                +0.594*(log(data0$RIDAGEYR)))
data0$METSVF.C<-categ(data0$METSVF, 7.18)

data0$t_SUA<-orderNorm(data0$LBXSUA)$x.t
data0$t_LUXSMED<-orderNorm(data0$LUXSMED)$x.t
data0$t_METSVF<-orderNorm(data0$METSVF)$x.t
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

data0$c_LUXCAPM<-NULL
data0$c_LUXCAPM[data0$LUXCAPM<=238]<-0
data0$c_LUXCAPM[data0$LUXCAPM>238&data0$LUXCAPM<=260]<-1
data0$c_LUXCAPM[data0$LUXCAPM>260&data0$LUXCAPM<=290]<-2
data0$c_LUXCAPM[data0$LUXCAPM>290]<-3

data0$c_LUXSMED<-NULL
data0$c_LUXSMED[data0$LUXSMED<=2]<-0
data0$c_LUXSMED[data0$LUXSMED>2&data0$LUXSMED<=7.5]<-1
data0$c_LUXSMED[data0$LUXSMED>7.5&data0$LUXSMED<=10]<-2
data0$c_LUXSMED[data0$LUXSMED>10&data0$LUXSMED<=14]<-3
data0$c_LUXSMED[data0$LUXSMED>14]<-4

set.seed(123)
md.pattern(data0, rotate.names=TRUE)

aggr_plot<-aggr(data0, col=c('navyblue','red'), 
                numbers=TRUE, sortVars=TRUE, labels=names(data0), 
                cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
print(aggr_plot)

imp1<-mice(data=data0, m=5, maxit=5, printFlag = FALSE, set.seed=123)
data1<-complete(imp1, "long")

####CORRELACIÓN####
df<-cbind(data0$t_SUA, data0$t_LUXSMED, data0$t_METSVF, data0$LUXCAPM)
res<-cor.mtest(df, conf.level = .95)
colnames(df)<-c("Uric\nAcid", "Fibrosis", "Visceral\nfat", "Steatosis")
cor<-cor(df, method="pearson", use = "complete.obs")
corrplot.mixed(cor,lower.col = "black", number.cex = 1, upper = "circle",
               p.mat=res$p, sig.level=.05)

df1<-cbind(data0$LBXSUA, data0$t_DR2TCARB, data0$t_DR2TKCAL, data0$t_DR2TCHOL, data0$t_DR2TFIBE,
           data0$t_DR2TSUGR, data0$t_DR2TTFAT, data0$t_DR2TSFAT, data0$t_DR2TMFAT,
           data0$t_DR2TPFAT, data0$t_DR2TPROT, data0$t_LUXSMED, data0$LUXCAPM)
res1<-cor.mtest(df1, conf.level = .95)
colnames(df1)<-c("Uric\nacid","Carbs", "Calories", "Cholesterol", "Fiber", "Sugar", "Total\nfat","Sat\nfat", "Mono\nfat", 
                 "Poly\nfat", "Protein", "Fibrosis", "Steatosis")
cor1<-cor(df1, method="spearman", use = "complete.obs")
corrplot.mixed(cor1,lower.col = "black", number.cex = 1, upper = "circle",
               p.mat=res1$p, sig.level=.05)

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

m1 <- mami(imp1, model.family="poisson", outcome=c("c_LUXSMED"), method="LASSO", alpha=0.05, kfold=10, 
           report.exp=FALSE, missing.data = "imputed", var.remove=c("SEQN","BMXWT","BMXHT","BMXWAIST","BMXHIP",
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
                                                                    "c_LUXCAPM", "LUXSMED"))
summary(m1)

m2 <- mami(data3, model.family="poisson", ycol=c("c_LUXCAPM"), method="LASSO", alpha=1, kfold=10, 
           report.exp=FALSE, missing.data = "imputed", var.remove=c("SEQN","BMXWT","BMXHT","BMXWAIST","BMXHIP",
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
                                                                    "c_LUXSMED"))
summary(m2)


mod <- glmnet(x=matrix(c(c(d1$t_DR2TPROT, d1$t_DR2TCARB, d1$t_DR2TSUGR, d1$t_DR2TFIBE, d1$t_DR2TTFAT, 
                                   d1$t_DR2TSFAT, d1$t_DR2TMFAT, d1$t_DR2TPFAT, d1$t_DR2TCHOL)), ncol=9), 
                      y=matrix(c(d1$c_LUXSMED), ncol=1), 
                      family="poisson",alpha = 1, k=5)

mod_cva <- cv.glmnet(x=matrix(c(c(d1$t_DR2TPROT, d1$t_DR2TCARB, d1$t_DR2TSUGR, d1$t_DR2TFIBE, d1$t_DR2TTFAT, 
                                 d1$t_DR2TSFAT, d1$t_DR2TMFAT, d1$t_DR2TPFAT, d1$t_DR2TCHOL)), ncol=9), 
                      y=matrix(c(d1$c_LUXSMED), ncol=1), 
                      family="poisson",alpha = 1, k=5)

coef <- coef(mod_cva, s = "lambda.min")
n3<-data.frame(name = coef@Dimnames[[1]][coef@i + 1], coefficient = coef@x)
n3

mod2 <- glmnet(x=matrix(c(c(d1$t_DR2TPROT, d1$t_DR2TCARB, d1$t_DR2TSUGR, d1$t_DR2TFIBE, d1$t_DR2TTFAT, 
                           d1$t_DR2TSFAT, d1$t_DR2TMFAT, d1$t_DR2TPFAT, d1$t_DR2TCHOL)), ncol=9), 
              y=matrix(c(d1$c_LUXCAPM), ncol=1), 
              family="poisson",alpha = 1, k=5)

mod_cva2 <- cv.glmnet(x=matrix(c(c(d1$t_DR2TPROT, d1$t_DR2TCARB, d1$t_DR2TSUGR, d1$t_DR2TFIBE, d1$t_DR2TTFAT, 
                                  d1$t_DR2TSFAT, d1$t_DR2TMFAT, d1$t_DR2TPFAT, d1$t_DR2TCHOL)), ncol=9), 
                     y=matrix(c(d1$c_LUXCAPM), ncol=1), 
                     family="poisson",alpha = 1, k=5)

coef <- coef(mod_cva2, s = "lambda.min")
n4<-data.frame(name = coef@Dimnames[[1]][coef@i + 1], coefficient = coef@x)
n4


alpha0<-NULL;mod_cva<-NULL;alpha<-NULL
set.seed(123)
for (i in 1:5) {
  ALPHA <- seq(0,1, by=0.01)
  mod_cva[[i]] <- cva.glmnet(x[[i]], y[[i]],family="poisson",alpha = ALPHA, k=5)
  alpha[[i]] <- ALPHA[which.min(sapply(mod_cva[[i]]$modlist, function(mod) min(mod$cvm)))];alpha
}

alpha0<-mean(alpha)

####PCA/CLUSTERING####
data2<-complete(imp1, 1)
df_diet<-scale(data2%>%dplyr::select(DR2TPROT, DR2TCARB, DR2TSUGR, DR2TFIBE, DR2TTFAT, 
                                     DR2TSFAT, DR2TMFAT, DR2TPFAT, DR2TCHOL))

clusters1<-NbClust(data = df_diet, method = "kmeans", distance = "euclidean", min.nc = 2, max.nc = 15, 
                   index = "all", alphaBeale = 0.1)

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

require(rgl);library(rgl)

plot3d(df_diet[,c(2,3,5)], col=data2$clusters1)

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

my_comparisons<-list(c("1", "2"))

plot2.0<-ggplot(data=data2,aes(x=factor(clusters1), y=DR2TTFAT, fill=factor(clusters1)))+
  geom_boxplot()+
  stat_compare_means(label.y=15)+
  stat_compare_means(comparisons = my_comparisons)+
  xlab("Diet cluster")+
  ylab("Total fat")+
  facet_wrap(~GENDER)+
  theme_minimal()

print(plot2.0)

plot2.1<-ggplot(data=data2,aes(x=factor(clusters1), y=DR2TSUGR, fill=factor(clusters1)))+
  geom_boxplot()+
  stat_compare_means(label.y=15)+
  stat_compare_means(comparisons = my_comparisons)+
  xlab("Diet cluster")+
  ylab("Total sugar")+
  facet_wrap(~GENDER)+
  theme_minimal()

print(plot2.1)

plot2.2<-ggplot(data=data2,aes(x=factor(clusters1), y=DR2TCARB, fill=factor(clusters1)))+
  geom_boxplot()+
  stat_compare_means(label.y=15)+
  stat_compare_means(comparisons = my_comparisons)+
  xlab("Diet cluster")+
  ylab("Total carbohydrates")+
  facet_wrap(~GENDER)+
  theme_minimal()

print(plot2.2)


plot2.3<-ggplot(data=data2,aes(x=factor(clusters1), y=DR2TMFAT, fill=factor(clusters1)))+
  geom_boxplot()+
  stat_compare_means(label.y=15)+
  stat_compare_means(comparisons = my_comparisons)+
  xlab("Diet cluster")+
  ylab("Monoinsaturated fat")+
  facet_wrap(~GENDER)+
  theme_minimal()

print(plot2.3)

plot2<-ggarrange(plot2.0, plot2.3, plot2.1, plot2.2, ncol=2, nrow=2, common.legend = T)
print(plot2)

#Clusters-> alto y bajo 

plot3<-ggplot(data=data2,aes(x=factor(clusters1), y=LBXSUA, fill=factor(c_LUXCAPM)))+
  geom_boxplot()+
  stat_compare_means(label.y=15)+
  stat_compare_means(comparisons = my_comparisons)+
  xlab("Diet cluster")+
  ylab("Serum uric acid [mg/dL]")+
  facet_wrap(~METSVF.C)+
  scale_fill_discrete(name="Degree of\nsteatosis" ,labels=c("S0","S1", "S2", "S3"))+
  theme_minimal()

print(plot3)

plot4<-ggplot(data=data2,aes(x=factor(clusters1), y=LBXSUA, fill=factor(c_LUXSMED)))+
  geom_boxplot()+
  stat_compare_means(label.y=15)+
  stat_compare_means(comparisons = my_comparisons)+
  xlab("Diet cluster")+
  ylab("Serum uric acid [mg/dL]")+
  scale_fill_discrete(name="Degree of\nfibrosis" ,labels=c("F0", "F1", "F2", "F3", "F4"))+
  facet_wrap(~GENDER)+
  theme_minimal()

print(plot4)

####MODELOS DE REGRESIÓN####

#Hacer un modelo por sexo?

#Esteatosis

r1<-with(data=imp1, exp=glm(c_LUXCAPM~METSVF+t_SUA+t_DR2TSUGR, family="gaussian"))
summary(pool(r1), conf.int = TRUE, exponentiate = FALSE)

r2<-with(data=imp1, exp=glm(c_LUXCAPM~METSVF+t_SUA, family="gaussian"))
summary(pool(r2), conf.int = TRUE, exponentiate = FALSE)

#Qué componentes dietéticos-> elevación de ácido úrico es lo q causa esteatosis
#1. Asociación a aumento de grasa visceral
#2. Asociación ajustada x ácido úrico 
#3. Esteatosis 
#-> ordinal logistic regression

#Fibrosis

r2<-with(data=imp1, exp=glm(LBXSUA~factor(c_LUXSMED)+GENDER+RIDAGEYR+t_DR2TSFAT+t_DR2TPFAT, 
                            family="poisson"))
summary(pool(r2), conf.int = TRUE, exponentiate = TRUE)



r3<-with(data=imp1, exp=glm(log(WHR)~t_DR2TTFAT+t_SUA+RIDAGEYR, family="gaussian"))
summary(pool(r3), conf.int = TRUE, exponentiate = FALSE)

#Modelos diferentes p hombres y mujeres
#AST/ALT
#Fatty liver index
