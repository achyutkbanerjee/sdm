Species distribution modelling-Part2
================
Achyut K Banerjee
September 10, 2021

**Premise1:** Niche characterization - **analysis of niche
conservatism**: *single species, two ranges*. This code was used to
characterize niche of *Butomus umbellatus* across its native and
invasive ranges.

**Data availability:** The associated data is available in
\[C:/Users/Achyut/Desktop/web/sdm2\]

*Install and load package*

``` r
source("niche.overlap.functions.R") #keep this file in the working directory#
source("occ.prep.functions.R") #keep this file in the working directory#
library(ecospat)
library(biomod2)
library(ade4)
library(adehabitatHR)
library(sp)
library(gam)
library(MASS)
library(mvtnorm)
library(gbm)
library(Rcpp)
library(dismo)
library(bindrcpp)
library(rlang)
library(factoextra)
library(FactoMineR)
library(corrplot)
```

*Preparing data*

``` r
#Note: The bioclimatic variables selected from the PCA (see script_sdm1) are considered here. The variables may change depending on the context.#

clim1<-na.exclude(read.delim("native background climate.txt",h=T,sep="\t"))
clim2<-na.exclude(read.delim("invasive background climate.txt",h=T,sep="\t"))
clim12<-rbind(clim1,clim2)
occ.sp.aggr<-na.exclude(read.delim("occ points.txt",h=T,sep="\t"))
occ.sp<-occ.desaggragation(df=occ.sp.aggr,colxy=1:2,min.dist=0,plot=F)
occ.sp1<-na.exclude(sample.sp.globvar(dfsp=occ.sp,colspxy=1:2,colspkept=NULL,dfvar=clim1,colvarxy=1:2,colvar="all",resolution=0.16666))
occ.sp2<-na.exclude(sample.sp.globvar(dfsp=occ.sp,colspxy=1:2,colspkept=NULL,dfvar=clim2,colvarxy=1:2,colvar="all",resolution=0.16666))
row.pa1<-sample.sp.globvar(dfsp=clim1,colspxy=1:2,colspkept=NULL,dfvar=occ.sp1,colvarxy=1:2,colvar=3,resolution=0)
pa<-data.frame((!is.na(row.pa1))*1);names(pa)<-"pa" #create 01 column
pa1<-cbind(clim1,pa)
row.pa2<-sample.sp.globvar(dfsp=clim2,colspxy=1:2,colspkept=NULL,dfvar=occ.sp2,colvarxy=1:2,colvar=3,resolution=0)
pa<-data.frame((!is.na(row.pa2))*1);names(pa)<-"pa" #create 01 column
pa2<-cbind(clim2,pa)
```

*Doing PCA*

``` r
PROJ = F
names(clim12)
Xvar<-c(3:8)
nvar<-length(Xvar)
iterations<-100
R=100

row.w.1.occ<-1-(nrow(occ.sp1)/nrow(rbind(occ.sp1,occ.sp2))) # prevalence of occ1
row.w.2.occ<-1-(nrow(occ.sp2)/nrow(rbind(occ.sp1,occ.sp2))) # prevalence of occ2
row.w.occ<-c(rep(0, nrow(clim1)),rep(0, nrow(clim2)),rep(row.w.1.occ, nrow(occ.sp1)),rep(row.w.2.occ, nrow(occ.sp2)))

row.w.1.env<-1-(nrow(clim1)/nrow(clim12))  # prevalence of clim1
row.w.2.env<-1-(nrow(clim2)/nrow(clim12))  # prevalence of clim2
row.w.env<-c(rep(row.w.1.env, nrow(clim1)),rep(row.w.2.env, nrow(clim2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))

fac<-as.factor(c(rep(1, nrow(clim1)),rep(2, nrow(clim2)),rep(1, nrow(occ.sp1)),rep(2, nrow(occ.sp2))))

data.env.occ<-rbind(clim1,clim2,occ.sp1,occ.sp2)[Xvar]
row.clim1<-1:nrow(clim1)
row.clim2<-(nrow(clim1)+1):(nrow(clim1)+nrow(clim2))
row.clim12<-1:(nrow(clim1)+nrow(clim2))
row.sp1<-(nrow(clim1)+nrow(clim2)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1))
row.sp2<-(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+nrow(occ.sp2))

#fitting occurrences from both ranges#
if(PROJ == F){          
  pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2)
}

scores.clim12<- pca.cal$li[row.clim12,]
scores.clim1<- pca.cal$li[row.clim1,]
scores.clim2<- pca.cal$li[row.clim2,]
scores.sp1<- pca.cal$li[row.sp1,]
scores.sp2<- pca.cal$li[row.sp2,]
```

*Niche visualization*

``` r
z1<-ecospat.grid.clim.dyn(scores.clim12,scores.clim1,scores.sp1,R=100)
z2<-ecospat.grid.clim.dyn(scores.clim12,scores.clim2,scores.sp2,R=100)

# Test of niche equivalency and similarity according to Warren et al. 2008#
a<-ecospat.niche.equivalency.test(z1,z2,rep=1)
b<-ecospat.niche.similarity.test(z1,z2,rep=100)
b2<-ecospat.niche.similarity.test(z2,z1,rep=100)
a
b
b2

#Plotting#
x11(); layout(matrix(c(1,1,2,2,1,1,2,2,3,3,4,5,3,3,6,7), 4, 4, byrow = TRUE)) #Required only to display all plots together#

ecospat.plot.niche(z1,title="PCA-env - native niche",name.axis1="PC1",name.axis2="PC2")
ecospat.plot.niche(z2,title="PCA-env - invasive niche",name.axis1="PC1",name.axis2="PC2")
ecospat.plot.contrib(pca.cal$co,pca.cal$eig)

plot.new(); text(0.5,0.5,paste("niche overlap:","\n","D=",round(as.numeric(ecospat.niche.overlap(z1,z2,cor=T)[1]),3)))

plot.overlap.test(a,"D","Equivalency")
plot.overlap.test(b,"D","Similarity 2->1")
plot.overlap.test(b2,"D","Similarity 1->2")
ecospat.plot.niche.dyn(z1,z2,title="niche occupancy between native and invasive region",quant=0.75)
ecospat.shift.centroids(scores.sp1,scores.sp2,scores.clim1,scores.clim2,"red")
R=10
test<-ecospat.niche.dyn.index(z1,z2,intersection=NA)
test
```

**END**

**Premise2:** Niche characterization - **analysis of niche
conservatism**: *two species, one range*. This code was used to
characterize niches of two mangroves Aegiceras corniculatum and Acanthus
ilicifolius in the Indo-West Pacific.

**Data availability:** The associated data is available in
\[C:/Users/Achyut/Desktop/web/sdm2\]

*Preparing data*

``` r
#Load background data#
clim_all<-read.csv("clim_all.csv") #Environmental variables (bioclimatic and demographic) in the IWP#
clim_all<-na.omit(as.data.frame(clim_all)) #Convert it to data frame#
Xvar <- c("bdod", "bio_1","bio_2", "bio_3","bio_4", "bio_5", "bio_6","bio_7","bio_8","bio_9",
          "bio_10","bio_11","bio_12","bio_13","bio_14","bio_15","bio_16","bio_17","bio_18","bio_19",
          "cec","clay","elev","hf","nitrogen","ph","sand","silt","soc") #Environmental variables#
nvar <- length(Xvar)
#Load species data#
ac<-read.csv("ac.csv") #Environmental variables (same set as the background) for species1#
ai<-read.csv("ai.csv") #Environmental variables (same set as the background) for species2#
```

*Doing PCA*

``` r
pca.cal <- dudi.pca(clim_all[, Xvar], center = TRUE,
                    scale = TRUE, scannf = FALSE, nf = 2)
scores.clim <- pca.cal$li
ac.scores <- suprow(pca.cal, ac[, Xvar])$li
ai.scores <- suprow(pca.cal, ai[, Xvar])$li
```

*Plotting niche*

``` r
tiff("ac.tiff", units="in", width=5, height=5, res=300,compression = "none") #AC/AI#
p<-plot(scores.clim, pch = 16, asp = 1,
     col = adjustcolor(1, alpha.f = 0.2), cex = 1,
     xlab = "PC1", ylab = "PC2")
p
p+points(ac.scores, pch = 1, col = "yellow", cex = 1) #AC/AI#
dev.off()
```

*Niche visualization*

``` r
z1 <- ecospat.grid.clim.dyn(scores.clim, scores.clim,
                            ac.scores, R = 100)
z2 <- ecospat.grid.clim.dyn(scores.clim, scores.clim,
                            ai.scores, R = 100)
plot(z1$z.uncor, legend = FALSE, axes = TRUE,
     main = "AC")
plot(z2$z.uncor, legend = FALSE, axes = TRUE,
     main = "AI")

#Niche similarity and equivalency tests#
a<-ecospat.niche.overlap (z1, z2, cor=TRUE) #estimate the D and I values#
a
b<-ecospat.niche.similarity.test (z1, z2, rep=100, alternative = "greater",
                               rand.type = 1, ncores= 1) #significance testing of D and I values#
b
c<-ecospat.niche.equivalency.test (z1, z2, rep=10, 
                                   alternative="greater", ncores = 1) #significance testing of D and I values#
c

ecospat.plot.overlap.test (a, "D", "Niche overlap")
ecospat.plot.overlap.test (b, "D", "Niche similarity")
```

**END**

**References**

-   \[ecospat package manual\]
    (<https://cran.r-project.org/web/packages/ecospat/ecospat.pdf>)
-   \[website\] (<https://plantarum.ca/notebooks/ecospat/>)
-   \[Paper - Hydrobiologia\]
    (<https://doi.org/10.1007/s10750-020-04205-1>)
