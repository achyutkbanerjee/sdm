Species distribution modelling-Part1
================
Achyut K Banerjee
September 10, 2021

**Premise:** To download bioclimatic variables, process occurrence and
climate data and **select variables through PCA or VIF analysis**. This
code was used to map potential distribution of an invasive plant
species, *Butomus umbellatus*.

**Data availability:** The associated data is available in
\[C:/Users/Achyut/Desktop/web/sdm1\]

*Install and load package*

``` r
library(biomod2)
library(ggplot2)
library(gridExtra)
library(rgdal)
library(raster)
library(ade4)
```

*Download bioclimatic variables*

``` r
#set working directory#
dir.create("Worldclim_data",showWarnings = F) #Create a directory for the data#

#Note: Data downloaded here are of 10 arc minute resolution. Change to 5/2.5 arc minute or 30 arc second where necessary#

#Current climate condition#
download.file(url="http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio_10m_esri.zip",destfile = "Worldclim_data/current_bioclim_10min.zip",method = "auto")
#Future climate condition#
download.file("http://biogeo.ucdavis.edu/data/climate/cmip5/10m/bc45bi50.zip",destfile="Worldclim_data/2050_BC_45_bioclim_10min.zip",method="auto") #GCM->BCC-CSM1-1,year->2050,RCP->4.5#
download.file("http://biogeo.ucdavis.edu/data/climate/cmip5/10m/bc45bi70.zip",destfile="Worldclim_data/2070_BC_45_bioclim_10min.zip",method="auto") #GCM->BCC-CSM1-1,year->2070,RCP->4.5#
#Unzip files#
unzip(zipfile = "Worldclim_data/current_bioclim_10min.zip",exdir = "Worldclim_data/current",overwrite = T)
list.files("Worldclim_data/current/bio/") #To visualize the unzipped files#
unzip(zipfile = "Worldclim_data/2050_BC_45_bioclim_10min.zip",exdir = "Worldclim_data/2050/BC_45",overwrite = T)
list.files("Worldclim_data/2050/BC_45")
unzip(zipfile = "Worldclim_data/2070_BC_45_bioclim_10min.zip",exdir = "Worldclim_data/2070/BC_45",overwrite = T)
list.files("Worldclim_data/2070/BC_45")
```

*Data preparation*

``` r
#Making stack of the bioclimatic variables#
bioclim_world<-stack(list.files("Worldclim_data/current/bio/",pattern="bio_",full.names = T),RAT=FALSE)

#Creating masks for native and invasive ranges#
mask_native<-shapefile("shapefiles/native_coun.shp")
mask_invasive<-shapefile("shapefiles/invasive_coun.shp")
bioclim_native<-mask(bioclim_world,mask_native)
bioclim_invasive<-mask(bioclim_world,mask_invasive)

#Loading occurrence data#
data<-read.csv("occ_data/native_occ.csv")#spatially rarefied occurrences#
data1<-read.csv("occ_data/invasive_occ.csv")
native_occ<-data[data$Species=="rush_native",]
invasive_occ<-data1[data1$Species=="rush_invasive",]
```

*Doing PCA*

``` r
#To obtain the identifiers where the species occurs#
points_native<-data.frame(native_occ[1:3641,c("Long","Lat")])
native_cell_id<-cellFromXY(subset(bioclim_native,1),points_native)
points_invasive<-data.frame(invasive_occ[1:1753,c("Long","Lat")])
invasive_cell_id<-cellFromXY(subset(bioclim_invasive,1),points_invasive)
bioclim_nativeocc<-extract(bioclim_native,points_native)
bioclim_nativeocc_df<-na.omit(as.data.frame(bioclim_nativeocc))
bioclim_invasive_occ<-extract(bioclim_world,points_invasive)
bioclim_invasiveocc_df<-na.omit(as.data.frame(bioclim_invasive_occ))

#Converting raster to data frame#
bioclim_native_df<-na.omit(as.data.frame(bioclim_native))
head(bioclim_native_df)
write.csv(bioclim_native_df,file = "native_back.csv")
bioclim_invasive_df<-na.omit(as.data.frame(bioclim_invasive))
head(bioclim_invasive_df)
write.csv(bioclim_invasive_df,file = "invasive_back.csv")

#Doing PCA and plotting#
pca_ZA<-dudi.pca(bioclim_native_df,scannf = F,nf=2)
plot(pca_ZA$li[,1:2])
par(mfrow=c(1,2))
s.class(pca_ZA$li[,1:2],fac=factor(rownames(bioclim_native_df)%in%native_cell_id,levels=c("FALSE","TRUE"),labels=c("back","occ")),col=c("red","blue"),csta=0,cellipse=4,cpoint=.3,pch=16)
mtext("(a)",side=3,line=3,adj = 0)
s.corcircle(pca_ZA$co,clabel = .5)
mtext("(b)",side = 3,line = 3,adj = 0)

#Subset of selected variables#
bioclim_native_sub<-stack(subset(bioclim_native,c("bio_1","bio_5","bio_7","bio_12","bio_14","bio_15")))
bioclim_invasive_sub<-stack(subset(bioclim_invasive,c("bio_1","bio_5","bio_7","bio_12","bio_14","bio_15")))
```

*Doing PCA - alternative approach*

``` r
#Install and load libraries#
library(bindrcpp)
library(factoextra)
library(rlang)
library(FactoMineR)
library(corrplot)

#Doing PCA and plotting#
res.pca<- PCA(bioclim_native_df, scale.unit = TRUE, ncp = 2, graph = TRUE)

eig.val <- get_eigenvalue(res.pca)
eig.val
var <- get_pca_var(res.pca)
var
print(var)
head(var$coord)
head(var$cos2)
head(var$contrib,19)
head(var$coord, 19)
contrib<-data.frame(var$contrib)
write.csv(contrib,file = "contrib.csv")
coord<-data.frame(var$coord)
write.csv(coord,file = "coord.csv")
fviz_pca_var(res.pca, col.var = "black")

corrplot(var$coord, is.corr=FALSE)  
fviz_contrib(res.pca, choice = "var", axes = 1, top = 3)
fviz_contrib(res.pca, choice = "var", axes = 2, top = 3)
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)
res.desc <- dimdesc(res.pca, axes = c(1,2), proba = 0.05)
res.desc$Dim.1
res.desc$Dim.2
```

*Variable selection through VIF analysis*

``` r
#Install and load library#
library(usdm)
vif(bioclim_native_df)
vifstep(bioclim_native_df)
vifcor(bioclim_native_df,th=.7) #Threshold can be modified#

#VIF on occurrence data#
bioclim_nativeocc<-extract(bioclim_native,points_native)
bioclim_nativeocc_df<-na.omit(as.data.frame(bioclim_nativeocc))
vif(bioclim_nativeocc_df)
vifstep(bioclim_nativeocc_df)
vifcor(bioclim_nativeocc_df,th=.7)

#Subset of selected variables#
bioclim_native_sub<-stack(subset(bioclim_native,c("bio_2","bio_7","bio_8",
                                                  "bio_9","bio_13","bio_15","bio_18","bio_19")))
bioclim_invasive_sub<-stack(subset(bioclim_invasive,c("bio_2","bio_7","bio_8",
                                                  "bio_9","bio_13","bio_15","bio_18","bio_19")))
```

**END**

**References**

-   \[Paper - Hydrobiologia\]
    (<https://doi.org/10.1007/s10750-020-04205-1>)
