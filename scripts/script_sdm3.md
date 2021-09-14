Species distribution modelling-Part3
================
Achyut K Banerjee
September 10, 2021

**Premise:** Testing of model transferability - niche analogy through
**MESS analysis**.This code was used to check the environmental analogy
between native and invasive ranges of *Mikania micrantha*.

**Data availability:** The associated data is available in
\[C:/Users/Achyut/Desktop/web/sdm3\]

*Install and load package*

``` r
library(dismo)
```

*Data preparation*

``` r
#Note: The bioclimatic variables are selected through PCA. The details of PCA are available in "script_sdm1", whereas the bioclimatic variables are introduced in the paper (see Reference list at the end). The variables may change depending on the context. The variables are in ASCII format.#

files <- list.files("current",pattern='asc', full.names=TRUE )
files
predictors<-stack(files) #Making stack#

mask_native<-shapefile("shp/native_range.shp")
mask_invasive<-shapefile("shp/invasive_range.shp")
mask_combined<-shapefile("shp/native_invasive.shp")
plot(mask_invasive)

stack_invasive<-mask(predictors,mask_invasive)
stack_native<-mask(predictors,mask_native)
stack_combined<-mask(predictors,mask_combined)
plot(stack_invasive)
plot(stack_native)
plot(stack_combined)

data_invasive<-read.csv("occ/invasive.csv")
invasive_occ<-data_invasive[data_invasive$Species=="Mikania_invasive",]
points_invasive<-data.frame(invasive_occ[1:43,c("Long","Lat")])

data_native<-read.csv("occ/native.csv")
native_occ<-data_native[data_native$Species=="Mikania_native",]
points_native<-data.frame(native_occ[1:1201,c("Long","Lat")])

data_combined<-read.csv("occ/combined.csv")
combined_occ<-data_combined[data_combined$Species=="Mikania_combined",]
points_combined<-data.frame(combined_occ[1:3426,c("Long","Lat")])

reference_points<-extract(stack_native,points_native)
reference_points1<-extract(stack_invasive,points_invasive)
reference_points2<-extract(stack_combined,points_combined)
```

*MESS analysis*

``` r
mss<-mess(x=stack_invasive,v=reference_points,full=TRUE)
plot(mss$rmess)

plot(mss$rmess,xlim=c(50,150),ylim=c(-20,60))
data_invasive<-extract(mss,points_invasive)
write.csv(data_invasive,file = "invasive_mess.csv")

#Plotting occurrences#
points_i<-data.frame(invasive_occ[1:43,c("Long","Lat")])
points_i <- SpatialPoints(coords = points_i, proj4string = CRS("+proj=longlat +datum=WGS84"))
fun <- function() {
  plot(points_i, add = TRUE, col = "red", pch = 20)
}
plot(mss$rmess,xlim=c(60,180),ylim=c(-40,60), addfun = fun)
```

**END**

**References** \* \[dismo package manual\]
(<https://mran.microsoft.com/snapshot/2014-08-18_0233/web/packages/dismo/dismo.pdf>)
\* \[Paper - Annals in Botany\] (<https://doi.org/10.1093/aob/mcaa044>)
