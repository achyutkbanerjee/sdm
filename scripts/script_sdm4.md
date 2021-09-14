Species distribution modelling-Part4
================
Achyut K Banerjee
September 13, 2021

**Premise:** To model distribution of *Butomus umbellatus*, an aquatic
invasive plant species, in current and future climate conditions.

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

*Data formatting*

``` r
#For loading occurrence data, refer to script_sdm1#
#Making subsets of bioclimatic variables for both native and invasive ranges has also been explained in script_sdm1#

native_data<-BIOMOD_FormatingData(resp.var = rep(1,nrow(native_occ)),
                                  expl.var = bioclim_native_sub,resp.xy = native_occ[,c('Long','Lat')],
                                  resp.name = "rush_native",PA.nb.rep = 3,PA.nb.absences = 500,
                                  PA.strategy = 'random')#PA points should be equal to number of occurrences#
native_data
plot(native_data)

invasive_data<-BIOMOD_FormatingData(resp.var = rep(1,nrow(invasive_occ)),
                                    expl.var = bioclim_invasive_sub,resp.xy = invasive_occ[,c('Long','Lat')],
                                    resp.name = "rush_invasive")
invasive_data
plot(invasive_data)
#native and invasive data can be jointly processed in case we want projection based on combined occurrences of both native and invasive ranges# 
```

*Biomod modeling*

``` r
#Create a working directory here. Put the maxent.jar file in there. The file can be downloaded from (http://www.cs.princeton.edu/~schapire/maxent)#
setwd("C:/Users/Achyut/Desktop/web/sdm/sdm4/model.run")
native_opt<-BIOMOD_ModelingOptions(GLM = list(type='quadratic',interaction.level=1),
                                   GBM=list(n.trees=1000),GAM=list(algo='GAM_mgcv'),
                                   MAXENT.Phillips = list(path_to_maxent.jar=getwd(),
                                                          maximumiterations=200,
                                                          visible=FALSE,
                                                          linear=TRUE,
                                                          quadratic=TRUE,
                                                          product=TRUE,
                                                          threhold=TRUE,
                                                          hinge=TRUE,
                                                          lq2lqptthreshold = 80,
                                                          l2lqthreshold = 10,
                                                          hingethreshold = 15,
                                                          beta_threshold = -1,
                                                          beta_categorical = -1,
                                                          beta_lqp = -1,
                                                          beta_hinge = -1,
                                                          defaultprevalence = 0.5))

#data splitting into 70% training and 30% testing (intrinsic evaluation)#
native_models<-BIOMOD_Modeling(data=native_data,models = c("GLM","GBM","RF","GAM","MAXENT.Phillips"),
                               models.options = native_opt,
                               NbRunEval = 4,DataSplit = 70,VarImport = 3,do.full.models = F,
                               modeling.id = "fr1") 
```

*Getting model evaluation scores*

``` r
native_models_scores<-get_evaluations(native_models)
dim(native_models_scores)
dimnames(native_models_scores)
nms_df<-data.frame(native_models_scores)

models_scores_graph(native_models,by="models",metrics = c("ROC","TSS"),
                    xlim=c(0.5,1),ylim=c(0.5,1))
models_scores_graph(native_models,by="cv_run",metrics = c("ROC","TSS"),
                    xlim=c(0.5,1),ylim=c(0.5,1))
models_scores_graph(native_models,by="data_set",metrics = c("ROC","TSS"),
                    xlim=c(0.5,1),ylim=c(0.5,1))
```

*Measure variable importance*

``` r
(native_models_var_import<-get_variables_importance(native_models))
apply(native_models_var_import,c(1,2),mean)
```

*Current projection*

``` r
native_models_proj_current<-BIOMOD_Projection(modeling.output = native_models,
                                              new.env = bioclim_invasive_sub,
                                              proj.name = "current",selected.models='all',
                                                binary.meth = "TSS",
                                              output_format=".img",do.stack=FALSE)
```

*External validation*

``` r
#External validation using invasive range occurrence data#
library(ecospat)
raster<-raster("rush.native/proj_current/individual_projections/proj_current_rush.native_PA1_RUN1_RF.grd")
ecospat.boyce (raster,
               points_invasive, 
               nclass=0, window.w="default", res=100, PEplot = TRUE)
#Check the BI values for each of the individual projections. Select the best model based on the average BI values. Develop the models with the selected model algorithm(s)#
```

*Selected model(s)*

``` r
setwd("C:/Users/Achyut/Desktop/web/sdm/sdm4/model.run.select")
native_opt1<-BIOMOD_ModelingOptions()
native_models1<-BIOMOD_Modeling(data=native_data,models = c("RF"),models.options = native_opt1,
                               NbRunEval = 4,DataSplit = 70,VarImport = 3,do.full.models = F,
                               modeling.id = "nk")
#For final model run, use DataSplit=100; i.e., allow the model to run with full occurrence dataset#
native_models_scores1<-get_evaluations(native_models1)
dim(native_models_scores1)
dimnames(native_models_scores1)
nms_df<-data.frame(native_models_scores1)
native_models_proj_current1<-BIOMOD_Projection(modeling.output = native_models1,
                                              new.env = bioclim_invasive_sub,
                                              proj.name = "current",
                                              binary.meth = "TSS",
                                              output_format=".img",do.stack=FALSE)
```

*Ensemble modeling*

``` r
get_built_models(native_models1)
native_ensemble_models<-BIOMOD_EnsembleModeling(modeling.output = native_models1,
                                                chosen.models ='all',
                                                em.by = 'all',eval.metric = 'all',
                                                models.eval.meth = c('KAPPA','TSS','ROC'),
                                                prob.mean = TRUE,
                                                prob.cv = TRUE,committee.averaging = TRUE,
                                                prob.mean.weight = TRUE,VarImport = 0)
#Getting ensemble model evaluation score#
(native_ensemble_models_scores<-get_evaluations(native_ensemble_models))
```

*Ensemble projection*

``` r
native_ensemble_models_proj_current<-BIOMOD_EnsembleForecasting(EM.output = native_ensemble_models,
                                                                projection.output = native_models_proj_current1,
                                                                binary.meth = "TSS",output_format=".grd",
                                                                do.stack=FALSE)
```

*Boyce index for extrinsic evaluation*

``` r
library(ecospat)
raster<-raster("rush.native/proj_current/individual_projections/rush.native_EMwmeanByROC_mergedAlgo_mergedRun_mergedData.grd")
ecospat.boyce (raster,
               points_invasive, 
               nclass=0, window.w="default", res=100, PEplot = TRUE)
```

*Future projection*

``` r
#load 2050 bioclim variables#
bioclim_world_2050_BC45<-stack(c(bio_1="Worldclim_data/2050/BC_45/bc45bi501.tif",
                                 bio_5="Worldclim_data/2050/BC_45/bc45bi505.tif",
                                 bio_7="Worldclim_data/2050/BC_45/bc45bi507.tif",
                                 bio_12="Worldclim_data/2050/BC_45/bc45bi5012.tif",
                                 bio_14="Worldclim_data/2050/BC_45/bc45bi5014.tif",
                                 bio_15="Worldclim_data/2050/BC_45/bc45bi5015.tif"))
#load 2070 bioclim variables#
bioclim_world_2070_BC45<-stack(c(bio_1="Worldclim_data/2070/BC_45/bc45bi701.tif",
                                 bio_5="Worldclim_data/2070/BC_45/bc45bi705.tif",
                                 bio_7="Worldclim_data/2070/BC_45/bc45bi707.tif",
                                 bio_12="Worldclim_data/2070/BC_45/bc45bi7012.tif",
                                 bio_14="Worldclim_data/2070/BC_45/bc45bi7014.tif",
                                 bio_15="Worldclim_data/2070/BC_45/bc45bi7015.tif"))
#Crop#
bioclim_invasive_2050_BC45<-crop(bioclim_world_2050_BC45,mask_invasive)
bioclim_invasive_2050_BC45<-mask(bioclim_invasive_2050_BC45,mask_invasive)
bioclim_invasive_2050_BC45<-stack(bioclim_invasive_2050_BC45)
plot(bioclim_invasive_2050_BC45)

bioclim_invasive_2070_BC45<-crop(bioclim_world_2070_BC45,mask_invasive)
bioclim_invasive_2070_BC45<-mask(bioclim_invasive_2070_BC45,mask_invasive)
bioclim_invasive_2070_BC45<-stack(bioclim_invasive_2070_BC45)
plot(bioclim_invasive_2070_BC45)

#Projection_2050#
rush_models_proj_2050_BC45<-BIOMOD_Projection(modeling.output = native_models1,
                                              new.env = bioclim_invasive_2050_BC45,
                                              proj.name = "2050_BC45",
                                              binary.meth = "TSS",
                                              output.format=".grd",
                                              do.stack=FALSE)
rush_ensemble_models_proj_2050_BC45<-BIOMOD_EnsembleForecasting(EM.output = native_ensemble_models,
                                                                projection.output = rush_models_proj_2050_BC45,
                                                                binary.meth = "TSS",
                                                                output.format=".grd",
                                                                do.stack=FALSE)

#Projection_2070#
rush_models_proj_2070_BC45<-BIOMOD_Projection(modeling.output = native_models1,
                                              new.env = bioclim_invasive_2070_BC45,
                                              proj.name = "2070_BC45",
                                              binary.meth = "TSS",
                                              output.format=".grd",
                                              do.stack=FALSE)
rush_ensemble_models_proj_2070_BC45<-BIOMOD_EnsembleForecasting(EM.output = native_ensemble_models,
                                                                projection.output = rush_models_proj_2070_BC45,
                                                                binary.meth = "TSS",
                                                                output.format=".grd",
                                                                do.stack=FALSE)
```

**END**

**References**

-   \[Paper - Hydrobiologia\]
    (<https://doi.org/10.1007/s10750-020-04205-1>)
-   \[Book - Habitat Suitability and Distribution Models - with
    Applications in R, by ANTOINE GUISAN et al.2017\]
