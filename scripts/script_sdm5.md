Species distribution modelling-Part5
================
Achyut K Banerjee
September 13, 2021

**Premise:** To model distribution of mangrove species *Aegiceras
corniculatum* in current and future climate conditions.

**Data availability:** The associated data is available in
\[C:/Users/Achyut/Desktop/web/sdm5\]

*Install and load package*

``` r
library(biomod2)
library(ggplot2)
library(gridExtra)
library(rgdal)
library(raster)
library(ade4)
```

*Processing data*

``` r
#Loading occurrence data#
occ_data<-read.csv("occ/ac_rar.csv") #rarefied occurrence data#
occ_data<-occ_data[occ_data$Species=="AC",]

#Loading environmental variables#
list.files("current") #30 environmental variables#
bioclim_world<-stack(list.files("current",full.names = T),RAT=FALSE)
plot(bioclim_world)

#Processing environmental variables#
#Reference layer is Koppen-Geiger climate classes having at least one occurrence#
mask_iwp<-shapefile("koppen/iwp_kop_cut_dissolve.shp")
bioclim_iwp<-mask(bioclim_world,mask_iwp)
bioclim_iwp<-crop(bioclim_iwp,mask_iwp)
plot(bioclim_iwp)
#Projection layer is the entire Indo-West Pacific# 
mask_iwp_project<-shapefile("shp/iwp1.shp")
bioclim_iwp_project<-mask(bioclim_world,mask_iwp_project)
bioclim_iwp_project<-crop(bioclim_iwp_project,mask_iwp_project)
plot(bioclim_iwp_project)

#To obtain the identifiers where the species occurs#
points_occ<-data.frame(occ_data[1:808,c("Long","Lat")])
occ_cell_id<-cellFromXY(subset(bioclim_world,1),points_occ)
bioclim_occ<-extract(bioclim_world,points_occ)
bioclim_occ_df<-na.omit(as.data.frame(bioclim_occ))
#Selection of variables to be included in the modeling framework, VIF analysis#
library(usdm)
vif(bioclim_occ_df)
vifstep(bioclim_occ_df)
vifcor(bioclim_occ_df,th=.7)

#Making subset of the selected variables#
bioclim_iwp_sub<-stack(subset(bioclim_iwp,c("bio_2","bio_3","bio_8","bio_13","bio_18","bio_19",
                                                    "elev","ndvi","hf","silt","sand","soc","bdod")))
bioclim_iwp_project_sub<-stack(subset(bioclim_iwp_project,c("bio_2","bio_3","bio_8","bio_13","bio_18","bio_19",
                                                    "elev","ndvi","hf","silt","sand","soc","bdod")))
```

*Processing data for past and future projections*

``` r
#Note: past projections were based on bioclimatic variables only, and so did future projections#
bioclim_world_lgm<-stack(list.files("Worldclim_data/lgm_ccsm",full.names = T),RAT=FALSE)
list.files("Worldclim_data/lgm_ccsm")

bioclim_iwp_lgm<-mask(bioclim_world_lgm,mask_iwp)#mask already loaded#
bioclim_iwp_lgm<-crop(bioclim_iwp_lgm,mask_iwp)
bioclim_iwp_project_lgm<-mask(bioclim_world_lgm,mask_iwp_project)#mask already loaded#
bioclim_iwp_project_lgm<-crop(bioclim_iwp_project_lgm,mask_iwp_project)

bioclim_world_cc45bi50<-stack(list.files("Worldclim_data/cc45bi50",full.names = T),RAT=FALSE)
list.files("Worldclim_data/cc45bi50")

bioclim_iwp_cc45bi50<-mask(bioclim_world_cc45bi50,mask_iwp)#mask already loaded#
bioclim_iwp_cc45bi50<-crop(bioclim_iwp_cc45bi50,mask_iwp)
bioclim_iwp_project_cc45bi50<-mask(bioclim_world_cc45bi50,mask_iwp_project)
bioclim_iwp_project_cc45bi50<-crop(bioclim_iwp_project_cc45bi50,mask_iwp_project)#mask already loaded#

#Note: since other environmental variable data were not available for past climate conditions, VIF analysis was done considering only the bioclimatic variables (current climate condition). The variables thus selected were used to subset the bioclimatic variables for past climate conditions. The same set of selected variables would be applicable for future projection also, since no other variable information except climate was available for future too.#
#For the additional VIF analysis (consider changing working directory)#
occ_data<-read.csv("occ/ac_rar.csv")
occ_data<-occ_data[occ_data$Species=="AC",]
points_occ<-data.frame(occ_data[1:808,c("Long","Lat")])

list.files("current_bio") #19 bioclimatic variables#
bioclim_world_new<-stack(list.files("current_bio",full.names = T),RAT=FALSE)
plot(bioclim_world_new)
occ_cell_id<-cellFromXY(subset(bioclim_world_new,1),points_occ)

bioclim_occ_new<-extract(bioclim_world_new,points_occ)
bioclim_occ_df_new<-na.omit(as.data.frame(bioclim_occ_new))
library(usdm)
vif(bioclim_occ_df_new)
vifstep(bioclim_occ_df_new)
vifcor(bioclim_occ_df_new,th=.7)

bioclim_iwp_sub_lgm<-stack(subset(bioclim_iwp_lgm,c("bio_2","bio_3","bio_8","bio_13","bio_18","bio_19")))
bioclim_iwp_project_sub_lgm<-stack(subset(bioclim_iwp_project_lgm,c("bio_2","bio_3","bio_8","bio_13","bio_18","bio_19")))

bioclim_iwp_sub_cc45bi50<-stack(subset(bioclim_iwp_cc45bi50,c("bio_2","bio_3","bio_8","bio_13","bio_18","bio_19")))
bioclim_iwp_project_sub_cc45bi50<-stack(subset(bioclim_iwp_project_cc45bi50,c("bio_2","bio_3","bio_8","bio_13","bio_18","bio_19")))
```

*Model evaluation*

``` r
#Data splitting into train and test#
data<-read.csv("occ/ac_rar.csv")
sample_size <- nrow(data)
set_proportions <- c(Training = 0.8, Test = 0.2) #20% data kept aside for extrinsic evaluation#
set_frequencies <- diff(floor(sample_size * cumsum(c(0, set_proportions))))
data$set <- sample(rep(names(set_proportions), times = set_frequencies))
data <- split(data, data$set)
write.csv(data$Training,"occ/train_ac.csv")
write.csv(data$Test,"occ/test_ac.csv")

#Loading train data#
occ_data<-read.csv("occ/train_ac.csv") 
occ_data<-occ_data[occ_data$Species=="AC",]

#Data formatting#
occ<-BIOMOD_FormatingData(resp.var = rep(1,nrow(occ_data)),expl.var = bioclim_iwp_sub,
                          resp.xy = occ_data[,c('Long','Lat')],resp.name = "AC",
                          PA.nb.rep = 3,PA.nb.absences =400,PA.strategy = 'random')
occ
plot(occ)

#Biomod modeling#
#Create one directory here (AC)#
ac_opt<-BIOMOD_ModelingOptions(GLM = list(type='quadratic',interaction.level=1),
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

ac_models<-BIOMOD_Modeling(data=occ,models = c("GLM","GBM","RF","GAM","CTA","FDA","SRE","MAXENT.Phillips"),
                           models.options = ac_opt,
                           NbRunEval = 4,DataSplit = 70,VarImport = 3,do.full.models = FALSE,
                           models.eval.meth = c('KAPPA','TSS','ROC'),
                           modeling.id = "ac1")

#Model evaluation and variable importance#
ac_models_scores<-get_evaluations(ac_models)
dim(ac_models_scores)
dimnames(ac_models_scores)
acms_df<-data.frame(ac_models_scores)
write.csv(acms_df,"ac_eval.csv")
models_scores_graph(ac_models,by="models",metrics = c("ROC","TSS"),
                    xlim=c(0.5,1),ylim=c(0.4,1))
models_scores_graph(ac_models,by="cv_run",metrics = c("ROC","TSS"),
                    xlim=c(0.5,1),ylim=c(0.5,1))
models_scores_graph(ac_models,by="data_set",metrics = c("ROC","TSS"),
                    xlim=c(0.5,1),ylim=c(0.5,1))
(ac_models_var_import<-get_variables_importance(ac_models))
apply(ac_models_var_import,c(1,2),mean)

#Current projection#
ac_models_proj_current<-BIOMOD_Projection(modeling.output = ac_models,
                                          new.env = bioclim_iwp_project_sub,
                                          proj.name = "current",selected.models='all',
                                          binary.meth = "TSS",
                                          output_format=".img",do.stack=FALSE)

#Ensemble modeling#
ac_ensemble_models<-BIOMOD_EnsembleModeling(modeling.output = ac_models,
                                            em.by = 'all',eval.metric = 'TSS',
                                            eval.metric.quality.threshold = 0.75,
                                            models.eval.meth = c('KAPPA','TSS','ROC'),
                                            prob.mean = FALSE,
                                            prob.cv = TRUE,committee.averaging = TRUE,
                                            prob.mean.weight = TRUE,VarImport = 0)
#Ensemble model evaluation scores#
(ac_ensemble_models_scores<-get_evaluations(ac_ensemble_models))

#Projection#
ac_ensemble_models_proj_current<-BIOMOD_EnsembleForecasting(EM.output = ac_ensemble_models,
                                                            projection.output = ac_models_proj_current,
                                                            binary.meth = "TSS",output_format=".img",
                                                            do.stack=FALSE)

#External validation of the ensemble modeled outputs#
#Load one ensemble raster#
ensemble<-raster("AC/proj_current/individual_projections/AC_EMcaByTSS_mergedAlgo_mergedRun_mergedData.grd")
#Load the test data (remember to change the working directory)#
occ_test<-read.csv("occ/test_ac.csv")
occ_test<-occ_test[occ_test$Species=="AC",]
points_test<-data.frame(occ_test[1:162,c("Long","Lat")])
library(ecospat)
ecospat.boyce (fit = ensemble , points_test, nclass=0, 
               window.w="default", res=100, PEplot = TRUE)
#Check the BI values for each of the individual projections#
```

*Full model run (with complete occurrence data)*

``` r
#Loading occurrence data#
occ_data<-read.csv("occ/ac_rar.csv")
occ_data<-occ_data[occ_data$Species=="AC",]

#To obtain the identifiers where the species occurs#
points_occ<-data.frame(occ_data[1:808,c("Long","Lat")])
occ_cell_id<-cellFromXY(subset(bioclim_iwp,1),points_occ)
bioclim_occ<-extract(bioclim_iwp,points_occ)
bioclim_occ_df<-na.omit(as.data.frame(bioclim_occ))

#Making subset of the selected variables#
bioclim_iwp_sub<-stack(subset(bioclim_iwp,c("bio_2","bio_3","bio_8","bio_13","bio_18","bio_19",
                                                    "elev","ndvi","hf","silt","sand","soc","bdod")))
bioclim_iwp_project_sub<-stack(subset(bioclim_iwp_project,c("bio_2","bio_3","bio_8","bio_13","bio_18","bio_19",
                                                    "elev","ndvi","hf","silt","sand","soc","bdod")))

#Data formatting#
occ<-BIOMOD_FormatingData(resp.var = rep(1,nrow(occ_data)),expl.var = bioclim_iwp_sub,
                          resp.xy = occ_data[,c('Long','Lat')],resp.name = "AC",
                          PA.nb.rep = 1,PA.nb.absences =450,PA.strategy = 'random') 
#equal number of PA points (check after omitting NAs in the occurrence data frame)#
occ
plot(occ)

#Biomod modeling#
#Create one directory here (AC_full)#
ac_opt<-BIOMOD_ModelingOptions(GLM = list(type='quadratic',interaction.level=1),
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

ac_models<-BIOMOD_Modeling(data=occ,models = c("GLM","GBM","RF","GAM","CTA","FDA","SRE","MAXENT.Phillips"),
                           models.options = ac_opt,
                           NbRunEval = 3,DataSplit = 100,VarImport = 3,do.full.models = TRUE,
                           models.eval.meth = c('KAPPA','TSS','ROC'),
                           modeling.id = "ac1") #Note: No data splitting#

#get model evaluation scores#
ac_models_scores<-get_evaluations(ac_models)
dim(ac_models_scores)
dimnames(ac_models_scores)
acms_df<-data.frame(ac_models_scores)
models_scores_graph(ac_models,by="models",metrics = c("ROC","TSS"),
                    xlim=c(0.5,1),ylim=c(0.4,1))
models_scores_graph(ac_models,by="cv_run",metrics = c("ROC","TSS"),
                    xlim=c(0.5,1),ylim=c(0.5,1))
models_scores_graph(ac_models,by="data_set",metrics = c("ROC","TSS"),
                    xlim=c(0.5,1),ylim=c(0.5,1))

#Measure variable importance#
(ac_models_var_import<-get_variables_importance(ac_models))
apply(ac_models_var_import,c(1,2),mean)

#Current projection#
ac_models_proj_current<-BIOMOD_Projection(modeling.output = ac_models,
                                          new.env = bioclim_iwp_project_sub,
                                          proj.name = "current",selected.models='all',
                                          binary.meth = "TSS",
                                          output_format=".img",do.stack=FALSE)

#Ensemble modeling#
ac_ensemble_models<-BIOMOD_EnsembleModeling(modeling.output = ac_models,
                                            em.by = 'all',eval.metric = 'TSS',
                                            eval.metric.quality.threshold = 0.75,
                                            models.eval.meth = c('KAPPA','TSS','ROC'),
                                            prob.mean = FALSE,
                                            prob.cv = TRUE,committee.averaging = TRUE,
                                            prob.mean.weight = TRUE,VarImport = 0)
#OR (depending upon the model run, choose models#
ac_ensemble_models<-BIOMOD_EnsembleModeling(modeling.output = ac_models,
                                            chosen.models = c('AC_PA1_Full_GLM','AC_PA1_Full_GBM','AC_PA1_Full_RF',
                                                            'AC_PA1_Full_GAM','AC_PA1_Full_FDA','AC_PA1_Full_CTA'),
                                            em.by = 'all',eval.metric = 'TSS',
                                            eval.metric.quality.threshold = 0.75,
                                            models.eval.meth = c('KAPPA','TSS','ROC'),
                                            prob.mean = FALSE,
                                            prob.cv = TRUE,committee.averaging = TRUE,
                                            prob.mean.weight = TRUE,VarImport = 0)

#get ENSEMBLE model evaluation scores#
(ac_ensemble_models_scores<-get_evaluations(ac_ensemble_models))

####Ensemble projection####
ac_ensemble_models_proj_current<-BIOMOD_EnsembleForecasting(EM.output = ac_ensemble_models,
                                                            projection.output = ac_models_proj_current,
                                                            binary.meth = "TSS",output_format=".img",
                                                            do.stack=FALSE)
```

*Past projection (with complete occurrence data)*

``` r
#Projection_LGM-CCSM#
ac_models_proj_lgm<-BIOMOD_Projection(modeling.output = ac_models,
                                      new.env = bioclim_iwp_project_sub_lgm,
                                      proj.name = "iwp_lgm",
                                      binary.meth = "TSS",
                                      output.format=".img",
                                      do.stack=FALSE)
ac_ensemble_models_proj_lgm<-BIOMOD_EnsembleForecasting(EM.output = ac_ensemble_models,
                                                        projection.output = ac_models_proj_lgm,
                                                        binary.meth = "TSS",
                                                        output.format=".grd",
                                                        do.stack=FALSE)
```

*Future projection (with complete occurrence data)*

``` r
#Projection_2050_CCSM_RCP45#
ac_models_proj_cc45bi50<-BIOMOD_Projection(modeling.output = ac_models,
                                           new.env = bioclim_iwp_project_sub_cc45bi50,
                                           proj.name = "iwp_cc45bi50",
                                           binary.meth = "TSS",
                                           output.format=".img",
                                           do.stack=FALSE)
ac_ensemble_models_proj_cc45bi50<-BIOMOD_EnsembleForecasting(EM.output = ac_ensemble_models,
                                                             projection.output = ac_models_proj_cc45bi50,
                                                             binary.meth = "TSS",
                                                             output.format=".grd",
                                                             do.stack=FALSE)
```

**END**

**References**

-   \[Book - Habitat Suitability and Distribution Models - with
    Applications in R, by ANTOINE GUISAN et al.2017\]
