#-----------------------------------------------------------------------------
#################MrIML analysis#################
#---------------------------------------------------------------------------
#The aim of this analysis is to select only the most important environmental variables for the MRF pipeline.

#This aslo provides some idea of the importance of environmental predictors in shaping microbial communities/

#can't compare directly back to mrf/crf findings    
library(mrIML)
#other package we need:
library(vip); library(tidymodels); library(randomForest);  library(caret); library(gbm);
library(tidyverse);library(parallel); library(doParallel); library(themis); library(viridis); library(janitor); library(hrbrthemes); library(xgboost); library(vegan);library(flashlight);
library(ggrepel); library(iml); library(plyr);
library(future.apply)

#load dat from previous step

tick_microASV<- readRDS('asvFilter_230asv')
env_data_withHost<- readRDS('env_data_withHost')

randomEff<- readRDS('randonEff')
#Host and enviro data
X <-  cbind(env_data_withHost, year=randomEff$year) # site couldnt be included  - too many levels #no need for plate. No [alte contamination was detected]

Xdummy <- fastDummies::dummy_cols(X, remove_first_dummy = TRUE, remove_selected_columns = TRUE) #interactions doesnt like factors
#ASV data
Y <-  tick_microASV

cl <- parallel::makeCluster(5) #increase for bigger runs
plan(cluster, workers=cl)

model_rf <- 
  rand_forest(trees = 100, mode = "classification", mtry = tune(), min_n = tune()) %>% #100 trees are set for brevity. Aim to start with 1000
  set_engine("randomForest")

#not ranger has a problem with importance measures

model_glm <- #model used to generate yhat
  logistic_reg() %>%
  set_engine("glm") %>%
  set_mode("classification") 

model_xgb <-  boost_tree(
  trees = 1000, 
  tree_depth = tune(), min_n = tune(), 
  loss_reduction = tune(),                     ## first three: model complexity
  sample_size = tune(), mtry = tune(),         ## randomness
  learn_rate = tune(),                         ## step size
) %>% 
  set_engine("xgboost") %>% 
  set_mode("classification")

#run the models

yhats_rf <- mrIMLpredicts(X=X,Y=Y, Model=model_rf, balance_data='no', mode='classification',
                          tune_grid_size=5, seed = 123 ) #for reproducibility
# save the model
saveRDS(yhats_rf, file='rf_model_updated')

#need to dummify sex for XGB and GLM



yhats_glm <- mrIMLpredicts(X=Xdummy,Y=Y, Model=model_glm, balance_data='no', mode='classification',   seed = sample.int(1e8, 1) ) ## in MrTidymodels. Balanced data= up upsamples and down downsampled to create a balanced set

# save the model
save(yhats_glm , file='glm_model_updated')

yhats_glm <- load('yhats_glm')

yhats_xgb <- mrIMLpredicts(X=Xdummy,Y=Y, Model=model_xgb, balance_data='no', mode='classification',   seed = sample.int(1e8, 1) ) ## in MrTidymodels. Balanced data= up upsamples and down downsampled to create a balanced set

# save the model
saveRDS(yhats_xgb , file='xgb_model_updated')
yhats_xgb <- readRDS(file = "xgb_model_updated'")

ModelPerf<- mrIMLperformance2(yhats_rf, Model=model_rf, Y=Y, mode='classification') #updated performance function
ModelPerf[2]

a <- ModelPerf[[1]]

yhats_rf <- readRDS('rf_model_updated')
#make it comparable to mrf
ModelPerf[[1]][is.na(ModelPerf[[1]])] <- 0#NAs = 0 here
ModelPerf[[1]][ModelPerf[[1]]==NaN] <- 0#NAs = 0 here (no true negatives detected for some taxa)

perf_data <- data.frame(ppv=as.numeric(ModelPerf[[1]]$ppv), spec = as.numeric(ModelPerf[[1]]$specificity), sens = as.numeric(ModelPerf[[1]]$sensitivity)) #make numeric

mean(perf_data$ppv)
mean(perf_data$sens)
mean(perf_data$spec)

ModelPerf_glm<- mrIMLperformance2(yhats_glm, Model=model_glm, Y=Y, mode='classification') #updated perfromance function

b <- ModelPerf_glm[[1]]

ModelPerf_glm[[1]][is.na(ModelPerf_glm[[1]])] <- 0#NAs = 0 here
ModelPerf_glm[[1]][ModelPerf_glm[[1]]==NaN] <- 0#NAs = 0 here (no true negatives detected for some taxa)

perf_data_glm <- data.frame(ppv=as.numeric(ModelPerf_glm[[1]]$ppv), spec = as.numeric(ModelPerf_glm[[1]]$specificity), sens = as.numeric(ModelPerf_glm[[1]]$sensitivity)) #make numeric

mean(perf_data_glm$ppv)
mean(perf_data_glm$sens)
mean(perf_data_glm$spec)

ModelPerf_xgb<- mrIMLperformance2(yhats_xgb, Model=model_xgb, Y=Y, mode='classification') #updated performance function

#make it comparable to mrf
ModelPerf_xgb[[1]][is.na(ModelPerf_xgb[[1]])] <- 0#NAs = 0 here
ModelPerf_xgb[[1]][ModelPerf_xgb[[1]]==NaN] <- 0#NAs = 0 here (no true negatives detected for some taxa)

perf_data_xgb <- data.frame(ppv=as.numeric(ModelPerf_xgb[[1]]$ppv), spec = as.numeric(ModelPerf_xgb[[1]]$specificity), sens = as.numeric(ModelPerf_xgb[[1]]$sensitivity)) #make numeric
mean(perf_data_xgb $ppv)
mean(perf_data_xgb $sens)
mean(perf_data_xgb $spec)


#plot custom boxplot

p1 <- ggplot(perf_data, aes(x='ppv', y= ppv))
p1 + geom_point(alpha=0.5, aes(middle = mean(ppv)))+
  stat_summary(fun= "mean", geom= "point", colour= "blue", position = position_dodge(width = 0.75))+
  theme_bw()

p2 <- ggplot(perf_data, aes(x='sensitivity', y= sens))
p2 + geom_boxplot(alpha=0.5, aes(middle = mean(sens)))+
  stat_summary(fun= "mean", geom= "point", colour= "blue", position = position_dodge(width = 0.75))+
  theme_bw()

glimpse(ppv_data)
#wroking on routine that we will refit models with different seeds

#for GLM - not a stocjastic algorithm
VI <- mrVip(yhats_glm, X) 
plot_vi(VI, Xdummy, Y)


#get summary of individual ASVs
#standard perf (no ppv)

ModelPerf_a<- mrIMLperformance(yhats_rf, Model=model_rf, Y=Y, mode='classification') #updated performance function

#plot it
VI <- mrVip(yhats_rf, X) 
cluster_pca-scores <- interpret_Mrvi(VI=VI,  X=X,Y=Y, modelPerf=ModelPerf_a,  cutoff= 0.7, mode='classification') 


#draft new MrIML function
VI_u <-mrVip_uncertainty(ity = 5, X=X, Y=Y, Model=model_rf, mode='classification', tune_grid_size=5, 
                        seed = sample.int(1e8, 1))

colmeans_function <- function( i ){
        colMeans(matrix(i,nrow=length(Y)))
      }

      vimp_summary <-  VI_u %>%  apply(2, colmeans_function) %>%
        as.data.frame()

      data_vi_global <-  vimp_summary  %>%
        gather(names(X), key = Xvar, value = importance)

#boxplot
  data_vi_global  %>%
  ggplot(aes(x=reorder(Xvar,importance), y=importance)) +
  geom_boxplot(fill = 'lightgrey')+
  coord_flip()+
  theme_bw()+
  labs(x="Features", y="Importance")
 
 #individual responses (example)
 data_vi<- VI_u 
 ity <-  5
 
 namesY <- rep(names(Y), ity)

 df_namesUpdated  <-  cbind (Ynames=namesY, data_vi )

 Y39 <-  df_namesUpdated %>% filter (Ynames == 'Borreliella')

 data_vi39 <-  Y39  %>%
   gather(names(X), key = Xvar, value = importance)

 data_vi39 %>%
   ggplot(aes(x=reorder(Xvar,importance), y=importance)) +
   geom_boxplot(fill = 'lightgrey')+
   coord_flip()+
   theme_bw()+
   labs(x="Features", y="Importance (Borreliella)")


#erorr in second plot

mfl <- mrFlashlight(yhats_rf, X=X, Y=Y, response = "multi", mode='classification') #needs a model perf argument 

profileData_ale <- light_profile(mfl, v = "prev_year_avg_tmean", type = "ale")
p1 <- mrProfileplot(profileData_ale , sdthresh =0.1) #plotting categorical predictors not working

profileData_ale <- light_profile(mfl, v = "max_snow_depth", type = "ale")
p2 <- mrProfileplot(profileData_ale , sdthresh =0.08)


profileData_ale <- light_profile(mfl, v = "prev_5years_avg_tmean", type = "ale")
p3 <- mrProfileplot(profileData_ale , sdthresh =0.05) #plotting categorical predictors not working

profileData_ale <- light_profile(mfl, v = "ndvi_point", type = "ale")
p3 <- mrProfileplot(profileData_ale , sdthresh =0.08) #sdthresh removes responses from the first plot that do not vary with the feature

profileData_ale <- light_profile(mfl, v = "max_snow_depth", type = "ale")
p4 <- mrProfileplot(profileData_ale , sdthresh =0.08)

profileData_ale <- light_profile(mfl, v = "sex", type = "ale")
p5 <- mrProfileplot(profileData_ale , sdthresh =0.05)#plotting categorical predictors not working

grid.arrange(p1[1], p2[1], p3[1], p4[1], p5[1])

#################individual models#################       

#Borreliella
Y[39]
VI[39,]

#Variable importance

#pvi <- plot_vi(VI=VI[39,],  X=X,Y=Y[39], modelPerf=ModelPerf_rf, cutoff= 0.05, mode='classification') 

sfl <- mrFlashlight(yhats_rf, X=X, Y=Y, response = "single", mode='classification', index = 39)


#ALE plots 

p1 <- light_ice(sfl, v = "prev_5years_avg_tmean", n_max = 200) %>%
  plot(alpha = 0.05, color = "darkred") +
  labs(title = "Centered ICE plot", y = "probability of occurence")+
  theme_bw()
p1


p2 <- light_ice(sfl, v = "ndvi_point", n_max = 200,) %>%
  plot(alpha = 0.05, color = "darkred") +
  labs(title = "Centered ICE plot", y = "probability of occurence")+
  theme_bw()
p2

p3 <- light_ale(sfl, v = "vpdmin", n_max = 200,) %>%
  plot(alpha = 0.05, color = "darkred") +
  labs(title = "Centered ICE plot", y = "probability of occurence")+
  theme_bw()

#################Interactions#################  
# try iml package

#No real strong interactions detected

miml <- MrIMLconverts(yhats_rf, X=Xdummy, mode='classification')

miml[[39]]

ia_all <- Interaction$new(miml[[39]])

ia_all$results %>%
  ggplot(aes(x = reorder(.feature, .interaction), y = .interaction, fill = .class)) +
  facet_wrap(~ .class, ncol = 2) +
  geom_bar(stat = "identity") +
  #scale_fill_tableau() +
  coord_flip() +
  guides(fill = FALSE)+
  theme_bw()

#pairwise
ia_vpd <- Interaction$new(miml[[39]], feature='vpdmin')

ia_all$results %>%
  ggplot(aes(x = reorder(.feature, .interaction), y = .interaction, fill = .class)) +
  facet_wrap(~ .class, ncol = 2) +
  geom_bar(stat = "identity", alpha = 0.8) +
  #scale_fill_tableau() +
  coord_flip() +
  guides(fill = FALSE)


#plot top interaction
ale <- FeatureEffect$predict(miml[[39]], feature = c("vpdmin", 'tmean'))

 plot(ale) +theme_bw()

 #vivid code - uses vInt as the basis. Much more stable that Friedmans
 
 #also seems to selecct the '0' class
 
 set.seed(113)
 viFit <- vivi(
   fit = extract_fit_parsnip(yhats_rf[[39]]$mod1_k),
   data = data_df,
   response = 'class',
   gridSize = 10,
  # importanceType= 'agnostic', #doesn't work
   normalized=FALSE
   
   )#,
 
 #
 viviHeatmap(mat = viFit) + ggtitle("rf heatmap")
 viviNetwork(mat = viFit)
 
 str(X)
 set.seed(1701)
 #this is definetly predicting the worng class
 
 pdpPairs(data = data_df, fit = fit, response = "class", nmax = 50,  gridSize = 30, #$class=1,
          convexHull = T, probability=T, vars = c('max_snow_depth'  , 'ndvi_point', 'prev_year_avg_tmean', 'vpdmin'))
 
 
 vintTidy(feature_names = n_features, train = dataAll[[i]], parallel = FALSE, mode=mode)
#Vint pdp version
 int.i <- vint(
   object = yhats_rf[[39]]$mod1_k,                    # fitted model object
   feature_names = paste0(names(Xdummy)),  # features for which to compute pairwise interactions statistics
  # n.trees = best.iter,                 # needed if object is of class "gbm"
   parallel = FALSE
 )
