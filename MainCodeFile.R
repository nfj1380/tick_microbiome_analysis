

#' Read in the raw input file. These can be BIOM file from QIIME and convert to a matrix

library(OTUtable)
library(tidyverse)
library(biomformat)


########################################################################################################################
#---------------------Data preparation--------------------------
########################################################################################################################

#import data
envdat <- read.csv('tick_env_final_data_namesFixedMar2021.csv', head=T)
envdatSort <- envdat[order(envdat$Site_Code),]

asvdatFiltered <- as.data.frame(read.csv('ASVtable.csv', head=T)) %>% 
   filter(!Site.Code.Name...Year == "BRACHYSPRIA") %>% 
   filter(!Site.Code.Name...Year == "BLANK") %>%  #blanks
   filter(!Site.Code.Name...Year == "WL19") %>% #site not sampled in the end
   filter(!Site.Code.Name...Year == "WT17")  

    
  #remove the soil as we dont have matching environmental data
  asvdataCleanNoSoil <-  filter(asvdatFiltered  , !Tick.Sex.Sample.Type == "S") 
  
  finalASVdat <- asvdataCleanNoSoil  [order(asvdataCleanNoSoil $Site.Code.Name...Year),]

    countSIte  <-finalASVdat %>% 
      group_by(Site.Code.Name...Year) %>%
      summarise(no_rows = length(Site.Code.Name...Year))
    
#check matching
    match(countSIte$Site.Code.Name...Year, envdatSort$Site_Code)
    
#repeat rows of environmental data to match ASV data
    

    envData_reps <- envdatSort[rep(seq_len(nrow(envdatSort)), countSIte$no_rows),]
    
    #add sequence id to make sure everything is the same
    envData_repsID <- cbind(finalASVdat$Sequncing.File.Name, envData_reps)
    
    #add other relvant data from the ASV data fram to the environmental data

    envData_final <- cbind(envData_repsID, sex=finalASVdat$Tick.Sex.Sample.Type)
    
#--------------------------------------------------------------------------------------------------
    #create random effect data (individual, plate, collection year and site and spatial data)
    indiv <- c(1:355)
    # as matrix
    randonEff <- as.data.frame(cbind(indiv, plate=finalASVdat$Plate.Code, year=finalASVdat$Sample.Collection.Year, site=finalASVdat$Site.Code.Only, latitude= envData_repsID$latitude, longitude=envData_repsID$longitude, season=envData_final$season ))
    #as factor dataframe
    randonEff_f <- data.frame(indiv, plate=finalASVdat$Plate.Code, year=finalASVdat$Sample.Collection.Year, site=finalASVdat$Site.Code.Only, latitude= envData_repsID$latitude, longitude=envData_repsID$longitude, season=envData_final$season)
#--------------------------------------------------------------------------------------------------
    
  #streamline data

    #no need for long/lat in main predictor set
   # envData_final$latitude <- NULL
    #envData_final$longitude <- NULL
    #envData_final$season <- NULL
    
  #reduce to predictor set and test for correlations
    envData_final_predictors <- envData_final[11:101]
    
    
#--------------------------------------------------------------------------------------------------    
#Predictor summary and correlations
    
    library(DataExplorer)
    create_report(envData_final_predictors ) #produces an html report

    #not too much missing data - lets impute!
  library(missForest)
  
    envData_final_noNA <- missForest(envData_final_predictors, variablewise=T) #default values have worked fine previously
    
    envData_final_noNA$OOBerror #all really low mse
    
    envData_final_noNA <- envData_final_noNA$ximp
    

#remove correlated predictors

  #tidymodel way  -adapted! removes less predictors than caret version. 16 predictors at 0.7 cor coefficent threshold
    
    library(tidymodels)
    
    rec <- recipe( tdmean  ~ .,
                  data = envData_final_noNA, retain = TRUE)
    
    corr_filter <- rec %>%
      step_corr(all_predictors(), threshold = .7)
    
    filter_obj <- prep(corr_filter, training = envData_final_noNA)
    
   envData_uncor<- bake(filter_obj, envData_final_noNA) #extract final data
   
   envData_uncor_WC <-cbind(envData_uncor, sex=envData_final$sex)
    
 #' Gather metadata
    ## ------------------------------------------------------------------------
    taxonomy <- read.table('combinded_table_taxa.txt') #got sites not included.
    
    #' 
    #' Remove duplicates at the same level. Level options are from column names from taxonomic dataset. No need to do this
    ## ------------------------------------------------------------------------
    #combine_otus('Genus', table, taxonomy)
  
    #asv data without the metadata and with NAs converted to 0s
    finalASVdataRefined <- finalASVdat[10:903]
    
    finalASVdataRefined[is.na( finalASVdataRefined)] = 0 # Converting NAs to 0s
    
    #' 
#To make a table containing only phyla with at least 10% abundance in any one sample and were observed at any abundance in at least 10% of samples. Convert the table to presence / absence and remove OTUs that are too common for analysis
    ## ------------------------------------------------------------------------
    asvFilter <- as.data.frame(t(filter_taxa(as.data.frame(t(finalASVdataRefined)), abundance=0.05, persistence =5)))
    
    #make presence/absence
    asvFilter[asvFilter>0] <-1  
    
    #remove these
    asvFilter$Total.Reads <- NULL #not sure how these survived the filter
    asvFilter$NA. <- NULL
    
    #riketsia is in every tick
    asvFilter$Rickettsia <- NULL
    
    #not sure about removing taxa yet....
    #sum rows
    #asvFilter$new <- rowSums(t(asvFilter )))
    #remove all common OTUs
    #asvFilterNoCommon <-subset(otuFilter, new < 45)
    #remove 'new' sort column
    #otuFilterNoCommon$new <- NULL

#-----------------------------------------------------------------------------
#################CMRF analysis#################
#---------------------------------------------------------------------------
library(MRFcov)

#' 
#' Extract coordinates columns to use for incorporating spatial splines as covariates in the model (for now, these need to be named `Latitude` and `Longitude`)
## ----message=FALSE, warning=FALSE----------------------------------------
Latitude <- envData_final$latitude
    envData_final$latitude <- NULL
Longitude <- envData_final$longitude
envData_final$longitude <- NULL
coords <- (data.frame(Latitude = Latitude,
                      Longitude = Longitude))


#' Prep the remaining covariates for CRF analysis. Here, we change categorical covariates to model matrix format and scale numeric covariates to unit variance
## ----message=FALSE, warning=FALSE----------------------------------------
envData_uncor_WC$sex <- as.factor(envData_uncor_WC$sex)

str(envData_uncor_WC)

envData_uncor_WC %>%
  cbind(.,data.frame(model.matrix(~.[,'sex'],.)[,-1])) %>%
  dplyr::select(-sex) %>%
  dplyr::rename_all(funs(gsub("\\.|model.matrix", "", .))) %>%
  dplyr::mutate_at(vars(vpdmax:tdmean),
                   funs(as.vector(scale(.)))) -> analysis.data

analysis.data$sexNA <- NULL
analysis.data$sexS <- NULL
#' 
#' We also create a dataset for non-spatial analysis for comparisons. For this, we add the `Latitude` and `Longitude` columns back in and scale them to unit variance
## ----message=FALSE, warning=FALSE----------------------------------------
analysis.data %>%
  dplyr::bind_cols(data.frame(Latitude = Latitude,
                              Longitude = Longitude)) %>%
  dplyr::mutate_at(vars(Latitude, Longitude),
                   funs(as.vector(scale(.)))) -> analysis.data.nonspatial


## ------------------------------------------------------------------------
# Generate a one-hot encoded site design matrix
## ------------------------------------------------------------------------

site_re <- mgcv::smooth.construct2(object = mgcv::s(Site_Code,
                                                    bs = "re"),
                                   data = envData_final, knots = NULL)
site_designmat <- data.frame(site_re$X)
colnames(site_designmat) <- paste0('site', 1:ncol(site_designmat))
head(site_designmat)

# Prep the data for CRF by cross-multiplication
nodes = 4 #number of species

data_prepped <- MRFcov::prep_MRF_covariates(data = asvFilter[,1:nodes], n_nodes = nodes) 
head(data_prepped)


# Add back in the site design matrix. This ensures the CRF coefficients are estimated conditionally on the 
# site effect


data_prepped_covar <- cbind(data_prepped, t=analysis.data[,1])

data_prepped_withSite <- cbind(data_prepped_covar, site_designmat)

#' 

#' 
#' Now we generate MRF and CRF models to determine which model fits the data most appropriately. First, the nonspatial MRF (without covariates)
## ------------------------------------------------------------------------
tick.crf <- MRFcov(data = data_prepped_covar, 
                    n_nodes = nodes, family = 'binomial',
                    n_cores = 3, n_covariates = 17, prep_covariates = FALSE
                    )
plotMRF_hm(tick.mrf,  plot_observed_vals = TRUE, data = data_prepped)

# Fit the model. Be sure to explicitly state how many covariates there were prior to prepping (1 in this case)
crf <- MRFcov(data = data_prepped, 
              n_covariates = 1, family = 'binomial')
tick.mrf$direct_coefs
tick.mrf$indirect_coefs

#' 
#' Use this model to generate predictions
## ------------------------------------------------------------------------
mrf.predictions <- predict_MRF(data = data_prepped[,1:nodes],
                               MRF_mod = tick.mrf)

#' 

#' 
#' Test how well the MRF predicts the data by comparing the predicted to the observed values. Here, we split the data into five folds and test for model specificity, sensitivity, positive predictive value and proportion of true predictions using the `cv_MRF_diag` function
## ------------------------------------------------------------------------
mrf.cv <- lapply(seq_len(100), function(x){
  cv_MRF_diag(data = analysis.data.nonspatial[,1:nodes], n_nodes = nodes,
              n_folds = 5,
              n_cores = 1, family = 'binomial',
              compare_null = FALSE, plot = FALSE,
              cached_model = moose.mrf,
              cached_predictions = list(predictions = mrf.predictions),
              sample_seed = 10)
})
mrf.cv <- do.call(rbind, mrf.cv)

#' 
#' Next, we repeat the above for a nonspatial CRF (including environmental covariates, but without spatial splines)
## ------------------------------------------------------------------------
moose.crf <- MRFcov(data = analysis.data.nonspatial, 
                    n_nodes = nodes, family = 'binomial',
                    n_cores = 3)
crf.predictions <- predict_MRF(data = analysis.data.nonspatial[,1:nodes],
                               MRF_mod = moose.crf)
crf.cv <- lapply(seq_len(100), function(x){
  cv_MRF_diag(data = analysis.data.nonspatial, n_nodes = nodes,
              n_folds = 5,
              n_cores = 1, family = 'binomial',
              compare_null = FALSE, plot = FALSE,
              cached_model = moose.crf,
              cached_predictions = list(predictions = crf.predictions),
              sample_seed = 10)
})
crf.cv <- do.call(rbind, crf.cv)

#' 
#' Next, we fit the spatial CRF. Here, the coordinates are used to produce spatial regression splines with the call 
#' `mgcv::smooth.construct2(object = mgcv::s(Latitude, Longitude, bs = "gp"), data = coords, knots = NULL)`. Splines are a series of different equations pieced together, 
#' which tends to increase the accuracy of interpolated values at the cost of  the ability to project outside the data range. These are highly appropriate here, as we are not interested in using them for prediction but more to account for non-independence in our observations. 
## ------------------------------------------------------------------------
moose.crf.spatial <- MRFcov_spatial(data = analysis.data, 
                                    n_nodes = nodes, family = 'binomial',
                                    coords = coords, n_cores = 3)
spatial.predictions <- predict_MRF(data = moose.crf.spatial$mrf_data,
                                   prep_covariates = F,
                                   MRF_mod = moose.crf.spatial)
spatial.cv <- lapply(seq_len(100), function(x){
  cv_MRF_diag(data = analysis.data, n_nodes = nodes,
              n_folds = 5,
              n_cores = 2, family = 'binomial',
              compare_null = FALSE, plot = FALSE,
              cached_model = moose.crf.spatial,
              cached_predictions = list(predictions = spatial.predictions),
              sample_seed = 10)
})
spatial.cv <- do.call(rbind, spatial.cv)

#' 
#' \pagebreak
#' 
#' Now that we have run all of the models and generated predictions, we can examine the proportion of unique observations that were correctly predicted by each of the different models
## ------------------------------------------------------------------------
quantile(mrf.cv$mean_tot_pred, probs = c(0.025, 0.5, 0.975))
quantile(crf.cv$mean_tot_pred, probs = c(0.025, 0.5, 0.975))
quantile(spatial.cv$mean_tot_pred, probs = c(0.025, 0.5, 0.975))








#-----------------------------------------------------------------------------
#################HMSC Analysis s#################
#-----------------------------------------------------------------------------

library(Hmsc)
set.seed(29698)

#Turn into matrices 
YMat<-as.matrix(asvFilter)
XMat<-model.matrix(~.,data=envData_uncor)


str(randonEff_f)

xycoords <- as.matrix(cbind(as.numeric(randonEff$latitude),as.numeric(randonEff$longitude ))) #needs to be a matrix

#drop unused levels
randonEff$site <- droplevels(randonEff$site)

row.names(xycoords) <- randonEff$site
#get data into right format

#add some random noise to make each tick have unique coordinates (0.0001)

p <- runif(355, 0.00001, 0.00002)
xycoords_n <- xycoords[,]+p

rL <-  HmscRandomLevel (sData = xycoords_n) #adding spatial component
studyDesign <- data.frame(sample = as.factor(randonEff$site)) #adding site random effect

m <-  Hmsc(Y = YMat, XData = as.data.frame(XMat), XFormula = ~XMat, distr="probit",
           studyDesign = studyDesign, ranLevels=list("sample"=rL) ) #probit model for pa data

#mcmc characteristics
nChains = 2
thin = 5
samples = 1000
transient = 500*thin
verbose = 500*thin

#run model
m1 <- sampleMcmc(m, thin = thin, samples = samples, transient = transient,
                 nChains = nChains, verbose = verbose)


saveRDS(m1, file = "m1.rds")


preds.spatial = computePredictedValues(m1)
MF.spatial = evaluateModelFit(hM=m1, predY=preds.spatial)
MF.spatial

#problem is sites have the same spatial coordinates
#check mcmc diagnostics
mpost <- convertToCodaObject(m1)
summary(mpost$Beta)
effectiveSize(mpost$Beta)
gelman.diag(mpost$Beta,multivariate=FALSE)$psrf


# check R2 etc
preds <- computePredictedValues(m1)
evaluateModelFit(hM=m1, predY=preds)

#plot paramters
plot(mpost$Beta)

#plot betas
postBeta = getPostEstimate(m1, parName = "Beta")
par(mar=c(5,11,2.5,0))
plotBeta(m1,
         post = postBeta, 
         plotTree = F,
         spNamesNumbers = c(T,F))
#only signif parameters

plotBeta(m1, 
         post = postBeta,
         param = "Mean",
         plotTree = F,  
         spNamesNumbers = c(T,F))

VP <-  computeVariancePartitioning(m1)
par(mar=c(4,4,4,4))
plotVariancePartitioning(m1, VP = VP,
                         las = 2, horiz=F)

### Mixing object
mixing <- as.mcmc(m1, parameters = "paramX")
### Draw trace and density plots for all combination of parameters
plot(mixing)

#formprior <- as.HMSCprior(formdata1)
#formparam <- as.HMSCparam(formdata1, formprior)      