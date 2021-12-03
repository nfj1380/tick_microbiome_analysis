

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
   filter(!Site.Code.Name...Year == "BRACHYSPRIA") %>% #splatter control %>% 
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
    #as factor data frame
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
    
   # library(DataExplorer)
   # create_report(envData_final_predictors ) #produces an html report

    #not too much missing data - lets impute!
  library(missForest)
  
    envData_final_noNA <- missForest(envData_final_predictors, variablewise=T) #default values have worked fine previously
    
    envData_final_noNA$OOBerror #all really low mse
    
    envData_final_noNA <- envData_final_noNA$ximp
    

#remove correlated predictors

  #tidymodel way  -adapted! removes less predictors than caret version. 16 predictors at 0.7 cor coefficent threshold
    
    library(GGally)
    
    ggpairs(envData_final_noNA)
    ggcorr(envData_final_noNA)
    
    library(tidymodels)
    
    rec <- recipe( vpdmax   ~ .,
                  data = envData_final_noNA, retain = TRUE)
    
    corr_filter <- rec %>%
      step_corr(all_predictors(), threshold = .8)
    
    filter_obj <- prep(corr_filter, training = envData_final_noNA)
    
   envData_uncor<- bake(filter_obj, envData_final_noNA) #extract final data
   
   corCheck <- cor(envData_uncor)
   
   envData_uncor_simp <- envData_uncor
   envData_uncor_simp[13:19] <- NULL
   envData_uncor_WC <-cbind(envData_uncor_simp, sex=envData_final$sex)
   
   
   #want to include tmax too (need to remove)
    
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
    asvFilter_230asv <- as.data.frame(t(filter_taxa(as.data.frame(t(finalASVdataRefined)), abundance=0.01, persistence =2.5)))
  #this threshold works best for mrf
    #riketsia is in every female tick and quite common in males too. 600 read threshold
    asvFilter_230asv$Rickettsia[asvFilter_230asv$Rickettsia<600] <- 0
    asvFilter_230asv$Rickettsia[asvFilter_230asv$Rickettsia>600] <- 1
    
    #remove the positive control for splatter from the ASV table too.
    asvFilter_230asv$Brachyspira <- NULL
    
    #make presence/absence
    asvFilter_230asv[asvFilter_230asv>0] <- 1  
    
    #remove these columns
    asvFilter_230asv$Total.Reads <- NULL
    asvFilter_230asv$NA. <- NULL
  

#-----------------------------------------------------------------------------
#################MRF analysis#################
#---------------------------------------------------------------------------
library(MRFcov)

#' 
#' Extract coordinates columns to use for incorporating spatial splines as covariates in the model (for now, these need to be named `Latitude` and `Longitude`)
## ----message=FALSE, warning=FALSE----------------------------------------
Latitude <- envData_final$latitude
    envData_final$latitude <- NULL #fix this
Longitude <- envData_final$longitude
envData_final$longitude <- NULL
coords <- (data.frame(Latitude = Latitude,
                      Longitude = Longitude))


#' Prep the remaining covariates for CRF analysis. Here, we change categorical covariates to model matrix format and scale numeric covariates to unit variance
## ----message=FALSE, warning=FALSE----------------------------------------
envData_uncor_WC$sex <- as.factor(envData_uncor_WC$sex)


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
                   funs(as.vector(scale(.)))) -> analysis.data.spatial




analysis.data.combined  <- cbind( asvFilter_230asv, analysis.data ) 

# Prep the data for CRF by cross-multiplication
nodes = 231 #number of species

## ------------------------------------------------------------------------
# Generate a one-hot encoded site design matrix
## ------------------------------------------------------------------------

site_re <- mgcv::smooth.construct(object = mgcv::s(Site_Code, bs = "re"),
                                   data = envData_final, knots = NULL)


site_designmat <- data.frame(site_re$X)
colnames(site_designmat) <- paste0('site', 1:ncol(site_designmat))

# How many covariates are there, assuming there are four nodes?
n_covariates <- ncol(analysis.data)

## ------------------------------------------------------------------------
# Non-spatial MRF
## ------------------------------------------------------------------------
library(MRFcov)

# Prep the data

# Put site variables back in with the prepped data
data_for_mrf_nonSpatial <- cbind(asvFilter_230asv, data.frame(site_designmat))

# Fit the model
tick_mrf_nonSpatial <- MRFcov(data = data_for_mrf_nonSpatial, n_nodes = nodes, family = 'binomial',
              n_cores = 4)

saveRDS(tick_mrf_nonSpatial, "mrf_model")

# Inspect outputs
tick_mrf_nonSpatial $intercepts
plotMRF_hm(tick_mrf_nonSpatial )

library(igraph)
net <- graph.adjacency(tick_mrf_nonSpatial$graph, weighted = T, mode = "undirected")

net_simplified_pos <- delete_edges(net, which(E(net)$weight< 0.15))

plot.igraph(net_simplified_pos, layout = igraph::layout_with_mds(net_simplified_pos, dim=2),
                    edge.width = abs(igraph::E(net_simplified_pos)$weight*4),
                    vertex.size = 2,
                    vertex.color = 'black',
                    vertex.label.cex = 0.6,
                    vertex.label.font=1,
                    vertex.label.color='black',
                    edge.color = ifelse(igraph::E(net_simplified_pos)$weight < 0, 
                                        'blue',
                                        'red'))

net_simplified_neg <- delete_edges(net, which(E(net)$weight> -0.1))

igraph::plot.igraph(net_simplified_neg , layout = igraph::layout_with_mds(net_simplified_neg , dim=2),
                    edge.width = abs(igraph::E(net_simplified_neg )$weight*4),
                    vertex.size = 2,
                    vertex.color = 'black',
                    vertex.label.cex = 0.6,
                    vertex.label.font=3,
                    vertex.label.color='black',
                    edge.color = ifelse(igraph::E(net_simplified_neg )$weight < 0, 
                                        'blue',
                                        'red'))


tick_mrf_nonSpatial $key_coefs$Borreliella #one group of interest

tick_mrf_nonSpatial $key_coefs$Wolbachia #one group of interest - only shaped by  Streptococcus 

# Add back in the site design matrix. This ensures the CRF coefficients are estimated conditionally on the 
# site effect

write.csv(data_prepped_withSite,"data_prepped_withSite.csv", row.names = TRUE)
 
#' Now we generate MRF and CRF models to determine which model fits the data most appropriately. First, the nonspatial MRF (without covariates)
## ------------------------------------------------------------------------
#' 
#' Use this model to generate predictions
## ------------------------------------------------------------------------
mrf.predictions <- predict_MRF(data = data_for_mrf_nonSpatial[1:nodes],
                               MRF_mod = tick_mrf_nonSpatial, prep_covariates = FALSE)

#' 
#' Test how well the MRF predicts the data by comparing the predicted to the observed values. Here, we split the data into five folds and test for model specificity, sensitivity, positive predictive value and proportion of true predictions using the `cv_MRF_diag` function
## ------------------------------------------------------------------------
mrf.cv <- lapply(seq_len(100), function(x){
  cv_MRF_diag(data = data_for_mrf_nonSpatial[,1:nodes], n_nodes = nodes,
              n_folds = 5,
              n_cores = 2, family = 'binomial',
              compare_null = FALSE, plot = FALSE,
              cached_model = tick_mrf_nonSpatial,
              cached_predictions = list(predictions = mrf.predictions),
              sample_seed = 10)
})
mrf.cv <- do.call(rbind, mrf.cv)

#' 
#' Next, we repeat the above for a nonspatial CRF (including environmental covariates, but without spatial splines)
## ------------------------------------------------------------------------
# Prep the data

data_nosite_prepped <- prep_MRF_covariates(data = analysis.data.combined , n_nodes = nodes)

# Put site variables back in with the prepped data
data_for_crf_non_spatial <- cbind(data_nosite_prepped, data.frame(site_designmat))

# Fit the model
tick_crf_non_spatial <- MRFcov(data =data_for_crf_non_spatial, n_nodes = nodes, family = 'binomial',
                   prep_covariates = FALSE, n_covariates = n_covariates,
                   n_cores = 4)

saveRDS(tick_crf_non_spatial, 'tick_crf_nonspatial') 
# Inspect outputs
crf$intercepts
plotMRF_hm(tick_crf )

crf.predictions <- predict_MRF(data = data_for_crf_non_spatial,
                               prep_covariates = F,
                               MRF_mod = tick_crf_non_spatial)

crf.cv <- lapply(seq_len(100), function(x){
  cv_MRF_diag(data = data_for_crf_non_spatial, n_nodes = nodes,
              n_folds = 5,
              n_cores = 3, family = 'binomial',
              compare_null = FALSE, plot = FALSE,
              cached_model = tick_crf_non_spatial,
              cached_predictions = list(predictions = crf.predictions),
              sample_seed = 10)
})
crf.cv <- do.call(rbind, crf.cv)

# Check how important each covariate is for predicting changing interactions
# mean of covariate absolute effect sizes

cov.imp.mean <- lapply(seq_along(tick_crf_non_spatial$indirect_coefs), function(x){
  graph <- tick_crf_non_spatial$indirect_coefs[[x]][[1]]
  coefs <- graph[upper.tri(graph)]
  round(mean(abs(coefs[which(coefs != 0 )])), 4)
})

impData <- cbind(cov.imp.mean, names(analysis.data))
#' 
#' Next, we fit the spatial CRF. Here, the coordinates are used to produce spatial regression splines with the call 
#' `mgcv::smooth.construct2(object = mgcv::s(Latitude, Longitude, bs = "gp"), data = coords, knots = NULL)`. Splines are a series of different equations pieced together, 
#' which tends to increase the accuracy of interpolated values at the cost of  the ability to project outside the data range. These are highly appropriate here, as we are not interested in using them for prediction but more to account for non-independence in our observations. 
## ------------------------------------------------------------------------
tick.mrf.spatial <- MRFcov_spatial(data = data_for_mrf_nonSpatial, 
                                    n_nodes = nodes, family = 'binomial',
 
                                   coords = coords, n_cores = 4)

saveRDS(tick.mrf.spatial, "FINAL_model")

FINAL_model <- readRDS('FINAL_model')

library(igraph)

net <- graph.adjacency(tick.mrf.spatial $graph, weighted = T, mode = "undirected")

net_simplified_pos <- delete_edges(net, which(E(net)$weight< 0.3))

igraph::plot.igraph(net_simplified_pos, layout = igraph::layout_with_mds(net_simplified_pos, dim=2),
                    edge.width = abs(igraph::E(net_simplified_pos)$weight*2),
                    vertex.size = 1,
                    vertex.color = 'black',
                    vertex.label.cex = 0.6,
                    vertex.label.font=3,
                    vertex.label.color='black',
                    edge.color = ifelse(igraph::E(net_simplified_pos)$weight < 0, 
                                        'blue',
                                        'red'))

net_simplified_neg <- delete_edges(net, which(E(net)$weight> -0.1))

igraph::plot.igraph(net_simplified_neg , layout = igraph::layout_with_mds(net_simplified_neg , dim=2),
                    edge.width = abs(igraph::E(net_simplified_neg )$weight*2),
                    vertex.size = 2,
                    vertex.color = 'black',
                    vertex.label.cex = 0.6,
                    vertex.label.font=3,
                    vertex.label.color='black',
                    edge.color = ifelse(igraph::E(net_simplified_neg )$weight < 0, 
                                        'blue',
                                        'red'))


spatial.predictions <- predict_MRF(data = tick.mrf.spatial$mrf_data,
                                   prep_covariates = F,
                                   MRF_mod = tick.mrf.spatial)

spatial.cv <- lapply(seq_len(100), function(x){
  cv_MRF_diag(data = data_for_mrf_nonSpatial, n_nodes = nodes,
              n_folds = 5,
              n_cores = 4, family = 'binomial',
              compare_null = FALSE, plot = FALSE,
              cached_model = tick.mrf.spatial ,
              cached_predictions = list(predictions = spatial.predictions),
              sample_seed = 10)
})
spatial.cv <- do.call(rbind, spatial.cv)

#' Next, we repeat the above for a nspatial CRF (including environmental covariates with splines)
## ------------------------------------------------------------------------
# Prep the data


# Fit the model
tick_crf_spatial <- MRFcov_spatial(data =data_for_crf_non_spatial, n_nodes = nodes, family = 'binomial',
                               prep_covariates = FALSE, n_covariates = n_covariates, coords=coords,
                                                              n_cores = 4)
saveRDS(tick_crf_spatial, 'tick_crf_spatial') 

tick_crf_spatial <- readRDS('tick_crf_spatial')
# Inspect outputs
crf$intercepts
plotMRF_hm(tick_crf_spatial )

#taxa of interest
tick_crf_spatial$key_coefs$Borreliella 
tick_crf_spatial$key_coefs$Borrelia 
tick_crf_spatial$key_coefs$Anaplasma
tick_crf_spatial$key_coefs$Ehrlichia 
tick_crf_spatial$key_coefs$Francisella 
tick_crf_spatial$key_coefs$Lawsonella
tick_crf_spatial$key_coefs$Rickettsia 

tick_crf_spatial$key_coefs$Wolbachia
#secondary importance
tick_crf_spatial$key_coefs$Bacillus
tick_crf_spatial$key_coefs$Corynebacterium
tick_crf_spatial$key_coefs$Legionella 
tick_crf_spatial$key_coefs$Mycobacterium
tick_crf_spatial$key_coefs$Streptococcus 
tick_crf_spatial$key_coefs$Pseudomonas
tick_crf_spatial$key_coefs$Ralstonia
tick_crf_spatial$key_coefs$Rhodococcus
tick_crf_spatial$key_coefs$Sphingomonas
tick_crf_spatial$key_coefs$Staphylococcus
tick_crf_spatial$key_coefs$Streptomyces

crf.spatial.predictions <- predict_MRF(data = tick_crf_spatial$mrf_data,
                                       prep_covariates = F,
                               MRF_mod = tick_crf_spatial)

crf.cv.spatial <- lapply(seq_len(100), function(x){
  cv_MRF_diag(data = tick_crf_spatial$mrf_data, n_nodes = nodes,
              n_folds = 5,
              n_cores = 4, family = 'binomial',
              compare_null = FALSE, plot = FALSE,
              cached_model = tick_crf_spatial,
              cached_predictions = list(predictions = crf.spatial.predictions),
              sample_seed = 10)
})
crf.cv.spatial <- do.call(rbind, crf.cv.spatial)

cov.imp.mean <- lapply(seq_along(tick_crf_spatial$indirect_coefs), function(x){
  graph <- tick_crf_spatial$indirect_coefs[[x]][[1]]
  coefs <- graph[upper.tri(graph)]
  round(mean(abs(coefs[which(coefs != 0 )])), 4)
})
impData <- cbind(cov.imp.mean, names(analysis.data))


#' Now that we have run all of the models and generated predictions, we can examine the proportion of unique observations that were correctly predicted by each of the different models
## ------------------------------------------------------------------------
quantile(mrf.cv$mean_tot_pred, probs = c(0.025, 0.5, 0.975))
quantile(crf.cv$mean_tot_pred, probs = c(0.025, 0.5, 0.975))
quantile(spatial.cv$mean_tot_pred, probs = c(0.025, 0.5, 0.975))
quantile(crf.cv.spatial$mean_tot_pred, probs = c(0.025, 0.5, 0.975))

quantile(crf.cv$Wolbachia, probs = c(0.025, 0.5, 0.975))

mrf.cv$model <- "MRF"
crf.cv$model <- "CRF"
spatial.cv$model <- "Spatial MRF"
crf.cv.spatial$model <- "Spatial CRF"
plot.dat <- rbind(mrf.cv, crf.cv, spatial.cv,
                  crf.cv.spatial)
plot.dat$model <- factor(plot.dat$model,
                         levels = c("MRF", "Spatial MRF", "CRF",
                                    "Spatial CRF"))
scaleFUN <- function(x) sprintf("%.2f", x)
preds.plot <- ggplot(plot.dat, aes(model,
                                   mean_tot_pred, colour = model)) + geom_boxplot() +
  labs(y = "True predictions", x = "") +
  scale_y_continuous(labels = scaleFUN,
                     breaks = seq(min(plot.dat$mean_tot_pred),
                                  max(plot.dat$mean_tot_pred),
                                  length.out = 5)) + theme_classic() +
  theme(legend.position = "")
spec.plot <- ggplot(plot.dat, aes(model,
                                  mean_specificity, colour = model)) +
  geom_boxplot() + labs(y = "Specificity",
                        x = "") + scale_y_continuous(labels = scaleFUN,
                                                     breaks = seq(min(plot.dat$mean_specificity),
                                                                  max(plot.dat$mean_specificity), length.out = 5)) +
  theme_classic() + theme(legend.position = "")
sens.plot <- ggplot(plot.dat, aes(model,
                                  mean_sensitivity, colour = model)) +
  geom_boxplot() + labs(y = "Sensitivity",
                        x = "Model") + scale_y_continuous(labels = scaleFUN,
                                                          breaks = seq(min(plot.dat$mean_sensitivity),
                                                                       max(plot.dat$mean_sensitivity), length.out = 5)) +
  theme_classic() + theme(legend.position = "")
ppv.plot <- ggplot(plot.dat, aes(model, mean_pos_pred,
                                 colour = model)) + geom_boxplot() + labs(y = "PPV",
                                                                          x = "") + theme_classic() + scale_y_continuous(labels = scaleFUN,
                                                                                                                         breaks = seq(min(plot.dat$mean_pos_pred),
                                                                                                                                      max(plot.dat$mean_pos_pred), length.out = 5)) +
  theme(legend.position = "")

gridExtra:: grid.arrange(ppv.plot, spec.plot,sens.plot, ncol = 3)


library(igraph)
net <- graph.adjacency(tick_crf_spatial$graph, weighted = T, mode = "undirected")

net_simplified_pos <- delete_edges(net, which(E(net)$weight< 0.2))

plot.igraph(net_simplified_pos, layout = igraph::layout_with_mds(net_simplified_pos, dim=2),
            edge.width = abs(igraph::E(net_simplified_pos)$weight*4),
            vertex.size = 2,
            vertex.color = 'black',
            vertex.label.cex = 0.4,
            vertex.label.font=1,
            vertex.label.color='black',
            edge.color = ifelse(igraph::E(net_simplified_pos)$weight < 0, 
                                'blue',
                                'red'))

net_simplified_neg <- delete_edges(net, which(E(net)$weight> -0.05))

igraph::plot.igraph(net_simplified_neg , layout = igraph::layout_with_mds(net_simplified_neg , dim=2),
                    edge.width = abs(igraph::E(net_simplified_neg )$weight*4),
                    vertex.size = 2,
                    vertex.color = 'black',
                    vertex.label.cex = 0.4,
                    vertex.label.font=3,
                    vertex.label.color='black',
                    edge.color = ifelse(igraph::E(net_simplified_neg )$weight < 0, 
                                        'blue',
                                        'red'))


#based on these results we'll choose the spatial MRF as it has higher performance
tick.bootstrap.spatial <- bootstrap_MRF(data = data_for_mrf_nonSpatial,
                                         n_nodes = nodes, family = 'binomial',
                                         spatial = TRUE,
                                         coords = coords, 
                                         n_cores = 5,
                                         n_bootstraps = 100)


#for CRF
tick.bootstrap.spatial.crf <- bootstrap_MRF(data = data_for_crf_non_spatial,
                                        n_nodes = nodes, family = 'binomial',
                                        n_covariates = n_covariates,
                                        spatial = TRUE,
                                        coords = coords, 
                                        n_cores = 4,
                                        n_bootstraps = 50)

saveRDS(tick.bootstrap.spatial.crf, 'tick.bootstrap.spatial.crf' )

tick.bootstrap.spatial <- readRDS('tick.bootstrap.spatial')
#things to do. 
#Add more species
#work out a way to plot networks 
tick.bootstrap.spatial$indirect_coef_mean
plotMRF_hm(MRF_mod = tick.bootstrap.spatial )
save(tick.bootstrap.spatial ,file="tick.bootstrap.spatial.RData")

#look at some ASVs
tick.bootstrap.spatial$mean_key_coefs$Borreliella #one group of interest

tick.bootstrap.spatial$mean_key_coefs$Wolbachia #one group of interest - only shaped by  Streptococcus
tick.bootstrap.spatial$mean_key_coefs$Borrelia #Massilia is by far the most important predictor

tick.bootstrap.spatial$mean_key_coefs$Anaplasma

str(tick.bootstrap.spatial$direct_coef_upper90)
str(tick.bootstrap.spatial)

summary(tick.bootstrap.spatial)
net <- igraph::graph.adjacency(tick.bootstrap.spatial$graph, weighted = T, mode = "undirected")
igraph::plot.igraph(net, layout = igraph::layout_with_mds(net, dim=2),
                    edge.width = abs(igraph::E(net)$weight),
                    vertex.size = 2,
                    vertex.color = 'black',
                    vertex.label.cex = 0.7,
                    edge.color = ifelse(igraph::E(net)$weight < 0, 
                                        'blue',
                                        'red'))

adj_mats <- predict_MRFnetworks(data = data_for_mrf_nonSpatial[1:nodes, ],
                                MRF_mod = tick.bootstrap.spatial,
                                metric = 'degree',
                                cutoff = 0.33,
                                prep_covariates=F
)
colnames(adj_mats) <- colnames(data_for_mrf_nonSpatial[, 1:nodes])
apply(adj_mats, 2, summary)

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