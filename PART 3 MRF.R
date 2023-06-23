#-----------------------------------------------------------------------------
#################MRF analysis#################
#---------------------------------------------------------------------------
library(MRFcov)
library(tidyverse)

asvFilter_230asv <- readRDS('asvFilter_230asv')
env_data_withHost<- readRDS('env_data_withHost')
complete_envdata <- readRDS('complete_data')

coords <- readRDS('coords')
# env_data_withHost_scaled <- env_data_withHost_scaled %>% 
#   mutate(sex = recode(sex, 
#                     "M" = "0", 
#                     "F" = "1"))
#' Prep the remaining covariates for CRF analysis. Here, we change categorical covariates to model matrix format and scale numeric covariates to unit variance
## ----message=FALSE, warning=FALSE----------------------------------------
env_data_withHost$sex <- as.factor(env_data_withHost$sex)

#scale
env_data_withHost %>%
  cbind(.,data.frame(model.matrix(~.[,'sex'],.)[,-1])) %>%
  dplyr::select(-sex) %>%
  dplyr::rename_all(funs(gsub("\\.|model.matrix", "", .))) %>%
  dplyr::mutate_at(vars(tmean:ndvi_point),
                   funs(as.vector(scale(.)))) -> analysis.data

analysis.data$sexNA <- NULL
analysis.data$sexS <- NULL
#' 
#' We also create a dataset for non-spatial analysis for comparisons. For this, we add the Latitude and Longitude` columns back in and scale them to unit variance
## ----message=FALSE, warning=FALSE----------------------------------------
analysis.data %>%
  dplyr::bind_cols(data.frame(Latitude = coords$Latitude,
                              Longitude = coords$Longitude)) %>%
  dplyr::mutate_at(vars(Latitude, Longitude),
                   funs(as.vector(scale(.)))) -> analysis.data.spatial

analysis.data.combined  <- cbind( asvFilter_230asv, analysis.data ) 

# Prep the data for CRF by cross-multiplication
nodes = 231 #number of species

## ------------------------------------------------------------------------
# Fix site codes - remove year and combine sites close to each other to reduce down
## ------------------------------------------------------------------------

library(stringr)

#remove year
complete_envdata$Site_Code <- str_sub(complete_envdata$Site_Code, end=-3)

complete_envdata$Site_Code[complete_envdata$Site_Code=='KF']= 'CA' #merge close sites
complete_envdata$Site_Code[complete_envdata$Site_Code=='GMW']= 'CA' 
complete_envdata$Site_Code[complete_envdata$Site_Code=='SU']= 'CNF' 
## ------------------------------------------------------------------------
# Generate a one-hot encoded site design matrix
## ------------------------------------------------------------------------

site_re <- mgcv::smooth.construct(object = mgcv::s(Site_Code, bs = "re"),
                                  data = complete_envdata, knots = NULL)


site_designmat <- data.frame(site_re$X)
colnames(site_designmat) <- paste0('site', 1:ncol(site_designmat))

# How many covariates are there, assuming there are four nodes?
n_covariates <- ncol(analysis.data)

## ------------------------------------------------------------------------
# Non-spatial MRF
## ------------------------------------------------------------------------

# Prep the data

# Put site variables back in with the prepped data
data_for_mrf_nonSpatial <- cbind(asvFilter_230asv, site_designmat)

glimpse(data_for_mrf_nonSpatial)
# Fit the model

if(T){
  
tick_mrf_nonSpatial <- MRFcov(data = data_for_mrf_nonSpatial, n_nodes = nodes, family = 'binomial',
                              n_cores = 5)

saveRDS(tick_mrf_nonSpatial, "mrf_model")

tick_mrf_nonSpatial <- readRDS('mrf_model')

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

degree(net_simplified_pos)
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

tick_mrf_nonSpatial $key_coefs$Wolbachia
tick_mrf_nonSpatial $key_coefs$Rickettsia #one group of interest - only shaped by  Streptococcus 

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

}

if(F){

#bootstrap
bootedMRF <- bootstrap_MRF(data = data_for_mrf_nonSpatial,
                               n_nodes = nodes,
                               family = 'binomial',
                               sample_prop = 0.7,
                               n_bootstraps = 5,
                               n_cores= 10
                           )

bootedMRF$mean_key_coefs$Borreliella
bootedMRF$mean_key_coefs$Anaplasma
}

#' 
#' Next, we repeat the above for a nonspatial CRF (including environmental covariates, but without spatial splines)
## ------------------------------------------------------------------------
# Prep the data

data_nosite_prepped <- prep_MRF_covariates(data = analysis.data.combined , n_nodes = nodes)

# Put site variables back in with the prepped data
data_for_crf_non_spatial <- cbind(data_nosite_prepped, site_designmat)
glimpse(data_for_crf_non_spatial )

if(F){
  
# Fit the model
tick_crf_non_spatial <- MRFcov(data =data_for_crf_non_spatial, n_nodes = nodes, family = 'binomial',
                               prep_covariates = FALSE, n_covariates = n_covariates,
                               n_cores = 5)

saveRDS(tick_crf_non_spatial, 'tick_crf_nonspatial') 

tick_crf_nonspatial<- readRDS('tick_crf_nonspatial')
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

}

#' 
#' Next, we fit the spatial CRF. Here, the coordinates are used to produce spatial regression splines with the call 
#' `mgcv::smooth.construct2(object = mgcv::s(Latitude, Longitude, bs = "gp"), data = coords, knots = NULL)`. Splines are a series of different equations pieced together, 
#' which tends to increase the accuracy of interpolated values at the cost of  the ability to project outside the data range. These are highly appropriate here, as we are not interested in using them for prediction but more to account for non-independence in our observations. 
## ------------------------------------------------------------------------

if(F){
  
tick.mrf.spatial <- MRFcov_spatial(data = data_for_mrf_nonSpatial, 
                                   n_nodes = nodes, family = 'binomial',
                                   coords = coords, n_cores = 5)

saveRDS(tick.mrf.spatial , "tick.mrf.spatial")

tick.mrf.spatial  <- readRDS('tick.mrf.spatial')

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

net <- graph.adjacency(tick.mrf.spatial$graph, weighted = T, mode = "undirected")

net_simplified_pos <- delete_edges(net, which(E(net)$weight< 0.4))

Isolated_pos <-  which(degree(net_simplified_pos)==0) #remove isolated vertices
G2 <-  delete.vertices(net_simplified_pos, Isolated_pos)

igraph::plot.igraph(G2, layout = igraph::layout_with_mds(G2, dim=2),
                    edge.width = abs(igraph::E(G2)$weight*2),
                    vertex.size = 1,
                    vertex.color = 'black',
                    vertex.label.cex = 0.6,
                    vertex.label.font=3,
                    vertex.label.color='black',
                    edge.color = ifelse(igraph::E(net_simplified_pos)$weight < 0, 
                                        'blue',
                                        'red'))


}

#' Next, we repeat the above for a nspatial CRF (including environmental covariates with splines)
## ------------------------------------------------------------------------
# Prep the data

#double check this. Cooords may not have worked
# Fit the model

if(F){
  
tick_crf_spatial <- MRFcov_spatial(data =data_for_crf_non_spatial, n_nodes = nodes, family = 'binomial',
                                   prep_covariates = FALSE, n_covariates = n_covariates, coords=coords,
                                   n_cores = 5)
saveRDS(tick_crf_spatial, 'tick_crf_spatial') 

tick_crf_spatial <- readRDS('tick_crf_spatial')
# Inspect outputs
crf$intercepts
plotMRF_hm(tick_crf_spatial )

#taxa of interest

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


}

if(F){
  
bootedSpatial <- bootstrap_MRF(data = data_for_crf_non_spatial, n_nodes = nodes,
                               family = 'binomial',
                               spatial = TRUE,
                               n_covariates = n_covariates,
                               coords = coords,
                               sample_prop = 0.7,
                               n_cores=4,
                               n_bootstraps = 4)
}

if(F){
  
library(igraph)

net <- igraph::graph.adjacency(tick_crf_spatial$graph, weighted = T, mode = "undirected")

net_simplified_pos <- delete_edges(net, which(E(net)$weight< 0.4))


Isolated_pos <-  which(degree(net_simplified_pos)==0) #remove isolated vertices
G2 <-  delete.vertices(net_simplified_pos, Isolated_pos)

igraph::plot.igraph(G2, layout = igraph::layout_with_mds(G2, dim=2),
                    edge.width = abs(igraph::E(G2)$weight*2),
                    vertex.size = 1,
                    vertex.color = 'black',
                    vertex.label.cex = 0.6,
                    vertex.label.font=3,
                    vertex.label.color='black',
                    edge.color = ifelse(igraph::E(net_simplified_pos)$weight < 0, 
                                        'blue',
                                        'red'))

#degree

degreeDist <- data.frame(deg=degree(G2))
degreeDist_ordered<- arrange(degreeDist, deg)

write.csv(degreeDist_ordered, 'degreeDist_ordered.csv')

net_simplified_neg <- delete_edges(net, which(E(net)$weight> -0.1))
#remove disconnected edged
Isolated_neg <-  which(degree(net_simplified_neg)==0)
G2 <-  delete.vertices(net_simplified_neg, Isolated_neg)

igraph::plot.igraph(G2, layout = igraph::layout_with_mds(G2 , dim=2),
                    edge.width = abs(igraph::E(G2)$weight*4),
                    vertex.size = 2,
                    vertex.color = 'black',
                    vertex.label.cex = 0.6,
                    vertex.label.font=3,
                    vertex.label.color='black',
                    edge.color = ifelse(igraph::E(net_simplified_neg )$weight < 0, 
                                        'blue',
                                        'red'))

#modularity - not included in the ms
wtc <- cluster_walktrap(G2)
modularity(wtc)
modularity(G2, membership(wtc))

str(wtc)
#' Now that we have run all of the models and generated predictions, we can examine the proportion of unique observations that were correctly predicted by each of the different models
## ------------------------------------------------------------------------
quantile(mrf.cv$mean_tot_pred, probs = c(0.025, 0.5, 0.975))
quantile(crf.cv$mean_tot_pred, probs = c(0.025, 0.5, 0.975))
quantile(spatial.cv$mean_tot_pred, probs = c(0.025, 0.5, 0.975))
quantile(crf.cv.spatial$mean_tot_pred, probs = c(0.025, 0.5, 0.975))
quantile(as.numeric(ModelPerf_rf[[1]]$ppv), probs = c(0.025, 0.5, 0.975), na.rm=T )
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

}
