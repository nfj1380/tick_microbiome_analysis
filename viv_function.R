if (!requireNamespace("graph", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("graph")
}
install.packages("zenplots")

install.packages("vivid")

library(vivid) # for visualisations 
library(tidymodels)


data_df <-  yhats_rf[[39]]$data
fit = extract_fit_parsnip(yhats_rf[[39]]$mod1_k)


#needs X var to be non catergorical

t <- yhats_rf[[39]]$data_testa
str(fit)
str(data_df)

library(randomForest)
set.seed(121)
viFit <- vivi(
  fit = extract_fit_parsnip(yhats_rf[[39]]$mod1_k),
  data = data_df,
  response = 'class',
  gridSize = 10)#,

viviHeatmap(mat = viFit) + ggtitle("rf heatmap")
viviNetwork(mat = viFit)

str(X)
set.seed(1701)
pdpPairs(data = data_df, fit = fit, response = "class", nmax = 50,  gridSize = 30, #$class=1,
         convexHull = T, probability=T, vars = c('max_snow_depth'  , 'ndvi_point', 'prev_year_avg_tmean', 'vpdmin'))

  #3  importanceType = 'agnostic')#,
 # predictF un = pred_fun )


