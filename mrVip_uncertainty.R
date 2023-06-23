#'Wrapper to estimate and plot model-agnostic variable importance with a measure of uncertainty for multi-response models. 
#'@param Y A \code{dataframe} is response variable data (species, OTUs, SNPs etc).
#'@param X A \code{dataframe} represents predictor or feature data.
#'@param dummy A \code{logical} 'TRUE or FALSE'. 
#'@param Model 1 A \code{list} can be any model from the tidy model package. See examples.
#'@param tune_grid_size A \code{numeric} sets the grid size for hyperparamter tuning. Larger grid sizes increase computational time.
#'@param k A \code{numeric} sets the number of folds in the 10-fold cross-validation. 10 is the default.
#'@seed A \code{numeric} as these models have a stochastic component, a seed is set to make to make the analysis reproducible. Defaults between 100 million and 1.
#'@param mode \code{character}'classification' or 'regression' i.e., is the generative model a regression or classification?

#'@details Calculates variable importance across multiple runs to quantify uncertainty in model estimates.
#'# this is particularly useful for smaller unbalanced data sets where the vip measure of variable
#'#importance can be a bit unstable.
#' @example 
#' # test <- mrVip_mutlirun(X=X,
#'                           Y=Y,
#'                           ity=4, #runs the model 4 times and summarizes. 
#'                          Model=model_rf,
#'                          mode='classification', 
#'                           tune_grid_size=5,
#'                          seed = sample.int(1e8, 1))#'@export 

mrVip_uncertainty<- function (X, Y, Model,
                            mode='classification',ity=5,
                            tune_grid_size= 10, #removed k fold cross validation
                            seed = sample.int(1e8, 1) ) { 

    
    internal_fit_function <- function( i ){
      
    yhats <- mrIMLpredicts(X,Y,
                              Model,
                              mode, 
                              tune_grid_size=tune_grid_size, balance_data='no',
                             seed = sample.int(1e8, 1)) #this will make sure dif seeds are used
    
    VI <- mrVip(yhats, X)   
    
    }
    
    im <- lapply(seq(1,ity), FUN=internal_fit_function)
    #doesn't need futures as each mrIMLpredicts function uses it.
    

    ImpGlobal <- as.data.frame(do.call(rbind, im))
    
  
# 
#     colmeans_function <- function( i ){
#       colMeans(matrix(i,nrow=length(Y)))
#     }
#     
#     
#     
#     vimp_summary <-  ImpGlobal %>%  apply(2, colmeans_function) %>% 
#       as.data.frame()
    
    

    #make wide to long
    # data_vi <-  vimp_summary  %>%
    #   gather(names(X), key = Xvar, value = importance)
    # 
    
    #boxplot
    # data_vi %>%
    #   ggplot(aes(x=reorder(Xvar,importance), y=importance)) +
    #   geom_boxplot(fill = 'lightgrey')+
    #   coord_flip()+
    #   theme_bw()+
    #   labs(x="Features", y="Importance")
    #
     # scale_fill_brewer(palette="BuPu")#from cbind

   #individual responses (example)

    # namesY <- rep(names(Y), ity)
    #
    # df_namesUpdated  <-  cbind (Ynames=namesY, ImpGlobal )
    #
    # Y39 <-  df_namesUpdated %>% filter (Ynames == 'Borreliella')
    #
    # data_vi39 <-  Y39  %>%
    #   gather(names(X), key = Xvar, value = importance)
    #
    # data_vi39 %>%
    #   ggplot(aes(x=reorder(Xvar,importance), y=importance)) +
    #   geom_boxplot(fill = 'lightgrey')+
    #   coord_flip()+
    #   theme_bw()+
    #   labs(x="Features", y="Importance (Borreliella)")


}