

#' Read in the raw input file. These can be BIOM file from QIIME and convert to a matrix

library(OTUtable)
library(tidyverse)


########################################################################################################################
#---------------------Data preparation--------------------------
########################################################################################################################

#import data
envdat <- read.csv('tick_env_final_data_namesFixedMar2021.csv', head=T)
envdatSort <- envdat[order(envdat$Site_Code),]

asvdatFiltered <- as.data.frame(read.csv('ASVtableSiteUpdated.csv', head=T)) 
glimpse(asvdatFiltered )

asvdatFiltered  <- asvdatFiltered %>% mutate(RicketsRatio =Rickettsia/Total.Reads)
#glimpse(contamIssue)

RicketsiaRatio <- asvdatFiltered %>% select(Sequncing.File.Name,Rickettsia, RicketsRatio)

ggplot(RicketsiaRatio, aes(x=RicketsRatio))+
  geom_histogram(bins = 30)+
 scale_x_continuous(breaks=seq(0,1,0.1))+
  theme_bw()

  
asvdatFiltered %>% 
  filter(!Site.Code.Name...Year == "BRACHYSPRIA") %>% #splatter control 
  filter(!Site.Code.Name...Year == "BLANK") %>%  #blanks
  filter(!Site.Code.Name...Year == "WL19") %>% #site not sampled in the end
  filter(!Site.Code.Name...Year == "WT17")  #singleton

glimpse(asvdatFiltered)

#remove ticks with < 1000 reads
asvdatFiltered <- asvdatFiltered %>%  filter(Total.Reads >1000)

#remove the soil as we dont have matching environmental data
asvdataCleanNoSoil <-  filter(asvdatFiltered  , !Tick.Sex.Sample.Type == "S") 

finalASVdat <- asvdataCleanNoSoil  [order(asvdataCleanNoSoil $Site.Code.Name...Year),]

countSIte  <-finalASVdat %>% 
  group_by(Site.Code.Name...Year) %>%
  summarise(no_rows = length(Site.Code.Name...Year))

countSIte  <-finalASVdat %>% 
  group_by(Site.Code.Name...Year) %>%
  summarise(no_rows = length(Site.Code.Name...Year))

#check matching
match(countSIte$Site.Code.Name...Year, envdatSort$Site_Code)

#repeat rows of environmental data to match ASV data

envData_reps <- envdatSort[rep(seq_len(nrow(envdatSort)), countSIte$no_rows),]

#add sequence id to make sure everything is the same
envData_repsID <- cbind(finalASVdat$Sequncing.File.Name, envData_reps)

#add other relevant data from the ASV data fram to the environmental data

envData_final <- data.frame(envData_repsID, sex=finalASVdat$Tick.Sex.Sample.Type)

glimpse(envData_final)
saveRDS(envData_final, 'complete_data')

#--------------------------------------------------------------------------------------------------
#create random effect data (individual, plate, collection year and site and spatial data)

#as factor data frame
randonEff_f <- data.frame( plate=finalASVdat$Plate.Code, year=finalASVdat$Sample.Collection.Year, 
                           site=finalASVdat$Site.Code.Only)

glimpse(randonEff_f )

saveRDS(randonEff_f
        , 'randonEff ')
#--------------------------------------------------------------------------------------------------

#streamline data

#reduce to predictor set and test for correlations
envData_final_predictors <- envData_final[11:128]


#--------------------------------------------------------------------------------------------------    
#Predictor summary and correlations
#-------------------------------------------------------------------------------------------------- 
# library(DataExplorer)
# create_report(envData_final_predictors ) #produces an html report

#not too much missing data - lets impute!
library(missForest)

envData_final_noNA <- missForest(envData_final_predictors, variablewise=T) #default values have worked fine previously

envData_final_noNA$OOBerror #all really low mse

envData_final_noNA <- envData_final_noNA$ximp

#remove spatial data
envData_final_noNA$latitude <- NULL
envData_final_noNA$longitude <- NULL

#remove correlated predictors

#tidymodel way  -adapted! removes less predictors than caret version. 16 predictors at 0.7 cor coefficent threshold

library(GGally)

ggpairs(envData_final_noNA)
cor <-  cor(envData_final_noNA)

library(tidymodels)

#approach one for removing correlated variables

# rec <- recipe( vpdmax   ~ .,
#                data = envData_final_noNA, retain = TRUE)
# 
# corr_filter <- rec %>%
#   step_corr(all_predictors(), threshold = .6)
# 
# filter_obj <- prep(corr_filter, training = envData_final_noNA)
# 
# envData_uncor<- bake(filter_obj, envData_final_noNA) #extract final data
# 
# corCheck <- cor(envData_uncor)
# 
# envData_uncor_simp <- envData_uncor
# envData_uncor_simp[13:19] <- NULL
# envData_uncor_WC <-cbind(envData_uncor_simp, sex=envData_final$sex)

# removing correlated variables  - biological focus with no seasonal variables. When correlated mean variables favoured

envData_uncor_simp <- envData_final_noNA %>%
  select(tmean, vpdmax, vpdmin, prev_30_day_avg_tmean , prev_30_day_avg_tmean,
         prev_30_day_avg_vpdmin, prev_year_avg_tmean, prev_year_avg_vpdmax,
         prev_year_avg_vpdmin, prev_2years_avg_tmean   , prev_2years_avg_vpdmax, prev_2years_avg_vpdmin,
         prev_5years_avg_tmean, prev_5years_avg_vpdmax, prev_5years_avg_vpdmin,
         max_snow_depth , ndvi_buff, ndvi_point
  )
corr_simp <- cor(envData_uncor_simp)

#VPD max strongly correlated with temp
# NDVi point strongly correlated with buffer. Buffer more correlated with other variables so dropped
#5 year previous vpd min correlated with 2 year average vpd min
#VPD min previos 2/5 years corr with temp so removed

#0.7 corr for temp across years but kept in to test for potential historical effects

envData_uncor_simp <- envData_final_noNA %>%
  select(tmean, vpdmin, #previous 30 day tmean correlated with tmean
         prev_year_avg_tmean, #previous 2 years tmean correlated with previous yera mean
         prev_5years_avg_tmean,
         max_snow_depth ,  ndvi_point
  )
corr_simp <- as.data.frame(cor(envData_uncor_simp))

env_data_withHost <-  cbind(envData_uncor_simp , sex= envData_final$sex) 


#' Extract coordinates columns to use for incorporating spatial splines as covariates in the model (for now, these need to be named `Latitude` and `Longitude`)
## ----message=FALSE, warning=FALSE----------------------------------------
Latitude <- envData_final$latitude
Longitude <- envData_final$longitude

coords <- (data.frame(Latitude = Latitude,
                      Longitude = Longitude))

saveRDS(coords , 'coords')

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

# 166 reads of Coryebacterium.1 occurred the extraction null - remove
asvFilter_230asv$Brachyspira <- NULL

#make presence/absence
asvFilter_230asv[asvFilter_230asv>0] <- 1  

#remove these columns
asvFilter_230asv$Total.Reads <- NULL
asvFilter_230asv$NA. <- NULL

#-------------------------------------------------------------------------------------------------- 
#ave files
saveRDS(asvFilter_230asv, 'asvFilter_230asv')
saveRDS(env_data_withHost , 'env_data_withHost')