####################################################################
### MAPC pilot Behavioural Preprocessing ###########################
####################################################################
### MAPC pilot Behavioural Analysis ################################
### Mixed-Effects Models ###########################################
####################################################################
# Re-edit by Mitchel Kappen 
# 27/07/2023

### MAPC pilot Behavioural Preprocessing ###########################


### LOAD LIBRARIES
library(tidyverse) 
library(plyr)
library(dplyr) 
# library(magrittr) 
# library(readxl)
library(lme4)
# library(gridExtra) 
# library(DHARMa)
library(broom.mixed)
# library(ggpirate) 
# library(gdata) 
library(car) 
# library(phia)
library(emmeans)
# library(simr) 
# library(pbkrtest)
library(lmerTest) 
# library(multcomp) 
# library(nlme) 
library(ggeffects)
# library(splines)
# library(stringr)
library(ggplot2)
library(writexl)

##### Set environment #####
rm(list = ls()) # Clear environment
cat("\014") # Clear console # # Or ctrl + l in VSCode
dev.off() # Clear plot window

# Set and Get directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #Set WD to script location
source("functions.R") # Load document where functions are stored
options(contrasts = c("contr.sum","contr.poly")) #use this for the p value of the t test
plotDirectory <- "../figures/"

nAGQ = 1
##### IMPORT DATA #####
# Load csv file
Data <- as.data.frame(read.csv("../data/MAPC_rawBEH_and_simulation_forR.csv", header = TRUE, dec = ",", sep = ";", fill = TRUE, na.strings = "."))
names(Data)[names(Data) == "ï..subject"] <- "subject" # something weird with loading in, so needs to be renamed

# Create variable for stimCondition condition (1= anodal, 2=sham)
Data$stimCondition <- NA  # Initialize a new column
for (i in 1:nrow(Data)) {
  if ((Data$group[i] == 1 & Data$session[i] == 1) | (Data$group[i] == 2 & Data$session[i] == 2)) {
    Data$stimCondition[i] <- 'sham'
  } else if ((Data$group[i] == 1 & Data$session[i] == 2) | (Data$group[i] == 2 & Data$session[i] == 1)) {
    Data$stimCondition[i] <- 'active'
  }
}

# Convert data types
factor_cols <- c("subject", "group", "session", "pictureCondition", "cancel", "stimCondition")
numeric_cols <- c("AQ", "daysInBetween", "ACC", "RT1", "RTTot", "delta_ACC", "delta_RT1", "delta_RTTot", 
                  "roi_efield", "P99_efield", "perc_overlap", "volume_roi", "volume_P99", "volume_overlap")

# Convert to factors
Data[factor_cols] <- lapply(Data[factor_cols], as.factor)

# Convert to numeric
Data[numeric_cols] <- lapply(Data[numeric_cols], as.numeric)

# Data clean up - create a minimalist dataframe
Data_backup <- Data # Backup
# List of variables not currently used
unused_vars <- c("delta_ACC", "delta_RT1", "delta_RTTot", "roi_efield", "P99_efield",
                 "perc_overlap", "volume_roi", "volume_P99", "volume_overlap")
# Create a new dataframe without the unused variables
Data <- Data[, !(names(Data) %in% unused_vars)]

##### PREPROCESSING #####
# Filter out cancelled trials
Data_filtered <- Data %>%
  subset(cancel == 0)

## General outlier 
#removing  general outliers --> not preferred due to differences in picture conditions
Data_filtered_out <- Data_filtered[DoubleMADsFromMedian(Data_filtered$RT1) <=3, ] #3 is how many SD from median, can be checked other values

# Create subsets by picture condition
Data_filtered_FB = Data_filtered%>%
  subset(pictureCondition == "3") # false belief
Data_filtered_TB = Data_filtered%>%
  subset(pictureCondition == "5") # true belief
Data_filtered_SS = Data_filtered%>%
  subset(pictureCondition == "2") # social script
Data_filtered_ME = Data_filtered%>%
  subset(pictureCondition == "1") # mechanical
# Outlier removal by picture condition --> preferred option
Data_filtered_FB_out <- Data_filtered_FB[DoubleMADsFromMedian(Data_filtered_FB$RT1) <=3, ] 
Data_filtered_TB_out <- Data_filtered_TB[DoubleMADsFromMedian(Data_filtered_TB$RT1) <=3, ] 
Data_filtered_SS_out <- Data_filtered_SS[DoubleMADsFromMedian(Data_filtered_SS$RT1) <=3, ] 
Data_filtered_ME_out <- Data_filtered_ME[DoubleMADsFromMedian(Data_filtered_ME$RT1) <=3, ] 
# Combine the data sets back together
Data_filtered_outBYpicture = rbind(Data_filtered_FB_out, Data_filtered_TB_out, Data_filtered_SS_out, Data_filtered_ME_out )

# For easier interpretability down the line, I am renaming the pictureCondition variable
# Define named vector of new levels
new_levels <- c("1" = "ME", "2" = "SS", "3" = "FB", "5" = "TB")

# Apply new levels
Data_filtered_outBYpicture$pictureCondition <- revalue(Data_filtered_outBYpicture$pictureCondition, new_levels)

# clean up envirnonment a bit
rm(Data_filtered_FB, Data_filtered_FB_out, Data_filtered_TB, Data_filtered_TB_out, 
   Data_filtered_SS, Data_filtered_SS_out, Data_filtered_ME, Data_filtered_ME_out)
### MAPC pilot Behavioural Analysis ################################
### LOAD LIBRARIES
# library(optimx) 
# library(tidyverse)
# library(lme4)
# library(mixedup)

##### GENERALIZED LINEAR MIXED MODELS ANALYSIS ##### 

# Input data -> choose one
# inputData = Data_filtered_out
inputData = Data_filtered_outBYpicture # Choose this one

# Formula -> Choose one
formula = "RT1 ~ session*stimCondition*pictureCondition + (1|subject)"
formula = "RT1 ~ stimCondition*pictureCondition + (1|subject)" #this one gives best results so far
formula = "RT1 ~ stimCondition*pictureCondition + (1|subject) + (1|session)"

formula = "RT1 ~ stimCondition*pictureCondition + (1+stimCondition*pictureCondition|subject)"
formula = "RT1 ~ stimCondition*pictureCondition + (1+stimCondition*pictureCondition*session|subject) "
formula = "RT1 ~ stimCondition*pictureCondition + (1+stimCondition+pictureCondition|subject)"

## Models --> Inverse Gaussian Identity link  gives the best results
### Mitch stuff
d0.1 <- lmer(formula,data=inputData)
d0.2 <- glmer(formula,data=inputData, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
# d0.3 <- glmer(formula,data=inputData, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
# d0.4 <- glmer(formula,data=inputData, family = gaussian(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)
d0.5 <- glmer(formula,data=inputData, family = Gamma(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = nAGQ)

# Model Selection
modelNames = c(d0.1,d0.2,d0.5)
tabel <- cbind(AIC(d0.1), AIC(d0.2), AIC(d0.5))
chosenModel = modelNames[which(tabel == min(tabel))] # Get model with lowest AIC

Anova(chosenModel[[1]], type = 'III')
plot(effect("stimCondition:pictureCondition", chosenModel[[1]]))

emmeans0.1 <- emmeans(chosenModel[[1]], pairwise ~ stimCondition | pictureCondition, adjust ="none", type = "response") #we don't adjust because we do this later
emmeans0.2 <- emmeans(chosenModel[[1]], pairwise ~ pictureCondition | stimCondition, adjust ="none", type = "response") #we don't adjust because we do this later
emm0.1 <- summary(emmeans0.1)$emmeans
emmeans0.1$contrasts

emmeans_interaction <- emmeans(chosenModel[[1]], ~ stimCondition:pictureCondition, type = "response")
pairwise_comparisons <- pairs(emmeans_interaction, adjust = "fdr")

## Plotting
# figure = behaviorplot(emm0.1, stimCondition, pictureCondition, "Stim X Pic") # Create plot
# figure = addpvalues(figure, emmeans0.1)
# figure = addpvaluesBetween(figure, emmeans0.2)
# savePlot(figure, "StimxPic") # Display and save plot

## plotting 2
inputData$interaction <- interaction(inputData$pictureCondition, inputData$stimCondition)

ggstatsplot::ggbetweenstats(data = inputData, 
                            x = interaction, 
                            y = RT1, 
                            title = "Interaction effect")

### Computation of effect sizes ####################################

## Fit the mixed-effects model
inputModel = chosenModel[[1]]

## extract the predicted values
predictedResults <- broom.mixed::augment(inputModel)
write_xlsx(predictedResults,"../output/predictedResults.xlsx")

## Filter data per picture condition
FB_predictedData <- predictedResults[predictedResults$pictureCondition == 'FB', ] 
TB_predictedData <- predictedResults[predictedResults$pictureCondition == 'TB', ]
SS_predictedData <- predictedResults[predictedResults$pictureCondition == 'SS', ] 
ME_predictedData <- predictedResults[predictedResults$pictureCondition == 'ME', ] 

### Effect of the stimCondition on each condition fitted values
c("1" = "ME", "2" = "SS", "3" = "FB", "5" = "TB")
## FALSE BELIEF
# Calculate the mean for stimCondition and sham
FBmeanFit_stim <- aggregate(.fitted ~ subject, data = FB_predictedData[FB_predictedData$stimCondition == 'active', ], FUN = mean)
FBmeanFit_sham <- aggregate(.fitted ~ subject, data = FB_predictedData[FB_predictedData$stimCondition == 'sham', ], FUN = mean)
# Merge the mean values for each stimCondition level
FBsubject_meanFit_diff <- merge(FBmeanFit_stim, FBmeanFit_sham, by = "subject", suffixes = c("_stim", "_sham"))
FBsubject_meanFit_diff$mean_difference <- FBsubject_meanFit_diff$.fitted_stim - FBsubject_meanFit_diff$.fitted_sham


## TRUE BELIEF

# Calculate the mean for stimCondition and sham
TBmeanFit_stim <- aggregate(.fitted ~ subject, data = TB_predictedData[TB_predictedData$stimCondition == 'active', ], FUN = mean)
TBmeanFit_sham <- aggregate(.fitted ~ subject, data = TB_predictedData[TB_predictedData$stimCondition == 'sham', ], FUN = mean)
# Merge the mean values for each stimCondition level
TBsubject_meanFit_diff <- merge(TBmeanFit_stim, TBmeanFit_sham, by = "subject", suffixes = c("_stim", "_sham"))
TBsubject_meanFit_diff$mean_difference <- TBsubject_meanFit_diff$.fitted_stim - TBsubject_meanFit_diff$.fitted_sham


## SOCIAL SCRIPT

# Calculate the mean for stimCondition and sham
SSmeanFit_stim <- aggregate(.fitted ~ subject, data = SS_predictedData[SS_predictedData$stimCondition == 'active', ], FUN = mean)
SSmeanFit_sham <- aggregate(.fitted ~ subject, data = SS_predictedData[SS_predictedData$stimCondition == 'sham', ], FUN = mean)
# Merge the mean values for each stimCondition level
SSsubject_meanFit_diff <- merge(SSmeanFit_stim, SSmeanFit_sham, by = "subject", suffixes = c("_stim", "_sham"))
SSsubject_meanFit_diff$mean_difference <- SSsubject_meanFit_diff$.fitted_stim - SSsubject_meanFit_diff$.fitted_sham


## MECHANICAL

# Calculate the mean for stimCondition and sham
MEmeanFit_stim <- aggregate(.fitted ~ subject, data = ME_predictedData[ME_predictedData$stimCondition == 'active', ], FUN = mean)
MEmeanFit_sham <- aggregate(.fitted ~ subject, data = ME_predictedData[ME_predictedData$stimCondition == 'sham', ], FUN = mean)
# Merge the mean values for each stimCondition level
MEsubject_meanFit_diff <- merge(MEmeanFit_stim, MEmeanFit_sham, by = "subject", suffixes = c("_stim", "_sham"))
MEsubject_meanFit_diff$mean_difference <- MEsubject_meanFit_diff$.fitted_stim - MEsubject_meanFit_diff$.fitted_sham

###############################################################################

### Effect of the stimCondition on each condition residuals

## FALSE BELIEF

# Calculate the mean for stimCondition and sham
FBmeanRes_stim <- aggregate(.resid ~ subject, data = FB_predictedData[FB_predictedData$stimCondition == 'active', ], FUN = mean)
FBmeanRes_sham <- aggregate(.resid ~ subject, data = FB_predictedData[FB_predictedData$stimCondition == 'sham', ], FUN = mean)
# Merge the mean values for each stimCondition level
FBsubject_meanRes_diff <- merge(FBmeanRes_stim, FBmeanRes_sham, by = "subject", suffixes = c("_stim", "_sham"))
FBsubject_meanRes_diff$mean_difference <- FBsubject_meanRes_diff$.resid_stim - FBsubject_meanRes_diff$.resid_sham


## TRUE BELIEF

# Calculate the mean for stimCondition and sham
TBmeanRes_stim <- aggregate(.resid ~ subject, data = TB_predictedData[TB_predictedData$stimCondition == 'active', ], FUN = mean)
TBmeanRes_sham <- aggregate(.resid ~ subject, data = TB_predictedData[TB_predictedData$stimCondition == 'sham', ], FUN = mean)

# Merge the mean values for each stimCondition level
TBsubject_meanRes_diff <- merge(TBmeanRes_stim, TBmeanRes_sham, by = "subject", suffixes = c("_stim", "_sham"))
TBsubject_meanRes_diff$mean_difference <- TBsubject_meanRes_diff$.resid_stim - TBsubject_meanRes_diff$.resid_sham


## SOCIAL SCRIPT

# Calculate the mean for stimCondition and sham
SSmeanRes_stim <- aggregate(.resid ~ subject, data = SS_predictedData[SS_predictedData$stimCondition == 'active', ], FUN = mean)
SSmeanRes_sham <- aggregate(.resid ~ subject, data = SS_predictedData[SS_predictedData$stimCondition == 'sham', ], FUN = mean)
# Merge the mean values for each stimCondition level
SSsubject_meanRes_diff <- merge(SSmeanRes_stim, SSmeanRes_sham, by = "subject", suffixes = c("_stim", "_sham"))
SSsubject_meanRes_diff$mean_difference <- SSsubject_meanRes_diff$.resid_stim - SSsubject_meanRes_diff$.resid_sham


## MECHANICAL

# Calculate the mean for stimCondition and sham
MEmeanRes_stim <- aggregate(.resid ~ subject, data = ME_predictedData[ME_predictedData$stimCondition == 'active', ], FUN = mean)
MEmeanRes_sham <- aggregate(.resid ~ subject, data = ME_predictedData[ME_predictedData$stimCondition == 'sham', ], FUN = mean)
# Merge the mean values for each stimCondition level
MEsubject_meanRes_diff <- merge(MEmeanRes_stim, MEmeanRes_sham, by = "subject", suffixes = c("_stim", "_sham"))
MEsubject_meanRes_diff$mean_difference <- MEsubject_meanRes_diff$.resid_stim - MEsubject_meanRes_diff$.resid_sham


