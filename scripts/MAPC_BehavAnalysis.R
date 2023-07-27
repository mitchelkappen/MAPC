####################################################################
### MAPC pilot Behavioural Preprocessing ###########################
####################################################################
### MAPC pilot Behavioural Analysis ################################
### Mixed-Effects Models ###########################################
####################################################################
# Code summary for Mitchel

### MAPC pilot Behavioural Preprocessing ###########################


### LOAD LIBRARIES
library(tidyverse) 
library(plyr)
library(dplyr) 
library(magrittr) 
library(readxl)
library(lme4) 
library(gridExtra) 
library(DHARMa)
library(broom.mixed)
library(ggpirate) 
library(gdata) 
library(car) 
library(phia)
library(emmeans)
library(simr) 
library(pbkrtest)
library(lmerTest) 
library(multcomp) 
library(nlme) 
library(ggeffects)
library(splines)
library(stringr)
library(ggplot2)



### IMPORT DATA
# Set directory
setwd ("C:/Users/Bea/Documents/R/MAPC_correlations simulations")
# Load csv file
Data <- read.csv(file = "MAPC_rawBEH_and_simulation_forR.csv", header = TRUE, dec = ",", sep = ";", fill = TRUE, na.strings = ".")

# Create variable for stimulation condition (1= anodal, 2=sham)
for (i in 1:nrow(Data)) {
  if ((Data$group[i] == 1 & Data$session[i] == 1) | (Data$group[i] == 2 & Data$session[i] == 2)) {
    stimCondition <- 2
  } else if ((Data$group[i] == 1 & Data$session[i] == 2) | (Data$group[i] == 2 & Data$session[i] == 1)) {
    stimCondition <- 1
  }
  Data[i, 20] <- stimCondition
}


# Convert data types
Data$subject = as.factor(Data$subject)
Data$AQ = as.numeric(Data$AQ) #autism quotient score
Data$group = as.factor(Data$group) #not relevant by itself, only relevant to know when they got sham and when stimulation
Data$session = as.factor(Data$session)
Data$daysInBetween = as.numeric(Data$daysInBetween) # this are the days that happened in between sessions
Data$pictureCondition = as.factor(Data$pictureCondition)
Data$cancel = as.factor(Data$cancel) # participants get the chance to cancel a trial is they made a mistake and they are aware of it
Data$ACC = as.numeric(Data$ACC) #accuracy
Data$RT1 = as.numeric(Data$RT1) #reaction time until 1st button press (ask me if you want to know more)
Data$RTTot = as.numeric(Data$RTTot) #reaction time until confirm press
Data$delta_ACC = as.numeric(Data$delta_ACC) #not being used currently
Data$delta_RT1 = as.numeric(Data$delta_RT1) #not being used currently
Data$delta_RTTot = as.numeric(Data$delta_RTTot) #not being used currently
Data$roi_efield = as.numeric(Data$roi_efield) #not being used currently
Data$P99_efield = as.numeric(Data$P99_efield) #not being used currently
Data$perc_overlap = as.numeric(Data$perc_overlap) #not being used currently
Data$volume_roi = as.numeric(Data$volume_roi) #not being used currently
Data$volume_P99 = as.numeric(Data$volume_P99) #not being used currently
Data$volume_overlap = as.numeric(Data$volume_overlap) #not being used currently
colnames(Data)[20] = "stimulation"
Data$stimulation = as.factor(Data$stimulation)


### PREPROCESSING


# Filter out cancelled trials
Data_filtered <- Data %>%
  subset(cancel == 0)


## General outlier 

# Outlier function for outlier removal (Lais told me to use this one)
DoubleMAD <- function(x, zero.mad.action="warn"){
  # The zero.mad.action determines the action in the event of an MAD of zero.
  # Possible values: "stop", "warn", "na" and "warn and na".
  x         <- x[!is.na(x)]
  m         <- median(x)
  abs.dev   <- abs(x - m)
  left.mad  <- median(abs.dev[x<=m])
  right.mad <- median(abs.dev[x>=m])
  if (left.mad == 0 || right.mad == 0){
    if (zero.mad.action == "stop") stop("MAD is 0")
    if (zero.mad.action %in% c("warn", "warn and na")) warning("MAD is 0")
    if (zero.mad.action %in% c(  "na", "warn and na")){
      if (left.mad  == 0) left.mad  <- NA
      if (right.mad == 0) right.mad <- NA
    }
  }
  return(c(left.mad, right.mad))
} ## outliers verwijderen niet obv std, verondersteld een normaalverdeling 
DoubleMADsFromMedian <- function(x, zero.mad.action="warn"){
  # The zero.mad.action determines the action in the event of an MAD of zero.
  # Possible values: "stop", "warn", "na" and "warn and na".
  two.sided.mad <- DoubleMAD(x, zero.mad.action)
  m <- median(x, na.rm=TRUE)
  x.mad <- rep(two.sided.mad[1], length(x))
  x.mad[x > m] <- two.sided.mad[2]
  mad.distance <- abs(x - m) / x.mad
  mad.distance[x==m] <- 0
  return(mad.distance)
}

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



### MAPC pilot Behavioural Analysis ################################


### LOAD LIBRARIES
library(optimx) 
library(tidyverse)
library(lme4)
library(mixedup)



### GENERALIZED LINEAR MIXED MODELS ANALYSIS

# Input data -> choose one
inputData = Data_filtered_out
inputData = Data_filtered_outBYpicture # Choose this one

# Formula -> Choose one
formula = "RT1 ~ session*stimulation*pictureCondition + (1|subject)"
formula = "RT1 ~ stimulation*pictureCondition + (1|subject)" #this one gives best results so far
formula = "RT1 ~ stimulation*pictureCondition + (1|subject) + (1|session)"

formula = "RT1 ~ stimulation*pictureCondition + (1+stimulation*pictureCondition|subject)"
formula = "RT1 ~ stimulation*pictureCondition + (1+stimulation*pictureCondition*session|subject) "
formula = "RT1 ~ stimulation*pictureCondition + (1+stimulation+pictureCondition|subject)"

## Models --> Inverse Gaussian Identity link  gives the best results

# Effects Model (Gamma family, log link) 
GLMER_model_gamma_log = glmer(formula,data=inputData, family = Gamma(link = "log"),glmerControl(optimizer= "bobyqa" , optCtrl = list(maxfun = 100000)),nAGQ = 1)
GLMER_model_gamma_log
car::Anova(GLMER_model_gamma_log, type=3)

# Generalized Linear Mixed-Effects Model (Gamma family, identity link)
GLMER_model_gamma_identity = glmer(formula,data=inputData, family = Gamma(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = 1)
GLMER_model_gamma_identity
car::Anova(GLMER_model_gamma_identity, type=3)

# Generalized Linear Mixed-Effects Model (Inverse Gaussian, identity link)
GLMER_model_invGaussian_identity = glmer(formula,data=inputData, family = inverse.gaussian(link = "identity"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = 1)
GLMER_model_invGaussian_identity
car::Anova(GLMER_model_invGaussian_identity, type=3)

#Generalized Linear Mixed-Effects Model (Gaussian, log link)
GLMER_model_gaussian_log = glmer(formula,data=inputData, family = gaussian(link = "log"),glmerControl(optimizer= "bobyqa", optCtrl = list(maxfun = 100000)),nAGQ = 1)
GLMER_model_gaussian_log
car::Anova(GLMER_model_gaussian_log, type=3)

# Linear mixed effects model
lmer_model = lmer(formula,data=inputData)
lmer_model
#car::Anova(lmer_model, type=3)



### Computation of effect sizes ####################################

## Fit the mixed-effects model
inputModel = GLMER_model_invGaussian_identity

## extract the predicted values
predictedResults <- broom.mixed::augment(inputModel)
library("writexl")
write_xlsx(predictedResults,"C:/Users/Bea/Desktop/predictedResults.xlsx")

## Filter data per picture condition
FB_predictedData <- predictedResults[predictedResults$pictureCondition == 3, ] 
TB_predictedData <- predictedResults[predictedResults$pictureCondition == 5, ]
SS_predictedData <- predictedResults[predictedResults$pictureCondition == 2, ] 
ME_predictedData <- predictedResults[predictedResults$pictureCondition == 1, ] 



### Effect of the stimulation on each condition fitted values

## FALSE BELIEF

# Calculate the mean for stimulation and sham
FBmeanFit_stim <- aggregate(.fitted ~ subject, data = FB_predictedData[FB_predictedData$stimulation == 1, ], FUN = mean)
FBmeanFit_sham <- aggregate(.fitted ~ subject, data = FB_predictedData[FB_predictedData$stimulation == 2, ], FUN = mean)
# Merge the mean values for each stimulation level
FBsubject_meanFit_diff <- merge(FBmeanFit_stim, FBmeanFit_sham, by = "subject", suffixes = c("_stim", "_sham"))
FBsubject_meanFit_diff$mean_difference <- FBsubject_meanFit_diff$.fitted_stim - FBsubject_meanFit_diff$.fitted_sham


## TRUE BELIEF

# Calculate the mean for stimulation and sham
TBmeanFit_stim <- aggregate(.fitted ~ subject, data = TB_predictedData[TB_predictedData$stimulation == 1, ], FUN = mean)
TBmeanFit_sham <- aggregate(.fitted ~ subject, data = TB_predictedData[TB_predictedData$stimulation == 2, ], FUN = mean)
# Merge the mean values for each stimulation level
TBsubject_meanFit_diff <- merge(TBmeanFit_stim, TBmeanFit_sham, by = "subject", suffixes = c("_stim", "_sham"))
TBsubject_meanFit_diff$mean_difference <- TBsubject_meanFit_diff$.fitted_stim - TBsubject_meanFit_diff$.fitted_sham


## SOCIAL SCRIPT

# Calculate the mean for stimulation and sham
SSmeanFit_stim <- aggregate(.fitted ~ subject, data = SS_predictedData[SS_predictedData$stimulation == 1, ], FUN = mean)
SSmeanFit_sham <- aggregate(.fitted ~ subject, data = SS_predictedData[SS_predictedData$stimulation == 2, ], FUN = mean)
# Merge the mean values for each stimulation level
SSsubject_meanFit_diff <- merge(SSmeanFit_stim, SSmeanFit_sham, by = "subject", suffixes = c("_stim", "_sham"))
SSsubject_meanFit_diff$mean_difference <- SSsubject_meanFit_diff$.fitted_stim - SSsubject_meanFit_diff$.fitted_sham


## MECHANICAL

# Calculate the mean for stimulation and sham
MEmeanFit_stim <- aggregate(.fitted ~ subject, data = ME_predictedData[ME_predictedData$stimulation == 1, ], FUN = mean)
MEmeanFit_sham <- aggregate(.fitted ~ subject, data = ME_predictedData[ME_predictedData$stimulation == 2, ], FUN = mean)
# Merge the mean values for each stimulation level
MEsubject_meanFit_diff <- merge(MEmeanFit_stim, MEmeanFit_sham, by = "subject", suffixes = c("_stim", "_sham"))
MEsubject_meanFit_diff$mean_difference <- MEsubject_meanFit_diff$.fitted_stim - MEsubject_meanFit_diff$.fitted_sham

###############################################################################

### Effect of the stimulation on each condition residuals

## FALSE BELIEF

# Calculate the mean for stimulation and sham
FBmeanRes_stim <- aggregate(.resid ~ subject, data = FB_predictedData[FB_predictedData$stimulation == 1, ], FUN = mean)
FBmeanRes_sham <- aggregate(.resid ~ subject, data = FB_predictedData[FB_predictedData$stimulation == 2, ], FUN = mean)
# Merge the mean values for each stimulation level
FBsubject_meanRes_diff <- merge(FBmeanRes_stim, FBmeanRes_sham, by = "subject", suffixes = c("_stim", "_sham"))
FBsubject_meanRes_diff$mean_difference <- FBsubject_meanRes_diff$.resid_stim - FBsubject_meanRes_diff$.resid_sham


## TRUE BELIEF

# Calculate the mean for stimulation and sham
TBmeanRes_stim <- aggregate(.resid ~ subject, data = TB_predictedData[TB_predictedData$stimulation == 1, ], FUN = mean)
TBmeanRes_sham <- aggregate(.resid ~ subject, data = TB_predictedData[TB_predictedData$stimulation == 2, ], FUN = mean)

# Merge the mean values for each stimulation level
TBsubject_meanRes_diff <- merge(TBmeanRes_stim, TBmeanRes_sham, by = "subject", suffixes = c("_stim", "_sham"))
TBsubject_meanRes_diff$mean_difference <- TBsubject_meanRes_diff$.resid_stim - TBsubject_meanRes_diff$.resid_sham


## SOCIAL SCRIPT

# Calculate the mean for stimulation and sham
SSmeanRes_stim <- aggregate(.resid ~ subject, data = SS_predictedData[SS_predictedData$stimulation == 1, ], FUN = mean)
SSmeanRes_sham <- aggregate(.resid ~ subject, data = SS_predictedData[SS_predictedData$stimulation == 2, ], FUN = mean)
# Merge the mean values for each stimulation level
SSsubject_meanRes_diff <- merge(SSmeanRes_stim, SSmeanRes_sham, by = "subject", suffixes = c("_stim", "_sham"))
SSsubject_meanRes_diff$mean_difference <- SSsubject_meanRes_diff$.resid_stim - SSsubject_meanRes_diff$.resid_sham


## MECHANICAL

# Calculate the mean for stimulation and sham
MEmeanRes_stim <- aggregate(.resid ~ subject, data = ME_predictedData[ME_predictedData$stimulation == 1, ], FUN = mean)
MEmeanRes_sham <- aggregate(.resid ~ subject, data = ME_predictedData[ME_predictedData$stimulation == 2, ], FUN = mean)
# Merge the mean values for each stimulation level
MEsubject_meanRes_diff <- merge(MEmeanRes_stim, MEmeanRes_sham, by = "subject", suffixes = c("_stim", "_sham"))
MEsubject_meanRes_diff$mean_difference <- MEsubject_meanRes_diff$.resid_stim - MEsubject_meanRes_diff$.resid_sham


