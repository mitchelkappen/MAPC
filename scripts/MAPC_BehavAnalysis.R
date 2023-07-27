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
library(lme4)
library(broom.mixed)
library(car) 
library(emmeans)
library(lmerTest) 
library(ggeffects)
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
pvalues = c() # Create a variable to store all p-values to correct later

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
##### GENERALIZED LINEAR MIXED MODELS ANALYSIS ##### 
# Input data for models
inputData = Data_filtered_outBYpicture # Choose this one

# Formula -> Choose one
# formula = "RT1 ~ session*stimCondition*pictureCondition + (1|subject)"
formula = "RT1 ~ stimCondition*pictureCondition + (1|subject)" #this one gives best results so far
# formula = "RT1 ~ stimCondition*pictureCondition + (1|subject) + (1|session)"
# 
# formula = "RT1 ~ stimCondition*pictureCondition + (1+stimCondition*pictureCondition|subject)"
# formula = "RT1 ~ stimCondition*pictureCondition + (1+stimCondition*pictureCondition*session|subject) "
# formula = "RT1 ~ stimCondition*pictureCondition + (1+stimCondition+pictureCondition|subject)"

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
pvalues  = append(pvalues ,summary(emmeans0.1$contrasts)$p.value) # Add pvalues to dataframe to later correct all

# emmeans_interaction <- emmeans(chosenModel[[1]], ~ stimCondition:pictureCondition, type = "response")
# pairwise_comparisons <- pairs(emmeans_interaction, adjust = "fdr")

## Plotting #####
## Get significant p-values
emmeans0.1$contrasts # just ME 1vs2

## Plot
pd = 0.2 #position dodge
p <- ggplot(emm_df, aes(x = pictureCondition, y = response, group = stimCondition)) +
  geom_point(aes(color = stimCondition), size = 3, position = position_dodge(pd)) +
  geom_errorbar(aes(ymin = response - SE, ymax = response + SE, color = stimCondition), width = 0.2, position = position_dodge(pd)) +
  labs(x = "Picture Condition", 
       y = "Estimated Marginal Means", 
       title = "Estimated Marginal Means for each Picture Condition by Stim Condition",
       color = "Stim Condition") +
  scale_colour_manual(values=cbPalette) +
  plot_theme_apa() +
  theme(legend.position = "bottom") #Change the position of the legend
p <- p + annotate(geom="text", x= 1.2, y= 14396, label= '*', size = 10) # Add the star to the plot## Significance

savePlot(p, "StimxPic") # Display and save plot

## Different plot #####
library(yarrr)

# Create the pirate plot
pirateplot(formula = RT1 ~ stimCondition * pictureCondition,
           data = inputData,
           pal = cbPalette,  # use custom colors
           xlab = "Stim Condition and Picture Condition",
           ylab = "RT1",
           main = "Comparison of RT1 across Stim Conditions and Picture Conditions",
           theme = 2  # try a different theme
)

##### Computation of effect sizes ######

##### Mitchel Effect Size #####
# Create empty forestplot df
forestdf <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("Outcome" = character(0), "D" = numeric(0), "Lower" = numeric(0), "Upper" = numeric(0), "Group" = character(0)))
forestdf = data.frame(Outcome=character(0), Group=character(0),  effectsize=numeric(0), Lower=numeric(0), Upper=numeric(0), SE=numeric(0))

# Get effect sizes
effSummary <- summary(eff_size(emmeans0.1, sigma=sigma(chosenModel[[1]]), edf=df.residual(chosenModel[[1]])))
# Cohen's D for Forest Plots
for(i in 1:length(effSummary$pictureCondition )){
  name = as.character(effSummary$pictureCondition[i])
  effectsize = effSummary$effect.size[i] * -1 # Inverted
  Upper = effSummary$asymp.LCL[i] * -1 # Inverted
  Lower = effSummary$asymp.UCL[i] * -1 # Inverted
  
  contrastdf = summary(emmeans0.1$contrasts) # get contrasts
  Beta = contrastdf$estimate[i] * -1 # Inverted
  SE = contrastdf$SE[i]
  t = contrastdf$t.ratio[i] * -1 # Inverted
  forestdf[nrow(forestdf) + 1,] = c(name, as.character(effSummary$pictureCondition[i]), effectsize, Lower, Upper, SE)
  print(forestdf)
}

# add p valuse
# Correct P values SPEECH ######
names = c("ME", "SS", "FB", "TB")
ps = list()
ps[names] = p.adjust(pvalues, method = "fdr", length(pvalues)) # Create list containing fdr corrected pvalues
collectedPvalues = ps

# Make forest
# you can do the factoring here
forestdf$Outcome = factor(forestdf$Outcome, levels = c("ME", "SS", "FB", "TB"))

# Make them numeric
forestdf$effectsize = round(as.numeric(forestdf$effectsize), digits = 2)
forestdf$Lower = round(as.numeric(forestdf$Lower), digits = 3)
forestdf$Upper = round(as.numeric(forestdf$Upper), digits = 3)
forestdf$SE = round(as.numeric(forestdf$SE), digits = 3)

forestdf$pvalues = round(as.numeric(collectedPvalues), digits = 3)

#define colours for dots and bars
# dotCOLS = c("#a6d8f0","#f9b282") # these are actually the bars
barCOLS = c("#56B4E9", "#E69F00", "#56B4E9", "#E69F00")
dotCOLS = c("#56B4E9", "#E69F00", "#56B4E9", "#E69F00")
boxlims = c(0.5, 6.5, 10.5, 12.5)

dodgevar = 0.5
##### Make Forest plot ####
forestplot <- ggplot(forestdf, aes(x=Outcome, y=effectsize, ymin=Upper, ymax=Lower,col=Group,fill=Group, group=Group)) + 
  # Draw forestplot
  geom_linerange(size=8,position=position_dodge(width = dodgevar), alpha = ifelse(forestdf$pvalues < .05, 1, 0.5)) +
  geom_hline(yintercept=0, lty=2, size = 1.5) + # Draw vertical 0 line
  
  # Create dots for effect sizes
  geom_point(size=4, shape=21, colour= ifelse(forestdf$pvalues < .05, "black", "white"), alpha = ifelse(forestdf$pvalues < .05, 1, 0.5),
             stroke = 1.4, position=position_dodge(width = dodgevar)) +
  
  # Set bar and dot colors
  scale_fill_manual(values=barCOLS)+
  scale_color_manual(values=dotCOLS)+
  
  # Set figure specifics
  scale_x_discrete(name="") +
  scale_y_continuous(limits = c(-.6, .8)) +
  
  # Set orientation and theme
  coord_flip()+
  theme_pubr() +
  plot_theme_apa()+
  ylab("Effect Size (Cohen's D) with 95% CI") + # Plot is flipped, this is actually the x-axis
  theme(legend.position = "bottom", legend.text = element_text(size = 18), legend.title = element_text(size = 18))

forestplot
savePlot(forestplot, "forestPlot") # Display and save plot

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

##### Effect of the stimCondition on each condition fitted values #####
## FALSE BELIEF ####
# Calculate the mean for stimCondition and sham
FBmeanFit_stim <- aggregate(.fitted ~ subject, data = FB_predictedData[FB_predictedData$stimCondition == 'active', ], FUN = mean)
FBmeanFit_sham <- aggregate(.fitted ~ subject, data = FB_predictedData[FB_predictedData$stimCondition == 'sham', ], FUN = mean)
# Merge the mean values for each stimCondition level
FBsubject_meanFit_diff <- merge(FBmeanFit_stim, FBmeanFit_sham, by = "subject", suffixes = c("_stim", "_sham"))
FBsubject_meanFit_diff$mean_difference <- FBsubject_meanFit_diff$.fitted_stim - FBsubject_meanFit_diff$.fitted_sham


## TRUE BELIEF ####

# Calculate the mean for stimCondition and sham
TBmeanFit_stim <- aggregate(.fitted ~ subject, data = TB_predictedData[TB_predictedData$stimCondition == 'active', ], FUN = mean)
TBmeanFit_sham <- aggregate(.fitted ~ subject, data = TB_predictedData[TB_predictedData$stimCondition == 'sham', ], FUN = mean)
# Merge the mean values for each stimCondition level
TBsubject_meanFit_diff <- merge(TBmeanFit_stim, TBmeanFit_sham, by = "subject", suffixes = c("_stim", "_sham"))
TBsubject_meanFit_diff$mean_difference <- TBsubject_meanFit_diff$.fitted_stim - TBsubject_meanFit_diff$.fitted_sham


## SOCIAL SCRIPT ####

# Calculate the mean for stimCondition and sham
SSmeanFit_stim <- aggregate(.fitted ~ subject, data = SS_predictedData[SS_predictedData$stimCondition == 'active', ], FUN = mean)
SSmeanFit_sham <- aggregate(.fitted ~ subject, data = SS_predictedData[SS_predictedData$stimCondition == 'sham', ], FUN = mean)
# Merge the mean values for each stimCondition level
SSsubject_meanFit_diff <- merge(SSmeanFit_stim, SSmeanFit_sham, by = "subject", suffixes = c("_stim", "_sham"))
SSsubject_meanFit_diff$mean_difference <- SSsubject_meanFit_diff$.fitted_stim - SSsubject_meanFit_diff$.fitted_sham


## MECHANICAL ####

# Calculate the mean for stimCondition and sham
MEmeanFit_stim <- aggregate(.fitted ~ subject, data = ME_predictedData[ME_predictedData$stimCondition == 'active', ], FUN = mean)
MEmeanFit_sham <- aggregate(.fitted ~ subject, data = ME_predictedData[ME_predictedData$stimCondition == 'sham', ], FUN = mean)
# Merge the mean values for each stimCondition level
MEsubject_meanFit_diff <- merge(MEmeanFit_stim, MEmeanFit_sham, by = "subject", suffixes = c("_stim", "_sham"))
MEsubject_meanFit_diff$mean_difference <- MEsubject_meanFit_diff$.fitted_stim - MEsubject_meanFit_diff$.fitted_sham

##### Effect of the stimCondition on each condition residuals #####

## FALSE BELIEF ####

# Calculate the mean for stimCondition and sham
FBmeanRes_stim <- aggregate(.resid ~ subject, data = FB_predictedData[FB_predictedData$stimCondition == 'active', ], FUN = mean)
FBmeanRes_sham <- aggregate(.resid ~ subject, data = FB_predictedData[FB_predictedData$stimCondition == 'sham', ], FUN = mean)
# Merge the mean values for each stimCondition level
FBsubject_meanRes_diff <- merge(FBmeanRes_stim, FBmeanRes_sham, by = "subject", suffixes = c("_stim", "_sham"))
FBsubject_meanRes_diff$mean_difference <- FBsubject_meanRes_diff$.resid_stim - FBsubject_meanRes_diff$.resid_sham


## TRUE BELIEF ####

# Calculate the mean for stimCondition and sham
TBmeanRes_stim <- aggregate(.resid ~ subject, data = TB_predictedData[TB_predictedData$stimCondition == 'active', ], FUN = mean)
TBmeanRes_sham <- aggregate(.resid ~ subject, data = TB_predictedData[TB_predictedData$stimCondition == 'sham', ], FUN = mean)

# Merge the mean values for each stimCondition level
TBsubject_meanRes_diff <- merge(TBmeanRes_stim, TBmeanRes_sham, by = "subject", suffixes = c("_stim", "_sham"))
TBsubject_meanRes_diff$mean_difference <- TBsubject_meanRes_diff$.resid_stim - TBsubject_meanRes_diff$.resid_sham


## SOCIAL SCRIPT ####

# Calculate the mean for stimCondition and sham
SSmeanRes_stim <- aggregate(.resid ~ subject, data = SS_predictedData[SS_predictedData$stimCondition == 'active', ], FUN = mean)
SSmeanRes_sham <- aggregate(.resid ~ subject, data = SS_predictedData[SS_predictedData$stimCondition == 'sham', ], FUN = mean)
# Merge the mean values for each stimCondition level
SSsubject_meanRes_diff <- merge(SSmeanRes_stim, SSmeanRes_sham, by = "subject", suffixes = c("_stim", "_sham"))
SSsubject_meanRes_diff$mean_difference <- SSsubject_meanRes_diff$.resid_stim - SSsubject_meanRes_diff$.resid_sham


## MECHANICAL ####

# Calculate the mean for stimCondition and sham
MEmeanRes_stim <- aggregate(.resid ~ subject, data = ME_predictedData[ME_predictedData$stimCondition == 'active', ], FUN = mean)
MEmeanRes_sham <- aggregate(.resid ~ subject, data = ME_predictedData[ME_predictedData$stimCondition == 'sham', ], FUN = mean)
# Merge the mean values for each stimCondition level
MEsubject_meanRes_diff <- merge(MEmeanRes_stim, MEmeanRes_sham, by = "subject", suffixes = c("_stim", "_sham"))
MEsubject_meanRes_diff$mean_difference <- MEsubject_meanRes_diff$.resid_stim - MEsubject_meanRes_diff$.resid_sham


