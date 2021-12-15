############

#### The Effects of Pulse Precipitation Reduction on Pinus edulis and Juniperus monosperma plant water status and leaf gas exchange ###

### Model assumptions, linear mixed models, ANOVA, and Post Hoc tests 

### R Script by: Jeremy Adkins, Taylor Brown, Mona Farnisa & Elise Pletcher

# Load packages
library(cowplot)
library(ggplot2)
library(jcolors)
library(RColorBrewer)
library(scale_fill_brewer)
library(tidyverse)
library(dplyr)
library(agricolae)
library(lme4)
library(lavaan)
library(effects)
library(car)
library(performance)
library(see)
library(ggeffects)
library(ggpubr)
library(multcompView)
library(emmeans)
library(multcomp)

#----------------Further data cleaning-----------------------------

# Read in CSV from cleanup script
df <- read.csv("pj_tidy_dataframe_v2.csv")

# convert treatment to categorical
df$Treatment <- as.factor(df$Treatment)

#pulse as a factor
df$pulse <- as.factor(df$pulse)

#convert species to categorical factor
df$Species <- as.factor(df$Species)

#convert TreeID / Aspect / Slope / Plot / week  to categorical factor
df$TreeID <- as.factor(df$TreeID)
df$Aspect_deg <- as.factor(df$Aspect_deg)
df$Slope_deg <- as.factor(df$Slope_deg)
df$Plot <- as.factor(df$Plot)
df$week <- as.factor(df$week)

### Summer statistics
summary(df)

# remove rows with TreeID = "124"  "128"  "153"  "158"  "164"  "1610" "1613" "1710" "1713"
df <- subset(df, TreeID != "124" & TreeID != "128" & TreeID != "153" & TreeID != "158" & TreeID != "164" & TreeID != "1610" & TreeID != "1613" & TreeID != "1710" & TreeID != "1713") # "1713"))

summary(df)

######################## Build LME Model #######################

#---------------------Midday Leaf Water Content-----by Taylor Brown --------#

#random mixed effects - final model
lmer_mixed_effects <-lmer(midday_leaf_wc ~ Species * pulse * Treatment +  (1|TreeID), data=df)
summary(lmer_mixed_effects)
Anova(lmer_mixed_effects)
hist(residuals(lmer_mixed_effects),breaks=10)
plot(lmer_mixed_effects)

###Diagnostic plots
qqnorm(resid(lmer_mixed_effects))
qqline(resid(lmer_mixed_effects))

####Tukey analyses
emm_wc <- emmeans (lmer_mixed_effects, ~ Species*Treatment | pulse)
emm_wc
plot(emm_wc)
pairs(emm_wc,adjust="tukey")

ltrs.log1 <- cld(emm_wc,
                alpha = 0.05, ###Use lower-case letters for .group
                Letters = letters, ##Report emmeans in orginal scale
                adjust = "tukey") ##Tukey adjustment for multiple comparison
ltrs.log1

#line graph patterns
plot(allEffects(lmer_mixed_effects))

#### final figure - box plots - midday leaf wc
plot(tukey_wc)
pulse.names <- as_labeller(
  c('0' = "Pre-pulse", '1' = "Post-pulse"))

(water_cont_box <- ggplot(df, aes(x=Treatment, y=midday_leaf_wc, fill=factor(Species)))+
    geom_boxplot() +
    labs(title="Mid-Day Leaf Water Content By Species, Treatment, and Pulse",
         x="Treatment (% Precipitation Reduction)",
         y= "Mid-Day Leaf Water Content") +
    facet_wrap(~pulse, labeller = pulse.names)+
    theme(axis.text.x = element_text(angle = 90) ))

#---------------- Water Potential--------by Jeremey Adkins-----------#

# Mixed-effects model
lmer2 <- lmer(water_pot ~ Species*Treatment*pulse + (1|TreeID), data=df)
summary(lmer2)

### check assumptions
check_model(lmer2)#, panel = FALSE)

plot(lmer2)

hist(residuals(lmer2),breaks=10)

qqnorm(resid(lmer2))
qqline(resid(lmer2))  # points don't quite line up with the line

### ANOVA
Anova(lmer2)

#(water_pot.pvals <- as.data.frame(Anova(lmer2)))
#write.csv(water_pot.pvals,"water_pot_pvals.csv")

####Tukey analyses
emm_wp <- emmeans (lmer2, ~ Species*Treatment | pulse)
emm_wp
plot(emm_wp)
pairs(emm_wp,adjust="tukey")

ltrs.log2 <- cld(emm_wp,
                 alpha = 0.05, ###Use lower-case letters for .group
                 Letters = letters, ##Report emmeans in orginal scale
                 adjust = "tukey") ##Tukey adjustment for multiple comparison
ltrs.log2


##### Boxplot for Water Potential

#Need to assign colors to each treatment
mypalette <- brewer.pal(4, "YlOrRd")
pulse.names <- as_labeller(
  c('0' = "Pre-pulse", '1' = "Post-pulse"))

(water_pot_box <- ggplot(df, aes(x=Treatment, y=water_pot, fill=factor(Species)))+
    geom_boxplot() +
    labs(title="Water Potential By Species, Treatment, and Pulse",
         x="Treatment (% Precipitation Reduction)",
         y= "Water Potential") +
    facet_wrap(~pulse, labeller = pulse.names)+
    theme(axis.text.x = element_text(angle = 90) ))


#-------------------- Max photosynthesis------by Elise Pletcher-------------------------#

#linear mixed effects model model
# linear with a normal distribution is probably fine, but include treeID and plot as random effects

# steps to improving model fit:
# log transform max photosynthesis data
# BUT first we need a constant, so: log(Y + a)

# define a constant: take the absolute value of the most negative value of max photsynthesis and add 0.02
(constant <- abs(min(df$max_photo, na.rm = T) - 0.02)) # 0.02 ensures distribution of values has no outliers

# now define new max_photo column
df$max_photo_log <- log(df$max_photo + constant)

# plot a histogram of new max photosynthesis model
hist(df$max_photo_log)

#specify non mixed model
lm.phot <- lm(max_photo_log ~ Species * Treatment * pulse, data = df)
# diagnostics
summary(lm.phot)
plot(lm.phot) # a couple extreme outliers but otherwise looks okay???
aov(lm.phot)

#specify mixed model
# this is the model we want, but we might need to simplify --it has convergence issues
# mixed.phot <- lmer(max_photo_log ~ Species + Treatment + pulse + Species:Treatment + Species:pulse + (1|TreeID), data = df)

df$Plot <- as.factor(df$Plot)

# This simplified model works!
mixed.phot <- lmer(max_photo_log ~ Species * pulse * Treatment + (1|Plot), data = df)

# diagnostic plots
plot(mixed.phot)
qqnorm(resid(mixed.phot))
qqline(resid(mixed.phot))  
hist(resid(mixed.phot))

#-------------------- Anova of LME ---------------------------
car::Anova(mixed.phot)

#(anova.pvals <- as.data.frame(car::Anova(mixed.phot))) 

#write.csv(anova.pvals, "anova_pvals_max_photosynthesis.csv")

#-------------------- Ad hoc tukey --------------------------

#specified by pre vs post pulse
emm_p <- emmeans(mixed.phot, ~ Species*Treatment | pulse)
emm_p

ltrs.log3 <- cld(emm_p,
                alpha   = 0.05,
                Letters = letters,    ### Use lower-case letters for .group
                type    = "response", ### Report emmeans in orginal scale
                adjust =  "tukey")    ### Tukey adjustment for multiple comparisons
ltrs.log3

#----------------- Plot A max box plots -----------------------------

(Amax_box <- ggplot(df, aes(x=Treatment, y=max_photo_log, fill=factor(Species)))+
   geom_boxplot() +
   labs(title="Maximum Photosynthesis by Treatment and pulse ocurrence",
        x="Treatment (% Precipitation Reduction)",
        y= "Maximum Photosynthesis (log transformed)") +
   facet_wrap(~pulse, labeller = pulse.names)+
   theme(axis.text.x = element_text(angle = 90) ))

#-------------------Stomatal Conductance ------by Mona Farnisa---------------------

df$stom_cond
hist(df$stom_cond)
sort(df$stom_cond)

### apply constant 
## log transform the data 
# run the model 
# check residuals and normality

# define constant: take the absolute value of the most negatove value of stom_cond and add 0.02
(constant <- abs(min(df$stom_cond, na.rm = T) - 0.01)) #0.01

# now define new stom_cond column
df$stom_cond_log <- log(df$stom_cond + constant)

# plot a histrogram of new stom_cond model 
hist(df$stom_cond_log)


mod1 <- lmer(stom_cond_log ~ Species*pulse*Treatment + (1|TreeID), data=df)
summary(mod1)

# check plots 
plot(mod1)
qqnorm(resid(mod1))
qqline(resid(mod1))
# residuals
hist(residuals(mod1), breaks=10)   # histogram of residuals

library(see)
library(performance)
### check assumptions
check_model(mod1)

car::Anova(mod1)

# post hoc
emmeans(mod1, specs=c("Species", "Treatment", "pulse"))

emm_sc <- emmeans (mod1, ~ Species*Treatment | pulse)
emm_sc

pwpp(emm_sc)

plot(emm_sc)

pairs(emm_sc, adjust="sidak")

contrast(emm_sc, method = "pairwise")

ltrs.log4 <- cld(emm_sc,
                alpha   = 0.05,
                Letters = letters,    ### Use lower-case letters for .group
                type    = "response", ### Report emmeans in orginal scale
                adjust =  "tukey")    ### Tukey adjustment for multiple comparisons
ltrs.log4

## stomatal boxplot figures 

mypalette <- brewer.pal(4, "YlOrRd")
pulse.names <- as_labeller(
  c('0' = "Pre-pulse", '1' = "Post-pulse"))

(water_pot_box <- ggplot(df, aes(x=Treatment, y=stom_cond_log, fill=factor(Species)))+
    geom_boxplot() +
    labs(title="Stomatal Conductance By Species, Treatment, and Pulse",
         x="Treatment (% Precipitation Reduction)",
         y= "Stomatal Conductance (log transformed)") +
    facet_wrap(~pulse, labeller = pulse.names)+
    theme(axis.text.x = element_text(angle = 90) ))
