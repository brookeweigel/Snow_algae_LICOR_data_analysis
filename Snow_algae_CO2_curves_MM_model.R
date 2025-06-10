
# Load libraries
library(plotrix)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

## Upload combined CO2 curve data
CO2_response <- read.csv("All_sites_CO2_curves.csv")
head(CO2_response)

## NOTE ABOUT UNIT CONVERSION:
# The rate of CO2_flux from the LI-COR is in units of µmol CO2 s⁻¹
# The units of Chl.a are µg ml⁻¹ (because you multiply by the chamber volume, 15 mL)
# The Chl.a normalized rate units are then: µmol CO2 µg⁻¹ Chl.a s⁻¹
# Convert to hourly rates (see below) so the final units will be: µmol CO2 µg⁻¹ Chl.a h⁻¹

## There are 3600 seconds in an hour, so convert from per second to per hour by x 3600
CO2_response <- CO2_response %>% mutate(CO2_Flux_hr = CO2_Flux*3600) %>%
  mutate(CO2_Flux_cell_hr = CO2_Flux_cell*3600) %>%
  mutate(CO2_Flux_chl.a_hr = CO2_Flux_chl.a*3600)

## Plot all CO2 curves separated per site (Y-axis fixed across sites)
CO2_curves_per_site_fixed <- ggplot(data = CO2_response, aes(x = CO2_sample, y = CO2_Flux_hr, fill = Patch)) +
  facet_wrap(~Site, nrow = 2, scales="fixed") +
  geom_point(aes(colour = Patch), size = 6) +
  geom_line(linetype = "dashed", color ="grey50") +
  theme_bw() + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25),
        legend.text = element_text(size = 25), legend.title = element_text(size = 25),
        strip.text.x = element_text(size = 25)) +
  ylab(bquote('CO2 Flux ('*mu*'mol' ~CO[2]~ h^-1*')')) +
  xlab ("CO2 (ppm)") +
  scale_x_continuous(breaks = c(0,200,400,600,800,1000,1200,1400,1600)) +
  scale_color_manual(values=c("slategray4", "skyblue", "steelblue", "royalblue", "royalblue3", "navy", "mediumpurple3", "plum3", "plum4")) +
  ## Note: strip.text.y changes the font size for the facet grid labels
  theme(strip.text.y = element_text(size=30), axis.text.x = element_text(angle = -45, hjust=0))
CO2_curves_per_site_fixed

# Save 8 x 16 PDF

######################################

## Plot all CO2 curves separated per site NORMALIZED BY CELL COUNTS

CO2_curves_cell_count_normalized <- ggplot(data = CO2_response, aes(x = CO2_sample, y = CO2_Flux_cell_hr, fill = Patch)) +
  facet_wrap(~Site, nrow = 2, scales="fixed") +
  geom_point(aes(colour = Patch), size = 6) +
  geom_line(linetype = "dashed", color ="grey50") +
  theme_bw() + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25),
        legend.text = element_text(size = 25), legend.title = element_text(size = 25),
        strip.text.x = element_text(size = 25)) +
  ylab(bquote('CO2 Flux ('*mu*'mol' ~CO[2]~ cell^-1* ~ h^-1*')')) +
  xlab ("CO2 (ppm)") +
  scale_x_continuous(breaks = c(0,200,400,600,800,1000,1200,1400,1600)) +
  scale_y_continuous(labels = scales::scientific) +
  scale_color_manual(values=c("slategray4", "skyblue", "steelblue", "royalblue", "royalblue3", "navy", "mediumpurple3", "plum3", "plum4")) +
  ## Note: strip.text.y changes the font size for the facet grid labels
  theme(strip.text.y = element_text(size=30), axis.text.x = element_text(angle = -45, hjust=0))
CO2_curves_cell_count_normalized

# Save 8 x 16 PDF

######################################

## Plot all CO2 curves separated per site NORMALIZED BY CHLOROPHYLL A

CO2_curves_Chl.a_normalized <- ggplot(data = CO2_response, aes(x = CO2_sample, y = CO2_Flux_chl.a_hr, fill = Patch)) +
  facet_wrap(~Site, nrow = 2, scales="fixed") +
  geom_point(aes(colour = Patch), size = 6) +
  geom_line(linetype = "dashed", color ="grey50") +
  theme_bw() + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25),
        legend.text = element_text(size = 25), legend.title = element_text(size = 25),
        strip.text.x = element_text(size = 25)) +
  ylab(bquote('CO2 flux ('*mu*'mol' ~CO[2] ~ mu*'g' ~Chl.a^-1* ~ h^-1*')')) +
  xlab ("CO2 (ppm)") +
  scale_x_continuous(breaks = c(0,200,400,600,800,1000,1200,1400,1600)) +
  scale_y_continuous(breaks = seq(-4, 10, 2)) +
  scale_color_manual(values=c("slategray4", "skyblue", "steelblue", "royalblue", "royalblue3", "navy", "mediumpurple3", "plum3", "plum4")) +
  ## Note: strip.text.y changes the font size for the facet grid labels
  theme(strip.text.y = element_text(size=30), axis.text.x = element_text(angle = -45, hjust=0))
CO2_curves_Chl.a_normalized

# Save 8 x 16 PDF

######################################
######################################

## Load packages for analyzing CO2 curves
library(drc) # for fitting Michaelis Menten model

## From Steensma et al. 2023: "Michaelis–Menten parameters (Km(CO2), Amax) were determined by using
## the R package “drc” (Version 3.0–1; Ritz et al. (2015)) to fit a two-parameter M-M equation to each replicate"

citation(package = "drc")

## Set grouping variable, filter one curve as an example
Single_curve <- CO2_response %>% group_by(LICOR_number) %>% filter(LICOR_number == "LI_10")

## Fit 3 parameter Michaelis–Menten model
model.drm <- drm (CO2_Flux ~ CO2_sample, data = Single_curve, fct = MM.3())

## This is the 3-parameter Michaelis–Menten model equation:
# f(x) = c + (d-c)/(1+(e/x))

# See: https://search.r-project.org/CRAN/refmans/drc/html/MM.html

# Dose (CO2 level) = x, 
# Lower limit c at dose 0 (x=0)
# Upper limit d
# Parameter e corresponds to the CO2 level yielding a response halfway between c and d (half saturation)

## Summarize output
summary(model.drm)

## Find P_max and Km model coefficients:
# Pmax = d:(Intercept)
# Km = e:(Intercept)
# Respiration = c:(Intercept)
coef(model.drm)

## Generate a model fit line based on the fitted parameters
Predicted_curve <- data.frame(CO2_level = seq(0, 1600, length.out = 100))
Predicted_curve$CO2_Flux <- predict(model.drm, newdata = Predicted_curve)

## Plot original data (points) with model fit (line):
ggplot(data = Predicted_curve, mapping = aes(CO2_level, CO2_Flux)) + geom_line() +
  geom_point(data = Single_curve, mapping = aes(CO2_sample, CO2_Flux)) +
  labs(x = "CO2 (ppm)", y = "CO2 Flux") + theme_bw()


####################################################
### Fit CO2 response curves for ALL replicates ###
####################################################

## Load libraries
library(purrr)
library(drc)

## Upload combined CO2 curve data
CO2_response <- read.csv("All_sites_CO2_curves.csv")
head(CO2_response)

## Fits all CO2 curves at once, grouped by "split" variable
## The function "map" applies the function using a "for loop" to each individual curve (grouped by split)
fits = CO2_response %>% filter (LICOR_number != "LI_15") %>% filter (LICOR_number != "LI_34") %>% filter (LICOR_number != "LI_45") %>%
  split(~ LICOR_number) %>% map(~ drm(CO2_Flux ~ CO2_sample, data = ., fct = MM.3()))

## Make a data frame of estimated parameters:
## Note, in the purr library, the function "map" makes a list
Model_fits <- fits %>% map(coef) %>% map(t) %>% map(as.data.frame) %>% imap_dfr(~ mutate(.x, LICOR_number = .y)) %>%
  rename(Respiration = "c:(Intercept)") %>% rename(Pmax = "d:(Intercept)") %>% rename(Km = "e:(Intercept)")
head(Model_fits)

## FOR LOOP to create model fit lines for EVERY LIGHT CURVE based on the fitted parameters:
# First, make a new data frame with CO2 levels from 0 to 1600
Predicted_model_curves <- data.frame(CO2 = seq(0, 1600, length.out = 100))

# Iterate over unique LICOR_number values
unique_LICOR_numbers <- unique(Model_fits$LICOR_number)

for (i in unique_LICOR_numbers) {
  # Subset data based on LICOR_number
  subset_data <- Model_fits[Model_fits$LICOR_number == i, ]
  # Extract required variables
  Respiration <- subset_data$Respiration
  Pmax <- subset_data$Pmax
  Km <- subset_data$Km
  # Fit the Michaelis–Menten model
  CO2_flux <- Respiration + (Pmax - Respiration) / (1 + (Km / Predicted_model_curves$CO2))
  # Add CO2_flux to the Predicted_model_curves data frame
  col_name <- paste(i, "_CO2_flux", sep = "")
  Predicted_model_curves[[col_name]] <- CO2_flux
}


### As a test, try to plot predicted model curves
All_model_curves_Bagley <- ggplot(Predicted_model_curves, aes(x = CO2_sample, y = CO2_flux)) + 
  geom_line(data = Predicted_model_curves, aes(x = CO2, y = LI_7_CO2_flux)) +
  geom_line(data = Predicted_model_curves, aes(x = CO2, y = LI_8_CO2_flux)) +
  geom_line(data = Predicted_model_curves, aes(x = CO2, y = LI_9_CO2_flux)) +
  geom_line(data = Predicted_model_curves, aes(x = CO2, y = LI_10_CO2_flux)) +
  geom_line(data = Predicted_model_curves, aes(x = CO2, y = LI_11_CO2_flux))
All_model_curves_Bagley

########################################################
### Add model fit to the original light curve graphs ###
########################################################

## First, I need to convert the format from multiple columns into one long row of model fits:
Predicted_model_curves_new <- data.frame(LICOR_number=unlist(Predicted_model_curves, use.names = TRUE))

## Save as excel file and fix names, add site and patch numbers
write.csv(Predicted_model_curves_new, "~/Documents/2_NSF_PRFB_Snow_algae/01_LI-COR-Experiments/2_Data_Analysis_R/Predicted_model_curves_CO2_new.csv", row.names=TRUE)

########################################################

## Upload final csv with re-formatted model fits, along with site and patch numbers
Model_fits_final <- read.csv("Predicted_model_curves_CO2_new.csv")
head(Model_fits_final)

## FINAL STEP: GRAPHING DATA WITH MODEL FITS
# Add model fitted lines to the original light curve graphs (add geom_line for each unique LICOR_number):

CO2_curves_with_model_fits <- ggplot(data = CO2_response, aes(x = CO2_sample, y = CO2_Flux, fill = Patch)) +
  geom_point(aes(colour = Patch), size = 6) +
  geom_line(data = Model_fits_final, aes(x = CO2, y = CO2_Flux, color = Patch), linetype = "solid") +
  facet_wrap(~Site, nrow = 2, scales="fixed") +
  theme_bw() + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25),
        legend.text = element_text(size = 25), legend.title = element_text(size = 25),
        strip.text.x = element_text(size = 25)) +
  ylab(bquote('CO2 Flux ('*mu~ 'mol' ~CO[2]~ s^-1*')')) +
  xlab ("CO2 (ppm)") +
  scale_y_continuous(breaks = seq(-0.0004, 0.0010, 0.0004)) +
  scale_x_continuous(breaks = c(0,200,400,600,800,1000,1200,1400,1600)) +
  scale_color_manual(values=c("slategray4", "skyblue", "steelblue", "royalblue", "royalblue3", "navy", "mediumpurple3", "plum3", "plum4")) +
  ## Note: strip.text.y changes the font size for the facet grid labels
  theme(strip.text.y = element_text(size=30), axis.text.x = element_text(angle = -45, hjust=0))
CO2_curves_with_model_fits

# Save 8 x 16 PDF

##############################################################################
### Fit CO2 response curves for ALL replicates NORMALIZED BY CHLOROPHYLL A ###
##############################################################################

## Load libraries
library(purrr)
library(drc)

## Upload combined CO2 curve data
CO2_response <- read.csv("All_sites_CO2_curves_400avg_R.csv")
head(CO2_response)

## There are 3600 seconds in an hour, so convert from per second to per hour by x 3600
CO2_response <- CO2_response %>% mutate(CO2_Flux_hr = CO2_Flux*3600) %>%
  mutate(CO2_Flux_cell_hr = CO2_Flux_cell*3600) %>%
  mutate(CO2_Flux_chl.a_hr = CO2_Flux_chl.a*3600)

## Fits all CO2 curves at once, grouped by "split" variable
## The function "map" applies the function using a "for loop" to each individual curve (grouped by split)
fits = CO2_response %>% filter (LICOR_number != "LI_15") %>% filter (LICOR_number != "LI_34") %>% filter (LICOR_number != "LI_45") %>%
  split(~ LICOR_number) %>% map(~ drm(CO2_Flux_chl.a_hr ~ CO2_sample, data = ., fct = MM.3()))

## Make a data frame of estimated parameters:
## Note, in the purr library, the function "map" makes a list
Model_fits <- fits %>% map(coef) %>% map(t) %>% map(as.data.frame) %>% imap_dfr(~ mutate(.x, LICOR_number = .y)) %>%
  rename(Respiration = "c:(Intercept)") %>% rename(Pmax = "d:(Intercept)") %>% rename(Km = "e:(Intercept)")
head(Model_fits)

### Save the Model Fits as excel file for graphing parameters later
write.csv(Model_fits, "~/Documents/2_NSF_PRFB_Snow_algae/01_LI-COR-Experiments/2_Data_Analysis_R/CO2_curves_Chl.a_normalized_parameters.csv", row.names=TRUE)

## FOR LOOP to create model fit lines for EVERY LIGHT CURVE based on the fitted parameters:
# First, make a new data frame with CO2 levels from 0 to 1600
Predicted_model_curves_Chl.a <- data.frame(CO2 = seq(0, 1600, length.out = 100))

# Iterate over unique LICOR_number values
unique_LICOR_numbers <- unique(Model_fits$LICOR_number)

for (i in unique_LICOR_numbers) {
  # Subset data based on LICOR_number
  subset_data <- Model_fits[Model_fits$LICOR_number == i, ]
  # Extract required variables
  Respiration <- subset_data$Respiration
  Pmax <- subset_data$Pmax
  Km <- subset_data$Km
  # Fit the Michaelis–Menten model
  CO2_flux <- Respiration + (Pmax - Respiration) / (1 + (Km / Predicted_model_curves_Chl.a$CO2))
  # Add CO2_flux to the Predicted_model_curves data frame
  col_name <- paste(i, "_CO2_flux", sep = "")
  Predicted_model_curves_Chl.a[[col_name]] <- CO2_flux
}


### As a test, try to plot predicted model curves
All_model_curves_Bagley <- ggplot(Predicted_model_curves_Chl.a, aes(x = CO2_sample, y = CO2_flux)) + 
  geom_line(data = Predicted_model_curves_Chl.a, aes(x = CO2, y = LI_7_CO2_flux)) +
  geom_line(data = Predicted_model_curves_Chl.a, aes(x = CO2, y = LI_8_CO2_flux)) +
  geom_line(data = Predicted_model_curves_Chl.a, aes(x = CO2, y = LI_9_CO2_flux)) +
  geom_line(data = Predicted_model_curves_Chl.a, aes(x = CO2, y = LI_10_CO2_flux)) +
  geom_line(data = Predicted_model_curves_Chl.a, aes(x = CO2, y = LI_11_CO2_flux))
All_model_curves_Bagley

########################################################
### Add model fit to the original light curve graphs ###
########################################################

## First, I need to convert the format from multiple columns into one long row of model fits:
Predicted_model_curves_Chl.a <- data.frame(LICOR_number=unlist(Predicted_model_curves_Chl.a, use.names = TRUE))

## Save as excel file and fix names, add site and patch numbers
write.csv(Predicted_model_curves_Chl.a, "~/Documents/2_NSF_PRFB_Snow_algae/01_LI-COR-Experiments/2_Data_Analysis_R/Predicted_model_curves_CO2_Chla.csv", row.names=TRUE)

########################################################

## Upload final csv with re-formatted model fits, along with site and patch numbers
Model_fits_final_Chla <- read.csv("Predicted_model_curves_CO2_Chl.a.csv")
head(Model_fits_final_Chla)

## FINAL STEP: GRAPHING DATA WITH MODEL FITS
# Add model fitted lines to the original light curve graphs (add geom_line for each unique LICOR_number):

CO2_curves_with_model_fits_Chla <- ggplot(data = CO2_response, aes(x = CO2_sample, y = CO2_Flux_chl.a_hr, fill = Patch)) +
  geom_point(aes(colour = Patch), size = 6) +
  geom_line(data = Model_fits_final_Chla, aes(x = CO2, y = CO2_Flux, color = Patch), linetype = "solid") +
  facet_wrap(~Site, nrow = 2, scales="fixed") +
  theme_bw() + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25),
        legend.text = element_text(size = 25), legend.title = element_text(size = 25),
        strip.text.x = element_text(size = 25)) +
  ylab(bquote('CO2 flux ('*mu*'mol' ~CO[2] ~ mu*'g' ~Chl.a^-1* ~ h^-1*')')) +
  xlab ("CO2 (ppm)") +
  scale_x_continuous(breaks = c(0,200,400,600,800,1000,1200,1400,1600)) +
  scale_y_continuous(breaks = seq(-4, 10, 2)) +
  scale_color_manual(values=c("slategray4", "skyblue", "steelblue", "royalblue", "royalblue3", "navy", "mediumpurple3", "plum3", "plum4")) +
  ## Note: strip.text.y changes the font size for the facet grid labels
  theme(strip.text.y = element_text(size=30), axis.text.x = element_text(angle = -45, hjust=0))
CO2_curves_with_model_fits_Chla

# Save 8 x 16 PDF

#######################################

