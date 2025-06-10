
# Load libraries
library(plotrix)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

## Upload combined light curve data
light_response <- read.csv("All_sites_light_curves.csv")
head(light_response)

## NOTE ABOUT UNIT CONVERSION:
# The rate of CO2_flux from the LI-COR is in units of µmol CO2 s⁻¹
# The units of Chl.a are µg ml⁻¹ (because you multiply by the chamber volume, 15 mL)
# The Chl.a normalized rate units are then: µmol CO2 µg⁻¹ Chl.a s⁻¹
# Convert to hourly rates (see below) so the final units will be: µmol CO2 µg⁻¹ Chl.a h⁻¹

## There are 3600 seconds in an hour, so convert from per second to per hour by x 3600
light_response <- light_response %>% mutate(CO2_Flux_hr = CO2_Flux*3600) %>%
  mutate(CO2_Flux_cell_hr = CO2_Flux_cell*3600) %>%
  mutate(CO2_Flux_chl.a_hr = CO2_Flux_chl.a*3600)

## Plot all light curves separated per site (Y-axis fixed across sites)
light_curves_per_site_fixed <- ggplot(data = light_response, aes(x = Light, y = CO2_Flux_hr, fill = Patch)) +
  facet_wrap(~Site, nrow = 2, scales="fixed") +
  geom_point(aes(colour = Patch), size = 6) +
  geom_line(linetype = "dashed", color ="grey50") +
  theme_bw() + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25),
        legend.text = element_text(size = 25), legend.title = element_text(size = 25),
        strip.text.x = element_text(size = 25)) +
  ylab(bquote('CO2 Flux ('*mu~ 'mol' ~CO[2]~ h^-1*')')) +
  xlab(bquote('Light level ('*mu*'mol' ~m^-2* ~ s^-1*')')) +
  scale_x_continuous(breaks = c(0,300,600,900,1200,1500,2000,2500,3000)) +
  scale_color_manual(values=c("slategray4", "skyblue", "steelblue", "royalblue", "royalblue3", "navy", "mediumpurple3", "plum3", "plum4")) +
  ## Note: strip.text.y changes the font size for the facet grid labels
  theme(strip.text.y = element_text(size=30), axis.text.x = element_text(angle = -45, hjust=0))
light_curves_per_site_fixed

# Save 8 x 16 PDF

######################################

## Plot all light curves separated per site NORMALIZED BY CELL COUNTS

light_curves_cell_count_normalized <- ggplot(data = light_response, aes(x = Light, y = CO2_Flux_cell_hr, fill = Patch)) +
  facet_wrap(~Site, nrow = 2, scales="fixed") +
  geom_point(aes(colour = Patch), size = 6) +
  geom_line(linetype = "dashed", color ="grey50") +
  theme_bw() + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25),
        legend.text = element_text(size = 25), legend.title = element_text(size = 25),
        strip.text.x = element_text(size = 25)) +
  ylab(bquote('CO2 Flux ('*mu*'mol'~CO[2]~ cell^-1* ~ h^-1*')')) +
  xlab(bquote('Light level ('*mu*'mol' ~m^-2* ~ s^-1*')')) +
  scale_x_continuous(breaks = c(0,300,600,900,1200,1500,2000,2500,3000)) +
  scale_y_continuous(labels = scales::scientific) +
  scale_color_manual(values=c("slategray4", "skyblue", "steelblue", "royalblue", "royalblue3", "navy", "mediumpurple3", "plum3", "plum4")) +
  ## Note: strip.text.y changes the font size for the facet grid labels
  theme(strip.text.y = element_text(size=30), axis.text.x = element_text(angle = -45, hjust=0))
light_curves_cell_count_normalized

# Save 8 x 16 PDF

######################################

## Plot all Light curves separated per site NORMALIZED BY CHLOROPHYLL A

## Note: switch to scales="free_y" to see variable y-axis per site

Light_curves_Chl.a_normalized <- ggplot(data = light_response, aes(x = Light, y = CO2_Flux_chl.a_hr, fill = Patch)) +
  facet_wrap(~Site, nrow = 2, scales="fixed") +
  geom_point(aes(colour = Patch), size = 6) +
  geom_line(linetype = "dashed", color ="grey50") +
  theme_bw() + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25),
        legend.text = element_text(size = 25), legend.title = element_text(size = 25),
        strip.text.x = element_text(size = 25)) +
  ylab(bquote('CO2 Flux ('*mu*'mol' ~CO[2] ~ mu*'g' ~Chl.a^-1* ~ h^-1*')')) +
  xlab(bquote('Light level ('*mu*'mol' ~m^-2* ~ s^-1*')')) +
  scale_x_continuous(breaks = c(0,300,600,900,1200,1500,2000,2500,3000)) +
  scale_color_manual(values=c("slategray4", "skyblue", "steelblue", "royalblue", "royalblue3", "navy", "mediumpurple3", "plum3", "plum4")) +
  ## Note: strip.text.y changes the font size for the facet grid labels
  theme(strip.text.y = element_text(size=30), axis.text.x = element_text(angle = -45, hjust=0))
Light_curves_Chl.a_normalized

# Save 8 x 16 PDF

######################################

## Load packages for analyzing light curves
library(broom)
library(dplyr)
library(tidyr)
library(photosynthesis)

## See: https://cran.r-project.org/web/packages/photosynthesis/readme/README.html

citation(package = 'photosynthesis')

# K_sat = maximum rate of net photosynthesis = light-saturated net CO2 assimilation rate
# phi_J = alpha = photochemical efficiency of photosynthesis at low light (or quantum yield of CO2 assimilation)
# theta_J = curvature of the light response
# Rd = rate of dark respiration
# LCP = light compensation point

## Upload combined light curve data
light_response <- read.csv("All_sites_light_curves.csv")
head(light_response)

## Set grouping variable, filter one curve as an example
Single_curve <- light_response %>% group_by(LICOR_number) %>% filter(LICOR_number == "LI_10")
  
## Fit one light-response curve
fit = fit_photosynthesis(.data = Single_curve, .photo_fun = "aq_response",
  .vars = list(.A = CO2_Flux, .Q = Light))

## Model summary:
summary(fit)

## This code extracts just the residual standard error from the model summary
summary(fit)$sigma

## Try to graph model fit with data

## Generate a model fit line based on the fitted parameters
b = coef(fit)
df_predict = data.frame(Qabs = seq(0, 3000, length.out = 100)) %>%
  mutate(A = marshall_biscoe_1980(Q_abs = Qabs, k_sat = b["k_sat"], b["phi_J"], b["theta_J"]) - b["Rd"])

## Plot original data (points) with model fit (line):
ggplot(data = df_predict, mapping = aes(Qabs, A)) + geom_line() +
  geom_point(data = Single_curve, mapping = aes(Light, CO2_Flux)) +
  labs(x = "Light", y = "CO2 Flux") + theme_bw()

####################################################
### Fit light response curves for ALL replicates ###
####################################################

## Load libraries
library(purrr)
library(photosynthesis)

## Upload combined light curve data
light_response <- read.csv("All_sites_light_curves.csv")
head(light_response)

## Fits all light curves at once, grouped by "split" variable
## The function "map" applies the function using a "for loop" to each individual curve (grouped by split)
fits = light_response %>% split(~ LICOR_number) %>%
  map(fit_photosynthesis, .photo_fun = "aq_response", .vars = list(.A = CO2_Flux, .Q = Light))

## Make a data frame of estimated parameters:
## Note, in the purr library, the function "map" makes a list
Model_fits <- fits %>% map(coef) %>% map(t) %>% map(as.data.frame) %>% imap_dfr(~ mutate(.x, LICOR_number = .y))

## FOR LOOP to create model fit lines for EVERY LIGHT CURVE based on the fitted parameters:

# First, make a new data frame with light levels from 0 to 3,000
Predicted_model_curves <- data.frame(Light = seq(0, 3000, length.out = 100))

# Iterate over unique LICOR_number values
unique_LICOR_numbers <- unique(Model_fits$LICOR_number)

for (i in unique_LICOR_numbers) {
  # Subset data based on LICOR_number
  subset_data <- Model_fits[Model_fits$LICOR_number == i, ]
  # Extract required variables
  k_sat <- subset_data$k_sat
  phi_J <- subset_data$phi_J
  theta_J <- subset_data$theta_J
  Rd <- subset_data$Rd
  # Calculate CO2_flux using the marshall_biscoe_1980 function
  CO2_flux <- marshall_biscoe_1980(Q_abs = Predicted_model_curves$Light, k_sat, phi_J, theta_J) - Rd 
  # Add CO2_flux to the Predicted_model_curves data frame
  col_name <- paste(i, "_CO2_flux", sep = "")
  Predicted_model_curves[[col_name]] <- CO2_flux
}

### As a test, try to plot predicted model curves
All_model_curves_Bagley <- ggplot(Predicted_model_curves, aes(x = Light, y = CO2_flux)) + 
  geom_line(data = Predicted_model_curves, aes(x = Light, y = LI_7_CO2_flux)) +
  geom_line(data = Predicted_model_curves, aes(x = Light, y = LI_8_CO2_flux)) +
  geom_line(data = Predicted_model_curves, aes(x = Light, y = LI_9_CO2_flux)) +
  geom_line(data = Predicted_model_curves, aes(x = Light, y = LI_10_CO2_flux)) +
  geom_line(data = Predicted_model_curves, aes(x = Light, y = LI_11_CO2_flux))
All_model_curves_Bagley

########################################################
### Add model fit to the original light curve graphs ###
########################################################

## First, I need to convert the format from multiple columns into one long row of model fits:
Predicted_model_curves_new <- data.frame(LICOR_number=unlist(Predicted_model_curves, use.names = TRUE))

## Save as excel file and fix names, add site and patch numbers
write.csv(Predicted_model_curves_new, "~/Documents/2_NSF_PRFB_Snow_algae/01_LI-COR-Experiments/2_Data_Analysis_R/Predicted_model_curves_light_new.csv", row.names=TRUE)

## Upload final csv with re-formatted model fits, along with site and patch numbers
Model_fits_final <- read.csv("Predicted_model_curves_light_new.csv")
head(Model_fits_final)

## FINAL STEP: GRAPHING DATA WITH MODEL FITS
# Add model fitted lines to the original light curve graphs (add geom_line for each unique LICOR_number):

light_curves_with_model_fits <- ggplot(data = light_response, aes(x = Light, y = CO2_Flux, fill = Patch)) +
  geom_point(aes(colour = Patch), size = 6) +
  geom_line(data = Model_fits_final, aes(x = Light, y = CO2_Flux, color = Patch), linetype = "solid") +
  facet_wrap(~Site, nrow = 2, scales="fixed") +
  theme_bw() + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25),
        legend.text = element_text(size = 25), legend.title = element_text(size = 25),
        strip.text.x = element_text(size = 25)) +
  ylab(bquote('CO2 Flux ('*mu~ 'mol' ~CO[2]~ s^-1*')')) +
  xlab ("Light level (PAR)") +
  scale_y_continuous(breaks = seq(-0.0003, 0.0010, 0.0003)) +
  scale_x_continuous(breaks = c(0,300,600,900,1200,1500,2000,2500,3000)) +
  scale_color_manual(values=c("slategray4", "skyblue", "steelblue", "royalblue", "royalblue3", "navy", "mediumpurple3", "plum3", "plum4")) +
  ## Note: strip.text.y changes the font size for the facet grid labels
  theme(strip.text.y = element_text(size=30), axis.text.x = element_text(angle = -45, hjust=0))
light_curves_with_model_fits

# Save 8 x 16 PDF
 
################################################################################
### Fit light response curves for ALL replicates NORMALIZED BY CHLOROPHYLL A ###
################################################################################

## Load libraries
library(purrr)
library(photosynthesis)
library(dplyr)

## Upload combined light curve data
light_response <- read.csv("All_sites_light_curves_R.csv")
head(light_response)

## There are 3600 seconds in an hour, so convert from per second to per hour by x 3600
light_response <- light_response %>% mutate(CO2_Flux_hr = CO2_Flux*3600) %>%
  mutate(CO2_Flux_cell_hr = CO2_Flux_cell*3600) %>%
  mutate(CO2_Flux_chl.a_hr = CO2_Flux_chl.a*3600)

## Fits all light curves at once, grouped by "split" variable
## The function "map" applies the function using a "for loop" to each individual curve (grouped by split)
fits = light_response %>% split(~ LICOR_number) %>%
  map(fit_photosynthesis, .photo_fun = "aq_response", .vars = list(.A = CO2_Flux_chl.a_hr, .Q = Light))

## Make a data frame of estimated parameters:
## Note, in the purr library, the function "map" makes a list
Model_fits <- fits %>% map(coef) %>% map(t) %>% map(as.data.frame) %>% imap_dfr(~ mutate(.x, LICOR_number = .y))

##  Solve for the half-saturation light level by solving for Q where k_sat / 2:
Model_fits <- Model_fits %>% mutate(Q_half = k_sat * (2 - theta_J) / (2 * phi_J))

### Save the Model Fits as excel file for graphing parameters later
write.csv(Model_fits, "~/Documents/2_NSF_PRFB_Snow_algae/01_LI-COR-Experiments/2_Data_Analysis_R/Light_curves_Chl.a_normalized_parameters.csv", row.names=TRUE)

## FOR LOOP to create model fit lines for EVERY LIGHT CURVE based on the fitted parameters:

# First, make a new data frame with light levels from 0 to 3,000
Predicted_model_curves_Chl.a <- data.frame(Light = seq(0, 3000, length.out = 100))

# Iterate over unique LICOR_number values
unique_LICOR_numbers <- unique(Model_fits$LICOR_number)

for (i in unique_LICOR_numbers) {
  # Subset data based on LICOR_number
  subset_data <- Model_fits[Model_fits$LICOR_number == i, ]
  # Extract required variables
  k_sat <- subset_data$k_sat
  phi_J <- subset_data$phi_J
  theta_J <- subset_data$theta_J
  Rd <- subset_data$Rd
  # Calculate CO2_flux using the marshall_biscoe_1980 function
  CO2_flux <- marshall_biscoe_1980(Q_abs = Predicted_model_curves_Chl.a$Light, k_sat, phi_J, theta_J) - Rd 
  # Add CO2_flux to the Predicted_model_curves data frame
  col_name <- paste(i, "_CO2_flux", sep = "")
  Predicted_model_curves_Chl.a[[col_name]] <- CO2_flux
}

### As a test, try to plot predicted model curves
All_model_curves_Bagley <- ggplot(Predicted_model_curves_Chl.a, aes(x = Light, y = CO2_flux)) + 
  geom_line(data = Predicted_model_curves_Chl.a, aes(x = Light, y = LI_7_CO2_flux)) +
  geom_line(data = Predicted_model_curves_Chl.a, aes(x = Light, y = LI_8_CO2_flux)) +
  geom_line(data = Predicted_model_curves_Chl.a, aes(x = Light, y = LI_9_CO2_flux)) +
  geom_line(data = Predicted_model_curves_Chl.a, aes(x = Light, y = LI_10_CO2_flux)) +
  geom_line(data = Predicted_model_curves_Chl.a, aes(x = Light, y = LI_11_CO2_flux))
All_model_curves_Bagley

########################################################
### Add model fit to the original light curve graphs ###
########################################################

## First, I need to convert the format from multiple columns into one long row of model fits:
Predicted_model_curves_Chl.a <- data.frame(LICOR_number=unlist(Predicted_model_curves_Chl.a, use.names = TRUE))

## Save as excel file and fix names, add site and patch numbers
# write.csv(Predicted_model_curves_Chl.a, "~/Documents/2_NSF_PRFB_Snow_algae/01_LI-COR-Experiments/2_Data_Analysis_R/Predicted_model_curves_light_Chl.a.csv", row.names=TRUE)

######### FINAL GRAPH HERE ############
## Upload final csv with re-formatted model fits, along with site and patch numbers
Model_fits_final_Chl.a <- read.csv("Predicted_model_curves_light_Chl.a.csv")
head(Model_fits_final_Chl.a)

## FINAL STEP: GRAPHING DATA WITH MODEL FITS
# Add model fitted lines to the original light curve graphs (add geom_line for each unique LICOR_number):

light_curves_with_model_fits_Chl.a <- ggplot(data = light_response, aes(x = Light, y = CO2_Flux_chl.a_hr, fill = Patch)) +
  geom_point(aes(colour = Patch), size = 6) +
  geom_line(data = Model_fits_final_Chl.a, aes(x = Light, y = CO2_Flux, color = Patch), linetype = "solid") +
  facet_wrap(~Site, nrow = 2, scales="fixed") +
  theme_bw() + theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), 
        axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25),
        legend.text = element_text(size = 25), legend.title = element_text(size = 25),
        strip.text.x = element_text(size = 25)) +
  ylab(bquote('CO2 flux ('*mu*'mol' ~CO[2] ~ mu*'g' ~Chl.a^-1* ~ h^-1*')')) +
  xlab(bquote('Light level ('*mu*'mol' ~m^-2* ~ s^-1*')')) +
  scale_x_continuous(breaks = c(0,300,600,900,1200,1500,2000,2500,3000)) +
  scale_color_manual(values=c("slategray4", "skyblue", "steelblue", "royalblue", "royalblue3", "navy", "mediumpurple3", "plum3", "plum4")) +
  ## Note: strip.text.y changes the font size for the facet grid labels
  theme(strip.text.y = element_text(size=30), axis.text.x = element_text(angle = -45, hjust=0))
light_curves_with_model_fits_Chl.a

# Save 8 x 16 PDF


