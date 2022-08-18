# ========================================================================================
# Project:  Occupation_downscaling
# Subject:  Script to prepare figures for si
# Author:   Michiel van Dijk
# Contact:  michiel.vandijk@wur.nl
# ========================================================================================

# ========================================================================================
# SETUP ----------------------------------------------------------------------------------
# ========================================================================================

# Load pacman for p_load
if(!require(pacman)) install.packages("pacman")
library(pacman)

# Load key packages
p_load(here, tidyverse, readxl, stringr, scales, glue)

# Load additional packages
p_load(tidymodels, stacks, sf, raster, viridis, maps, DALEXtra, tictoc, egg, exactextractr)
p_load(viridis, cowplot, ggforce, maps, raster, sf, lessR, patchwork, ggpubr, DALEX)

# Set path
source(here("scripts/set_path.r"))

# R options
options(scipen = 999)
options(digits = 4)


# ========================================================================================
# SET ISO3C ------------------------------------------------------------------------------
# ========================================================================================

iso3c_sel <- "VNM"

# ========================================================================================
# SET DATASET ----------------------------------------------------------------------------
# ========================================================================================

# Change to the date of your results folder
results_date <- "2022-02-12"


# ========================================================================================
# LOAD DATA ------------------------------------------------------------------------------
# ========================================================================================

# model input
source(here("scripts/prepare_input_data.r"))

# grid level predictors
predictors_raw <- readRDS(file.path(proc_path, "model_input/VNM_predictors.rds"))

# adm
adm <- readRDS(file.path(proc_path, glue("adm/{iso3c_sel}_adm.rds")))
adm_r <- raster(file.path(proc_path, glue("adm/{iso3c_sel}_adm_r_30sec.tif")))

# Working age population
wp_raw <- raster(file.path(proc_path, glue("worldpop/{iso3c_sel}_working_population_2009.tif")))

# Model results
db_pred <- readRDS(file.path(proc_path, glue("results/{results_date}/db_pred.rds")))
db_pred_error <- readRDS(file.path(proc_path, glue("results/{results_date}/db_pred_error.rds")))
db_test_pred <- readRDS(file.path(proc_path, glue("results/{results_date}/db_test_pred.rds")))
db_test_stat <- readRDS(file.path(proc_path, glue("results/{results_date}/db_test_stat.rds")))
db_ensemble_tab <- readRDS(file.path(proc_path, glue("results/{results_date}/db_ensemble_tab.rds")))
db_ale <- readRDS(file.path(proc_path, glue("results/{results_date}/db_ale.rds")))

# List of predictors and groupings
pred_info <- read_excel(file.path(proc_path, "predictors/tables.xlsx"), sheet = "predictors")

# We load pre-processed figures created in the post-processing script creating them takes a long time
fig_vim_variable <- readRDS(file.path(proc_path, glue("results/{results_date}/fig_vim_variable.rds")))


# ========================================================================================
# PREPARE --------------------------------------------------------------------------------
# ========================================================================================

# rename to value for easy of using general functions below
db_pred <- db_pred %>%
  rename(value = .pred)

# colorblind palette
col   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# List of variables
var_list <- c("off_mgr_pros", "tech_aspros", "service_shop", "agric", "othlowsk", "lfpr", "elementary")

# Create lookup table
lookup <- data.frame(
  variable = var_list,
  var_name = c("Managers and professionals", "Technicians and associate professionals", "Clerks and service workers",
               "Agricultural workers","Craft workers and operators", "Labour force participation rate", "Elementary occupations"),
  tag = c("a", "b", "c", "d", "e", "g", "f")
)


# ========================================================================================
# NATIONAL OCCUPATION SHARES AND LABOUR FORCE PARTICIPATION RATE -------------------------
# ========================================================================================

# Total working population
wp_sum <- cellStats(wp_raw, sum)

# Working population and lfpr adm level
wp_adm <- adm %>%
  st_drop_geometry() %>%
  mutate(
    wp = exact_extract(wp_raw, adm, fun = "sum")) %>%
  dplyr::select(adm_code, wp) %>%
  left_join(
    db_adm %>%
        filter(variable == "lfpr") %>%
        dplyr::select(adm_code, lfpr = value)
  )

# National occupation shares
occ_nat <- db_adm %>%
  filter(variable != "lfpr") %>%
  dplyr::select(adm_code, variable, value) %>%
  left_join(wp_adm) %>%
  group_by(variable) %>%
  summarize(value = sum(value * lfpr * wp, value, na.rm = TRUE)/1000000,
            .groups = "drop") %>%
  mutate(share = value/sum(value, na.rm = TRUE)) %>%
  left_join(lookup)

fig_occ_nat <- ggplot(data = occ_nat, aes(x = "", y = value, fill = var_name, label = scales::percent(share, accuracy = 1))) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  labs(x = NULL, y = "Number of workers (million)", fill = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), breaks = pretty_breaks()) +
  geom_text(position=position_stack(vjust=0.5), size = 3) +
  scale_fill_manual(values = col[-c(1,5)]) +
  theme_classic() +
  rremove("x.ticks")


# ========================================================================================
# INPUT DATA -----------------------------------------------------------------------------
# ========================================================================================

# ADM OCCUPATION SHARES ------------------------------------------------------------------

# Occupation rates
fig_occ <- db_adm %>%
  filter(variable != "lfpr") %>%
  left_join(adm,.) %>%
  left_join(lookup) %>%
  mutate(value = cut(value, breaks = c(0, 0.02, 0.10, 0.5, 1),
                     labels = c("0-2", "2-10", "10-50", "50-100"))) %>%
  ggplot() +
  geom_sf(color = "transparent", aes(fill = value)) +
  labs(fill = "Occupation share (%)") +
  scale_fill_viridis_d(direction = 1) +
  facet_wrap(~tag, nrow = 2) +
  theme_void() +
  theme(strip.text = element_text(hjust = 0))


# ADM LABOUR FORCE PARTICIPATION SHARES --------------------------------------------------

fig_lfpr <- db_adm %>%
  filter(variable == "lfpr") %>%
  left_join(adm,.) %>%
  left_join(lookup) %>%
  mutate(value = cut(value, breaks = c(0, 0.65, 0.75, 0.85, 1),
                     labels = c("<65", "65-75", "75-85", "85-100"))) %>%
  ggplot() +
  geom_sf(color = "transparent", aes(fill = value)) +
  scale_fill_viridis_d() +
  theme_void() +
  labs(fill = "Labour force\nparticipation rate (%)") +
  facet_wrap(~tag) +
  theme(strip.text = element_text(hjust = 0))


# WORKING POPULATION ---------------------------------------------------------------------

wp_s <- sampleRegular(wp_raw, size=250000, asRaster = TRUE)
names(wp_s) <- "value"
wp_df <- as.data.frame(wp_s, xy = T) %>%
  na.omit() %>%
  mutate(value = cut(value, breaks = c(-1, 100, 500, 5000, +Inf), labels = c("0-100", "100-500", "500-5.000", ">5,000")),
         tag = "h")

fig_wp <- ggplot() +
  geom_tile(data = wp_df, aes(x = x, y = y, fill = value)) +
  geom_sf(data = adm, fill = "transparent", color = "transparent") +
  scale_fill_viridis_d() +
  theme_void() +
  labs(fill = "Working age\npopulation (persons)") +
  facet_wrap(~tag) +
  theme(strip.text = element_text(hjust = 0))


# COMBINE ADM FIGURES -------------------------------------------------------------------

fig_adm <- fig_occ /(fig_lfpr | fig_wp) +
  plot_layout(heights = c(2, 1))


# ========================================================================================
# PREDICTORS -----------------------------------------------------------------------------
# ========================================================================================

# CREATE RASTER --------------------------------------------------------------------------
# Process predictors for plotting. We winsorize all values > 3 sd to 3sd to improve visualization
predictors <- predictors_raw %>%
  pivot_longer(-c(x, y, ID, adm0_code, adm_code), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  mutate(value = as.vector(scale(value))) %>%
  ungroup %>%
  left_join(pred_info, by = c('variable' = 'Predictor variable')) %>%
  dplyr::select(x, y, variable, Topic, value, value) %>%
  mutate(value = case_when(
    value > 3 ~ 3,
    value < -3 ~ -3,
    TRUE ~ value)
  )

# Function to create predictor figure by type
# Raster is reduced to reduce file size
# https://stackoverflow.com/questions/60345163/how-to-generate-high-resolution-temperature-map-using-unevenly-spaced-coordinate
plot_predictor <- function(t, n = 5, s){
  cat("Plotting", t, "\n")
  p <- predictors %>%
    filter(Topic == t) %>%
    mutate(
      x = plyr::round_any(x, 0.01),
      y = plyr::round_any(y, 0.01)
    ) %>%
    unique() %>%
    ggplot() +
    geom_raster(aes(x =x, y = y, fill = value)) +
    coord_sf() +
    labs(fill = "Scaled predictor value", title = t) +
    scale_fill_viridis(direction = -1, na.value = "transparent") +
    facet_wrap(~variable, ncol = n) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "bottom")
  #plot(p)
  return(p)
}

# Plot for each predictor type
pred_type <- unique(predictors$Topic)

fig_trans <- plot_predictor(pred_type[1], 3)
fig_clim <- plot_predictor(pred_type[2], 3)
fig_econ <- plot_predictor(pred_type[3])
fig_night <- plot_predictor(pred_type[4])
fig_urban <- plot_predictor(pred_type[5], 3)
fig_land <- plot_predictor(pred_type[6], 3)
fig_topo <- plot_predictor(pred_type[7])



# ========================================================================================
# PREDICTIONS ----------------------------------------------------------------------------
# ========================================================================================

# OCC ------------------------------------------------------------------------------------
fig_pred_occ <- db_pred %>%
  filter(variable != "lfpr") %>%
  mutate(value = cut(value, breaks = c(-Inf, 0.02, 0.10, 0.5, Inf),
                     labels = c("0-2", "2-10", "10-50", "50-100"))) %>%
  left_join(lookup) %>%
  ggplot() +
    geom_tile(aes(x = x, y = y, fill = value)) +
  scale_fill_viridis_d(direction = 1) +
  labs(x = NULL, y = NULL, fill = "Occupation\nshare (%)") +
  coord_sf() +
  theme_void() +
  facet_wrap(~tag, nrow = 2) +
  theme(strip.text = element_text(hjust = 0))


# LFPR -----------------------------------------------------------------------------------
fig_pred_lfpr <- db_pred %>%
  filter(variable == "lfpr") %>%
  mutate(value = cut(value, breaks = c(0, 0.65, 0.75, 0.85, 1),
                     labels = c("<65", "65-75", "75-85", "85-100"))) %>%
  left_join(lookup) %>%
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = value)) +
  scale_fill_viridis_d(direction = 1) +
  labs(x = NULL, y = NULL, fill = "Labour force\nparticipation rate (%)") +
  coord_sf() +
  theme_void() +
  facet_wrap(~tag, nrow = 1) +
  theme(strip.text = element_text(hjust = 0))


# COMBINE FIGURES ------------------------------------------------------------------------
fig_pred <- plot_grid(fig_pred_occ, fig_pred_lfpr, rel_heights = c(2, 1), ncol = 1)


# ========================================================================================
# PREDICTION ERRORS ----------------------------------------------------------------------
# ========================================================================================

# OCC ------------------------------------------------------------------------------------
fig_pred_error <- db_pred_error %>%
  left_join(lookup) %>%
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = rmspe)) +
  scale_fill_viridis(direction = 1) +
  labs(x = NULL, y = NULL, fill = "Prediction error") +
  coord_sf() +
  theme_void() +
  facet_wrap(~tag, nrow = 2) +
  theme(strip.text = element_text(hjust = 0))


# ========================================================================================
# HYPERPARAMETERS ------------------------------------------------------------------------
# ========================================================================================

# Function to extract hyperparameters for ensemble members
extract_parameter <- function(var, mt){
  cat("Loading", var, "file\n")
  df <- readRDS(file.path(proc_path, glue("results/{results_date}/{var}_ensemble_model.rds")))
  em_list <- db_ensemble_tab %>%
    filter(variable == var,
           member != "ensemble model")

  select_member <- function(mt){
    l <- collect_parameters(df, mt) %>%
      filter(member %in% em_list$member) %>%
      mutate(variable = var) %>%
      rename(weight = coef) %>%
      dplyr::select(variable, everything())
    return(l)
  }

  cat("Selecting hyperparameters\n")
  l <- lapply(mt, select_member)
  l <- l[sapply(l, function(x) dim(x)[1]) > 0]
  l <- bind_rows(l) %>%
    mutate(across(is.numeric, round, digits=3))
  l <- em_list %>%
    filter(transformation == "logit") %>%
    dplyr::select(-transformation, -mae, -mape, -mpe, -weight) %>%
    left_join(l, by = c("member", "variable")) %>%
    arrange(rmse) %>%
    mutate(n = row_number()) %>%
    dplyr::select(-variable) %>%
    dplyr::select(n, member, rmse, rsq, weight, everything())

  return(l)
}

# For each variable, create a list with hyperparameters per model type
model_list <- c("ml_neural_network", "ml_svm_radial", "ml_svm_poly", "ml_random_forest", "ml_xgboost", "ml_glmnet")
tab_hp <- lapply(var_list, extract_parameter, mt = model_list)
names(tab_hp) <- lookup$var_name


# ========================================================================================
# ACCUMULATED EFFECTS PLOTS --------------------------------------------------------------
# ========================================================================================

# Function to make ale plot
plot_ale <- function(var){
  p <- db_ale %>%
    filter(variable == var) %>%
    ggplot() +
    geom_line(aes(x = `_x_`, y = `_yhat_`, color = `_label_`), size = 1) +
    labs(x = NULL, y = "Average prediction", color = NULL) +
    facet_wrap(~`_vname_`, ncol = 4, scales = "free") +
    theme_bw() +
    labs(subtitle = NULL, title = NULL) +
    theme(legend.position = "bottom") +
    guides(color= guide_legend(nrow=2 ,byrow=TRUE)) +
    scale_color_manual(values = rev(col))
  print(p)
}

# Create plots
fig_ale <- lapply(var_list, plot_ale)

# ========================================================================================
# AGGREGATE COMPARISON -------------------------------------------------------------------
# ========================================================================================

# function to create rasters
create_raster <- function(var, df_in, lbl = ""){
  cat("Processing", var ,"\n")
  df <- df_in %>%
    filter(variable == var) %>%
    dplyr::select(x, y, value)
  r <- rasterFromXYZ(df, crs = CRS('+init=EPSG:4326'))
  plot(r, main = var)
  names(r) <- var
  return(r)
}

# Create rasters
pred_stack <- stack(lapply(var_list, create_raster, df_in = db_pred, lbl = "share"))

# Harmonize extent of working population and add to stack
wp <-resample(wp_raw, pred_stack, method = 'bilinear')
names(wp) <- "working_population"
pred_stack <- stack(pred_stack, wp)

# Add ADM info
# Harmonize extent of adm and add to stack
adm_r <- resample(adm_r, pred_stack, method = 'ngb')
names(adm_r) <- "ID"
pred_stack <- stack(pred_stack, adm_r)

# Aggregate predictions
adm_pred <- as.data.frame(pred_stack, xy = TRUE) %>%
  pivot_longer(-c(x, y, working_population, lfpr, ID), names_to = "variable", values_to = "value") %>%
  na.omit %>%
  left_join(adm %>%
              st_drop_geometry()) %>%
  group_by(variable, adm_code, adm_name) %>%
  summarize(
    value = sum(value * working_population * lfpr, na.r = TRUE)/1000,
    .groups = "drop"
  ) %>%
  mutate(source = "ml")

# Calculate adm workers per occupation
adm_ag <- db_adm %>%
  filter(variable != "lfpr") %>%
  dplyr::select(adm_code, adm_name, variable, value) %>%
  left_join(wp_adm) %>%
  group_by(variable, adm_code, adm_name) %>%
  summarize(
    value = sum(value * lfpr * wp, value, na.rm = TRUE)/1000,
    .groups = "drop") %>%
  mutate(source = "census")

# Combine
adm_comp <- bind_rows(
  adm_pred,
  adm_ag
) %>%
  pivot_wider(names_from = source, values_from = value) %>%
  na.omit

# Plot
fig_adm_comp <- adm_comp %>%
  left_join(lookup) %>%
  ggplot(aes(y = ml, x = census)) +
  geom_point(alpha = 0.5) +
  geom_abline(linetype = "dashed", color = col[6], size = 1) +
  labs(x = "Observations (1000 workers)", y = "Predictions (1000 workers)") +
  geom_smooth(method = "lm", color = col[6]) +
  stat_cor(
    aes(label = paste(..rr.label.., sep = "")),
    label.y= Inf, label.x = Inf, vjust = 1, hjust = 2, size = 4
  ) +
  theme(plot.title = element_text(hjust = 1)) +
  facet_wrap(~tag, scales = "free") +
  theme(aspect.ratio = 1) +
  theme_classic2() +
  rremove("grid") +
  theme(strip.text = element_text(hjust = 0),
        strip.background = element_rect(colour="transparent", fill="transparent"))


