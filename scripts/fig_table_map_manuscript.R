# ========================================================================================
# Project:  Occupation_downscaling
# Subject:  Script to prepare figures for manuscript
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
p_load(tidymodels, stacks, sf, raster, viridis, maps, DALEXtra, tictoc, exactextractr)
p_load(viridis, cowplot, ggforce, maps, raster, sf, lessR, patchwork, ggpubr, DALEX, broom,
       stars, egg)

#devtools::install_github("wilkelab/ungeviz")
library(ungeviz)

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
db_test_pred <- readRDS(file.path(proc_path, glue("results/{results_date}/db_test_pred.rds")))
db_test_stat <- readRDS(file.path(proc_path, glue("results/{results_date}/db_test_stat.rds")))
db_test_sum_stat <- readRDS(file.path(proc_path, glue("results/{results_date}/db_test_sum_stat.rds")))
db_ensemble_tab <- readRDS(file.path(proc_path, glue("results/{results_date}/db_ensemble_tab.rds")))


# ========================================================================================
# PREPARE --------------------------------------------------------------------------------
# ========================================================================================

# rename to value for easy of using general functions below
db_pred <- db_pred %>%
  dplyr::rename(value = .pred)

# colorblind palette
col   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#show_col(col)

# List of variables
var_list <- c("off_mgr_pros", "tech_aspros", "service_shop", "agric", "othlowsk", "lfpr", "elementary")

# Get cities in Vietnam and select Hanoi and Ho Chi Minh
data(world.cities)
city <- world.cities %>%
  filter(country.etc == "Vietnam") %>%
  slice_max(pop, n = 2)

# city coordinates for plotting
city_coord <- city %>%
  st_as_sf(coords = c("long", "lat"), crs = raster::crs(adm))

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

# Share of low-skilled occupations
lowskill_sh <- db_adm %>%
  dplyr::filter(variable != "lfpr") %>%
  dplyr::select(adm_code, variable, value) %>%
  left_join(wp_adm) %>%
  group_by(variable) %>%
  summarize(value = sum(value * lfpr * wp, value, na.rm = TRUE)/1000000,
            .groups = "drop") %>%
  mutate(share = value/sum(value, na.rm = TRUE)) %>%
  filter(variable %in% c("agric", "othlowsk", "service_shop", "elementary")) %>%
  summarize(share = sum(share, na.rm = TRUE)) %>%
  pull



# ========================================================================================
# OCC PERSON PREDICTIONS -----------------------------------------------------------------
# ========================================================================================

# LINK ADM AND PREDICTIONS ---------------------------------------------------------------
# function to create rasters
create_raster <- function(var, df_in){
  cat("Processing", var ,"\n")
  df <- df_in %>%
    ungroup() %>%
    filter(variable == var) %>%
    dplyr::select(x, y, value)
  r <- rasterFromXYZ(df, crs = crs(wp_raw))
  #plot(r, main = var)
  names(r) <- var
  return(r)
}

# Create rasters
pred_stack <- stack(lapply(var_list, create_raster, df_in = db_pred))

# Harmonize extent of working population and add to stack
wp <-resample(wp_raw, pred_stack, method = 'bilinear')
names(wp) <- "working_population"
pred_stack <- stack(pred_stack, wp)

# Harmonize extent of adm and add to stack
adm_r <-resample(adm_r, pred_stack, method = 'ngb')
names(adm_r) <- "ID"
pred_stack <- stack(pred_stack, adm_r)

# Calculate occupation in persons by multiplying working population, occ share and lfpr
occ_person <- as.data.frame(pred_stack, xy = TRUE) %>%
  pivot_longer(-c(x, y, working_population, lfpr, ID), names_to = "variable", values_to = "occ_share") %>%
  mutate(value = occ_share * lfpr * working_population) %>%
  na.omit

adm_pred <- as.data.frame(pred_stack, xy = TRUE) %>%
  pivot_longer(-c(x, y, working_population, lfpr, ID), names_to = "variable", values_to = "value") %>%
  na.omit %>%
  left_join(adm %>%
              st_drop_geometry()) %>%
  group_by(variable, adm_code, adm_name) %>%
  summarize(
    value = sum(value * working_population * lfpr, na.r = TRUE),
    .groups = "drop"
  ) %>%
  mutate(source = "ml")


# PLOT -----------------------------------------------------------------------------------

# Function to create plot per variable
plot_variable <- function(var){
  h_w <- 1
  city_rec <- city %>%
    mutate(xmin = long-h_w,
           xmax = long+h_w,
           ymin = lat-h_w,
           ymax = lat+h_w)
  df <- occ_person %>%
    filter(variable == var)

  p_vnm <- df %>%
    mutate(value = cut(value, breaks = c(-Inf, 20, 100, 500, 1000, Inf),
                       labels = c("0-20", "20-100", "100-500", "500-1,000", ">1,000"))) %>%
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = value)) +
    geom_rect(data = city_rec, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color = "red", fill = "transparent") +
    scale_fill_viridis_d(direction = 1) +
    labs(x = NULL, y = NULL, fill = "Workers per grid cell") +
    coord_sf() +
    theme_void()+
    theme(legend.position = "none")


  p_hanoi <- df %>%
    mutate(value = cut(value, breaks = c(-Inf, 20, 100, 500, 1000, Inf),
                       labels = c("0-20", "20-100", "100-500", "500-1,000", ">1,000"))) %>%
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = value)) +
    scale_fill_viridis_d(direction = 1) +
    coord_sf(xlim = c(city$long[2]-h_w, city$long[2]+h_w),
             ylim = c(city$lat[2]-h_w, city$lat[2]+h_w),
             expand = FALSE) +
    theme_void() +
    theme(legend.position = "none") +
    theme(plot.margin = margin(0, 0.1, 0.1, 0, "cm"))

  p_hcmc <- df %>%
    mutate(value = cut(value, breaks = c(-Inf, 20, 100, 500, 1000, Inf),
                       labels = c("0-20", "20-100", "100-500", "500-1,000", ">1,000"))) %>%
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = value)) +
    scale_fill_viridis_d(direction = 1) +
    coord_sf(xlim = c(city$long[1]-h_w, city$long[1]+h_w),
             ylim = c(city$lat[1]-h_w, city$lat[1]+h_w),
             expand = FALSE) +
    theme_void() +
    theme(legend.position = "none") +
    theme(plot.margin = margin(0.1, 0.1, 0, 0, "cm"))
  p <- plot_grid(p_vnm, p_hanoi, p_hcmc, ncol = 1, rel_heights = c(1, 0.5, 0.5))
  p
}

# Compose plot
p_list <- lapply(var_list[var_list != "lfpr"], plot_variable)
p <- (p_list[[1]] | p_list[[2]] | p_list[[3]] | p_list[[4]] | p_list[[5]] | p_list[[6]]) +
  plot_annotation(tag_levels = 'a')

# Create figure from which legend can be stripped
fig_legend <- occ_person %>%
  filter(variable == "agric") %>%
  mutate(value = cut(value, breaks = c(-Inf, 20, 100, 500, 1000, Inf),
                     labels = c("0-20", "20-100", "100-500", "500-1,000", ">1,000"))) %>%
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = value)) +
  scale_fill_viridis_d(direction = 1) +
  labs(x = NULL, y = NULL, fill = "Workers per grid cell") +
  coord_sf() +
  theme(legend.position = "bottom")

legend <- get_legend(
  # create some space to the left of the legend
  fig_legend +  guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom") +
    theme(legend.box.margin = margin(0, 0, 0, 12))
)

# Create composed plot with legend
fig_occ_person <- plot_grid(p, legend, ncol = 1, rel_heights = c(1, .1))


# ========================================================================================
# SUMMARY OF MODEL PERFORMANCE -----------------------------------------------------------
# ========================================================================================

# PLOT -----------------------------------------------------------------------------------
# Function to create rmse plot
plot_rmse <- function(var,  trans = "logit"){
  df_rmse <- db_ensemble_tab %>%
    pivot_longer(-c(member, transformation, variable, weight), names_to = "metric", values_to = "value") %>%
    filter(transformation == trans, metric %in% c("rmse"), variable == var) %>%
    mutate(metric = "")

  p_rmse <- ggplot(data = df_rmse, aes(x = metric, y = value)) +
    geom_hpline(color = "grey50", size = 1, alpha = 0.5, width = 0.9) +
    geom_hpline(data = filter(df_rmse, member == "ensemble model"),
                color = "red", size = 1.5,  width = 1) +
    labs(x = "RMSE", y = NULL) +
    theme_pubr() +
    rremove("grid") +
    rremove("x.ticks")
  return(p_rmse)
}

# Function to create rsq plot
plot_rsq <- function(var, trans = "logit"){
  df_rsq <- db_ensemble_tab %>%
    pivot_longer(-c(member, transformation, variable, weight), names_to = "metric", values_to = "value") %>%
    filter(transformation == trans, metric %in% c("rsq"), variable == var) %>%
    mutate(metric = "")

  p_rsq <- ggplot(data = df_rsq, aes(x = metric, y = value)) +
    geom_hpline(color = "grey50", size = 1, alpha = 0.5, width = 0.9) +
    geom_hpline(data = filter(df_rsq, member == "ensemble model"),
                color = "red", size = 1.5,  width = 1) +
    labs(x = "R2", y = NULL) +
    theme_pubr() +
    rremove("grid") +
    rremove("x.ticks")
  return(p_rsq)
}

# Function to plot compare observations versus predictions
plot_compare <- function(var, min, max,  trans = "logit"){
  df_bias <- bind_rows(db_test_pred) %>%
    dplyr::select(variable, value, .pred) %>%
    left_join(
      member_range <- bind_rows(db_test_pred) %>%
        pivot_longer(-c(variable, .pred, value), names_to = "member", values_to = "member.pred") %>%
        na.omit %>%
        group_by(variable, value, .pred) %>%
        summarize(
          lower = min(member.pred,na.rm = TRUE),
          upper = max(member.pred, na.rm = TRUE),
          n = n(),
          .groups = "drop"
        )
    ) %>%
    filter(variable == var)

  p_compare <- ggplot(data = df_bias, aes(x = value, y = .pred)) +
    geom_point(alpha = 0.5) +
    geom_abline(linetype = "dashed", color = col[6], size = 1) +
    labs(x = "Observations", y = "Prediction") +
    scale_x_continuous(limits = c(min, max)) +
    scale_y_continuous(limits = c(min, max)) +
    # geom_smooth(method = "lm", color = col[6]) +
    # stat_cor(
    #   aes(label = paste(..rr.label.., sep = "")),
    #   label.y= Inf, label.x = Inf, vjust = 1, hjust = 2.5, size = 4
    # ) +
    theme_pubr() +
    rremove("grid") +
    theme(plot.title = element_text(hjust = 1))

  return(p_compare)
}

# Create rmse and rsq plots
p_rmse <- lapply(var_list, plot_rmse)
p_rsq <- lapply(var_list, plot_rsq)

# Create comparison plots and set range
p_agric <- plot_compare(var_list[4], -7.5, 2.5)
p_off_mgr_pros <- plot_compare(var_list[1], -4.5, -1)
p_othlowsk <- plot_compare(var_list[5], -5, 0.5)
p_service_shop <- plot_compare(var_list[3], -4, 0)
p_tech_aspros <- plot_compare(var_list[2], -4.5, -2.5)
p_lfpr <- plot_compare(var_list[6], 0, 3)
p_elementary <- plot_compare(var_list[7], -6.5, -2.5)

# Combine into one plot
fig_performance <- plot_grid(p_off_mgr_pros, p_rmse[[1]], p_rsq[[1]], p_tech_aspros, p_rmse[[2]], p_rsq[[2]],
          p_service_shop, p_rmse[[3]], p_rsq[[3]], p_agric,  p_rmse[[4]], p_rsq[[4]],
          p_othlowsk, p_rmse[[5]], p_rsq[[5]], p_elementary, p_rmse[[7]], p_rsq[[7]], p_lfpr, p_rmse[[6]], p_rsq[[6]],
          labels = "auto", align = "hv", axis = "tblr", ncol = 6, rel_widths = c(1, 0.3, 0.3, 1, 0.3, 0.3))

# Calculate mean error
mean_error <- bind_rows(db_test_pred) %>%
  dplyr::select(variable, value, .pred) %>%
  group_by(variable) %>%
  summarize(e = mean(value-.pred, na.rm = TRUE))

# Check super learner performance in comparison to other models using RMSE
sl_performance <- db_ensemble_tab %>%
  filter(transformation == "logit") %>%
  group_by(variable) %>%
  slice_min(rmse, n = 3)

# Which are best performing ml models per variable
db_ensemble_tab %>%
  filter(transformation == "logit", member != "ensemble model") %>%
  group_by(variable) %>%
  slice_min(rmse, n = 1)


# ========================================================================================
# ACCURACY MEASURES ENSEMBLE MODEL -------------------------------------------------------
# ========================================================================================

# Combine summary statistics and accuracy measures
tab_ac_stat <- db_test_sum_stat %>%
 filter(transformation == "logit") %>%
  left_join(
    db_test_stat %>%
      dplyr::filter(member == ".pred",
                    indicator %in% c("mae", "rmse", "rsq"),
                    transformation == "logit") %>%
      pivot_wider(
        names_from = c(indicator),
        values_from = .estimate)
  ) %>%
  dplyr::select(-c(transformation, .estimator, member)) %>%
  mutate(rmse = round(rmse,2),
         rsq = round(rsq, 2))


# ========================================================================================
# COMPARE WITH WEALTH INDEX --------------------------------------------------------------
# ========================================================================================

# OCC SHARES -----------------------------------------------------------------------------

# Wealth index
wi_raw <- readRDS(file.path(proc_path, "wealth_index/wealth_index.rds"))

# Stack occupation share maps
files <- list.files(file.path(proc_path, glue("results/{results_date}/maps")), pattern = "share\\.tif$", full.names = TRUE)
occ_share <- read_stars(files) %>%
  merge

# Warp occ_share to wi_raw, split stack again to add the wealth index layer
occ_share_w <- st_warp(occ_share, wi_raw) %>%
  split
names(occ_share_w) <- gsub("_share.tif$", "", names(occ_share_w))

# Add wealth_index and convert to data.frame where rwi is a dependent variable
wi_s <- c(occ_share_w, wi_raw) %>%
  merge %>%
  as.data.frame %>%
  na.omit %>%
  pivot_wider(names_from = attributes, values_from = X) %>%
  na.omit %>%
  pivot_longer(-c(x, y, rwi), names_to = "variable", values_to = "value")

# Regress rwi on occupation shares
wi_reg <- wi_s %>%
  nest(data = -variable) %>%
  mutate(
    fit = purrr::map(data, ~ lm(rwi ~ value, data = .x)),
    #tidied = purrr::map(fit, tidy),
    glanced = purrr::map(fit, glance),
  ) %>%
  dplyr::select(-data, -fit) %>%
  unnest(glanced)

# Adj R2 for low-skilled occupations, which have best fit
r2_ls <- wi_reg %>%
  filter(variable %in% c("agric", "othlowsk", "service_shop", "elementary")) %>%
  dplyr::select(adj.r.squared) %>%
  mutate(adj.r.squared = adj.r.squared*100) %>%
  pull()


# Function to plot wi
plot_wi <- function(var, x_pos){
  cat("Plotting", var, "\n")
  df <- wi_s %>%
    filter(variable == var)
  p <- ggplot(data = df, aes(x = value, y = rwi)) +
      #geom_point(alpha = 0.5) +
      geom_hex(color = "transparent", bins = 100) +
      scale_fill_viridis() +
      stat_cor(
        aes(label = paste(..rr.label.., sep = "")),
        label.y= Inf, label.x = Inf, vjust = 1, hjust = 3, size = 4
      ) +
      scale_y_continuous(limits = c(-2, 2.5)) +
      labs(x = "Occupation share", y = "Relative wealth index") +
      geom_smooth(method = lm, se = TRUE, color = "grey70") +
      theme_classic2() +
      rremove("grid") +
      theme(legend.position = c(x_pos, 0.2)) +
    theme(legend.key.height= unit(0.35, 'cm'),
          legend.key.width= unit(0.35, 'cm'))
  return(p)
}

# Plot
fig_wi <- map2(var_list[var_list != "lfpr"], c(0.8, 0.8, 0.8, 0.2, 0.8, 0.8), plot_wi)
fig_wi <- plot_grid(plotlist = fig_wi, nrow = 2, labels = "auto")



