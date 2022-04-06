# ========================================================================================
# Project:  occupation downscaling
# Subject:  Post process machine learning model results
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
p_load(tidymodels, stacks, sf, raster, viridis, maps, DALEXtra, tictoc, cowplot, ranger,
       forestError)

# Set path
source(here("scripts/set_path.r"))

# R options
options(scipen = 999)
options(digits = 4)


# ========================================================================================
# SET DATE -------------------------------------------------------------------------------
# ========================================================================================

# Ensemble model statistics for all occupation
results_date <- "2022-02-12"


# ========================================================================================
# LOAD DATA ------------------------------------------------------------------------------
# ========================================================================================

# Model dependent and predictor data
source(here("scripts/prepare_input_data.r"))

# Grid level predictors
predictors_raw <- readRDS(file.path(proc_path, "model_input/VNM_predictors.rds")) %>%
  ungroup

# Adm
adm_raw <- readRDS(file.path(proc_path, glue("adm/VNM_adm.rds")))

# Load working age population
wp_raw <- raster(file.path(proc_path, glue("worldpop/VNM_working_population_2009.tif")))

# List of predictors and groupings
pred_info <- read_excel(file.path(proc_path, glue("predictors/tables.xlsx")), sheet = "predictors")


# ========================================================================================
# PREPARE --------------------------------------------------------------------------------
# ========================================================================================

# Set variables
var_list <- c("off_mgr_pros", "tech_aspros", "service_shop", "agric", "othlowsk", "lfpr", "elementary")

# Prepare predictor file for predictions
predictors <- predictors_raw %>%
  left_join(adm_raw %>%
              st_drop_geometry())

# Create lookup table
lookup <- data.frame(
  variable = var_list,
  var_name = c("Managers and professionals", "Technicians and associate professionals", "Clerks and service workers",
               "Agricultural workers","Craft workers and operators", "Labour force participation rate", "Elementary occupations"),
  tag = c("a", "b", "c", "d", "e", "g", "f")
)


# ========================================================================================
# CREATE RASTER FROM PREDICTORS AND SAVE AS TIF ------------------------------------------
# ========================================================================================

# function to create rasters and save
create_raster <- function(var, df_in, lbl = ""){
  cat("Processing", var ,"\n")
  df <- df_in %>%
    filter(variable == var) %>%
    dplyr::select(x, y, value)
  r <- rasterFromXYZ(df, crs = CRS('+init=EPSG:4326'))
  plot(r, main = var)
  names(r) <- var
  writeRaster(r, file.path(proc_path, glue("results/{results_date}/maps/{var}_{lbl}.tif")), overwrite = TRUE)
  return(r)
}

# Process predictors for plotting. We winsorize all values > 3 sd to 3sd to improve visualization
db_predictors <- predictors_raw %>%
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

# Create rasters and save
walk(unique(db_predictors$variable), create_raster, df_in = db_predictors)


# ========================================================================================
# SPATIAL PREDICTIONS --------------------------------------------------------------------
# ========================================================================================

# CREATE PREDICTIONS ---------------------------------------------------------------------
# Function to create predictions
predict_var <- function(var){
  cat("Loading", var, "file\n")
  ensemble_model <- readRDS(file.path(proc_path, glue("results/{results_date}/{var}_ensemble_model.rds")))
  cat("Predicting", var, "\n")
  predictors <- predictors %>%
    mutate(variable = var)
  pr <- predict(ensemble_model, predictors)
  pr <- bind_cols(pr, predictors) %>%
    dplyr::select(adm_code, adm_name, x, y, .pred) %>%
    mutate(variable = var)
  return(pr)
}

# Create predictions database for all variables
tic()
db_pred <- map_df(var_list, predict_var)
toc()

# Apply inverse logit
db_pred <- db_pred %>%
  mutate(.pred = exp(.pred)/(1+exp(.pred)))

# Check if occ shares sum to one
occ_share_check <- db_pred %>%
  filter(variable != "lfpr") %>%
  group_by(x, y, adm_code, adm_name) %>%
  summarize(value = sum(.pred, na.rm = TRUE),
            .groups = "drop")
summary(occ_share_check$value)

# Rebalance so shares add up to 1
db_pred_s <- db_pred %>%
  filter(variable != "lfpr") %>%
  group_by(x, y, adm_code, adm_name) %>%
  mutate(.pred = .pred/sum(.pred, na.rm = TRUE)) %>%
  ungroup()

# Add lfpr
db_pred_s <- bind_rows(
  db_pred_s,
  db_pred %>%
    filter(variable == "lfpr")
)

# Save
saveRDS(db_pred_s, file.path(proc_path, glue("results/{results_date}/db_pred.rds")))


# CREATE RASTER FROM PREDICTIONS ---------------------------------------------------------
# Create maps folder
dir.create(file.path(proc_path, glue("results/{results_date}/maps")), showWarnings = FALSE, recursive = TRUE)

# rename to value for easy of using general functions below
db_pred_s <- db_pred_s %>%
  rename(value = .pred)

# Create rasters
pred_stack <- stack(lapply(var_list, create_raster, df_in = db_pred_s, lbl = "share"))

# Harmonize extent of working population, save and add to stack
wp <-projectRaster(wp_raw, pred_stack, method = 'bilinear')
names(wp) <- "working_population"
writeRaster(wp, file.path(proc_path, glue("results/{results_date}/maps/working_population.tif")), overwrite = TRUE)
pred_stack <- stack(pred_stack, wp)

# Calculate occupation in persons by multiplying working population, occ share and lfpr
occ_person <- as.data.frame(pred_stack, xy = TRUE) %>%
  pivot_longer(-c(x, y, working_population, lfpr), names_to = "variable", values_to = "occ_share") %>%
  mutate(value = occ_share * lfpr * working_population) %>%
  na.omit

# Create occ_person rasters
occ_person_stack <- lapply(var_list[var_list != "lfpr"], create_raster, df_in = occ_person, lbl = "persons")

# Clean up
rm(wp, wp_raw, pred_stack, occ_person, occ_person_stack)


# ========================================================================================
# PREDICTION INTERVAL --------------------------------------------------------------------
# ========================================================================================

# CREATE INTERVALS -----------------------------------------------------------------------
pred_interval <- function(var){
  cat("Loading", var, "file\n")
  ensemble_model <- readRDS(file.path(proc_path, glue("results/{results_date}/{var}_ensemble_model.rds")))
  em <- ensemble_model$data_stack
  w <- stacks:::top_coefs(ensemble_model, n = Inf) %>%
    dplyr::select(-type)

  cat("Predicting", var, "\n")
  predictors <- predictors %>%
    mutate(variable = var)
  pr <- predict(ensemble_model, predictors, member = TRUE) %>%
    dplyr::select(-.pred)

  train <- em %>%
    dplyr::select(c("value", w$member))
  eml <- ranger::ranger(value ~ ., train, num.trees = 85, importance = "impurity",
                        quantreg = TRUE, keep.inbag = TRUE)

  quantiles = c((1-.682)/2, 1-(1-.682)/2) # note that we construct 1 sd intervals (=68%)
  n.cores = parallel::detectCores()
  xtr <- train %>%
    dplyr::select(-value)
  ytr <- train %>%
    dplyr::select(value) %>% pull
  cat("Calculate mean square prediction error", var, "\n")
  pred_e = forestError::quantForestError(eml,
                                         X.train = xtr,
                                         X.test = pr,
                                         Y.train = ytr,
                                         alpha = (1-(quantiles[2]-quantiles[1])), n.cores=n.cores)
  error <- predictors %>%
    mutate(rmspe =  sqrt(pred_e$estimates$mspe),
           l =  pred_e$estimates$lower_0.318,
           u =  pred_e$estimates$upper_0.318,
           interval = u-l,
           variable = var) %>%
    dplyr::select(variable, adm_code, ID, x, y, rmspe, l, u, interval)
  return(error)
}

db_pred_error <- map_df(var_list, pred_interval)

# Save
saveRDS(db_pred_error, file.path(proc_path, glue("results/{results_date}/db_pred_error.rds")))


# CREATE RASTER FROM ERRORS ---------------------------------------------------------
# select rmspe
db_pred_error <- db_pred_error %>%
  rename(value = rmspe)

# Create rasters
walk(unique(db_pred_error$variable), create_raster, df_in = db_pred_error, lbl = "error")


# ========================================================================================
# MODEL FIT ------------------------------------------------------------------------------
# ========================================================================================

# CREATE DATABASE WITH PREDICTIONS FOR TEST SET ------------------------------------------
# Function to create predictions database
create_pred <- function(var){
  cat("Loading", var, "file\n")
  ensemble_model <- readRDS(file.path(proc_path, glue("results/{results_date}/{var}_ensemble_model.rds")))
  db_test <- readRDS(file.path(proc_path, glue("results/{results_date}/{var}_db_test.rds")))
  cat("Predicting", var, "file\n")
  df <- db_test %>%
    dplyr::select(value) %>%
    bind_cols(predict(ensemble_model, db_test, members = TRUE)) %>%
    mutate(variable = var) %>%
    dplyr::select(variable, everything())
  return(df)
}

# Create database for all variables
tic()
db_test_pred <- lapply(var_list, create_pred)
toc()

# Save
saveRDS(db_test_pred, file.path(proc_path, glue("results/{results_date}/db_test_pred.rds")))


# CREATE DATABASE WITH ACCURACY STATISTICS ------------------------------------------------

# Function to create accuracy statistics
calc_fit <- function(df){
  var <- unique(df$variable)
  df <- df %>%
    dplyr::select(-variable)
  df_fit <-  bind_rows(
    map_dfr(df, rmse, truth = value, data = df) %>%
      mutate(member = colnames(df)),
    map_dfr(df, mae, truth = value, data = df) %>%
      mutate(member = colnames(df)),
    map_dfr(df, mpe, truth = value, data = df) %>%
      mutate(member = colnames(df)),
    map_dfr(df, mape, truth = value, data = df) %>%
      mutate(member = colnames(df)),
    map_dfr(df, rsq, truth = value, data = df) %>%
      mutate(member = colnames(df))
  ) %>%
    rename(indicator = .metric) %>%
    arrange(indicator, .estimate) %>%
    mutate(variable = var)
  return(df_fit)
}

# Create database with statistics using logit transformed variables
db_test_stat_l <- map_df(db_test_pred, calc_fit)

# Function to apply inverse logit to db_test_pred
calc_il <- function(df){
  df <- df %>%
    mutate(.pred = exp(.pred)/(1+exp(.pred)),
           value = exp(value)/(1+exp(value)))
  return(df)
}

# Apply inverse logit to create 'normal' predictions
db_test_pred_n <- lapply(db_test_pred, calc_il)

# Create database with statistics for all variables after inverse logit
db_test_stat_n <- map_df(db_test_pred_n, calc_fit)

# Combine test stat
db_test_stat <- bind_rows(
  db_test_stat_l %>%
    mutate(transformation = "logit"),
  db_test_stat_n %>%
    mutate(transformation = "none")
)

# Save
saveRDS(db_test_stat, file.path(proc_path, glue("results/{results_date}/db_test_stat.rds")))

# CREATE DATABASE WITH TEST SET SUMMARY STATISTICS ---------------------------------------

# Function to create summary statistics
create_sum_stat <- function(var){
  db_test <- readRDS(file.path(proc_path, glue("results/{results_date}/{var}_db_test.rds")))
  df <- bind_rows(
    db_test %>%
      group_by(variable) %>%
      summarize(
        n = n(),
        min = min(value, na.rm = TRUE),
        mean = mean(value, na.rm = TRUE),
        max = max(value, na.rm = TRUE),
        .groups = "drop") %>%
      mutate(transformation = "logit"),
    db_test %>%
      mutate(value = exp(value)/(1+exp(value))) %>%
      group_by(variable) %>%
      summarize(
        n = n(),
        min = min(value, na.rm = TRUE),
        mean = mean(value, na.rm = TRUE),
        max = max(value, na.rm = TRUE),
        .groups = "drop") %>%
      mutate(transformation = "none")
  )
    return(df)
}

# Create database for all variables
db_test_sum_stat <- map_df(var_list, create_sum_stat)

# Save
saveRDS(db_test_sum_stat, file.path(proc_path, glue("results/{results_date}/db_test_sum_stat.rds")))


# ========================================================================================
# ENSEMBLE MODEL SPECIFICATIONS ----------------------------------------------------------
# ========================================================================================

# Function to create table with ensemble model and member statistics
create_ensemble_table <- function(var){
  cat("Loading", var, "file\n")
  ensemble_model <- readRDS(file.path(proc_path, glue("results/{results_date}/{var}_ensemble_model.rds")))
  w <- stacks:::top_coefs(ensemble_model, n = Inf) %>%
    dplyr::select(-type)
  tab <- db_test_stat %>%
    filter(variable == var) %>%
    dplyr::filter(member != "value") %>%
    mutate(member = ifelse(member == ".pred", "ensemble model", member)) %>%
    dplyr::select(-.estimator) %>%
    pivot_wider(names_from = indicator, values_from = .estimate) %>%
    left_join(w) %>%
    arrange(rmse)
  return(tab)
}

# Create database for all variables
tic()
db_ensemble_tab <- map_df(var_list, create_ensemble_table)
toc()

# Save
saveRDS(db_ensemble_tab, file.path(proc_path, glue("results/{results_date}/db_ensemble_tab.rds")))


# ========================================================================================
# VARIABLE IMPORTANCE PLOTS --------------------------------------------------------------
# ========================================================================================

# https://www.tmwr.org/explain.html
# https://ema.drwhy.ai/featureImportance.html
# https://christophm.github.io/interpretable-ml-book/feature-importance.html
# Variable importance plot when applied in a model agnostic way is referred to as feature importance

# Function to create variable importance measures
create_vim <- function(var, n_model, n_predictors){
  cat("Loading", var, "file\n")
  ensemble_model <- readRDS(file.path(proc_path, glue("results/{results_date}/{var}_ensemble_model.rds")))
  db_test <- readRDS(file.path(proc_path, glue("results/{results_date}/{var}_db_test.rds")))
  w <- stacks:::top_coefs(ensemble_model, n = Inf) %>%
    arrange(desc(weight)) %>%
    slice_head(n = n_model)

  create_vim_member <- function(member_name, var, n_predictors){
    cat("Create variable importance measures", member_name, "\n")
    member_model <- ensemble_model$member_fits[names(ensemble_model$member_fits) %in% member_name][[1]]
    explainer <-
      explain_tidymodels(
        member_model,
        data = db_test %>%
          dplyr::select(-value),
        y = db_test$value,
        label = member_name,
        verbose = FALSE)
    df <- model_parts(explainer, loss_function = loss_root_mean_square)
    df <- plot(df, max_vars = n_predictors)
    return(df)
  }
  db_vim <- w %>%
    group_by(member) %>%
    mutate(vim = purrr::map(member, create_vim_member, var = var, n_predictors = n_predictors),
           variable = var)
  return(db_vim)
}

# Create database for all variables
tic()
db_vim <- map_df(var_list, create_vim, n_model = 5, n_predictors = 10)
toc()

# Save
# db_vim is extremely large so we do not save and create relevant figures here
#saveRDS(db_vim, file.path(proc_path, glue("results/{results_date}/db_vim.rds")))

# Function to create top n vim figures per var ranked by ensemble model weight
create_vim <- function(var, n){
  df <- db_vim %>%
    filter(variable == var) %>%
    arrange(desc(weight))

  plot_vim <- function(x){
    weight <- df$weight[x]
    p <- df$vim[[x]] +
      theme_bw()  +
      labs(title = NULL, subtitle = glue("weight is {round(weight, 3)}"), NULL, y = "RMSE loss") +
      guides(color = "none")
  }

  p <- lapply(c(1:n), plot_vim)

  title <- ggdraw() +
    draw_label(
      var,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 250)
    )

  p <- cowplot::plot_grid(plotlist = p, ncol = 2, axis = "lr", align = "hv")
  plot_grid(p, rel_heights = c(1), ncol = 1)
}

# Create plots for each variable
fig_vim_variable <- lapply(var_list, create_vim, 5)
names(fig_vim_variable) <- lookup$var_name

# Combine vim plots for ensemble model members with highest weight across variables
p_vim_hw <- db_vim %>%
  group_by(variable) %>%
  slice_head(n= 1)

# Function to take first vim for each model
plot_vim <- function(var){
  df <- p_vim_hw %>%
    filter(variable == var)
  p <- df$vim[[1]] +
    theme_bw()  +
    labs(title = var, subtitle = NULL, y = "RMSE loss") +
    guides(color = "none")
  p
}

# Plot
fig_vim <- lapply(var_list, plot_vim)
names(fig_vim) <- lookup$var_name
fig_vim <- cowplot::plot_grid(plotlist = fig_vim, ncol = 3, axis = "lr", align = "hv")

# Save
saveRDS(fig_vim, file.path(proc_path, glue("results/{results_date}/fig_vim.rds")))
saveRDS(fig_vim_variable, file.path(proc_path, glue("results/{results_date}/fig_vim_variable.rds")))
rm(fig_vim, fig_vim_variable)


# ========================================================================================
# ACCUMULATED LOCAL EFFECTS PLOTS---------------------------------------------------------
# ========================================================================================

# Function to create local accumulation plots
create_ale <- function(var, n_model = 5){
  cat("Loading", var, "file\n")
  ensemble_model <- readRDS(file.path(proc_path, glue("results/{results_date}/{var}_ensemble_model.rds")))
  db_test <- readRDS(file.path(proc_path, glue("results/{results_date}/{var}_db_test.rds")))
  w <- stacks:::top_coefs(ensemble_model, n = Inf) %>%
    arrange(desc(weight)) %>%
    slice_head(n = n_model)

  create_ale_member <- function(member_name, var){
    cat("Create accumulated local effects plot", member_name, "\n")
    member_model <- ensemble_model$member_fits[names(ensemble_model$member_fits) %in% member_name][[1]]
    explainer <-
      explain_tidymodels(
        member_model,
        data = db_test %>%
          dplyr::select(-value),
        y = db_test$value,
        label = member_name,
        verbose = FALSE)
    df <- model_profile(explainer, type = "accumulated")$agr_profiles %>%
      mutate(variable = var)
    return(df)
  }
  db_ale <- map_df(w$member, create_ale_member, var = var)
  return(db_ale)
}

tic()
db_ale <- map_df(var_list, create_ale, n_model = 5)
toc()

# Save
saveRDS(db_ale, file.path(proc_path, glue("results/{results_date}/db_ale.rds")))
