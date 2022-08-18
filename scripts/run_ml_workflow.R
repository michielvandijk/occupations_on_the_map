# ========================================================================================
# Project:  KB project, people on the map
# Subject:  Machine learning model
# Author:   Michiel van Dijk; Thijs de Lange
# Contact:  michiel.vandijk@wur.nl; thijs.delange@wur.nl
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
p_load(tidymodels, stacks, tictoc)

# Set path
source(here("set_path.r"))

# R options
options(scipen = 999)
options(digits = 4)

# Set tidymodels as standard to avoid package conflicts
tidymodels_prefer()


# ========================================================================================
# SOURCES --------------------------------------------------------------------------------
# ========================================================================================

# https://ema.drwhy.ai/breakDown.html
# https://www.tmwr.org/
# https://koalaverse.github.io/machine-learning-in-R/
# https://medium.com/nerd-for-tech/extrapolation-is-tough-for-trees-tree-based-learners-combining-learners-of-different-type-makes-659187a6f58d
# https://cran.r-project.org/web/packages/stacks/vignettes/basics.html


# ========================================================================================
# LOAD DATA ------------------------------------------------------------------------------
# ========================================================================================

# Model dependent and predictor data
source(here("scripts/prepare_input_data.r"))


# ========================================================================================
# PREPARE DATA ---------------------------------------------------------------------------
# ========================================================================================

# Remove adm identifiers
db <- db_adm %>%
  dplyr::select(!c(
    "adm1_code", "adm1_name", "adm0_code"
  ))

# Logit transformation
summary(db$value)
db <- db %>%
  mutate(value = log(value/(1-value)))
summary(db$value)


# ========================================================================================
# DEFINE MODELS --------------------------------------------------------------------------
# ========================================================================================

# RANDOM FOREST -------------------------------------------------------------------------
rf_spec <-
  rand_forest(mtry = tune(), min_n = tune(), trees = 1000) %>%
  set_mode("regression") %>%
  set_engine("ranger", regularization.factor = tune("regularization"))

# The mtry parameter does not have an upper bound as this depends on the number of predictors
# We set it manually to the number of columns in the database - 3, two adm variables and value
parameters(rf_spec) %>%
  pull_dials_object("mtry")

rf_param <-
  rf_spec %>%
  parameters() %>%
  update(mtry = mtry(c(1, ncol(db)-3)))


# XGBOOST --------------------------------------------------------------------------------
xgb_spec <-
  boost_tree(trees = tune(), min_n = tune(), tree_depth = tune(), learn_rate = tune(),
             loss_reduction = tune(), sample_size = tune()) %>%
  set_mode("regression") %>%
  set_engine("xgboost")


# GLMNET ---------------------------------------------------------------------------------
glmnet_spec <-
  linear_reg(penalty = tune(), mixture = tune()) %>%
  set_mode("regression") %>%
  set_engine("glmnet")


# NEURAL NETWORK -------------------------------------------------------------------------
nnet_spec <-
  mlp(hidden_units = tune(), penalty = tune(), epochs = tune()) %>%
  set_mode("regression") %>%
  set_engine("nnet")

# The analysis in Kuhn and Johnson (2013) specifies that the neural network should have up
# to 27 hidden units in the layer. The parameters() function creates a parameter set that
# we modify to have the correct parameter range.
nnet_param <-
  nnet_spec %>%
  parameters() %>%
  update(hidden_units = hidden_units(c(1, 27)))


# RADIAL BASIS SUPPORT VECTOR MACHINE ----------------------------------------------------
svm_rbf_spec <-
  svm_rbf(cost = tune(), rbf_sigma = tune()) %>%
  set_mode("regression") %>%
  set_engine("kernlab")


# POLYNOMINAL SUPPORT VECTOR MACHINE -----------------------------------------------------
svm_poly_spec <-
  svm_poly(cost = tune(), degree = tune()) %>%
  set_mode("regression") %>%
  set_engine("kernlab")


# ========================================================================================
# RUN ML MODELS --------------------------------------------------------------------------
# ========================================================================================

# Function to tune range of machine learning models
run_ml <- function(occ_s, n_grid, k_fold, n_rep, prop_value){
  cat("running ml", occ_s, "\n")

  library(doParallel)
  all_cores <- parallel::detectCores(logical = FALSE)
  cl <- makePSOCKcluster(all_cores-1)
  registerDoParallel(cl)

  set.seed(100)
  db <- filter(db, variable == occ_s)
  db_split <- initial_split(db, strata = value, prop = prop_value)
  db_train <- training(db_split)
  db_test <- testing(db_split)
  saveRDS(db_test, file.path(result_path, glue("{occ_s}_db_test.rds")))
  saveRDS(db_train, file.path(result_path, glue("{occ_s}_db_train.rds")))

  folds <- vfold_cv(db_train, v = k_fold, repeats = n_rep, strata = value)

  main_recipe <- recipe(formula = value ~ ., data = db_train) %>%
    update_role(adm_code, new_role = "adm_code") %>%
    update_role(adm_name, new_role = "adm_name") %>%
    update_role(variable, new_role = "variable") %>%
    step_YeoJohnson(all_numeric_predictors()) %>%
    step_normalize(all_numeric_predictors())  %>%
    step_corr(all_predictors(), threshold = .7) # Set correlation to 0.7

  all_workflows <-
    workflow_set(
      preproc = list(ml = main_recipe),
      models = list(svm_radial = svm_rbf_spec, svm_poly = svm_poly_spec,
                    neural_network = nnet_spec, random_forest = rf_spec,
                    xgboost = xgb_spec, glmnet = glmnet_spec)
    ) %>%
    option_add(param_info = nnet_param, id = "ml_neural_network") %>%
    option_add(param_info = rf_param, id = "ml_random_forest")

  grid_ctrl <-
    control_grid(
      save_pred = TRUE,
      parallel_over = "everything",
      save_workflow = TRUE
    )

  full_results_time <-
    system.time(
      tuned_models <-
        all_workflows %>%
        workflow_map(
          seed = 1003, resamples = folds, grid = n_grid,
          control = grid_ctrl, verbose = TRUE)
    )
  saveRDS(tuned_models, file.path(result_path, glue("{occ_s}_tuned_models.rds")))
  stopImplicitCluster()
  cat("saved tuned models", occ_s, "\n")
}


# Function to run super learner
run_superlearner <- function(occ_s, sd_min){
  library(doParallel)
  all_cores <- parallel::detectCores(logical = FALSE)
  cl <- makePSOCKcluster(all_cores-1)
  registerDoParallel(cl)
  cat("running superlearner", occ_s, "\n")

  tuned_models <- readRDS(file.path(result_path, glue("{occ_s}_tuned_models.rds")))

  model_stack <-
    stacks() %>%
    add_candidates(tuned_models)

  sd_calc <-  model_stack %>%
    dplyr::select(-value) %>%
    pivot_longer(everything(), names_to = "member", values_to = "value") %>%
    group_by(member) %>%
    summarize(sd = sd(value)) %>%
    filter(sd <= sd_min)

  model_stack <- model_stack[!names(model_stack) %in% sd_calc$member]

  ensemble_model <- blend_predictions(model_stack, metric = metric_set(rmse)) %>%
    fit_members
  print(ensemble_model)
  saveRDS(ensemble_model, file.path(result_path, glue("{occ_s}_ensemble_model.rds")))
  stopImplicitCluster()
  cat("saved superlearner results", occ_s, "\n")
}


# RUN MODELS -----------------------------------------------------------------------------

# set results path
result_path <- file.path(proc_path, glue("results/{Sys.Date()}"))
dir.create(result_path, showWarnings = FALSE, recursive = TRUE)

# Tune models for all occ
tic()
walk(unique(db$variable), run_ml, n_grid = 30, k_fold = 10, n_rep = 5, prop_value = 0.80)


# Run superlearner
walk(unique(db$variable), run_superlearner, sd_min = 0.001)
toc()
warnings()

