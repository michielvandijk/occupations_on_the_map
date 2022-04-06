# ========================================================================================
# Project:  KB project, people on the map
# Subject:  Prepare input data for machine learning model
# Author:   Michiel van Dijk; Thijs de Lange
# Contact:  michiel.vandijk@wur.nl; thijs.delange@wur.nl
# ========================================================================================

# ========================================================================================
# SETUP ----------------------------------------------------------------------------------
# ========================================================================================

# Load pacman for p_load
if (!require(pacman)) install.packages("pacman")
library(pacman)

# Load key packages
p_load(here, tidyverse, readxl, stringr, scales, glue)

# Load additional packages
p_load(naniar, DescTools, DataExplorer)

# Set path
source(here("scripts/set_path.r"))

# R options
options(scipen = 999)
options(digits = 4)


# ========================================================================================
# LOAD DATA ------------------------------------------------------------------------------
# ========================================================================================

lfpr_raw <- readRDS(file.path(proc_path, "model_input/VNM_lfpr.rds")) %>%
  ungroup
occ_raw <- readRDS(file.path(proc_path, "model_input/VNM_occ.rds")) %>%
  ungroup
predictors_raw <- readRDS(file.path(proc_path, "model_input/VNM_predictors.rds")) %>%
  ungroup


# ========================================================================================
# PREPARE ADM-LEVEL DATA -----------------------------------------------------------------
# ========================================================================================

# Create database with predictors at adm-level to estimate model
db_adm <- predictors_raw %>%
  group_by(adm_code) %>%
  summarize(across(c(airports:viirs),
                   mean,
                   na.rm = TRUE),
            .groups = "drop"
            )
summary(db_adm)
vis_miss(db_adm)

# Link lfpr and occ with pred and combine
db_adm <- bind_rows(
  lfpr_raw %>%
    dplyr::select(-working_population) %>%
    dplyr::rename(value = lfpr) %>%
    full_join(db_adm) %>%
    mutate(variable = "lfpr"),
  occ_raw %>%
    full_join(db_adm) %>%
    dplyr::select(-value) %>%
    dplyr::rename(value = share)
)


# ========================================================================================
# CLEAN UP -------------------------------------------------------------------------------
# ========================================================================================

rm(lfpr_raw, predictors_raw)

