# occupations_on_the_map

This repository contains the scripts to reproduce the analysis presented in Van Dijk et al. (2022), Occupations on the map: Using a super learner algorithm to downscale labor statistics:

__set_path.r:__ Contains the path where the input is saved. At the moment this is `c:/temp/occupations_on_the_map` but can be replaced by any folder name. This script is sourced by all the other scripts so input data can be found and further processed.

__run_ml_workflow.r:__ Loads the data and runs six different machine learning models and the superlearner for each of the six occupations for which maps are created as well as the labour force participation rate.

__prepare_input_data.r:__ Processes the input data (occupation shares and labour force participation rate at the subnational level and database with geospatial predictors). It is sourced by `run_ml_workflow.r`.

__post_process_model_results.r:__ Post-processes the machine learning output in order to create maps and figures and conduct additional analysis as described in the paper.

__fig_table_map_manuscript.r:__ Creates all the figures in the main manuscript.

__fig_table_map_si.r:__ Creates all the figures in the Supplementary Information.



