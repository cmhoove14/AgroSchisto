# Source functions and models for agrochemical review models and response functions  

## `litchfield_wilcoxon_get_b1_from_slope`  
Brainstorming/process used to arrive at method for estimating the slope of an underlying linear model fit to probit transformed outcome and log10 transformed dose data (i.e. via the [litchfield and wilcoxon method](http://jpet.aspetjournals.org/content/96/2/99)) when only the slp and lc50 parameters are reported  
+ .Rmd file contains notes and entire thought process used to produce the .pdf file
+ .R file contains the function loaded in response function fitting scripts  

## `d-r_function_compare_methods.R`  
Exploring different ways of reproducing dose response function from reported slope and lc50 parameters reported from studies using litchfield and wilcoxon 1949 method for estimating dose response relationships. Mainly motivated by an inability to reproduce a d-r function that makes sense for [Abdel Ghaffar et al 2016](http://www.bioone.org/doi/abs/10.4002/040.059.0201). This script demonstrates that fitting a linear function to probit transformed mortality and log10 LC_ values produces a d-r function similar to the litchfield and wilcoxon d-r function at lower concentrations where we find the EEC (i.e. the range of concentrations we're actually ultimately modeling). Furthermore, the variance of these functions is higher than that of the reproduced LW1949 function, implying that this methods also has the desirable trait of being more conservative given it is less certain that the directly reproduced function.

## `r0_functions.R`  
R0 functions for agrochemical model, base model (no agrochemicals, no predators) and predator only model

## `r0_explore_parameter_influences.R`  
Exploration of influence of different agrochemical-sensitive parameters on R0  

## `models.R`  
Functions to simulate full model under difference agrochemical exposure scenarios  

## `model_helper_functions.R`
Small functions used in the R0 and dynamic model code to make things a bit more readable

## `mod_q_test.R`  
Exploration of full model    