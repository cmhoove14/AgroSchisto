# Response_Fxs README  
Scripts for fitting dose-response functions relating agrochemical concentrations to components of schistosomiasis transmission.

## Plots  
Plots of all data, fitted dose response functions, and results of sampling the d-r function with uncertainty (separated by study) 
## Data  
Folder containing data extracted from studies identified in review  
## Scripts  

### `bakry2011_malathion_deltamethrin_snails_fit.R`  
Data from [Bakry et al 2011](https://www.sciencedirect.com/science/article/pii/S0048357511001283) investigating effect of **malathion** and **deltamethrin** on *Helisoma duryi* snail mortality and reproduction  
#### Data reported  
+ Table 1: LC50, confidence limits of LC50, slope function, LC0, LC10, LC25, LC90 corresponding to 24-hr mortality and estimated using litchfield and wilcoxon method for both malathion and deltamethrin
+ Figure 2: Longitudinal survival of control snail cohort and snail cohorts exposed to LC10 of malathion and deltamethrin
+ Table 2: Longitudinal egg production of control snail cohort and snail cohorts exposed to LC10 of malathion and deltamethrin
+ Tables 3,4,5: Physiological effects of exposure on enzymatic activity and biochemical measures  
#### Response functions fit  
+ Daily mortality rate for *H. duryi* exposed to malathion, produced from reported parameters derived using Litchfield and Wilcoxon method  
  + Function: `muNq_mal_Bakry11_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Bakry_2011/bakry2011_malathion_snail_mortality.png)  
+ Daily mortality rate for *H. duryi* exposed to deltamethrin, produced from reported parameters derived using Litchfield and Wilcoxon method  
  + Function: `muNq_del_Bakry11_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Bakry_2011/bakry2011_deltamethrin_snail_mortality.png)  
#### Notes
+ Longitudinal survival could be used to investigate effects of cumulative exposure with survival analysis
+ Snail reproduction (measured as eggs/snail at discrete times over a four week period) sharply reduced in snails exposed to LC10 of malathion and deltamethrin, but only able to compare control group to each exposed group, therefore unable to fit dose-response function [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Bakry_2011/bakry2011_snail_reproduction.png)  

### `bakry2012_atrazine_glyphosate_snails_fit.R`  
Data from [Bakry et al 2012](https://www.researchgate.net/publication/256457311_Influence_of_Atrazine_and_Roundup_pesticides_on_biochemical_and_molecular_aspects_of_Biomphalariaalexandrina_snails) investigating effect of **atrazine** and **glyphosate** on *Biomphalaria alexandrina* snail mortality and reproduction  
#### Data reported  
+ Table 3: LC50, confidence limits of LC50, slope function, LC10, LC90 corresponding to 24-hr mortality and estimated using litchfield and wilcoxon method for both atrazine and glyphosate
+ Table 4: Longitudinal egg production, egg abnormalities, and hatching of control snail cohort and snail cohorts exposed to LC10 of atrazine and glyphosate
+ Table 5 & Figs2&5: Physiological effects of exposure on enzymatic activity and biochemical measures  
#### Response functions fit  
+ Daily mortality for *B. alexandrina* exposed to atrazine, produced from reported parameters derived using Litchfield and Wilcoxon method  
  + Function: `muNq_atr_Bakry12_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Bakry_2012/bakry2012_atrazine_snail_mortality.png)  
+ Daily mortality for *B. alexandrina* exposed to glyphosate, produced from reported parameters derived using Litchfield and Wilcoxon method  
  + Function: `muNq_gly_Bakry12_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Bakry_2012/bakry2012_glyphosate_snail_mortality.png)  
#### Notes
+ Longitudinal survival could be used to investigate effects of cumulative exposure with survival analysis
+ Snail reproduction (measured as eggs/snail in addition to hatchlings/snail at discrete times over a four week period) sharply reduced in snails exposed to LC10 of atrazine and glyphosate, but only able to compare control group to each exposed group, therefore unable to fit dose-response function [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Bakry_2012/bakry2012_snail_reproduction.png)  

### `bakry2016_paraquat_snails_fit.R`
Data from [Bakry et al 2016](https://www.ncbi.nlm.nih.gov/pubmed/24081640) investigating effect of **paraquat** on *Lymnaea natalensis* snails
#### Data reported  
+ Table 3: LC25, LC50, LC90 and slope function for 24-hr paraquat exposure, produced from reported parameters derived using Litchfield and Wilcoxon method  
#### Response functions fit  
+ Daily mortality for *L. natalensis* exposed to paraquat  
  + Function fit: `muNq_prq_Bakry16_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Bakry_2016/bakry2016_paraquat_snail_mortality.png) 

### `Baxter_Rohr2011_reanalysis_atrazine_snail_carrying_capacity_fit.R`  
Data from [Baxter et al 2011](http://onlinelibrary.wiley.com/doi/10.1002/etc.552/abstract) and reanalysis of it by [Rohr et al 2012](https://www.researchgate.net/publication/224284702_The_herbicide_atrazine_algae_and_snail_populations?enrichId=rgreq-687b4b78-6cbe-4c48-9c04-d495370159ae&enrichSource=Y292ZXJQYWdlOzIyNDI4NDcwMjtBUzoxMDEyNDM0Njg5MTA1OTlAMTQwMTE0OTczMjY1Mg%3D%3D&el=1_x_3) investigating effect of **atrazine** on algal dynamics and ultimate effects on *Physella* spp and *Stagnicola elodes* snail population dynamics  
#### Data reported  
+ Fig 5: Snail population over time, reanalyzed by Rohr to estimate peak growth rate prior to population collapse due to resource exhaustion  
+ Figs 4&6: Snail direct effects (reproduction, size, and mortality)  
#### Response functions fit  
+ Snail population carrying capacity as a function of atrazine concentration, produced via a log-linear glm relating log+1 atrazine to peak growth rate  
  + Function: `phi_Nq_atr_baxrohr.no30`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Baxter_Rohr_2011/Baxter_rohr_atrazine_carrying_capacity_sim.png)
+ Snail population carrying capacity as a function of any herbicide concentration, assuming response is proportional to its EEC value and is the same as the atrazine response `phi_Nq_rel_baxrohr.no30`  
  + no30 refers to exclusion of data point at 30ppb atrazine due to large influence of spatial block on results. See [Rohr et al 2012](https://www.researchgate.net/publication/224284702_The_herbicide_atrazine_algae_and_snail_populations?enrichId=rgreq-687b4b78-6cbe-4c48-9c04-d495370159ae&enrichSource=Y292ZXJQYWdlOzIyNDI4NDcwMjtBUzoxMDEyNDM0Njg5MTA1OTlAMTQwMTE0OTczMjY1Mg%3D%3D&el=1_x_3) for details  
  
### `Ghaffar2016_butralin_glyphosate_pendimethalin_cercariae_fit.R`  
Data from [Abdel Ghaffar et al 2016](http://www.bioone.org/doi/abs/10.4002/040.059.0201) investigating effect of **butralin**, **glyphosate**, and **pendimethalin** on *Schistosoma mansoni* cercariae  
#### Data reported  
+ Table 6: Longitudinal survival of cercariae exposed to 5 concentrations of each chemical
#### Response functions fit  
+ Cercarial survival as a function of butralin concentration, produced via estimation of LC50 and slope parameters as a function of butralin concentration incorporated into log-logistic d-r curve. See [INSERT LINK TO LARVAL SURVIVAL DYNAMICS EXPLANATION HERE](link.linkylink.com)  
  + Function: `piC.ghaf_butr.exp_unc`
  + Sampling and data: [Raw data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_raw_cercariae_mortality_data_butralin.png), [Parameter models](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_fitted_d-r_parameters_and_functions_cercariae_butralin.png), [Sampling and observed data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_function_simulate_cercariae_butralin.png)  
+ Cercarial survival as a function of glyphosate concentration, produced via estimation of LC50 and slope parameters as a function of butralin concentration incorporated into log-logistic d-r curve. See [INSERT LINK TO LARVAL SURVIVAL DYNAMICS EXPLANATION HERE](link.linkylink.com)  
  + Function: `piC.ghaf_gly.exp_unc`
  + Sampling and data: [Raw data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_raw_cercariae_mortality_data_glyphosate.png), [Parameter models](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_fitted_d-r_parameters_and_functions_cercariae_glyphosate.png), [Sampling and observed data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_function_simulate_cercariae_glyphosate.png)  
+ Cercarial survival as a function of pendimethalin concentration, produced via estimation of LC50 and slope parameters as a function of butralin concentration incorporated into log-logistic d-r curve. See [INSERT LINK TO LARVAL SURVIVAL DYNAMICS EXPLANATION HERE](link.linkylink.com)  
  + Function: `piC.ghaf_pen.exp_unc`
  + Sampling and data: [Raw data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_raw_cercariae_mortality_data_pendimethalin.png), [Parameter models](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_fitted_d-r_parameters_and_functions_cercariae_pendimethalin.png), [Sampling and observed data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Cercariae/Ghaffar_2016_function_simulate_cercariae_pendimethalin.png)  

### `Ghaffar2016_butralin_glyphosate_pendimethalin_miracidia_fit.R`  
Data from [Abdel Ghaffar et al 2016](http://www.bioone.org/doi/abs/10.4002/040.059.0201) investigating effect of **butralin**, **glyphosate**, and **pendimethalin** on *Schistosoma mansoni* miracidia  
#### Data reported  
+ Table 6: Longitudinal survival of miracidia exposed to 5 concentrations of each chemical
#### Response functions fit  
+ miracidial survival as a function of butralin concentration, produced via estimation of LC50 and slope parameters as a function of butralin concentration incorporated into log-logistic d-r curve. See [INSERT LINK TO LARVAL SURVIVAL DYNAMICS EXPLANATION HERE](link.linkylink.com)  
  + Function: `piM.ghaf_butr.exp_unc`
  + Sampling and data: [Raw data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_raw_miracidia_mortality_data_butralin.png), [Parameter models](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_fitted_d-r_parameters_and_functions_miracidia_butralin.png), [Sampling and observed data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_function_simulate_miracidia_butralin.png)  
+ Cercarial survival as a function of glyphosate concentration, produced via estimation of LC50 and slope parameters as a function of butralin concentration incorporated into log-logistic d-r curve. See [INSERT LINK TO LARVAL SURVIVAL DYNAMICS EXPLANATION HERE](link.linkylink.com)  
  + Function: `piC.ghaf_gly.exp_unc`
  + Sampling and data: [Raw data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_raw_miracidia_mortality_data_glyphosate.png), [Parameter models](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_fitted_d-r_parameters_and_functions_miracidia_glyphosate.png), [Sampling and observed data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_function_simulate_miracidia_glyphosate.png)  
+ Cercarial survival as a function of pendimethalin concentration, produced via estimation of LC50 and slope parameters as a function of butralin concentration incorporated into log-logistic d-r curve. See [INSERT LINK TO LARVAL SURVIVAL DYNAMICS EXPLANATION HERE](link.linkylink.com)  
  + Function: `piC.ghaf_pen.exp_unc`
  + Sampling and data: [Raw data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_raw_miracidia_mortality_data_pendimethalin.png), [Parameter models](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_fitted_d-r_parameters_and_functions_miracidia_pendimethalin.png), [Sampling and observed data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Miracidia/Ghaffar_2016_function_simulate_miracidia_pendimethalin.png)  