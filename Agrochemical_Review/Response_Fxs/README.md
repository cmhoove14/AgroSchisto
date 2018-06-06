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

### `barbieri2016_carbofuran_predators_fit.R`  
Data from [Barbieri et al 2016](https://www.ncbi.nlm.nih.gov/pubmed/23847016) investigating effect of **carbofuran** on *Macrobrachium olfersii* prawns  
#### Data reported  
+Table 1: 24-hr mortality rates of cohorts of 135 juvenile prawns exposed to 8 concentrations of pesticide  
#### Response functions fit  
+ Daily mortality of *Macrobrachium olfersii* using both DRC log-logistic fit and LW1949 methods   
  + Function fit: `barb_carbofuran_muPq_uncertainty` (LW1949 function - used for comparability to other functions)  
  + Function fit: `muPq_carb_barb16_uncertainty` (DRC function)  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Barbieri2016/carbofuran_mupq_sim_observed.png)  

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
  
### `Benli2007_2-4D_predators_plot.R`  
Data from [Benli et al 2007](https://www.sciencedirect.com/science/article/pii/S0048357507000053#bib17) investigating effects of the herbicide **2,4D** on survival of adult *Astacus leptodactylus* crayfish  
#### Data reported  
+ Table 1 and text: lc1, lc5, lc10, lc15, lc50, lc85, lc90, lc95, lc99, lc50CIs, slope
#### Response functions fit  
+ `muPq_24d_benli07_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Benli2007/benli2007_2-4D_predator_mortality.png)
  
### `Bhavan2010_carbaryl_predator_feed_rate_fit.R`  
Data from [Bhavan et al 2010](https://www.ncbi.nlm.nih.gov/pubmed/22167951) investigating effect of **Carbaryl** on *Macrobrachium malcomsonii* food consumption  
#### Data reported  
+ Table 1 - food consumed in grams per day across three concentrations and a control group  
#### Response functions fit  
+ Predator feeding rate as a function of carbaryl concentration  
  + Function: `psi_q_carb_bhavan10_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Bhavan2010/bhavan2010_carbaryl_feed_rate_data_sim.png)  
  
### `Browne&Moore2014_2-4D_pred_consumption_fit.R`   
Data from [Browne and Moore 2014](https://www.ncbi.nlm.nih.gov/pubmed/24799048) investigating effects of **2,4-D** on feeding behavior of *Orconectes rusticus* crayfish  
#### Data reported  
+ Figure 6: consumption rate of crayfish measured as weight of gelatin fish food eaten across dose groups  
#### Response functions fit  
+ `psi_q_24D_browne14_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Browne2014/pred_consumption_2-4D_data_sim.png)
  
###`Fornstrom1997_terbufos_predators_fit.R`  
Data from [Fornstrom et al 1997](https://setac.onlinelibrary.wiley.com/doi/abs/10.1002/etc.5620161212) investigating effects of **Terbufos** on survival of *Procambarus clarkii*  
#### Data reported  
+ Figure 1: Longitudinal mortality rate of test organisms in each concentration group; data extracted at 24-hr mark for 24-hr D-R analysis  
#### Response functions fit  
+ `muPq_terb_Fornstrom_uncertainty`  
  + Sampling and Data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Fornstrom1997/Fornstrom1997_pred_mortality_data&sim.png)
  
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
  
### `Ghaffar2016_butralin_glyphosate_pendimethalin_snails_fit.R`  
Data from [Abdel Ghaffar et al 2016](http://www.bioone.org/doi/abs/10.4002/040.059.0201) investigating effect of **butralin**, **glyphosate**, and **pendimethalin** on *Biomphalaria alexandrina* snails  
#### Data reported  
+ Table 1: 24-hr LC0, LC10, LC25, LC50, LC50 95%CI, LC90, LC90 95%CI, and "Slope of refression mortality on conc."  
+ This data is somewhat questionable and reproducing dose-response relationships from it relies on estimating d-r based on the reported LC values rather than LC and slope as for other studies (see `Models/d-r_function_compare_methods.R`. Cited method is Litchfield and Wilcoxon, but methods to reproduce d-r function from reported Litchfield and Wilcoxon parameters in other studies from the same authors are not working here. 
#### Response functions fit  
+ Butralin  
  + `mu_Nq_butr_gaf16_uncertainty` snail mortality  
    + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Snails/butralin_muN_simulate)  
  + `fN.butr.fx.uncertainty` snail reproduction  
    + Sampling and data [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Snails/butralin_fN_simulate.png)
+ Glyphosate    
  + `mu_Nq_gly_gaf16_uncertainty` snail mortality  
    + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Snails/glyphosate_muN_simulate)  
  + `fN.gly.fx.uncertainty` snail reproduction  
    + Sampling and data [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Snails/glyphosate_fN_simulate.png)
+Pendimethalin      
  + `mu_Nq_pen_gaf16_uncertainty` snail mortality  
    + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Snails/pendimethalin_muN_simulate)  
  + `fN.pen.fx.uncertainty` snail reproduction  
    + Sampling and data [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ghaffar_2016/Snails/pendimethalin_fN_simulate.png)

### `griggs08_atrazine_metolachlor_cercariae_fit.R`  
Data from [Griggs et al 2008](https://www.ncbi.nlm.nih.gov/pubmed/17763881) investigating effects of **atrazine** and **metloachlor** on *Echinistoma trivolvis* cercariae  
#### Data reported  
+ Fig 1: cercarial die off over time in water control, solvent control, low dose group, and high dose group. Low dose group contains 10ppb metolachlor **and** 15ppb atrazine; high dose group contains 85ppb metolachlor **and** 100ppb atrazine therefore other studies that estimate response to these chemicals alone and in schistosome (rather than *Echistoma*) will be prioritized  
#### Response functions fit  
+ cercarial survival as a function of metolachlor/atrazine concentration, produced via estimation of LC50 and slope parameters as a function of metolachlor/atrazine concentration incorporated into log-logistic d-r curve. See [INSERT LINK TO LARVAL SURVIVAL DYNAMICS EXPLANATION HERE](link.linkylink.com)  
  + Function: `piC.grg08_atr_unc2`
  + Sampling and data: [Raw data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Griggs2008/Griggs2008_raw_cercariae_mortality_data_atrazine.png), [Parameter models](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Griggs2008/Griggs2008_fitted_d-r_parameters_and_functions_cercariae_atrazine.png), [Sampling and observed data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Griggs2008/Griggs2008_function_simulate_cercariae_atrazine.png)  

### `Gustafson2016_atrazine_theta_fit.R`  
Data from [Gustafson et al 2016](https://www.ncbi.nlm.nih.gov/pubmed/26762862) investigating effects of **atrazine** on parasite host nteractions in an amphibian trematode system. data of use is cercarie produced per day by *Physa acuta* snails infected with *Halipegus eccentricus*  
#### Data reported  
+ Table 1/ Figure 2B: cercarie production per snail in control and two atrazine dose groups  
#### Response functions fit  
+ Relative cercariae produced per snail
  + Response function: `gust_theta_uncertainty`
  + sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Gustafson2016/atrazine_theta_data_sim.png)

### `Halstead2015_insecticides_predators_fit.R`  
Data from [Halstead et al 2015](https://www.sciencedirect.com/science/article/pii/S0045653515003410?via%3Dihub) investigating effects of **malathion**, **chlorpyrifos**, **terbufos**, **esfenvalerate**, **lambda-cyhalothrin**, and **permethrin** on *Procambarus clarkii* (crayfish) adults  
#### Data reported  
+ Received the raw data from Neal Halstead (lead author), but could also extract the data from Fig 1  
  + Data for 4-day and 10-day mortality endpoints  

#### Response functions fit  
+ For each insecticide, 4-day mortality as a function of insecticide concentration  
  + malathion  
    + function: `muPq_mal_Halstead_uncertainty`  
    + sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Halstead2015/halstead2015_malathion_predator_mortality.png)  
  + chlorpyrifos  
    + function: `muPq_chlor_Halstead_uncertainty`  
    + sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Halstead2015/halstead2015_chlorpyrifos_predator_mortality.png)  
  + terbufos  
    + function: `muPq_terb_Halstead_uncertainty`  
    + sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Halstead2015/halstead2015_terbufos_predator_mortality.png)  
  + esfenvalerate  
    + function: `muPq_esfen_Halstead_uncertainty`  
    + sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Halstead2015/halstead2015_esfenvalerate_predator_mortality.png)  
  + lambda-cyhalothrin  
    + function: `muPq_lamcy_Halstead_uncertainty`  
    + sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Halstead2015/halstead2015_lamcy_predator_mortality.png)  
  + permethrin  
    + function: `muPq_perm_Halstead_uncertainty`  
    + sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Halstead2015/halstead2015_permethrin_predator_mortality.png)  
    
### `Halstead_meso2017_fit.R`   
Data from [Halstead et al 2018](https://www.nature.com/articles/s41467-018-03189-w) investigating effects of **Atrazine**, **Chlorpyrifos**, and **Fertilizer** on snail populations and infection rates in a mesocosm setting  
#### Data reported  
+ Received raw data on snail and predator populations in different mesocosm treatments from Neal Halstead (lead author)

#### Response functions fit  
+ `halstead17_phiN_fe_uncertainty`
+ `halstead17_phiN_at_uncertainty`

### `Hasheesh2011_chlorpyrifos_profenofos_larvae_fit.R`  
Data from [Hasheesh & Mohamed 2011](https://www.sciencedirect.com/science/article/pii/S0048357511000186) investigating effects of **Chlorpyrifos** and **Profenofos** on survival of *S. haematobium* miracidia and cercarie and *Bulinus truncatus* snails (see below)  
#### Data reported  
+ Table 5: 8 hr LC50. LC50 95%CI, LC90, and slp parameters for miracidia and cercariae from Litchfield and Wilcoson method  
#### Response functions fit  
+ Cercariae
  + Chlorpyrifos 
    + Function: `piC_ch_Hash11_uncertainty`
    + sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Hasheesh2011/Larvae/hasheesh2011_chlorpyrifos_cercariae_mortality.png)
  + Profenofos: `piC_prof_Hash11_uncertainty`
    + sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Hasheesh2011/Larvae/hasheesh2011_profenofos_cercariae_mortality.png)
+ Miracidia
  + Chlorpyrifos 
    + Function: `piM_ch_Hash11_uncertainty`
    + sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Hasheesh2011/Larvae/hasheesh2011_chlorpyrifos_miracidia_mortality.png)
  + Profenofos: `piM_prof_Hash11_uncertainty`
    + sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Hasheesh2011/Larvae/hasheesh2011_profenofos_miracidia_mortality.png)  
    
### `Hasheesh2011_chlorpyrifos_profenofos_snails_fit.R`    
Data from [Hasheesh & Mohamed 2011](https://www.sciencedirect.com/science/article/pii/S0048357511000186) investigating effects of **Chlorpyrifos** and **Profenofos** on survival of *Bulinus truncatus* snails  
#### Data reported  
+ Table 1: 24 hr LC25, LC50. LC50 95%CI, LC90, and slp parameters for *Bulinus truncatus* adult snails from Litchfield and Wilcoson method  
+ Table 2: Longitudinal survival and reproduction of snails exposed to LC25 of each chemical **NOTE** Only two data points (control and single dose group), so not possible to fit a dose response, but will consider how effects at this single concentration affects transmission  
+ Table 3: Growth rate of snails exposed LC25 of each chemical **NOTE** No growth in our model, so don't fit a function to this data  
+ Table 4: Effect of LC25 of Chlorpyrifos and Profenophos pesticides on the infectivity of Schistosoma haematobium miracidia to Bulinus truncatus snails. **NOTE** We assume that reduction in infectivity captured here are mostly due to increased mortality of miracidia therefore don't attempt to separately fit a response to this data 
#### Response functions fit  
+ Chlorpyrifos  
  + Function: `muNq_ch_hash11_uncertainty`  
    + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Hasheesh2011/Snails/hasheesh2011_chlorpyrifos_snail_mortality.png)  
  + Function: `fN.hash.chlor.uncertainty` (samples measurement of eggs/snail/day in the chlorpyrifos group and in the control group and returns the ratio which is interpreted as relative change in snail reproduction when exposed to LC25 of chlorpyrifos)  
    + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Hasheesh2011/Snails/hasheesh2011_snail_reproduction_compare.png)  
+ Profenofos  
  + Function: `muNq_prof_hash11_uncertainty`  
    + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Hasheesh2011/Snails/hasheesh2011_profenofos_snail_mortality.png)  
  + Function: `fN.hash.prof.uncertainty` (samples measurement of eggs/snail/day in the profenofos group and in the control group and returns the ratio which is interpreted as relative change in snail reproduction when exposed to LC25 of profenofos)  
    + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Hasheesh2011/Snails/hasheesh2011_snail_reproduction_compare.png)     

### `hua_insecticides_cercariae_fit.R`  
Data from [Hua et al 2016](doi:10.1017/S0031182015001894) investigating effects of **carbaryl**, **malathion**, **cypermethrin**, **permethrin**, **imidacloprid**, and **thiamethoxam** on *Echinoparyphium spp*  
#### Data reported  
+ Figure 1: mortality-dose plots for all six insecticides  
#### Response functions fit  
+ NONE RELIABLE: They tested insecticide toxicity to two populations: "PWA" is from a pond far from any agricultural influence and "ICP" is from a pond next to ag sources. We'll use PWA since we don't want to consider the role of adaptation in any quantitative analysis here. Also have a feeling that within the small concentration they investigated, they may be picking up more natural variability and/or error in their actual concentrations than any sort of direct effect of the insecticides themselves. In addition, they should really test a range of concentrations that spans from no effect to complete effect (e.g. 100% mortality within 6 hours) in order to accurately estimate a dose-response relationship. Attempted to proceed regardless of these concerns, but was unable to estimate reliable d-r functions for most insecticides, so excluding this study. See R script for details

### `Hussein2016_fertilizer_snails_fit.R`  
Data from [Hussein et al 2016](https://www.researchgate.net/publication/305149067_Effects_of_Three_Inorganic_Fertilizers_on_the_Biology_and_Histopathology_of_infected_Biomphalaria_alexandrina_snails?enrichId=rgreq-cf1a38509b7c5460cbc30f698ea58594-XXX&enrichSource=Y292ZXJQYWdlOzMwNTE0OTA2NztBUzozODI2MjU5NTgxMjE0NzJAMTQ2ODIzNjU0NTE4NQ%3D%3D&el=1_x_3&_esc=publicationCoverPdf) investigating effects of three fertilizers on *Biomphalaria alexandrina* survival  
#### Data reported  
+ Table 1: LC50, LC50 CIs, slope function, LC90 from LW1949 method  
#### Response functions fit  
+ balanced fertilizer  
  + Function: `muNq_balanced_hussein16_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Hussein2016/balanced_fertilizer_snail_mortality.png)  
+ High phosphorous fertilizer  
  + Function: `muNq_highp_hussein16_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Hussein2016/highp_fertilizer_snail_mortality.png)  
+ High nitrogen fertilizer     
  + Function: `muNq_highn_hussein16_uncertainty`   
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Hussein2016/highn_fertilizer_snail_mortality.png)   
    
### `Ibrahim1992_chlorpyrifos_snails_fit.R`    
Data from [Ibrahim et al 1992](https://www.ncbi.nlm.nih.gov/pubmed/1379273) investigating effects of **chlorpyrifos** on survival and reproduction of *Biomphalaria alexandrina*  
#### Data reported   
+ Lots of different series tested based on when chlorpyrifos was added/removed from tanks containing snails, series 1 with data contained in Table 1 is most applicable  
#### Response functions fit  
+ Function: `mu_N_chlor_ibr92_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ibrahim1992/ibrahim1992_chlorpyrifos_snail_mortality.png)
+ Function: `f_N_chlor_ibr92_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ibrahim1992/ibrahim1992_chlorpyrifos_snail_reproduction.png)  
  
### `Johnson2007_fertilizer_fit.R`  
Data from [Johnson et al 2007](https://doi.org/10.1073/pnas.0707763104) investigating effects of eutrophication from **fertilizer** on snail populations and cercarial output and infection rates in tadpoles for the *Ribeiroa ondatrae* system. No dose-response functions here because only one treatment group compared to a control group.
#### Data reported  
+ Cercariae shed per infected snail per day in control and fertilizer groups (Fig 3c)  
+ Snail eggs collected in samplers at 5 sampling points
+ Total snail biomass over the course of the experiment
#### Response functions fit  
+ Function: `johnson07_theta_uncertainty` (proportional increase in cercarial shedding in fertilizer group)
+ Function: `johnson07_fN_uncertainty` (proportional increase in snail eggs sampled in fertilizer group)
+ Function: `johnson07_phin_uncertainty` (proportional increase in snail biomass in fertilizer group)

### `koprivnikar2006_atrazine_cercariae_fit.R`  
Data from [Koprivnikar et al 2006](https://www.ncbi.nlm.nih.gov/pubmed/16729687) investigating effects of atrazine on survival of *Echinistoma trivolvis* cercariae  
#### Data reported  
+ Cercarial survival and standard error across three concentration levels  
#### Response functions fit  
+ Function: `piC_kop_atr_unc`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Koprivnikar2006/Koprivnikar_2006_function_simulate_piC_atrazine.png)

### `Kristoff2011_azinphos-methyl_snail_reproduction_fit.R`  
Data from [Kristoff et al 2011](https://www.ncbi.nlm.nih.gov/pubmed/21524429) investigating reprodcutive toxicity of **azinphos-methyl**, an organophosphate insecticide, to *Biomphalaria glabrata* snails  
#### Data reported  
+ Table 1: reproductive parameters of snail cohorts exposed to azinphos-methyl
#### Response functions fit  
+ Function: `fNq_azmeth_kristoff11_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Kristoff2011/kristoff2011_azinphos-methyl_snail_reproduction_data_sim.png)

### `Leung1980_fenitrothion_predators_fit.R`  
Data from [Leung et al 1980](https://www.ncbi.nlm.nih.gov/pubmed/7426798) investigating effects of **fenitrothion** on mortality of adult and juvenile *Procambarus clarkii*  
#### Data reported  
+ Tables 2 and 3: cohort size and mortality rate (24-hr mortality endpoint used)
#### Response functions fit  
+ `muPq_fenitrothion_Leung_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Leung1980/Leung1980_pred_mortality_data&sim.png)
+ `muPq_juv_fenitrothion_Leung_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Leung1980/Leung1980_juv_pred_mortality_data&sim.png)

### `Mohamed2012_profenofos_diazinon_snails_fit.R`  
Data from [Mohamed et al 2012](https://www.sciencedirect.com/science/article/pii/S0048357511001696) investigating effect of **Diazinon** and **profenofos** on mortality and reproduction of *Biomphalaria alexandrina*  
#### Data reported  
+ Table 1: LC0, LC10, LC25, LC50, LC90, and slp parameters from Litchfield and Wilcoxon d-r model. No reported of uncertainty, so borrowed a value of 0.12 which represents the upper range of uncertainty from similar studies  
#### Response functions fit  
+ Diazinon  
  + Function: `muNq_diaz_mohamed_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Mohamed2012/diazinon_snails_data_sim.png)
  + Function: `fNq_moh_diaz_moh12_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Mohamed2012/diazinon_snail_reproduction_data_sim.png)
+ profenofos  
  + Function: `muNq_prof_mohamed_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Mohamed2012/profenofos_snails_data_sim.png)
  + Function: `fNq_moh_prof_moh12_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Mohamed2012/profenofos_snail_reproduction_data_sim.png)

### `Monde2016_endosulfan_snail_fit.R`  
Data from [Monde et al 2016](https://link.springer.com/article/10.1007/s00244-016-0275-7) investigating effects of **endosulfan** on *Bulinus globusus* survival and predation by catfish  
#### Data reported  
+ Table 2: lc10, lc50, and lc90 parameters with confidence intervals  
#### Response functions fit  
+ `muNq_endo_monde16_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Monde2016/endosulfan_snails_sim.png)

### `Monte2016_glyphosate_larvae_fit.R`  
Data from [Monte et al 2016](https://www.ncbi.nlm.nih.gov/pubmed/27373431) investigating effects of roundup (active ingredient: 41% **glyphosate**) on survival of *Echinistoma paraensei* cercariae and miracidia
#### Data reported  
+ Table 2: longitudinal survival of miracidia in control and 5 dose groups  
+ Table 3: longitudinal survival of cercariae in control and 9 dose groups  
#### Response functions fit  
+ `piC.mont16_gly_unc`
  + Raw time series data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Monte2016/monte2016_glyphosate_cercariae_data.png)
  + Time d-r parameters over concentration: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Monte2016/monte2016_glyphosate_cercariae_dose_pars.png)
  + Time series validation: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Monte2016/monte2016_glyphosate_cercariae_validation.png)
  + Sampling and data of d-r function: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Monte2016/monte2016_glyphosate_cercariae_piC_validation.png)
+ `piM.mont16_gly_unc`
  + Raw time series data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Monte2016/monte2016_glyphosate_miracidia_data.png)
  + Time d-r parameters over concentration: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Monte2016/monte2016_glyphosate_miracidia_dose_pars.png)
  + Time series validation: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Monte2016/monte2016_glyphosate_miracidia_validation.png)
  + Sampling and data of d-r function: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Monte2016/monte2016_glyphosate_miracidia_piM_validation.png)

### `Naqvi1983_trifluralin_oryzalin_predators_fit.R`  
Data from [Naqvi et al 1983](https://www.ncbi.nlm.nih.gov/pubmed/6688742) investigating effects of  **trifluralin** and **oryzalin** on survival of *Procambarus clarkii*
#### Data reported  
+ Tables 1 & 3: surivial rate and cohort size over tested concentrations of each chemical
#### Response functions fit  
+ `muPq_trifluralin_Naqvi83_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Naqvi/Naqvi1983_trifluralin_p_clarkii.png)
+ `muPq_oryzalin_Naqvi83_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Naqvi/Naqvi1983_oryzalin_p_clarkii.png)

### `Naqvi1987_endosulfan_etc_predators_fit.R`  
Data from [Naqvi et al 1987](https://www.sciencedirect.com/science/article/pii/0269749187901084) investigating effects of Thiodan (active ingredient **endosulfan**), Treflan (active ingredient **trifluralin**), **MSMA**, and Oust (active ingredient **sulfometuron-methyl**) on survival of *Procambarus clarkii*. Other ingredient Cutrine-Plus not considered because it's an algicide
#### Data reported  
+ Table 2: LC50, LC50 95%CIs, and slope parameters derived from Probit regression
#### Response functions fit  
+ `muPq_endo_Naqvi87_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Naqvi/Naqvi1987_endosulfan_adult_p_clarkii.png)
+ `muPq_trif_Naqvi87_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Naqvi/Naqvi1987_treflan_adult_p_clarkii.png)
+ `muPq_msma_Naqvi87_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Naqvi/Naqvi1987_msma_adult_p_clarkii.png)
+ `muPq_oust_Naqvi87_uncertainty`  
  + Not found to be toxic to adults therefore function just returns 0

### `Oliveira_endosulfan_snail_fit.R`  
Data from [Oliveira et al 2009](https://www.ncbi.nlm.nih.gov/pubmed/19138786) investigating effects of **endosulfan** on reproduction of *Biomphlaria tenagophila* snails  
#### Data reported  
+ Figs 1&2: eggs per snail over eights weeks, mortality of eggs and hathability of eggs. Together provide good estimate of hatchlings/snail  
#### Response functions fit  
+ `fNq_endo_oliv09_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Oliveira2009/endosulfan_snail_reproduction_data_sim.png)

### `Omkar1985_endosulfan_phosphamidon_carbaryl_predators_fit.R`  
Data from [Omkar and Multi 1985](https://www.jstor.org/stable/20104053?seq=1#page_scan_tab_contents) investigating effects of **endosulfan**, **phosphamidon**, and **carbaryl** on adult (50-55mm) Macrobrachium dayanum  
#### Data reported  
+ Table 1: 24hr lc25, lc50, lc75 and CI of lc50  
#### Response functions fit  
+ `muPq_endo_omkar85_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Omkar1985/omkar1985_predator_mortality_endosulfan.png)
+ `muPq_phos_omkar85_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Omkar1985/omkar1985_predator_mortality_phosphamidon.png)
+ `muPq_carb_omkar85_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Omkar1985/omkar1985_predator_mortality_carbaryl.png)

### `Omran&Salama2013_atrazine_glyphosate_snails_fit.R`  
Data from [Omran and Salama 2013](https://www.ncbi.nlm.nih.gov/pubmed/24215068) investigating effects of atrazine on survival of *Biomphalaria alexandrina* snails  
#### Data reported  
+ Figure 1 log concentration vs probit mortality for 7 dose groups    
#### Response functions fit  
+ Function: `ons.munq.atr`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Omran&Salama2013/Omran&Salama2013_function_simulate_muN_atrazine.png)
+ Function: `ons.munq.gly`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Omran&Salama2013/Omran&Salama2013_function_simulate_muN_glyphosate.png)

### `ragab2006_fertilizers_snails_fit`  
Data from [Ragab and Shoukry 2006](https://www.ncbi.nlm.nih.gov/pubmed/17153706) investigating effects of fertilizers on *Biomphalaria alexandrina* snails  
#### Data reported  
+ Table 1 LC50, 95%CI of LC50, LC90, and slp parameters  
#### Response functions fit  
+ Function: `rag06_mun_amm`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ragab&Shoukry2006/Ragab&Shoukry2006_function_simulate_muN_ammonium_nitrate.png)
+ Function: `rag06_mun_pot`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ragab&Shoukry2006/Ragab&Shoukry2006_function_simulate_muN_potassium_sulphate.png)
+ Function: `rag06_mun_urea`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Ragab&Shoukry2006/Ragab&Shoukry2006_function_simulate_muN_urea.png)

### `rohr08_nature_atrazine_fit.R`   
Data from [Rohr et al 2008](https://www.ncbi.nlm.nih.gov/pubmed/18972018) investigating effects of **Atrazine** on snail populations and infection rates in a an amphibian trematode system    
#### Data reported  
+ Table 1 data on snail population dynamics in mesocosm setting comparing control to atrazine treatment group  
#### Response functions fit  
+ `rohr08_fN_uncertainty` - function using snail eggs data
+ `rohr08_fN_uncertainty2` - function using snail hatchlings data (preferred due to parameter interpretation)

### `rohr08_atrazine_cercariae_fit.R`   
Data from [Rohr et al 2008](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/07-1429.1) investigating effects of **Atrazine** on survival of  *Echinistoma trivolvis* cercariae    
#### Data reported  
+ Figure 3 %mortality in hours 14-18 of cercariae in 6 different dose groups  
#### Response functions fit  
+ Function: `piC.atr.rohr08.lin`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Rohr08/Rohr08_piC_atrazine_function_test.png)

### `Satapornvanit2009_insecticides_predators_fit.R`  
Data from [Satapornvanit et al 2009](https://www.ncbi.nlm.nih.gov/pubmed/19103457) investigating effects of **Zinc**, **Chlorpyrifos**, **Dimethoate**, and **profenofos** on mortality and feeding rate of post-larvae *Macrobrachium rosenbergii*  
#### Data reported  
+ Figure 1 contains data on mean mortality rate for each chemical and dose group as well as reduction in feeding rate  
#### Response functions fit  
*Functions were fit using both `drc` and `LW1949` with results basically the same*
+ Function: `muPq_chlor_satapornvanit09_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Satapornvanit2009/Satapornvanit2009_data_functions_chlorpyrifos_mupq.png)
+ Function: `muPq_dim_satapornvanit09_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Satapornvanit2009/Satapornvanit2009_data_functions_dimethoate_mupq.png)
+ Function: `muPq_prof_satapornvanit09_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Satapornvanit2009/Satapornvanit2009_data_functions_profenofos_mupq.png)
+ Function: `muPq_zinc_satapornvanit09_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Satapornvanit2009/Satapornvanit2009_data_functions_zinc_mupq.png)
+ Function: `psi_q_chlor_satapornvanit09_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Satapornvanit2009/Satapornvanit2009_data_functions_chlorpyrifos_psiq.png)
+ Function: `psi_q_zinc_satapornvanit09_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Satapornvanit2009/Satapornvanit2009_data_functions_zinc_psiq.png)

### `Shukla1985_quinalphos_dichlorvos_monocrotophos_carbaryl_predators_plot.R`  
Data from [Omran and Shukla 1985](https://www.jstor.org/stable/20104000) investigating effects of **quinalphos**, **dichlorvos**, **monocrotophos**, and **carbaryl** on *Macrobrachium lamarrei*  
#### Data reported  
+ Table 1: lc25, lc50, lc75 and CIs for lc50
#### Response functions fit  
+ `muPq_quin_shukla85_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/shukla1985/shukla1985_predator_mortality_quinalphos.png)
+ `muPq_dich_shukla85_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/shukla1985/shukla1985_predator_mortality_dichlorvos.png)
+ `muPq_mono_shukla85_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/shukla1985/shukla1985_predator_mortality_monocrotophos.png)
+ `muPq_carb_shukla85_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/shukla1985/shukla1985_predator_mortality_carbaryl.png)
  
### `Shukla1984_endosulfan_methyl-demeton_carbaryl_predators_fit.R`  
Data from [Omran and Shukla 1984](http://www.jstor.org/stable/20103942) investigating effects of **endosulfan**, **methyldemeton**, and **carbaryl** on *Macrobrachium lamarrei*  
#### Data reported  
+ Table 1: lc25, lc50, lc75 and CIs for lc50
#### Response functions fit  
+ `muPq_endo_shuk84_uncertainty`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/shukla1984/shukla1984_predator_mortality_endosulfan.png)
+ `muPq_mede_shuk84_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/shukla1984/shukla1984_predator_mortality_methyl-demeton.png)
+ `muPq_carb_shuk84_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/shukla1984/shukla1984_predator_mortality_carbaryl.png)

### `tantawy2002_butachlor_fpb_snails_fit.R`   
Data from [Tantawy 2002](https://www.ncbi.nlm.nih.gov/pubmed/12512816) investigating effects of **Butachlor**, **fluazifop-p-butyl** on mortality of *Biomphalaria alexandrina* snails and *Schistosoma mansoni* miracidia and cercariae  
#### Data reported  
+ **Snails -** Table 1: 24-hr LC50, 95%CI LC50, LC0, LC10, LC25, LC90, slp for *Biomphalaria alexandrina*  
#### Response functions fit  
+ `muN.tant.but_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Snails/tantawy2002_butachlor_muN_data_sim.png)
+ `muN.tant.fpb_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Snails/tantawy2002_fpb_muN_data_sim.png)

### `tantawy2002_butachlor_fpb_miracidia_fit.R`   
Data from [Tantawy 2002](https://www.ncbi.nlm.nih.gov/pubmed/12512816) investigating effects of **Butachlor**, **fluazifop-p-butyl** on mortality of *Biomphalaria alexandrina* snails and *Schistosoma mansoni* miracidia and cercariae  
#### Data reported  
+ **Miracidia -** Table 4&5: longitudinal survival of *Schistosoma mansoni* miracidia in control and 4 dose groups 
#### Response functions fit  
+ Miracidial survival as a function of butralin concentration, produced via estimation of LC50 and slope parameters as a function of butralin concentration incorporated into log-logistic d-r curve. See [INSERT LINK TO LARVAL SURVIVAL DYNAMICS EXPLANATION HERE](link.linkylink.com)  
+ Butachlor linear function: `piM.tant02_but.lin_unc` - linear relationship between time-LC50 and concentration  
+ Butachlor exponential function: `piM.tant02_but.exp_unc` - exponential relationship between time-LC50 and concentration  
  + Sampling and data for both: [Raw data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Miracidia/tantawy2002_butachlor_piM_data.png), [Parameter models](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Miracidia/tantawy2002_butachlor_piM_parameter_models.png), [Sampling and observed data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Miracidia/tantawy2002_butachlor_piM_data_sim.png)  
+ Fluazifop-p-butyl linear function: `piM.tant02_fpb.lin_unc` - linear relationship between time-LC50 and concentration  
+ Fluazifop-p-butyl exponential function: `piM.tant02_fpb.exp_unc` - exponential relationship between time-LC50 and concentration  
  + Sampling and data for both: [Raw data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Miracidia/tantawy2002_fpb_piM_data.png), [Parameter models](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Miracidia/tantawy2002_fpb_piM_parameter_models.png), [Sampling and observed data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Miracidia/tantawy2002_fpb_piM_data_sim.png)  

### `tantawy2002_butachlor_fpb_cercariae_fit.R`   
Data from [Tantawy 2002](https://www.ncbi.nlm.nih.gov/pubmed/12512816) investigating effects of **Butachlor**, **fluazifop-p-butyl** on mortality of *Biomphalaria alexandrina* snails and *Schistosoma mansoni* miracidia and cercariae  
#### Data reported  
+ **Cercariae -** Table 7&8: longitudinal survival of *Schistosoma mansoni* miracidia in control and 4 dose groups 
#### Response functions fit  
+ Cercarial survival as a function of butralin concentration, produced via estimation of LC50 and slope parameters as a function of butralin concentration incorporated into log-logistic d-r curve. See [INSERT LINK TO LARVAL SURVIVAL DYNAMICS EXPLANATION HERE](link.linkylink.com)  
+ Butachlor linear function: `piC.tant02_but.lin_unc` - linear relationship between time-LC50 and concentration  
+ Butachlor exponential function: `piC.tant02_but.exp_unc` - exponential relationship between time-LC50 and concentration  
  + Sampling and data for both: [Raw data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Cercariae/tantawy2002_butachlor_piC_data.png), [Parameter models](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Cercariae/tantawy2002_butachlor_piC_parameter_models.png), [Sampling and observed data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Cercariae/tantawy2002_butachlor_piC_data_sim.png)  
+ Fluazifop-p-butyl linear function: `piC.tant02_fpb.lin_unc` - linear relationship between time-LC50 and concentration  
+ Fluazifop-p-butyl exponential function: `piC.tant02_fpb.exp_unc` - exponential relationship between time-LC50 and concentration  
  + Sampling and data for both: [Raw data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Cercariae/tantawy2002_fpb_piC_data.png), [Parameter models](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Cercariae/tantawy2002_fpb_piC_parameter_models.png), [Sampling and observed data](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tantawy2002/Cercariae/tantawy2002_fpb_piC_data_sim.png)  

### `tchounwou91_malathion_snails_fit.R`  
Data from [Tchounwou et al](https://www.ncbi.nlm.nih.gov/pubmed/1953025) investigating effects of **malathion** on mortality of *Bulinus havanensis* and *Helisoma trivolvis* snails and eggs  
#### Data reported  
+ LC50, LC0, LC90 all with 95%CIs, and slope (not slope parameter, raw slope of probit-log10 transformed D-R function) for adult snails, juvenile snails, and snail eggs  
#### Response functions fit  
+ `muNq_mal_tch91_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Snails/tchounwou1991_malathion_snail_mortality.png)
+ `fNq_mal_tch91_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Snails/tchounwou1991_malathion_snail_fecundity.png)
+ `muNq_trivolvis_mal_tch91_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Snails/tchounwou1991_malathion_snail_trivolvis_mortality.png)
+ `fNq_trivolvis_mal_tch91_uncertainty`
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Snails/tchounwou1991_malathion_snail_trivolvis_fecundity.png)

### `tchounwou91_malathion_miracidia_fit.R`  
Data from [Tchounwou et al](https://www.ncbi.nlm.nih.gov/pubmed/2066529) investigating effects of **malathion** on mortality of *Schistosoma mansoni*  
#### Data reported  
+ Time series of miracidia mortality across 6 malathion concentrations. Also assessed bayluscide, but not interested in effects of molluscicides    
#### Response functions fit  
+ `piM.tch91_mal_unc`
  + Raw time series data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_malathion_miracidia_time_data.png)
  + Time d-r parameters over concentration: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_malathion_miracidia_dose_pars.png)
  + Time series validation: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_malathion_miracidia_validate.png)
  + Sampling and data of d-r function: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_malathion_piM_validate.png)

### `tchounwou91b_fertilizer_miracidia_fit.R`  
Data from [Tchounwou et al](https://www.ncbi.nlm.nih.gov/pubmed/1906907) investigating effects of **ammonium phosphate** and **urea** on mortality/survival of *Schistosoma mansoni*  
#### Data reported  
+ Time series of miracidia mortality/survival across 6 concentrations of urea and ammonium phosphate fertilizers    
#### Response functions fit  
+ `piM.tch91_amm_unc`
  + Raw time series data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_ammonium_phosphate_miracidia_time_data.png)
  + Time d-r parameters over concentration: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_ammonium_phosphate_miracidia_d-r_pars.png)
  + Time series validation: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_ammonium_phosphate_miracidia_validate.png)
  + Sampling and data of d-r function: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_ammonium_phosphate_miracidia_piM_validate.png)
+ `piM.tch91_ure_unc`
  + Raw time series data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_urea_miracidia_time_data.png)
  + Time d-r parameters over concentration: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_urea_miracidia_piM_dose_pars.png)
  + Time series validation: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_urea_miracidia_validate.png)
  + Sampling and data of d-r function: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Miracidia/tchounwou1991_urea_miracidia_piM_validate.png)

### `tchounwou91b_fertilizer_egg_viability_fit.R`  
Data from [Tchounwou et al](https://www.ncbi.nlm.nih.gov/pubmed/1906907) investigating effects of **ammonium phosphate** and **urea** on hatchability of *Schistosoma mansoni* eggs
#### Data reported  
+ Tables 1 & 2: proportion of eggs that hatch into miracidia across concentrations of fertilizers    
#### Response functions fit  
+ `tch91_amm_v_unc`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Eggs/amm_sulphate_egg_viability_data_sim.png)
+ `tch91_ure_v_unc`  
  + Sampling and data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1991/Eggs/urea_egg_viability_data_sim.png)
  
### `tchounwou92_malathion_cercariae_fit.R`  
Data from [Tchounwou et al (1992)](https://onlinelibrary.wiley.com/doi/abs/10.1002/tox.2530070202) investigating effects of **malathion** on mortality/survival of *Schistosoma mansoni* cercariae  
#### Data reported  
+ Figure 2: Time series of cercariae mortality/survival across 6 concentrations of malathion    
#### Response functions fit  
+ `piC.tch92_mal_unc`
  + Raw time series data: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1992/tchounwou92_malathion_cercariae_data.png)
  + Time d-r parameters over concentration: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1992/tchounwou92_malathion_cercariae_dose_pars.png)
  + Time series validation: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1992/tchounwou92_malathion_cercariae_validation.png)
  + Sampling and data of d-r function: [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Tchounwou1992/tchounwou92_malathion_piC_validation.png)
