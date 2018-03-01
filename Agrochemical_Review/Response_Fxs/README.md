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
+ Daily mortality rate for *H. duryi* exposed to malathion `muNq_mal_Bakry11_uncertainty` [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Bakry_2011/bakry2011_malathion_snail_mortality.png)
+ Daily mortality rate for *H. duryi* exposed to deltamethrin `muNq_del_Bakry11_uncertainty` [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Bakry_2011/bakry2011_deltamethrin_snail_mortality.png) 
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
+ Daily mortality for *B. alexandrina* exposed to atrazine `muNq_atr_Bakry12_uncertainty` [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Bakry_2012/bakry2012_atrazine_snail_mortality.png)
+ Daily mortality for *B. alexandrina* exposed to glyphosate `muNq_gly_Bakry12_uncertainty` [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Bakry_2012/bakry2012_glyphosate_snail_mortality.png)
#### Notes
+ Longitudinal survival could be used to investigate effects of cumulative exposure with survival analysis
+ Snail reproduction (measured as eggs/snail in addition to hatchlings/snail at discrete times over a four week period) sharply reduced in snails exposed to LC10 of atrazine and glyphosate, but only able to compare control group to each exposed group, therefore unable to fit dose-response function [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Bakry_2012/bakry2012_snail_reproduction.png)  

### `Baxter_Rohr2011_reanalysis_atrazine_snail_carrying_capacity_fit.R`  
Data from [Baxter etal 2011](http://onlinelibrary.wiley.com/doi/10.1002/etc.552/abstract) and reanalysis of it by [Rohr et al 2012](https://www.researchgate.net/publication/224284702_The_herbicide_atrazine_algae_and_snail_populations?enrichId=rgreq-687b4b78-6cbe-4c48-9c04-d495370159ae&enrichSource=Y292ZXJQYWdlOzIyNDI4NDcwMjtBUzoxMDEyNDM0Njg5MTA1OTlAMTQwMTE0OTczMjY1Mg%3D%3D&el=1_x_3) investigating effect of **atrazine** on algal dynamics and ultimate effects on *Physella* spp and *Stagnicola elodes* snail population dynamics  
### Data reported  
+ Fig 5: Snail population over time, reanalyzed by Rohr to estimate peak growth rate prior to population collapse due to resource exhaustion  
+ Figs 4&6: Snail direct effects (reproduction, size, and mortality)  
### Response functions fit  
+ Snail population carrying capacity as a function of atrazine concentration `phi_Nq_atr_baxrohr.no30` [Image](https://github.com/cmhoove14/AgroSchisto/blob/master/Agrochemical_Review/Response_Fxs/Plots/Baxter_Rohr_2011/Baxter_rohr_atrazine_carrying_capacity_sim.png)
+ Snail population carrying capacity as a function of any herbicide concentration, assuming response is proportional to its EEC value and is the same as the atrazine response `phi_Nq_rel_baxrohr.no30`  
  + no30 respods to exclusion of data point at 30ppb atrazine due to large influence of spatial block on results. See [Rohr et al 2012](https://www.researchgate.net/publication/224284702_The_herbicide_atrazine_algae_and_snail_populations?enrichId=rgreq-687b4b78-6cbe-4c48-9c04-d495370159ae&enrichSource=Y292ZXJQYWdlOzIyNDI4NDcwMjtBUzoxMDEyNDM0Njg5MTA1OTlAMTQwMTE0OTczMjY1Mg%3D%3D&el=1_x_3) for details  
  
###   
