# Agrochemical Review README
Code repository for review and synthesis of agrochemical effects on schistosomiasis transmission 

## `Response_Fxs`  
Folder containing all data, code, and resulting plots used to estimate and analyze agrochemical response functions relating a particular model parameter to agrochemical concentration

### `Response_Fxs/Plots`  
Plots of all data, fitted dose response functions, and results of sampling the d-r function with uncertainty  
### `Response_Fxs/Data`  
Folder containing data extracted from studies identified in review  
### `Response_Fxs/Summary`  
Folder containing summary csv of all response functions and some code to get summary statistics for all response functions   

## `Models`  
Mathematical models, R0(q) functions, parameter sets, and uncertainty analyses  

## `Sims`  
Code for simulating model and estimating R0 in the presence of agrochemical pollution  
### EEC
Simulations and results at each agrochemical's peak Expected Environmental Concentration (*EEC*)  
#### Parameter Sets  
Intermediate data matrices with parameter sets used for simulations at each agrochemical's EEC

### Range
Simulations and results across a range of agrochemical concentrations from 0 - 2*EEC  
#### Parameter Sets  
Intermediate data matrices with parameter sets used for simulations across relevant ranges of agrochemical concnetration  
