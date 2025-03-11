# Data and code associated with the paper "Nematomorph parasites potentially drive nutritional flow-from terrestrial to-aquatic ecosystems"
## Code files
> 1. SourceCode_EPAcont.R
> 2. SourceCode_EPA_per_Individual.R
> 3. SourceCode_DailyEPAintake.R

## Data files
> 1. SourceData_EPAcont_Aqua.csv
> 2. SourceData_EPAcont_Cricket.csv
> 3. SourceData_BodyMass.csv
> 4. SourceData_ConsumptionRate_SCW.csv
> 5. SourceData_ConsumptionRate_WT.csv
https://doi.org/10.5281/zenodo.15003386

## Description
### 1. SourceCode_EPAcont.R
This is the code for drawing Figure 1B, probability distributions of EPA content for camel crickets and aquatic invertebrates.
"SourceData_EPAcont_Aqua.csv" will be applied to a Bayesian random-effects meta-analysis to estimate posterior distributions of EPA content by synthesizing multiple studies.
"SourceData_EPAcont_Cricket.csv" will be applied to a Bayesian random sampling approach to estimate EPA content of camel cricket.

### 2. SourceCode_EPA_per_Individual.R
This is the code for drawing Figure 1C, the mean and 95% credible intervals of EPA content per individual prey ingested by fish consumers.
"SourceData_BodyMass.csv" will be applied to a Bayesian random sampling approach to estimate body mass of aquatic invertebrates and camel cricket.
EPA content per individual will be  calculated by multiplying 1000 randomly drawn samples from the posterior distributions of EPA content and body mass.

### 3. SourceCode_DailyEPAintake.R
This is the code for drawing Figure 1D, The daily area-based EPA intake of fish from camel crickets and aquatic invertebrates.
"SourceData_ConsumptionRate_SCW.csv" and "SourceData_ConsumptionRate_WT.csv" will be applied to the model of daily consumption rate.
The details are available on Supplemental information.



