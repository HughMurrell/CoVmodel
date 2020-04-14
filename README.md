# CoVmodel

## Simple SIR models for COVID-19 breakout

In this repository we present a Julia notebook for
simulating the COVID-19 breakout in South Africa and outher countries using 
a simple SIR model and case counts from open data sets.

### `SIRjulia`

This directory contains a Jupyter notebook that uses `julia` and
the `Optim` package to perform discrete simulation
and optimisation to fit a simple piecewise SIR model to
John Hopkins case count data. The piecewise nature of the model allows
the user to observe the effect of government interventions.

### `SIRjh`

This javascript code runs piecewise SIR simulations in the web browser. 
Data from John Hopkins, cleaned by the julia notebook
above is loaded by the script and then the user can experiment with
the SIR parameters and obtain simulated predictions of how the 
pandemic will play out. The link to the script is
[here](https://hughmurrell.github.io/CoVmodel/SIRjh/index.html)

### `SIRou`

This javascript code also runs SIR simulations in the web browser. 
Case counts and Stringency indices from Oxford University are used 
by the script and then the user can experiment with
the SIR parameters for simulated predictions. 
An attempt will be made to enable the app itself to optimise SIR parameters.
The link to the script is
[here](https://hughmurrell.github.io/CoVmodel/SIRou/index.html)

### Data Sets

`data_jh_clean.csv` cleaned John Hopkins case counts for javascript simulation runs.
`confirmedcases-Table 1.csv` case counts from Oxford
`stringencyindex-Table 1.csv` daily stringency index from Oxford
`wb_population.csv` population data from World Bank

### Visualisation

Below is a a visualisation of the early breakout of COVID-19 in South Africa.

![alt text](https://github.com/HughMurrell/CoVmodel/blob/master/SIRjulia/covid_plots/SIRfitSouthAfrica.png "COVID-19 South Africa")

This graphic can be generated using the `covid_sir.ipynb` notebook from the `SIRjulia`
directory.

