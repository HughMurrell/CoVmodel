# CoVmodel

## Simple SIR models for COVID-19 breakout

In this repository we present Julia notebooks for modeling the
COVID-19 breakout in South Africa and outher countries using 
case counts from open data sets.

### `SIRjulia`

This directory contains a Jupyter notebook that uses `julia` and
the `Optim` package to perform discrete simulation
and optimisation to fit a simple piecewise SIR model to
John Hopkins case count data. The piecewise nature of the model allows
the user to observe the effect of government interventions.

An accompanying javascript code, `SIRjh`,  runs piecewise SIR simulations 
in the web browser.  Data from John Hopkins, cleaned and augmented by 
the julia notebook above is loaded by the script and then the user can 
experiment with the SIR parameters and obtain simulated predictions of how the 
pandemic will play out. The link to the javascript is
[here](https://hughmurrell.github.io/CoVmodel/SIRjh/index.html)


Another javascript code,  `SIRou`, also runs SIR simulations in the web browser. 
Case counts and Stringency indices from Oxford University are used 
by the script and then the user can experiment with
the SIR parameters for simulated predictions. This app is independent of
julia as an attempt has been made to enable the app itself to optimise SIR parameters.
The link to the javascript is
[here](https://hughmurrell.github.io/CoVmodel/SIRou/index.html)

### `RtLive`

This directory contains a Jupyter notebook with Julia code for estimating
$R-T$ in real time. The code uses maximum likelyhood to estimate the
"most likely" $R$ value for a population at time $t$. The code is a
direct translation from the Python code discribed at [http://rt.live/](http://rt.live/)
and made available to the public as a Jupyter notebook with accompanying
explanations of the methods employed.

### Data Sets

`data_jh_clean.csv` cleaned John Hopkins case counts for javascript simulation runs.
`confirmedcases.csv` case counts from Oxford
`stringencyindex.csv` daily stringency index from Oxford
`wb_population.csv` population data from World Bank

### Visualisation

Below is a visualisation of the early breakout of COVID-19 in South Africa.

![alt text](https://github.com/HughMurrell/CoVmodel/blob/master/SIRjulia/covid_plots/SIRfitSouthAfrica.png "COVID-19 South Africa")

This graphic can be generated using the `covid_sir.ipynb` notebook from the `SIRjulia`
directory.

Next is a visualisation of $R_t$ over time for the COVID-19 South African case counts.

![alt text](https://github.com/HughMurrell/CoVmodel/blob/master/RtLive/plots/Rt_SouthAfrica.png "Rt Live, South Africa")

This graphic can be generated using the `RtLive.ipynb` notebook from the `RtLive`
directory.
