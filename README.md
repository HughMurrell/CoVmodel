# CoVmodel

## Simple SIR models for COVID-19 breakout

In this repository we present Julia notebooks for modeling the
COVID-19 breakout in South Africa and outher countries using 
case counts from open data sets.

### `RtLive`

This directory contains a selection of Jupyter notebooks with Julia code for estimating
Rt in real time. The code uses maximum likelyhood to estimate the
"most likely" R value for a population at time t. Version 1 of the code is a
direct translation from Kevin Systrom's Python code discribed at [http://rt.live/](http://rt.live/)
and made available to the public as a Jupyter notebook with accompanying
explanations of the methods employed. Version 2 and 3 are modifications
of that code with version 3 implementing a simple fast discrete estimation of
Rt from daily case counts.

### `RtPaper`

This directory contains latex source for an article describing how we perform simple
Bayes-free estimation of Rt from daily case counts.

### `SIRjulia`

This directory contains a Jupyter notebook that uses `julia` to perform discrete simulation
and optimisation to fit a simple piecewise SIR model to John Hopkins case count data. 
The piecewise nature of the model allows the user to observe the effect of government 
interventions. A spreadsheet is also provided to create simple SIR simulations.

### Visualisation

Below is a visualisation of the early breakout of COVID-19 in South Africa.

![alt text](https://github.com/HughMurrell/CoVmodel/blob/master/SIRjulia/covid_plots/SIRfitSouthAfrica.png "COVID-19 South Africa")

This graphic can be generated using the `covid_sir.ipynb` notebook from the `SIRjulia`
directory.

Next is a visualisation of Rt over time for the COVID-19 South African case counts.

![alt text](https://github.com/HughMurrell/CoVmodel/blob/master/RtLive/plots/Rt_SouthAfrica.png "Rt Live, South Africa")

This graphic can be generated using the `RtLive.ipynb` notebook from the `RtLive`
directory.
