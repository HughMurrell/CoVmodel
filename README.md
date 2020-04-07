# CoVmodel

## Simple SIR model for COVID-19 breakout from John Hopkins data

In this repository we present a Julia notebook and a JavaScript app for
simulating the COVID-19 breakout in South Africa using a simple SIR model
and the case counts from John Hopkins data set.

### `SIRjulia`

This directory contains a Jupyter notebook that uses discrete simulation
and QuadDirect optimisation to fit a simple piecewise SIR model to
John Hopkins case coun data. The piecewise nature of the model allows
the user to measure the effect of government interventions.

### `SIRjs`

This directory contains a javascript code for running SIR simulations
in the web browser. Data from John Hopkins, cleaned by the julia notebook
above is loaded by the script and then the user can experiment with
the SIR parameters and obtain simulated predictions of how the 
pandemic will play out.

### Data Sets

`data_jh_clean.csv` cleaned John Hopkins case counts for javascript simulation runs

### Visualisation

Below is a a visualisation of the early breakout of COVID-19 in South Africa.

![alt text](https://github.com/HughMurrell/CoVmodel/blob/master/SIRjulia/covid_plots/SIRfitSouthAfrica.png "COVID-19 South Africa")

This graphic can be generated using the `covid_sir` notebook.

