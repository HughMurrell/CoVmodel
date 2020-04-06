# CoVmodel

## Simple SIR model for COVID-19 breakout in South Africa

In this repository we present a Julia notebook and a JavaScript app for
simulating the COVID-19 breakout in South Africa using a simple SIR model
and the case counts from John Hopkins data set.

### `TempDataClean`

This notebook gives some idea on how the data sets from NOAA were cleaned
in order to generate the latest recorded average monthly temperatures for
weather stations around the globe.

It is not possible to run this notebook as is because the input data from
NOAA is too large storage on github repositories


### `TempInterp`

This notebook reads the cleaned version of NOAA monthly temperature data from weather stations
around the world and then iterpolates weekly temperatures for any point on the globe
from the nearest neighbor in the cleaned data set.

### Data Sets

`data_clean.csv` cleaned NOAA average monthly temperatures for stations on the globe.

`airports.csv` list of airports around the globe with their geo-coordinates

`airports_temp_augmented.csv` same as `airports` data but augmented by average weekly temperatures
obtained by running `TempInterp` notebook


### Visualisation

Below is a a visualisation of the early breakout of COVID-19 in South Africa.

![alt text](https://github.com/HughMurrell/CoVmodel/SIRjulia/covid_plots/SIRfitSouthAfrica.png "COVID-19 South Africa")

This graphic can be generated using the `covid_sir` notebook.

