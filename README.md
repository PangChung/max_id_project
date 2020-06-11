# Assessing-Non-Stationary-Heatwave-Hazard-with-Magnitude-Dependent-Spatial-Extremal-Dependence
It contains the code for a toy example and the data used in our paper https://arxiv.org/abs/2006.01569.

The original data are public available from https://www.ecad.eu//dailydata/predefinedseries.php, which contains daily maximum temperature recored over 20,000 monitor stations across Europe. We selected 44 of them. 

To reproduce the results, you can load the data.Rdata file and replace the corresponding variables in the example.R. It might take a fair long time to run. 


Varaiables in the data.Rdata file:

data.info: the station index and details extracted from the original dataset available from https://www.ecad.eu//dailydata/predefinedseries.php.

coord: distance matrix (unit: 1000km) of 44 stations.

data: aggregated yearly maximum temperature (unit: 0.1C) for 44 stations.

U: pesudo uniform data transfered using the fitted Generalized additive model with penalized cubic regression spline.

reg: spatial covariates linked to range parameter (intercept and altitiude (km)).

reg.t: temporal covariates (equally spaced between 0 and 1).
