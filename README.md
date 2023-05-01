# TPE-Sim1
Work by the EIWG Estimation sub-team on treatment policy estimation. Simulation study comparing MI and MLE estimation methods.

## Description of the work
To add.

Results are included in 7Zip archive file format ( https://www.7-zip.org/ )


### Data Generation Code

|File name | Description|
|------------|----------|
| dar_drop_dgm.R | Data simulation for scenario with discontinuation at random and immediate change to off treatment effects.|
| dnar_drop_dgm.R | Data simulation for scenario with discontinuation not at random and immediate change to off treatment effects.|
| dar_linearchg_dgm.R | Data simulation for scenario with discontinuation at random and linear change to off treatment effects.|
| dnar_linearchg_dgm.R | Data simulation for scenario with discontinuation not at random and linear change to off treatment effects.|
| dar_drop_true.R | Creation of large data to estimate the true value in scenario with discontinuation at random and immediate change to off treatment effects.|
| dnar_drop_true.R | Creation of large data to estimate the true value in scenario with discontinuation not at random and immediate change to off treatment effects.|
| dar_linearchg_true.R | Creation of large data to estimate the true value in scenario with discontinuation at random and linear change to off treatment effects.|
| dnar_linearchg_true.R | Creation of large data to estimate the true value in scenario with discontinuation not at random and linear change to off treatment effects.|


### Estimation Code


|File name | Description|
|------------|----------|
| dar_drop_mmrm1 | Simple MMRM estimation using all observed data and a basic MAR assumption.|
| dar_drop_mi1 | Simple MI estimation using all observed data and a basic MAR assumption|

... To be updated.
