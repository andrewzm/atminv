atminv
=====

This page hosts a supplement package to the paper 'Spatio-temporal bivariate
    statistical models for atmospheric trace-gas inversion by Zammit-Mangion
    et al. (2016, in press). The accompanying vignette in this package
    may be used to reproduce the results shown in Section 4.1 of the paper.

To download the vignette without reproducing it on your machine, please view it in the `vignettes` folder or by directly click on the following links:

[Vignette (Section 4.1)](https://github.com/andrewzm/atminv/blob/master/vignettes/chemo_sim_study.pdf?raw=true)


Presentations
-----

The presentation I gave for the 27th Goulburn Research Fellows Meeting (2015) that is mostly focused on bivariate conditional models but that also contains some slides on this topic is [here](https://github.com/andrewzm/bicon/blob/master/pres/2015_03_Goulburn.pdf?raw=true).

The presentation I gave in the Spatial Statistics Conference in Avignon (2015) that is centred around this work is [here](https://github.com/andrewzm/bicon/blob/master/pres/2015_06_Zammit.pdf?raw=true).



Reproducibility 
-------
    
If you wish to reproduce the results, you will need to install the package and its dependencies. One of its dependencies, `hmc`, needs to be installed manually. To do this please open an R consoles, make sure you have `devtools` installed, and then type

    library(devtools)
    install_github("andrewzm/hmc",,dependencies=T)

To install the `atminv` package and the other dependencies, please make sure you have `devtools` loaded and then type

    install_github("andrewzm/atminv",build_vignettes=T,dependencies=T)
  
If you do not wish to compile the vignette (which takes a while) please set `build_vignettes=F` above. Otherwise, to view the vignette please type

    library(atminv)
    vignette("chemo_sim_study")
    
and select the vignette under `atminv`.



Details
-------

Section 4.1 details three sub-experiments. The first is the vanilla case with no model misspecification. The second assumes model misspecification, specifically that the flux field is uncorrelated even though it is not. The third assumes no misspecification but uses 1000 observations and shows that the EM-Laplace algorithm can run even in this scenario. If you wish to run these three different scenarios please download the `.Rnw` file from the `vignettes`, load `RStudio` and set `Tools -> Global Options -> Sweave -> Weave Rnw files using knitr`, and then press `Compile PDF`. Before doing so:

- To run the first case make sure `model = "full"` and `misspecification = 0`.
- To run the second case make sure `model = "full"` and `misspecification = 1`.
- To run the third case make sure `model = "full_big"` and `misspecification = 0`.

Note that if `load_results <- 1` then the computationally intensive results are simply loaded from cache. To run from scratch please set `load_results <- 0`.

UK and Ireland case study
-------

The file `chemo_UK_study.R` contains the code for Example 4.2. Unfortunately the files are too large to upload on Github so they have been uploaded [here](http://hpc.niasra.uow.edu.au/ckan/dataset/example-dataset-for-atmospheric-trace-gas-inversion) instead. We only provide an `R` file for this example, which, however is self-contained. The `R` file is named `chemo_UK_study.R` and is found in the root folder on this page. It should run without any problems when all the datasets have been extracted and put in a `../data/` folder. Obviously this path can be changed within the code to suit one's needs.



Versions
--------

The results in the paper were carried out with 

    R version 3.2.0 (2015-04-16)
    Platform: x86_64-pc-linux-gnu (64-bit)
    Running under: Ubuntu 14.04.2 LTS

The `Matrix` package is v1.2-0. Results may change slightly across different `R` and package versions and across different platforms. To see the vignettes obtained under different platforms, see the pdf files within the folder `reproducibility`.

References
-----

Zammit-Mangion, A., Cressie, N., Ganesan, A.L., O'Doherty, S., Manning, A.J. (2015). Spatio-temporal bivariate statistical models for atmospheric trace-gas inversion. Chemometrics and Intelligent Laboratory Systems, doi: 10.1016/j.chemolab.2015.09.006, in press.