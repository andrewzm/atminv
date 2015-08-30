Reproducibility
===============

Here we store vignettes that were compiled under different systems. All the vignettes have slightly different results but not significantly so. The biggest gap appears to be between the vignettes compiled on Windows platforms and the Linux distros. In particular, the acceptance rate of the HMC sampler jumps to slightly above 60\% on a Windows platforms. These differences lie in the underlying linear algebraic roundoffs that vary across OSs and CPU types.

The filenames should be self-explanatory. The vignette name is followed by the Operating System (with `W7` implying Windows 7 and `W10` implying Windows 10). If the system had `OpenBLAS` running this is also indicated in the filaname. Otherwise, the standard R BLAS package is used.