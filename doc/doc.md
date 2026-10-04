# Documentation
This documentation serves to describe the usage of the two command line scripts, `do_fourier_calculation.py` and `plot_fourier_pattern.py`, which serve to carry out the numerical analysis and plot the results respectively.
## do_fourier_calculation.py
### Arguments
* `--file` image of material to evaluate
* `--intervall` Maximum range of the considered spatial frequencies, which all lie in the square region $\[-\mathrm{intervall},\mathrm{intervall}\]^2$
* `--threshold` required minimum amplitude $\left|\hat{\rho}(\mathbf{k})\right|$ for a fourier coefficient to be shown in the diffraction diagram
* `--thresholdgauge` required minimum amplitude $\left|\hat{\rho}(\mathbf{k})\right|$ for a fourier coefficient to be considered in the determination of symmetries
* `--multiples` All considered spatial frequencies lie on the orbit under the holohedry the set of all linear combinations of `multiples` basis vectors
* `--symmetries` One or more .csv files, linear isometries whose compatibility with the symmetry criterion is evaluated
### Input files
* `file` Grayscale .png image file to evaluate
### Ouput files
## plot_fourier_pattern.py
