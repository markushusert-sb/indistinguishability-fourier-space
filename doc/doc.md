# Documentation
This documentation serves to describe the usage of the two command line scripts, `do_fourier_calculation.py`, `identify_peaks.py` and `plot_fourier_pattern.py`, which serve to carry out the numerical analysis and plot the results respectively.\
The intened way to use these tools is the following:
1. Execute `do_fourier_calculation.py` and `plot_fourier_pattern.py` without specifying the basis vectors (and ignore the resulting error messages) in order to plot the diffraction diagram of the image
2. Identify the rough location of the fundamental frequencies in the diffraction diagram and execute `identify_peaks.py` in order to determine them precisely as the local maximum of the amplitude $\left|\hat{\rho}(\mathbf{k})\right|$.
3. Save the desired basis vectors in `basevectors.csv` and indicate the holohedry as well as the investigated symmetries with regards to these basis vectors. Rerun `do_fourier_calculation.py` and `plot_fourier_pattern.py`
## do_fourier_calculation.py
### Arguments
* `--file` image of material to evaluate
* `--intervall` Maximum range of the considered spatial frequencies, which all lie in the square region $\[-\mathrm{intervall},\mathrm{intervall}\]^2$
* `--threshold` required minimum amplitude $\left|\hat{\rho}(\mathbf{k})\right|$ for a fourier coefficient to be shown in the diffraction diagram
* `--thresholdgauge` required minimum amplitude $\left|\hat{\rho}(\mathbf{k})\right|$ for a fourier coefficient to be considered in the determination of symmetries
* `--multiples` All considered spatial frequencies lie on the orbit under the holohedry of the set of all linear combinations of `multiples` basis vectors
* `--symmetries` One or more .csv files, linear isometries whose compatibility with the symmetry criterion is evaluated
### Input files
* `file` Grayscale .png image file to evaluate
* ``
### Ouput files
## plot_fourier_pattern.py
