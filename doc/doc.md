# Documentation
This documentation serves to describe the usage of the two command line scripts, `do_fourier_calculation.py`, `identify_peaks.py` and `plot_fourier_pattern.py`, which serve to carry out the numerical analysis and plot the results respectively.\
The intened way to use these tools is the following:
1. Execute `do_fourier_calculation.py` and `plot_fourier_pattern.py` without specifying the basis vectors (and ignore the resulting error messages) in order to plot the diffraction diagram of the image
2. Identify the rough location of the fundamental frequencies $\mathbf{k}_i$ in the diffraction diagram and execute `identify_peaks.py` in order to determine them precisely as the local maximum of the amplitude $\left|\hat{\rho}(\mathbf{k})\right|$.
3. Save the desired basis vectors in `basevectors.csv` and indicate the holohedry as well as the investigated symmetries with regards to these basis vectors. Rerun `do_fourier_calculation.py` and `plot_fourier_pattern.py`
## do_fourier_calculation.py
This function carries out the numerical analysis of the image and produces several `.csv` as outputs which are documented below
### Arguments
* `--file` image of material to evaluate
* `--interval` Maximum range of the considered spatial frequencies, which all lie in the square region $\[-\mathrm{interval},\mathrm{interval}\]^2$
* `--threshold` required minimum amplitude $\left|\hat{\rho}(\mathbf{k})\right|$ for a fourier coefficient to be shown in the diffraction diagram
* `--thresholdgauge` required minimum amplitude $\left|\hat{\rho}(\mathbf{k})\right|$ for a fourier coefficient to be considered in the determination of symmetries
* `--multiples` All considered spatial frequencies lie on the orbit under the holohedry of the set of all linear combinations of `multiples` basis vectors
* `--symmetries` One or more .csv files, linear isometries whose compatibility with the symmetry criterion is evaluated
### Input files
* `file`: Grayscale .png image file to evaluate
* `basevectors.csv:` A .csv file with with two columns corresponding to x and y directions and $\mu$ rows corresponding to the $\mu$ different basis vectors
* .csv files indicated under `symmetries`: Transformations whose symmetric properties are to be evaluated. Each transformation $\mathbf{q}$ is indicated as an element of $\mathrm{GL}(\mathbb{N},\mu)$ by its action on the $\mathbf{k}_i$:
```math
\mathbf{q}:\mathcal{M}^{\mathrm{s}}\rightarrow\mathcal{M}^{\mathrm{s}}\quad \alpha_i\mathbf{k}_i\mapsto q_{ij}\alpha_j\mathbf{k}_i
```
### Ouput files
In the following, `<filename>` corresponds to the name of the investigated image file without the ending `.png`, and `<symmetryname>` corresponds to the name of a `.csv` file provided under `--symmetries` without the ending `.csv`:
* `diffractiondiagram_<filename>.csv` A .csv file where: with four columns describing 
	* Each row describes the coefficient $\hat{\rho}(\mathbf{k})$ of the discrete Fourier transform of the image with amplitude is above `threshold`.
	* The four columns correspond to: the horizontal component $k_x$, the vertical component $k_y$, the amplitude $\left|\hat{\rho}(\mathbf{k})\right|$ and the complex phase $\hat{\rho}(\mathbf{k})$
* `wavevectors_base_<file>.csv` and `wavevectors_other_<file>.csv` A .csv file with one row for each element of $\mathcal{M}^{\mathrm{s}}$ and four columns corresponding to:
	1. The horizontal component $k_x$
	2. The vertical component $k_y$
	3. The amplitude $\left|\hat{\rho}(\mathbf{k})\right|$
	4. The complex phase $\mathrm{arg}(\hat{\rho}(\mathbf{k}))$
The fourier coefficient $\hat{\rho}(\mathbf{k})$ is calculated as
```math
\hat{\rho}(\mathbf{k})=\frac{\sum_{i=1}^{m}\sum_{j=1}^{n}\rho_{ij} \mathrm{e}^{-2 \pi \mathrm{i} (k_x(j-n/2-1)+k_y(i-m/2-1))/\sqrt{mn}}}{nm}
```
where $\rho_{ij}$ denotes the grid of grayscale values of the image (using 0-based indexing).
The exponent $(k_xj+k_yi)/\sqrt{mn}$ is normed to $\sqrt{mn}$ and not the respective numbers of pixels $m$ or $n$ in each direction in order to maintain the orthonormality of the basis of the Fourier space in case of a rectangular image
* `gaugeerrors_<filename>_<symmetryname>.csv` A .csv file with one row for each element of $\mathcal{M}^{\mathrm{s}}$ and six columns corresponding to
	1. The horizontal component $k_x$
	2. The vertical component $k_y$
	3. The phase function $\Phi_{\mathbf{q}}(\mathbf{k})$
	4. Its associated error measure $\Delta\Phi_{\mathbf{q}}(\mathbf{k})$
	5. The amplitude $\left|\hat{\rho}(\mathbf{k})\right|$
	6. Its associated error measure $\Delta \hat{\rho}$
* `Indices_p4msquarerandomFFTbase_pixel_3201.csv` A .csv file with one row for each element of $\mathcal{M}^{\mathrm{s}}$ and $\mu+2$ columns corresponding to:
	* First $\mu$ columns: the indices $\alpha_i$ of the frequency $\mathbf{k}=\alpha_i\mathbf{k}_i$ 
	* Last two columns: the horizontal and vertical components $k_x$ and $k_y$
## identify_peaks.py
This script allows to quantify the position of a basis vector of $\mathcal{M}^{\mathrm{s}}$ as the local maximum of $\left|\hat{\rho}(\mathbf{k})\right|$ around an estimation $\mathbf{k}^{\mathrm{g}}$
### Arguments
* `basisvector` two positional arguments, components $\mathbf{k}^{\mathrm{g}}_x$ and $\mathbf{k}^{\mathrm{g}}_y$
* `--file` image of material to evaluate
* `--nraster` the number of points in each dimension constituting the raster of the region $\mathbf{k}^{\mathrm{g}}+\[-0.5,0.5\]^2$ on which the local maximum is determined.
## plot_fourier_pattern.py
This script plots the data exported in the `.csv` files of `do_fourier_calculation.py`
### Arguments
execute `plot_fourier_pattern.py -h` for a description of the arguments
### Input files
`.csv` files created by `do_fourier_calculation.py`
### Ouput files
In the following, `<filename>` corresponds to the name of the investigated image file without the ending `.png`, and `<symmetryname>` corresponds to the name of a `.csv` file provided under `--symmetries` without the ending `.csv`:
* `wavevectors_<filename>.png`, `wavevectors_<filename>amplitudes.png`: studied Fourier coefficients $\hat{\rho}(\mathcal{M}^{\mathrm{s}})$. The former represents the complex phase by color-coding, the latter does not.
* `diffraction_threshold_<filename>.png`, `diffraction_threshold_<filename>amplitudes.png`: Coeficcients of the discrete Fourier transform that lie within the region specified by `interval` and whose amplitudes lie above `threshold`. The former represents the complex phase by color-coding, the latter does not.
* `fourier_module_<filename>.png`: The set $\mathcal{M}^{\mathrm{s}}$ of studied frequencies.
* `overlay_considered_vectors_diffraction_<filename>.png`: The set $\mathcal{M}^{\mathrm{s}}$ of studied frequencies overlayed with the coefficients of the diffraction diagram plotted in `diffraction_threshold_<filename>.png`.
* `(amplitudeerrors|gaugeerrors|combinederrors)_<filename>_<symmetryname>.png` values of one of the three deviation measures on $\mathcal{M}^{\mathrm{s}}$
* `(amplitudeerrors|gaugeerrors|combinederrors)_<filename>_<symmetryname>.png` histogram of the values of one of the three deviation measures
