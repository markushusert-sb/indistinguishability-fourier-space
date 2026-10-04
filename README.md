# Detection of symmetries under the criterion of indistinguishability in Fourier space
## Setup
* Clone this repository to /path/to/repo
* Add /path/to/repo  to your $PATH
* Add /path/to/repo to your Wolfram init.m file
	This file can be located (normally ~/.Wolfram/Kernel/init.m) using
  ```bash
	wolframscript -code 'Print[FileNameJoin[{$UserBaseDirectory, "Kernel", "init.m"}]]'
  ```
 	and the path can be added by appending the line
	```
	AppendTo[$Path, "/path/to/repo"];
  ```
## Theory
This code evaluates the symmetries of a material under the weak symmetry criterion of indistinguishability in the Fourier transform of a provided grayscale image.
A full documentation of the theory and methodology can be found in our published [paper](https://www.sciencedirect.com/science/article/pii/S0022509626001201).

The material is described by a characteristic function $\rho:\mathbb{R}^2\rightarrow\[0,1\]$ which is given by the pixels of the grayscale image.
In real space, the criterion of indistinguishability between two functions $\rho$ and $\rho'$ amounts to the equality of the correlation functions
```math
C_{n}(\mathbf{x}_1,..,\mathbf{x}_{n-1})=\lim_{\Omega \to \infty} \frac{1}{\left|{\Omega}\right|} \int_{\Omega}\rho(\mathbf{x}_1)\rho(\mathbf{x}_1+\mathbf{x})..\rho(\mathbf{x}_{n-1}+\mathbf{x}) \mathrm{d} \mathbf{x}
```
for any $n \in \mathbb{N}$.\
However, this program operates on the Fourier coefficients $\hat{\rho}(\mathbf{k})$, where the criterion takes the equivalent form
```math
\hat{\rho}'(\mathbf{k})=\mathrm{e}^{2\pi\mathrm{i}\chi(\mathbf{k})}\hat{\rho}(\mathbf{k})
```
The function $\chi(\mathbf{k})$, called gauge function, is a linear function whose codomain is the unit torus $\mathbf{R}/\mathbf{Z}$.
It follows that a transformation $\mathbf{q}\in\mathrm{O}(2)$ is a symmetry if it can be associated to such a gauge function $\Phi_{\mathbf{q}}(\mathbf{k})$.
## Method
This code numerically evaluates the deviation from the criterion of indistinguishability for given transformations $\mathbf{q}\in\mathrm{O}(2)$ on a finite set $\mathcal{M}^\mathrm{s}$ of spatial frequencies.
This set is determined as follows:
* A set of $\mu$ basis vectors $\mathbf{k}_i$ is chosen;
* The set
```math
\mathcal{M}(m)=\left\{\sum_{i=1}^{\mu}\alpha_i \mathbf{k}_i \enspace|\enspace\alpha_i \in \mathbb{N} \text{ and } \sum_{i=1}^{\mu}\alpha_i \leq m\right\}
```
  of the linear combinations of up to $m=$`multiples` basis vectors is defined;
* The orbit $\mathcal{H}\bigstar\mathcal{M}(m)$ under the holohedry $\mathcal{H}$ is generated;
* The studied frequencies $\mathcal{M}^\mathrm{s}$ are selected as the orbit under $\mathcal{H}$ of those frequencies of $\mathcal{H}\bigstar\mathcal{M}(m)$ whose amplitudes $\left|\hat{\rho}(\mathbf{k})\right|$ are above `thresholdgauge`.
 
The deviation of the operation $\mathbf{q}$ from indistinguishability is evaluated numerically using three deviation measures:
1. A deviation measure on the amplitudes $\left|\hat{\rho}(\mathbf{k})\right|$, i.e. does the action of $\mathbf{q}$ amount to arbitrary phase shift $\Phi_{\mathbf{q}}$.
2. A deviation measure on the complex phases $\mathrm{arg}(\hat{\rho}(\mathbf{k}))$, i.e. does $\Phi_{\mathbf{q}}$ respect the gauge-linearity.
3. A resulting deviation measure that combines the two previous ones.
## Usage
This code may be used in two different ways:
* From the command line, using the script `run.sh` followed by `plot.sh`;
* From a mathematica notebook, using directly the relevant functions of `Fouriergroups.m`.

Both are examplified using a periodic and a quasiperiodic material in `exa/cmd_line/` and `exa/notebook/` respectively.
The documentation of the inputs and outputs of the command line scripts can be found [here](doc/doc.md).
