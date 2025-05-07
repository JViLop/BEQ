# Posterior Exploration and Predictive Anlysis of Bayesian finite-fault earhtqyahe models

## Requirements

* okada4py (https://github.com/jolivetr/okada4py.git)
* CSI (https://github.com/jolivetr/csi.git)
* EDKS (https://github.com/JViLop/EDKS_py.git)
* MudPy (https://github.com/UO-Geophysics/MudPy.git)
* Standard libraries (e.g., numpy, meson, ninja, ...)

## Initialize

* Change paths to python libraries in `init.bash` in own directories and `.edks_config` files

```
cd path/to/BEQ
source init.bash

```
* Remove `.gitignore` files in `path\to\BEQ\EQ\name\kinematic\original_model_file`

## Pre-process
Given earthquake called `name`, in terminal run `source ${name}.py`.Number of samples `nsamples` can be adjusted in python module.

