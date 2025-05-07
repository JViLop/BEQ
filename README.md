# Posterior Exploration and Predictive Anlysis of Bayesian finite-fault earthquake models

## Requirements

* okada4py (https://github.com/jolivetr/okada4py.git)
* CSI (https://github.com/jolivetr/csi.git)
* EDKS (https://github.com/JViLop/EDKS_py.git)
* MudPy (https://github.com/UO-Geophysics/MudPy.git)
* Standard libraries (e.g., numpy, meson, ninja, ...)

## Initialize

* Change paths to python libraries in `init.bash` in own directories and `/.edks_config` files

```
cd path/to/BEQ
source init.bash

```
* Remove `.gitignore` files in `path/to/BEQ/EQ/name/kinematic/original_model_file`

## Pre-process

Given earthquake called `Name`, in terminal run `source ${Name}.py` (e.g., `Iquique.py`). Other functionalities and parameters (such as `nsamples`) can be set up in python module. 
This creates formatted files for subsequent analysis/calculations.

## Posterior Analysis

* For individual event correlation, in terminal run `source corr_parameters.py`.

* For multi-event correlation comparison, run `source cross_event_corr.py`. 

## Run Okada

For all events,  in terminal run `source run_okada.py` for displacement and stress calculations. This produces stress change calculations.

## Plot stress/strain drop distributions

For all avents, in terminal run `source average_stress_strain_drop.py` for stress/strain drop histograms

 
## Run EDKS
Given earthquake name 

