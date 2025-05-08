# Posterior Exploration and Predictive Anlysis of Bayesian finite-fault earthquake models

## Requirements

* okada4py (https://github.com/jolivetr/okada4py.git)
* CSI (https://github.com/jolivetr/csi.git)
* EDKS (https://github.com/JViLop/EDKS_py.git)
* MudPy (https://github.com/UO-Geophysics/MudPy.git)
* Standard packages (in `requirements.txt`)

## Initialize

* Change paths to python libraries in `init.bash` and `/.edks_config` files according to own directory structure

```
cd path/to/BEQ
source init.bash
```
* Remove `.gitignore` files in `path/to/BEQ/EQ/name/kinematic/original_model_file` and place respective model file `.dat` or `.h5`
* Copy python modules in `mudpy_src` and place them in `MudPy/src/python`

## Pre-process

For event called `Name`, in terminal run `python ${Name}.py` (e.g., `Iquique.py`). Other functionalities and parameters (such as `nsamples`) can be set up in python module. 
This creates formatted files for subsequent analysis/calculations in `BEQ/INPUT`.
## Posterior Analysis

* For individual event correlation, run `python corr_parameters.py`.
* For multi-event correlation comparison, run `python cross_event_corr.py`. 

## Run Okada

For all events,  run `python run_okada.py` for displacement and stress calculations. This produces stress change calculations.
Note: Outputs are stored in `BEQ/OUTPUT`

## Plot stress/strain drop distributions

For all avents, run `python average_stress_strain_drop.py` for stress/strain drop histograms

 
## Run EDKS

For event called `Name`, run the following commands in sequence:
* `source build_edks_files.bash` to create edks files in layered media
* `source run_edks.bash` to calculate mean surface displacement and uncertainty.
 
Notes: Outputs are stored in `BEQ/OUTPUT`. Update number of processors accordingly in both `.bash` files. Check line 102 in `BEQ/utils/run_ensemble_edks_GFs.py` and update location of `BIN_EDKS` accordingly
 
## Compare curves of slip, (EDKS) displacement and (Okada) stress

For all events, run `python curves_all.py`

## Compare EDKS displacement correlation matrix with distance 

For all events, run `python corr_matrix.py`

## Run dynamic simulations

For event called `Name`, run:
*  `source run_FK.bash`to calculate and assemble displacement waveforms.
* `python save_FK.py` to store waveform ensemble in h5 files. 

Note: Outputs are stored in `BEQ/Dynamic_Simulations`


## Create dynamic displacement field snapshots

For all events, run:
* `python vector_field_snapshots.py`to create dynamic displacement field snapshots
* `python std_vector_field_snapshots.py` to create uncertainty snapshots (pending integration)

 
## Compare PGD with distance/azimuth

For all avents, run `python PGD_wrt_distance_and_azimuth.py` (pending integration)

## Create animations of dynamic displacement field

For all events, run `python animate.py` (pending integration)

## Plot slip and displacement with maps

* For event called `Name`, run `python ${name}_georef.py` to plot slip with pygmt (segmentation fault issues when run in linux)
* For all events, run `python deformation_map.py` to plot EDKS displacement with cartopy (segmentation fault issues when run in linux)






