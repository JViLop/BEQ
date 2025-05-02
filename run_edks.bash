
name=Iquique
python utils/run_ensemble_edks_GFs.py $name
mpirun -n 12 python utils/run_ensemble_edks_okada_parallel.py $name
python utils/run_edks_mean_errors_parallel.py $name
