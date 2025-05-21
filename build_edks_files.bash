name=Gorkha
nproc=10
mpirun -n ${nproc} python utils/MPI_EDKS.py edks_INPUT/${name}.edks_config
rm *.asc
