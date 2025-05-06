name=Tohoku
mpirun -n 4 python utils/MPI_EDKS.py edks_INPUT/${name}.edks_config
rm *.asc
