name=Tohoku
mpirun -n 4 python MPI_EDKS.py edks_INPUT/${name}.edks_config
rm *.asc
