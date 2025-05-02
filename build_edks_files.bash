name=Iquique
mpirun -n 10 python MPI_EDKS.py edks_INPUT/${name}.edks_config
rm *.asc
