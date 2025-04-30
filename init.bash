# activate virtual environment 
source ~/.venv/benchmark/bin/activate

# EDKS 
export EDKS_HOME=${HOME}/EDKS_py
export PYTHONPATH=${EDKS_HOME}/MPI_EDKS:$PYTHONPATH
export PYTHONPATH=${EDKS_HOME}/MPI_EDKS_Py3:$PYTHONPATH
export PATH=${EDKS_HOME}/MPI_EDKS:$PATH
export PATH=${EDKS_HOME}/src/sum_layered_sum:$PATH
export EDKS_BIN=${EDKS_HOME}/bin

# NUMBER OF CORES FOR OPENMPI
export OMP_NUM_THREADS=4

# MPI4PY 
mpi4py_HOME=/home/josevilo/.venv/benchmark/lib/python3.8/site-packages/mpi4py
export PYTHONPATH=${mpi4py_HOME}:${PYTHONPATH}

# CSI
export PYTHONPATH=$HOME/static:$PYTHONPATH 

# MudPy
export PYTHONPATH=$HOME/MudPy/src/python:$PYTHONPATH
export PATH=$HOME/MudPy/src/fk:$PATH
export MUD=$HOME/MudPy

