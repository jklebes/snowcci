ensemble.py uses compil.sh to compile the fortran code in src and runs an ensemble of simulation with inputs from met/.

Notes:

- install fortran netcdf, on ubuntu ``apt install libnetcdff-dev``
- Install a python ``<``3.11 environment with dask
- Currently fortran CPU and memory activity does not show up on dask taskboard
