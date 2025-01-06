ensemble.py uses compil.sh to compile the fortran code in src and runs an ensemble of simulation with inputs from met/.

Version using gnu-parallel bash command

Notes:
- install fortran netcdf, on ubuntu ``apt install libnetcdff-dev``
- install gnu-parallel with `conda install parallel -c conda-forge`
- edit nCPUs in script ensemble.py
