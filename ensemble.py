"""
Run model ensemble 
"""

# minimal change version using gnu-parallel bash command.
# To use: `conda install parallel -c conda-forge` or 
# distribution-specific gnu-parallel install

import numpy as np
import os

nCPUs = 2
Nens = 100    # ensemble size
Pstd = 0.63   # precipitation multiplier standard deviation
Tstd = 1.0    # temperature offset standard deviation

# compile fortran code
os.system('./compil.sh') 

# generate meteorology perturbations
Pmlt = np.random.lognormal(-0.5*Pstd**2,Pstd,Nens)
Tadd = np.random.normal(0.,Tstd,Nens)

# write the different input files nlst
for n in range(Nens):
    nlst = open('nlst_'+str(n),'w')
    nlst.write('&params \n') 
    nlst.write('  Pmlt = '+str(Pmlt[n])+' \n') 
    nlst.write('  Tadd = '+str(Tadd[n])+' \n') 
    nlst.write('/ \n') 
    nlst.write('&outputs \n') 
    nlst.write('  runid = '+str(n)+'_ \n') 
    nlst.write('/ \n') 
    nlst.close()

# process input files by launching parallel ./FSM2 commands 
os.system("ls nlst_* | parallel --jobs "+str(nCPUs)+" '(cat {} | ./FSM2 )'")

