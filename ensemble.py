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
print("Compiling")
os.system('./compil.sh') 

# generate meteorology perturbations
Pmlt = np.random.lognormal(-0.5*Pstd**2,Pstd,Nens)
Tadd = np.random.normal(0.,Tstd,Nens)

# write the different input files nlst
print("Writing nlst_<i> files")
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

# submit task array on forth 
# request for example 100 simulations with nCPUs running at once
submit_command = "sbatch --array=1-"+str(Nens)+"%"+str(nCPUs)+" submit.sh"
print("Submitting job: " , submit_command)
os.system(submit_command)

