"""
Run model ensemble
"""
import numpy as np
import os

Nens = 100    # ensemble size
Pstd = 0.63   # precipitation multiplier standard deviation
Tstd = 1.0    # temperature offset standard deviation

# compile fortran code
os.system('./compil.sh') 

# open loop simulation
nlst = open('nlst','w')
nlst.write('&params \n') 
nlst.write('/ \n') 
nlst.close()
os.system('./FSM2 < nlst')
os.system('mv FSM2out.nc FSM2out_opn.nc')

# generate meteorology perturbations
Pmlt = np.random.lognormal(-0.5*Pstd**2,Pstd,Nens)
Tadd = np.random.normal(0.,Tstd,Nens)

# run the ensemble and copy outputs
for n in range(Nens):
    nlst = open('nlst','w')
    nlst.write('&params \n') 
    nlst.write('  Pmlt = '+str(Pmlt[n])+' \n') 
    nlst.write('  Tadd = '+str(Tadd[n])+' \n') 
    nlst.write('/ \n') 
    nlst.close()
    os.system('./FSM2 < nlst')
    os.system('mv FSM2out.nc FSM2out_'+str(n).zfill(3)+'.nc')

