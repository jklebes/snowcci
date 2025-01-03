"""
Run model ensemble
"""
import numpy as np
import os
import dask.bag as db
from dask.distributed import Client
import dask.delayed
import time

# function which writes info from the tuple t=(label, Pmlt, Tadd) to 
# input file nlst_i and returns file name
def generateInputFile(t):
    filename = "nlst_"+t[0]
    nlst = open(filename,'w')
    print("Writing ", filename)
    nlst.write('&params \n') 
    nlst.write('  Pmlt = '+str(t[1])+' \n') 
    nlst.write('  Tadd = '+str(t[2])+' \n') 
    nlst.write('/ \n') 
    nlst.write('&outputs \n') 
    nlst.write("  runid = '"+t[0]+"_' \n") 
    nlst.write('/ \n') 
    nlst.close()
    #t1 =  t[1]
    #while t1>0:
    #    t1 -=.00000001
    return filename

def runFSM2(filename):   
    print("Running ./FSM2 <", filename)
    os.system('./FSM2 < '+ filename)
    # return exit code ?

if __name__=="__main__":
    # just for info, print link to show computation
    client = Client()
    print(client)
    print("Click here to see dask dashboard: ", client.dashboard_link)


    Nens = 100    # ensemble size
    Pstd = 0.63   # precipitation multiplier standard deviation
    Tstd = 1.0    # temperature offset standard deviation

    # compile fortran code
    os.system('./compil.sh') 

    # Start dask bags

    # 100 parallel dask tasks to generate 100 pairs of random numbers
    # as initital content of the bag
    get_Pmlt = lambda : np.random.lognormal(-0.5*Pstd**2, Pstd)
    get_Tadd = lambda : np.random.normal(0., Tstd )
    tasks = db.range(Nens, npartitions=4)
    def getRandom(i):
        return((str(i), get_Pmlt(), get_Tadd()))

    
    # process random number pair to argument string of FSM2
    tasks.map(getRandom).map(generateInputFile).map(runFSM2).compute()
    """
    nlst = open('nlst','w')
    nlst.write('&params \n') 
    nlst.write('/ \n') 
    nlst.close()
    os.system('./FSM2 < nlst')
    os.system('mv FSM2out.nc FSM2out_opn.nc')

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
    """
