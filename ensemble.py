"""
Run model ensemble
"""
import numpy as np
import os
import dask.bag as db
from dask.distributed import Client
import dask.delayed
import time

nCPUs = 4
Nens = 100    # ensemble size
Pstd = 0.63   # precipitation multiplier standard deviation
Tstd = 1.0    # temperature offset standard deviation

# function to transform int to tuple (label, Pmlt, Tadd)
def getRandom(i):
     return((str(i), np.random.lognormal(-0.5*Pstd**2, Pstd), np.random.normal(0., Tstd) ) )

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
    return filename

# run the simulation given file name
def runFSM2(filename):   
    print("Running ./FSM2 <", filename)
    os.system('./FSM2 < '+ filename)
    # return exit code ?


if __name__=="__main__":
    # just for info, print link to show computation
    client = Client()
    print(client)
    print("Click here to see dask dashboard: ", client.dashboard_link)

    # compile fortran code
    os.system('./compil.sh') 
    
    # create Bag containing ints 0..99
    tasks = db.range(Nens, npartitions=nCPUs)


    # process random number pair to argument string of FSM2
    tasks.map(getRandom).map(generateInputFile).map(runFSM2).compute()
