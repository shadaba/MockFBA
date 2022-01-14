#Author Shadab Alam, January 2022
#This manages and runs the inherent Julia package for the fibre-assignment on mocks

import os
from multiprocessing import Pool
#import yaml

def multi_run_wrapper(args):
   #return add(*args)
   return run_julia(*args)



def run_julia(tile_pass,ncpu,this_cpu):
    comm="julia --threads 1 MockFBA_Assign.jl config.yaml %d %d %d"%(tile_pass,ncpu,this_cpu)
    return os.system(comm)
    
    
if __name__ == "__main__":
    #from multiprocessing import Pool
    ncpu=6
    tile_pass=5

    pool = Pool(ncpu)

    part=[]
    for ii in range(0,ncpu):
        part.append((tile_pass,ncpu,ii))

    #results = pool.map(multi_run_wrapper,[(1,2),(2,3),(3,4)])
    results = pool.map(multi_run_wrapper,part)
    print results
