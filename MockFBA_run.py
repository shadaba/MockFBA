#Author Shadab Alam, January 2022
#This manages and runs the inherent Julia package for the fibre-assignment on mocks

import os
from multiprocessing import Pool
import yaml
from datetime import datetime
import glob 
import fitsio as F
import numpy as np

def utcnow():
    utc_time = datetime.utcnow()
    return utc_time.strftime('%Y-%m-%dT%H:%M:%S')

def combin_fits(file_list,outfile):
    '''combines the fits files assumin same columns in each one of them'''

    with F.FITS(file_list[0]) as fin:
        data_type=fin[1][0].dtype


    if(os.path.isfile(outfile)):
        print('File exists quiting: ',outfile)
        return

    fout=F.FITS(outfile,'rw')

    for ff,fname in enumerate(file_list):
        print('working on %d %s'%(ff,fname))

        with F.FITS(fname) as fin:
            nrow=fin[1][data_type.names[0]][:].size

            dataout=np.zeros(nrow,dtype=data_type)

            for cc,col in enumerate(data_type.names):
                dataout[col]=fin[1][col][:]

        if(ff==0):
            fout.write(dataout)
        else:
            fout[-1].append(dataout)

    fout[-1].write_checksum()

    fout.close()

    print('close fits: %s'%outfile)

    return


def run_julia_preprocess(config_file,tracer):
    comm="julia MockFBA_PreProcess.jl %s %s"%(config_file,tracer)
    return os.system(comm)

def wrapper_run_julia_preprcoess(args):
   return run_julia_preprocess(*args)


"""runs the pre-processing steps calling Julia
For each entry in the target runs this throuhg a different process"""
def run_pre_process(config,config_file,ncpu):
    print("%s Begin PreProcess"%(utcnow()))

    list_tracer=list(config["target"].keys())
    ntracer=len(list_tracer)

    pool = Pool(min(ntracer,ncpu))
    part=[]
    for ii in range(0,ntracer):
        tracer=list_tracer[ii]
        if(os.path.isfile(config["target"][tracer]["JLDfile"])):
            print("File Exists, No preprocessing for %s\n using: %s"%(tracer,config["target"][tracer]["JLDfile"]))
            continue
        part.append((config_file,list_tracer[ii]))
    results = pool.map(wrapper_run_julia_preprcoess,part)
    
    if(sum(results)!=0):
        print("PreProcess step failed with error")
        exit(1)

    print("%s End PreProcess"%(utcnow()))
    return 


"""
To run the Julia code for fiber-assignment
"""
def run_julia_assign(config_file,group,tile_pass,ncpu,this_cpu):
    comm="julia --threads 1 MockFBA_Assign.jl %s %s %d %d %d"%(config_file,group,tile_pass,ncpu,this_cpu)
    print(comm)
    return os.system(comm)
    
def wrapper_run_julia_assign(args):
   return run_julia_assign(*args)

    
"""Run each of the passes in parallel one by one
"""
def run_FBA_passes(config_file,group,npass,ncpu):

    for tile_pass in range(0,npass):
        print("Begin pass=%d"%(tile_pass))
        pool = Pool(ncpu)
        part=[]
        for ii in range(0,ncpu):
            part.append((config_file,group,tile_pass,ncpu,ii))
        results = pool.map(wrapper_run_julia_assign,part)
    
        if(sum(results)!=0):
            print("EXITING: FBA assignment failed for pass=%d"%(tile_pass))
            exit(1)

        print("End pass=%d"%(tile_pass))

    return


"""Runs the post-process steps
This involves two steps done for each tracer in the target
1)Julia in parallel for each zone: Combine the tiles FBA output with zone JLD2 file and write FITS for each zone
2)Python: Combines all the zone fits to single fits file for each tracer and group
"""
def run_post_process(config,config_file,ncpu,group):
    print("%s Begin PostProcess"%(utcnow()))

    list_tracer=list(config["target"].keys())
    ntracer=len(list_tracer)

    pool = Pool(ncpu)
    part=[]
    for ii in range(0,ntracer):
        tracer=list_tracer[ii]
        if(tracer!="ELG"):
            npart=1
        else:
            if(ncpu<=ntracer+2):
                npart=ncpu
            else:
                npart=ncpu-ntracer+1

        for mypart in range(0,npart):
            part.append((config_file,list_tracer[ii],group,npart,mypart))

    results = pool.map(wrapper_run_julia_postprcoess,part)

    if(sum(results)!=0):
        print("PostProcess step failed with error")
        exit(1)

    #finally combine them to one gits file
    pool = Pool(ntracer)
    part=[]
    for ii in range(0,ntracer):
        tracer=list_tracer[ii]
        part.append((config,group,list_tracer[ii]))

    results = pool.map(wrapper_combine_group_to_singlefits,part)

    print("%s End PostProcess"%(utcnow()))
    return

def run_julia_postprcoess(config_file,tracer,group,npart,mypart):
    comm="julia MockFBA_PostProcess.jl %s %s %s %d %d"%(config_file,tracer,group,npart,mypart)
    return os.system(comm)
    
def wrapper_run_julia_postprcoess(args):
    return run_julia_postprcoess(*args)

#combine all the zone_fits in single fits file
def combine_group_to_singlefits(config,group,tracer):
    
    #The directory to write the fits files
    outdir="%s%s/"%(config["OUTPUT"]["FITS_dir"],config["OUTPUT"]["FBA-tag"])
    #file names of the zone based fits file
    outfile_str="%s/tmp/%s_%s_zone_*.fits.gz"%(outdir,group,tracer)

    #print(outfile_str)
    zone_fits=glob.glob(outfile_str)
    outfile="%s%s_%s.fits.gz"%(outdir,tracer,group)
    combin_fits(zone_fits,outfile)
    print("written: ",outfile)

    return

def wrapper_combine_group_to_singlefits(args):
    return combine_group_to_singlefits(*args)

def proces_FBA(config,config_file,ncpu,step="assignment"):
    #Runs the actual fiber-assignment steps per tile
    if(config["SplitSky"]>0):
        if(config["SplitSky"]==1 or config["SplitSky"]==3):
            print("%s Begin %s for NGC"%(step,utcnow()))
            if(step=="assignment"):
                run_FBA_passes(config_file,"NGC",config["TILES"]["NumPass"],ncpu)
            else:
                run_post_process(config,config_file,ncpu,"NGC")

        if(config["SplitSky"]==2 or config["SplitSky"]==3):
            print("%s Begin %s for SGC"%(step,utcnow()))
            if(step=="assignment"):
                run_FBA_passes(config_file,"SGC",config["TILES"]["NumPass"],ncpu)
            else:
                run_post_process(config,config_file,ncpu,"SGC")
    else:
        print("%s Begin %s for NS-combined"%(step,utcnow()))
        if(step=="assignment"):
            run_FBA_passes(config_file,"NS",config["TILES"]["NumPass"],ncpu)
        else:
            run_post_process(config,config_file,ncpu,"NS")


    return


if __name__ == "__main__":
    #from multiprocessing import Pool
    #config_file="config.yaml"
    config_file="config_tmp.yaml"
    with open(config_file, 'r') as file:
        config=yaml.safe_load(file)

    ncpu=config["ncpu"]

    if("pre-process" in config["steps"]):
        #Runs the actual fiber-assignment steps per tile
        run_pre_process(config,config_file,ncpu)
    

    if("assignment" in config["steps"]):
        proces_FBA(config,config_file,ncpu,step="assignment")


    if("post-process" in config["steps"]):
        #post-process the assignment, this converts the tiles files to zone_fits files
        proces_FBA(config,config_file,ncpu,step="post-process")

    #completion message
    print("%s All Processing Finished"%(utcnow()))


