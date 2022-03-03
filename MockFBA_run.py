#Author Shadab Alam, January 2022
#This manages and runs the inherent Julia package for the fibre-assignment on mocks

import os
import yaml
import glob 
import fitsio as F
import numpy as np
import warnings
import subprocess
from multiprocessing import Pool
from datetime import datetime

import argparse


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='MockFBA: Runs the fibreassignment in mocks',
            formatter_class=argparse.RawTextHelpFormatter)
    #these should be input
    parser.add_argument('-config_file',default=None,help='Input can be provided using a config file.')
    parser.add_argument('-steps' ,nargs='+',type=str,default=["pre-process","assignment","post-process"], 
                   help='''The process have three steps and all of them can be run in one go or multiple call
You can provide a list of values from below 
   pre-process: This is the step when the fits file is converted to the efficient JLD2 format
   assignment: runs the actual assignment for all the tiles
   post-process: post-process all the tiles and generates easy to use fits file''')
    parser.add_argument('-ncpu',type=int, default=4,help='Number of cpu can be used to parallelize the process')

    args = parser.parse_args()


def julia_executable(step="pre-process",verbose=False):
    if("MOCKFBA_PATH" in os.environ.keys()):
        julia_root=os.environ["MOCKFBA_PATH"]
    else:
        julia_root=""
        if(verbose):
            warnings.warn("MOCKFBA_PATH is not set")

    if(step=="pre-process"):
        jl_exe="MockFBA_PreProcess.jl"
    elif(step=="post-process"):
        jl_exe="MockFBA_PostProcess.jl"
    elif(step=="FBA"):
        jl_exe="MockFBA_Assign.jl"

    return julia_root+jl_exe

def utcnow():
    utc_time = datetime.utcnow()
    return utc_time.strftime('%Y-%m-%dT%H:%M:%S')

def combin_fits(file_list,outfile,original_fits=None,Original_columns=None,index_column="FITS_index_julia",index_one=True):
    '''combines the fits files assumin same columns in each one of them
    additional pulls the additional columns from original_fits and add them based on
    index_col, assuming indexing from 1 you can set this to zero convetion by index_one=False
    '''

    #for the original fits file
    forg=None

    with F.FITS(file_list[0],"r") as fin:
        data_type=fin[1][0].dtype

    if(Original_columns is not None):
        if(not os.path.isfile(original_fits)):
            miss_col=""
            for col in Original_columns:
                if(col in data_type.names):
                    continue
                miss_col="%s %s,"%(miss_col,col)

            tmsg="Warning: original fits not found: %s"%(original_fits)
            tmsg=tmsg+"\n   Following original columns can't be added: %s"%(miss_col)
            tmsg=tmsg+"\n   to file %s"%(outfile)
            print(tmsg)
            out_data_type=data_type
        else:
            forg=F.FITS(original_fits,"r")
            org_dtype=forg[1][0].dtype

            new_list=[]
            for col in data_type.names:
                new_list.append((col,data_type[col]))

            for cc,col in enumerate(org_dtype.names):
                if(col in data_type.names):
                    continue
                elif(col in Original_columns):
                    #print(cc,col,org_dtype[col])
                    new_list.append((col,org_dtype[col]))

            #debug
            #if(True):
            #    new_list.append(("RA_org",">f4"))
            #    new_list.append(("DEC_org",">f4"))

            out_data_type=np.dtype(new_list)
    else:
        out_data_type=data_type


    #Now combine the FBA_BITS in one variable
    if(True):
        new_list=[]
        nFBA_BITS=0
        for col in data_type.names:
            if('FBA_BITS'== col[:8]):
                nFBA_BITS=nFBA_BITS+1
            else:
                new_list.append((col,data_type[col]))

        if(nFBA_BITS==1):
            new_list.append(('FBA_BITS',data_type['FBA_BITS1']))
        else:
            new_list.append(('FBA_BITS',data_type['FBA_BITS1'],(nFBA_BITS,)))
        out_data_type=np.dtype(new_list)




    if(os.path.isfile(outfile)):
        print('File exists quiting: ',outfile)
        return


    fout=F.FITS(outfile,'rw')


    for ff,fname in enumerate(file_list):
        print('Combining %s'%(fname))

        with F.FITS(fname) as fin:
            nrow=fin[1][data_type.names[0]][:].size

            dataout=np.zeros(nrow,dtype=out_data_type)

            for cc,col in enumerate(out_data_type.names):
                if('FBA_BITS'== col):
                    if(nFBA_BITS==1):
                        dataout['FBA_BITS']=fin[1]['FBA_BITS1'][:]
                    else:
                        #transfering 1d arrays to 2d arrays for FBA_BITS
                        for ifba in range(0,nFBA_BITS):
                            dataout['FBA_BITS'][:,ifba]=fin[1]['FBA_BITS%d'%(ifba+1)][:]

                elif col in data_type.names:
                    dataout[col]=fin[1][col][:]
                else:
                    index=fin[1][index_column][:]
                    #if(col in ["RA_org","DEC_org"]):
                    #    dataout[col]=forg[1][col[:-4]][index-1]
                    if(index_one):
                        dataout[col]=forg[1][col][index-1]
                    else:
                        dataout[col]=forg[1][col][index]

        if(ff==0):
            fout.write(dataout)
        else:
            fout[-1].append(dataout)

    fout[-1].write_checksum()

    fout.close()

    if(forg is not None):
        forg.close()

    print('%s Written Combined fits: %s'%(utcnow(),outfile))

    return



def run_julia_preprocess(config_file,tracer,focal_plane_preprocess):
    #comm="julia %s %s %s %d"%(julia_executable(step="pre-process"),config_file,tracer,focal_plane_preprocess)
    #return os.system(comm)
    comm=["julia",julia_executable(step="pre-process"),config_file,tracer,'%d'%focal_plane_preprocess]
    return subprocess.run(comm)

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

        #This is to also preprocess the focal plane file should be done only once
        if(ii==0):
            focal_plane_process=1
        else:
            focal_plane_process=0

        part.append((config_file,list_tracer[ii],focal_plane_process))
    results = pool.map(wrapper_run_julia_preprcoess,part)
   
    for res in results:
        if(res.returncode!=0):
            print("PreProcess step failed with error: ",res)
            exit(1)

    print("%s End PreProcess"%(utcnow()))
    return 


"""
To run the Julia code for fiber-assignment
"""
def run_julia_assign(config_file,group,tile_pass,ncpu,this_cpu):
    #comm="julia --threads 1 %s %s %s %d %d %d"%(julia_executable(step="FBA"),config_file,group,tile_pass,ncpu,this_cpu)
    #print(comm)
    #return os.system(comm)
    comm=["julia",julia_executable(step="FBA"),config_file,group,'%s'%tile_pass,'%d'%ncpu,'%d'%this_cpu]
    return subprocess.run(comm)

def wrapper_run_julia_assign(args):
   return run_julia_assign(*args)

    
"""Run each of the passes in parallel one by one
"""
def run_FBA_passes(config_file,group,npass,ncpu):

    for tile_pass in range(0,npass):
        print("%s Begin pass=%d"%(utcnow(),tile_pass))
        pool = Pool(ncpu)
        part=[]
        for ii in range(0,ncpu):
            part.append((config_file,group,tile_pass,ncpu,ii))
        results = pool.map(wrapper_run_julia_assign,part)
    
        for res in results:
            if(res.returncode!=0):
                print("EXITING: FBA assignment failed for pass=%d"%(tile_pass),res)
                exit(1)

        print("%s End pass=%d"%(utcnow(),tile_pass))

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

    for res in results:
        if(res.returncode!=0):
            print("PostProcess step failed with error: ",res)
            exit(1)

    #finally combine them to one fits file
    pool = Pool(ntracer)
    part=[]
    for ii in range(0,ntracer):
        tracer=list_tracer[ii]
        part.append((config,group,list_tracer[ii]))

    results = pool.map(wrapper_combine_group_to_singlefits,part)

    print("%s End PostProcess"%(utcnow()))
    return

def run_julia_postprcoess(config_file,tracer,group,npart,mypart):
    #comm="julia %s %s %s %s %d %d"%(julia_executable(step="post-process"),config_file,tracer,group,npart,mypart)
    #return os.system(comm)
    comm=["julia",julia_executable(step="post-process"),config_file,tracer,group,'%d'%npart,'%d'%mypart]
    return subprocess.run(comm)

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
    #combin_fits(zone_fits,outfile,config)
    combin_fits(zone_fits,outfile,original_fits=config["target"][tracer]["FITSfile"],
            Original_columns=config["target"][tracer]["carry_columns"],
            index_column="FITS_index_julia",index_one=True)
    #print("written: ",outfile)

    return

def wrapper_combine_group_to_singlefits(args):
    return combine_group_to_singlefits(*args)

def proces_FBA(config,config_file,ncpu,step="assignment"):
    #Runs the actual fiber-assignment steps per tile
    if(config["SplitSky"]>0):
        if(config["SplitSky"]==1 or config["SplitSky"]==3):
            print("%s Begin %s for NGC"%(utcnow(),step))
            if(step=="assignment"):
                run_FBA_passes(config_file,"NGC",config["TILES"]["NumPass"],ncpu)
            else:
                run_post_process(config,config_file,ncpu,"NGC")

        if(config["SplitSky"]==2 or config["SplitSky"]==3):
            print("%s Begin %s for SGC"%(utcnow(),step))
            if(step=="assignment"):
                run_FBA_passes(config_file,"SGC",config["TILES"]["NumPass"],ncpu)
            else:
                run_post_process(config,config_file,ncpu,"SGC")
    else:
        print("%s Begin %s for NS-combined"%(utcnow(),step))
        if(step=="assignment"):
            run_FBA_passes(config_file,"NS",config["TILES"]["NumPass"],ncpu)
        else:
            run_post_process(config,config_file,ncpu,"NS")


    return


"""validates the config file and write an extended config file for julia
"""
def validate_config(config_file):
    #check if MOCKFBA_PATH is setup
    for ss,step in enumerate(["pre-process","post-process","FBA"]):
        if(ss==0):
            verbose=True
        else:
            verbose=False
        
        julia_jl=julia_executable(step="pre-process",verbose=verbose)
        if( not os.path.isfile(julia_jl)):
            err_msg="Cannot find julia executeable on following path:\n  %s\n"%(julia_jl)
            err_msg=err_msg+"Either MOCKFBA_PATH is not set properly or developer not in the MockFBA directory\n"
            raise ValueError(err_msg)

    #load the yaml config file
    with open(config_file, 'r') as file:
        config=yaml.safe_load(file)

    config_file_julia="%s_julia.yaml"%(config_file[:-5])
    
    config=check_config_file(config)

    #write the updte config file for julia
    with open(config_file_julia,"w") as fout:
        yaml.dump(config, fout,sort_keys=True)

    return config,config_file_julia

def make_dir(dirin):
    '''checks and creates a directory'''
    if(not os.path.isdir(dirin)):
        os.mkdir(dirin)
        print('created: ',dirin)
    else:
        print('Already exists: ',dirin)
        
    return

def check_config_file(config):
    #create output directory
    make_dir(config["OUTPUT"]["OUT_dir"])
    for tdir in ["JLD2_dir","FITS_dir"]:
        if(tdir in config["OUTPUT"].keys()):
            this_dir=config["OUTPUT"][tdir]
        else:
            this_dir="%s%s-%s/"%(config["OUTPUT"]["OUT_dir"],tdir,config["OUTPUT"]["pre-process-tag"])
            config["OUTPUT"][tdir]=this_dir
        make_dir(this_dir)

        #To generate the directory for writing tiles as output
        if(tdir=="JLD2_dir"):
            tile_dir="%s%s"%( config["OUTPUT"][tdir],config["OUTPUT"]["FBA-tag"])
            make_dir(tile_dir)
        elif(tdir=="FITS_dir"):#Directory to write the fits file
            fits_dir="%s%s"%(config["OUTPUT"][tdir],config["OUTPUT"]["FBA-tag"])
            make_dir(fits_dir)
            #A temp directory to write zone fits file
            make_dir(fits_dir+"/tmp")


    #check for the focal_plane directory
    make_dir(config["focal_plane"]["focalplane_dir_jld2"])

    #check if focalplane directory exists
    if(not os.path.isdir(config["focal_plane"]["focalplane_dir"])):
        print("Warning the foalplane firctory not found: %s"%(config["focal_plane"]["focalplane_dir"]))

    #check if tile file exists
    if(not os.path.isfile(config["TILES"]["tile_file"])):
        print("Error: Tile file not found: %s"%(config["TILES"]["tile_file"]))
        raise
    else: #make sure the tile id is greate than zero
        ftile=config["TILES"]["tile_file"]
        if(np.min(F.FITS(ftile,"r")[1]["TILEID"][:])<1):
            print("Error: The TILEID in %s is less than 1, TILEID must be >=1"%(ftile))
            raise 

    #Now check if all the target files exists and add JLD file columns
    for tracer in list(config["target"].keys()):
        fitsfile=config["target"][tracer]["FITSfile"]
        if(not os.path.isfile(fitsfile)):
            print("Erorr: The %s FitsFile %s not found"%(tracer,fitsfile))
            raise
        elif("JLDfile" not in config["target"][tracer].keys()):
            fname=fitsfile.split("/")[-1]
            jldfile="%s%s.jld2"%(config["OUTPUT"]["JLD2_dir"],fname[:-5])
            config["target"][tracer]["JLDfile"]=jldfile

    return config



if __name__ == "__main__":
    #from multiprocessing import Pool
    config_file_in=args.config_file

    #config_file="config_tmp.yaml"
    
    config,config_file=validate_config(config_file_in)

    ncpu=args.ncpu

    if("pre-process" in args.steps):
        #Runs the actual fiber-assignment steps per tile
        run_pre_process(config,config_file,ncpu)
    

    if("assignment" in args.steps):
        proces_FBA(config,config_file,ncpu,step="assignment")


    if("post-process" in args.steps):
        #post-process the assignment, this converts the tiles files to zone_fits files
        proces_FBA(config,config_file,ncpu,step="post-process")

    #completion message
    print("%s All Processing Finished"%(utcnow()))


