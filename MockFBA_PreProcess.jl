#!/bin/bash
#=
exec julia --color=yes --startup-file=no "${BASH_SOURCE[0]}" "$@"
=#

using Pkg
Pkg.activate(".")
#include("srg/MockFBA.jl")
import MockFBA

using YAML


function main(ARGS)
    #load the config file 
    config=YAML.load_file(ARGS[1])
    tracer=ARGS[2]
    #1 to preprocess focalp plan
    focal_plane=parse(Int32,ARGS[3])

    
    #get the needed info from config
    fname=config["target"][tracer]["FITSfile"]
    priority_list=config["target"][tracer]["Priorities"]
    priority_frac_list=config["target"][tracer]["Priority_fraction"]
    outfile=config["target"][tracer]["JLDfile"]
    tile_file=config["TILES"]["tile_file"]
    tile_radius_indeg=config["TILES"]["tile_radius"]
    npass=config["TILES"]["NumPass"]
    program=config["TILES"]["PROGRAM"]

    #pre-process the tracer
    MockFBA.preprocess_tracer(fname,priority_list,priority_frac_list,outfile,
        tile_file,tile_radius_indeg,npass,program)

    #pre-process the focal plane
    if(focal_plane==1)
        fp_dic,exc_dic=MockFBA.load_hw_full_FocalPlane!(config["focal_plane"];date=config["focal_plane"]["date"])
    end


end

@show ARGS
main(ARGS)

