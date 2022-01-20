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
    group =ARGS[3]

    #how many partition we want this pass to run with
    npart_zone=parse(Int32,ARGS[4])
    #which partition should I execute now
    mypart_zone=parse(Int32,ARGS[5])

    #process the zones
    MockFBA.PostFBA_JLD2FITS(config,tracer,group,mypart_zone,npart_zone,;columns=["RA","DEC","FITS_index_julia"])
end

@show ARGS
main(ARGS)

