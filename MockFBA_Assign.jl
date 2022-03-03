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
    #An option to split the calculation in NGC/SGC or use NS for combined
    group=ARGS[2] 
    tile_pass=parse(Int32,ARGS[3])
    #how many partition we want this pass to run with
    npartition=parse(Int32,ARGS[4])
    #which partition should I execute now
    mypart=parse(Int32,ARGS[5])

    #load the tiles
    tiles_dic=MockFBA.load_DESI_tiles(config["TILES"]["tile_file"],config["TILES"]["PROGRAM"],
                  tile_pass,group)

    #Number of tiles for the current setting
    ntile=size(tiles_dic["TILEID"],1)


    if(ntile==0)
        println("No tile to assign in pass=$(tile_pass)\n")
        return
    elseif(mypart>=ntile)
	println("Nothing to be done, have only ntile=$(ntile), but this cpu=$(mypart)")
	return
    end
    
    #The indices of tiles I will execute
    beg_index,end_index= MockFBA.partition_indices(ntile,npartition,mypart)

    #println("my_indices: $(beg_index) $(end_index) , total_tile: $(ntile)")

    this_index=collect(beg_index:end_index)
    #tile_date="2019-09-16T00:00:00"
    tile_date=config["focal_plane"]["date"]
    MockFBA.Run_Many_Tile(config,this_index,tiles_dic;group=group,tile_date=tile_date, plate_scale=nothing,fp_dic=nothing,exc_dic=nothing,verbose=2)
end

#@show ARGS
main(ARGS)

