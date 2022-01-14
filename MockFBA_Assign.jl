#!/bin/bash
#=
exec julia --color=yes --startup-file=no "${BASH_SOURCE[0]}" "$@"
=#

using Pkg
Pkg.activate(".")
#include("srg/MockFBA.jl")
import MockFBA

using YAML

@show ARGS
#config=YAML.load_file("NoteBooks/test-output.yml")
config=YAML.load_file(ARGS[1])

tile_pass=parse(Int32,ARGS[2])
#how many partition we want this pass to run with
npartition=parse(Int32,ARGS[3])
#which partition should I execute now
mypart=parse(Int32,ARGS[4])

#load the tiles
tiles_dic=MockFBA.load_DESI_tiles(config["TILES"]["tile_file"],config["TILES"]["PROGRAM"],
    tile_pass,config["sky"])

#Number of tiles for the current setting
ntile=size(tiles_dic["TILEID"],1)
#How many should I execute
ntile_per_part= Int(ceil(ntile/npartition))

#The indices of tiles I will execute
beg_index=(mypart*ntile_per_part)+1
end_index=minimum([beg_index+ntile_per_part-1,ntile])









println("my_indices: $(beg_index) $(end_index) , total_tile: $(ntile)")



this_index=collect(beg_index:end_index)
tile_date="2019-09-16T00:00:00"
MockFBA.Run_Many_Tile(config,this_index,tiles_dic;tile_date=tile_date, plate_scale=nothing,fp_dic=nothing,exc_dic=nothing,verbose=2)

