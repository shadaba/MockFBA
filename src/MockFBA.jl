module MockFBA

#using LoopVectorization
using FITSIO
using YAML
using CSV
using JLD2
using Formatting
using Printf
using DelimitedFiles
using DataInterpolations
using StaticArrays
using LoopVectorization
using DataFrames
using CodecLz4

export preprocess_tracer

include("FBA_utility.jl")
include("Harware.jl")
include("Tiles.jl")
include("Geometry.jl")
include("Positioner.jl")
include("ZonesTargets.jl")

#Should be load for development but not deployment
include("PlotLibrary.jl")

"""This goes through the full proe-process step for a particular tracer
 Steps:
  1) call the pre_process_FITS2JLD_tracer which splits the tracer in NGC/SGC and 20x20 deg zones
  2) Assign TILEID per_pass for the objects by load the tile file
  3) Generate indexing fro RA,DEC 
  4) Generate the listing for TILEID across all passes

   # Arguments
   - `fname::String`: Name of the fits file containing tracer objects, Only needs RA,DEC columns
   - `priority_list::Array[Integer]`: The list of priority to be assigned randomly
   - `priority_frac_list::Array[Float]`: Fraction of objects in each priority
   - `outfile::String`: Filename for JLD2 file to be created for efficient I/O
   - `tile_file::String`: fits file containing list tiles

"""
function preprocess_tracer(fname,priority_list,priority_frac_list,outfile,
        tile_file,tile_radius_indeg,npass,program)

    #generates the JLD2 file
    pre_process_FITS2JLD_tracer(fname,priority_list,priority_frac_list,outfile)
    
    #Assign TILEID to the object and ad to JLD2 file for each pass
    overwrite=false #To not overwrite the tileid field if it alredy exists
    #For SGC first and then NGC
    for group in ["NGC","SGC"]
        Append_TileID(outfile,tile_file,group,program,npass,tile_radius_indeg,overwrite)
        MockFBA.Append_column(outfile,group,["XYZ_UNIT"]) #Adding XYZ_UNIT
    end
    
    
    #Add the indexing to the file for RA and DEC and listing for tileid
    tileid_cols=[]
    for pass in 0:npass-1
        append!(tileid_cols,["TILEID_PASS$(pass)"])
    end
    
    for group in ["NGC","SGC",""]
        Indexing_column(outfile,group,["RA","DEC"],nothing)
        udic=LISTING_column(outfile,group,"TILEID",tileid_cols)
    end
end

end #end of module
