
"""
Loads all the tiles as per config file and generate a tileid to pass map
"""
function map_TILEID2PASS(config,group)
    #load all the files load tiles from all the pass
    if(group=="NS" || group=="")
       sky=nothing
    else
       sky=group
    end

    tiles_dic=load_DESI_tiles(config["TILES"]["tile_file"],config["TILES"]["PROGRAM"],
                  nothing,sky)
    tiles_pass=Dict{Integer,Integer}()
    
    for ii in 1:size(tiles_dic["TILEID"],1)
        tiles_pass[tiles_dic["TILEID"][ii]]=tiles_dic["PASS"][ii]
    end
        
    return tiles_pass
end

"""For each zone and each and given tracer writes the fits file including FBA information
config: input config dictionary loaded from config file
tracer: One of the tracer in the targets categories
group: The overall group we are working with this should be consistent among all the steps
mypart: which partition of the zones needs to be executed
npart_zone: How many partition we want for zones, ech partition we work in parallell cpu
columns: The list of columns tow write apart from FBA related output read from tile files
nbit_pack: How many bits to pack together, the code works for 64 or 32 but FITSIO library can-not write
UInt64 and hence only 32 will work at the moment.
"""
function PostFBA_JLD2FITS(config,tracer,group,mypart_zone,npart_zone,;
        columns=["RA","DEC","FITS_index_julia"],nbit_pack=32)
    #= for the group in config file
    scan alls the zones and write them to fits including the FBA info from tile file
    =#
    
    @assert ((nbit_pack==32) || (nbit_pack==64)) "Invalid bitpacking size $(nbit_pack), must be 32 or 64"

    jldfile=config["target"][tracer]["JLDfile"]
    num_fba=config["NumFBARealization"]
    
    
    #file=jldopen(jldfile,"r")
    #@show keys(file[group]["zone_77"])
    #close(file)
    #return
    #load the zones in group
    jld_order=traverse_group(jldfile,group)
    
    nzone=size(jld_order["zone_order"],1)
    izone_min,izone_max=partition_indices(nzone,npart_zone,mypart_zone)
    
    #load the listing to identify the zones which has tile
    listing=load(jldfile,"$(group)/listing")

    #maximum number of objects among the zones of interest
    max_nobj=maximum(jld_order["nobj"][izone_min:izone_max])
    
    #declare array to store info to be written
    tracer_targets=Dict{String,Array}("RA"=>zeros(Float32,max_nobj),"DEC"=>zeros(Float32,max_nobj))
    for col in columns
        if(col in ["FITS_index_julia","TILEID"])
            tracer_targets[col]=zeros(Int64,max_nobj)
        else
            tracer_targets[col]=zeros(Float32,max_nobj)
        end
    end

    #FBA related quantities
    tracer_targets["FBA_PASS_R1"]=zeros(Int8,max_nobj) .-1
    tracer_targets["FBA_LOCATION_R1"]=zeros(Int16,max_nobj)
    
    #additionally we need to load the group index to match the object in FBA
    tracer_targets["index_group"]=zeros(Int64,max_nobj)
    #To keep a count of how many times and object is assigned for R1, THis is for sanity checks
    tracer_targets["Count_R1"]=zeros(Int32,max_nobj)
    
    #number of bitweights are needed given fba and packing rate
    nbitweight=Int(floor(num_fba/nbit_pack))
    if(nbitweight*nbit_pack<num_fba)
        println("PostP Warning: Want to pack $(num_fba) bits in $(nbit_pack), will use only first $(nbitweight*nbit_pack)")
        println("PostP Warning: Ignoring the last $(num_fba-(nbitweight*nbit_pack)) realizations")
    end
    
    for ii in 1:nbitweight
        if(nbit_pack==32)
            tracer_targets["FBA_BITS$(ii)"]=zeros(UInt32,max_nobj)
        else
            tracer_targets["FBA_BITS$(ii)"]=zeros(UInt64,max_nobj)
        end
    end
    
    #generate and tileif 2d passs maps
    tiles_pass=map_TILEID2PASS(config,group)
    
    @show izone_min,izone_max, max_nobj
    
    for iz in izone_min:izone_max
        fill!(tracer_targets["FBA_PASS_R1"],-1)
        fill!(tracer_targets["FBA_LOCATION_R1"],-1)
        
        #check the number of tiles in the zone 
        #if it is one with tileid=0 mean none of the tiles in zone and hence skip
        if(listing["TILEID"][iz]==[0])
            println("POSTP Warning: No tile in zone $(iz), skipping")
            continue
        end
        
        
        #Number of objects in this zone
        nobj_zone=jld_order["nobj"][iz]
        
        #poiting to the group in the file for this zone
        tg=jld_order["zone_order"][iz]

        zone_imin=jld_order["beg_index"][iz]
        zone_imax=zone_imin+jld_order["nobj"][iz]-1
        zone_nobj=jld_order["nobj"][iz]
        
        
        #indices of objects in the zone
        tracer_targets["index_group"][1:nobj_zone] .= collect(zone_imin:zone_imax)
        
        
        #collect different columns directly available in JLD2 zones file
        for col in columns
            #get the columns based on the need
            #tval=load(jldfile,"$(tg)$(col)")
            tracer_targets[col][1:nobj_zone] .= load(jldfile,"$(tg)$(col)")
            
            #if(col in keys(tracer_targets))
            #    tracer_targets[col]=[tracer_targets[col];tval]
            #else
            #    tracer_targets[col][1:nobj_zone]=tval
            #end
        end
        
        #collect the columns from FBA procedure by reading the tile files
        fill_FBA_results!(config,tiles_pass,tracer,listing["TILEID"][iz],tracer_targets,nobj_zone,
            num_fba,nbitweight,nbit_pack) 
        
        
        
        #columns to write
        col_write=[]
        array_write=[]
        tmp_dic=Dict{String,Array}()
        for col in keys(tracer_targets)
            if(col in ["index_group","Count_R1"]) #columns to skip
                continue
            end
            #append!(col_write,[col])
            #append!(array_write,view(tracer_targets[col],1:nobj_zone))
            #println(col,nobj_zone,size(tracer_targets[col]))
            tmp_dic[col]=view(tracer_targets[col],1:nobj_zone)
        end
        
        #fits file to write this zone 
        outfits=ZoneFITS_FileName(config,iz,tracer,group)
        FITS(outfits,"w") do fout
            #write(fout,tracer_targets)
            #write(fout,col_write,array_write)
            write(fout,tmp_dic) #Here we write to the fits file
        end
        println(now(UTC)," written: $(nobj_zone)",outfits)
        #break

    end #end of the iz loop
end
    
    
function fill_FBA_results!(config,tiles_pass,tracer,tile_list,tracer_targets,nobj_zone,
            num_fba,nbitweight,nbit_pack)    
    #index to sort the objects in tile by index field
    zone_isort=sortperm(tracer_targets["index_group"][1:nobj_zone])
    
    #To account for how many times an object is assigned for realization 1
    fill!(tracer_targets["Count_R1"],0)

    for tile_id in tile_list
        if(tile_id==0) #Not an actuall tile used to fill values
            continue
        end
        #FBA jldfile for the tile
        fba_jldfile=TileFBA_FileName(config,tile_id)

        index_pre_tile=load(fba_jldfile,"$(tracer)_index")
        bool_ass=load(fba_jldfile,"$(tracer)_Bool_ass")
        location_ass=load(fba_jldfile,"$(tracer)_location_ass")
        
        #make sure the tile was assigned with consistent num_fba
        @assert num_fba==size(bool_ass,2) "inconsistent num_fba \nfor tile: $(col_tile) in $(fba_jldfile)"

        #find the sorting permutation
        isort_perm=sortperm(index_pre_tile)

        #rearrange with sorted index for this tile
        index_pre_tile=index_pre_tile[isort_perm]
        bool_ass=bool_ass[isort_perm,:]
        location_ass=location_ass[isort_perm,:]
        
        #find the pass for this tile
        #@show tile_id, tiles_pass[tile_id]
        tile_pass=tiles_pass[tile_id]

        #append the info from FBA #follow the logic of sorted indices
        ctile=1
        for ti in 1:size(index_pre_tile,1)
            while(index_pre_tile[ti]>tracer_targets["index_group"][zone_isort[ctile]])
                ctile +=1
                if(ctile>nobj_zone)
                    break
                end
            end

            if(ctile>nobj_zone)
                break
            end

            if(index_pre_tile[ti]==tracer_targets["index_group"][zone_isort[ctile]])
                ind_zone=zone_isort[ctile]
                
                #Now fill the needed information for first realization
                if(bool_ass[ti,1])
                    #track the counts for first realization
                    #println(typeof(ind_zone),typeof(tracer_targets["index_group"]))
                    tracer_targets["Count_R1"][ind_zone] +=1
                    tracer_targets["FBA_PASS_R1"][ind_zone]= tile_pass
                    tracer_targets["FBA_LOCATION_R1"][ind_zone]=location_ass[ti,1]
                end
                #Now fill the FBA_BITS
                PackBITS!(nbitweight,nbit_pack,bool_ass[ti,:],tracer_targets,ind_zone)
            end
        end
    
    end
    #Sanity checks
    if(maximum(tracer_targets["Count_R1"]) > config["target"][tracer]["Num_obs"])
        nmore=size(findall(tracer_targets["Count_R1"] .> config["target"][tracer]["Num_obs"]),1)
	msg_warn="$(tracer): $(nmore) objects have more than $(config["target"][tracer]["Num_obs"]) assignment for realization 1"
	println("PostP Warning: $(msg_warn)")
    end
    println("counts (min,max):",minimum(tracer_targets["Count_R1"]),' ',maximum(tracer_targets["Count_R1"]))
end

"""Converts and array of boolean to UInt64 representation
example: value_out =UInt64(0)
bool_in=ones(Bool,32)
bool_in[1]=false #[false,false,false]#tar[1,1:4]
#@show value_out
value_out=Bool2UInt(bool_in,value_out)
@show bool_in
@show value_out,typeof(value_out)
@show bitstring(value_out)
"""
function Bool2UInt(bool_in,value_out::UInt64)
    #@show value_out
    value_out=UInt64(0)
    for ii in 1:size(bool_in,1)
        if(bool_in[ii])
            value_out += UInt64(2)^(ii-1)
            #println(ii,' ',value_out,' ',typeof(value_out))
        end
    end
    return value_out
end

"""Converts and array of boolean to UInt32 representation
example:
value_out =UInt32(0)
bool_in=ones(Bool,32)
bool_in[1]=false #[false,false,false]#tar[1,1:4]
#@show value_out
value_out=Bool2UInt(bool_in,value_out)
@show bool_in
@show value_out,typeof(value_out)
@show bitstring(value_out)
"""
function Bool2UInt(bool_in,value_out::UInt32)
    value_out=UInt32(0)
    #@show value_out
    for ii in 1:size(bool_in,1)
        if(bool_in[ii])
            value_out += UInt32(2)^(ii-1)
            #println(ii,' ',value_out,' ',typeof(value_out))
        end
    end
    return value_out
end

"""
Transfers the bollean array to array of UInt32 or UInt64
nbitweight: Number of bitweights columns
nbit_pack: 32 or 64 the unitsize of packing
bool_in: Input boolean array of atleast nbitweight*nbit_pack size
pack_out: output array of size nbitweight and type UInt32 or UInt64 based on nbit_pack
"""
function PackBITS!(nbitweight,nbit_pack,bool_in,tracer_targets,ind_zone)
    for ii in 1:nbitweight
        imin=((ii-1)*nbit_pack)+1
        imax=imin+nbit_pack-1
        tracer_targets["FBA_BITS$(ii)"][ind_zone]=Bool2UInt(bool_in[imin:imax],
            tracer_targets["FBA_BITS$(ii)"][ind_zone])
        #println(ii,' ',bool_in[imin:imax])
        #println(bitstring(pack_out[ii]))
    end
end

