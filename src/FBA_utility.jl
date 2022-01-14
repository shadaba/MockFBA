function costheta_uv(u,varr)
    #=Evaluates the cosine of theta between a fixe vector u and and array of vector v
    Returns the dot product of the u with each vector in v
    =#
    a = fill(zero(eltype(varr)),size(varr,1) )
    for k in 1:size(varr, 2)
        for j in 1:size(varr, 1)
            a[j] += u[k] * varr[j, k]
        end
    end
    return a
end

function norm(varr)
    a = fill(zero(eltype(varr)),size(varr,1) )
    for k in 1:size(varr, 2)
        for j in 1:size(varr, 1)
            a[j] += varr[j,k] * varr[j, k]
        end
    end
    return a
end

function assign_sky_zones(ra,dec,dRAzone,dDECzone)
    #assigns the zone for each object
    RA_ed=collect(0:dRAzone:360)
    DEC_ed=collect(-90:dDECzone:90)
    nRA=size(RA_ed,1)
    nDEC=size(DEC_ed,1)
    
    zone_num=zeros(Int32,size(ra,1))
    
    for ii in 1:nRA-1
        ibin=findall( (ra.> RA_ed[ii]) .& (ra.<=RA_ed[ii+1]))
        zone_num[ibin].=ii
    end
        
    for jj in 1:nDEC-1
        ibin=findall((dec.> DEC_ed[jj]) .& (dec.<=DEC_ed[jj+1]))
        zone_num[ibin].= zone_num[ibin] .+ (nRA*jj)
    end
        
    
    return zone_num
end

function assign_priority(nobj,priority,priority_frac)
    #= Generates uniform random numbers
    split from priority fraction
    assign priority
    =#
    
    @assert sum(priority_frac)==1
    
    obj_priority=zeros(Int32,nobj)
    rand_num=rand(nobj)
    #println(minimum(rand_num),' ',maximum(rand_num))
    n_class=size(priority,1)
    
    
    max_lim=0;
    #count_ass=zeros(Int32,nobj)
    for ii in 1:n_class
        min_lim=max_lim
        max_lim=min_lim+priority_frac[ii]
        ind_this=findall( (rand_num.>= min_lim) .& (rand_num.<=max_lim))
        obj_priority[ind_this] .= priority[ii]
        #count_ass[ind_this] .= count_ass[ind_this].+1
    end    
    return obj_priority
end

function sort_by_zone_priority(indsel,zone_num,obj_priority)
    #= first sort by zone and then sort by priority for each zone
    This construct a new variable to sort finally
     first sort by main variable, find minimum difference
    take second variable and set in the minimum difference
    add the scale second variable to first one and sort by this new quantity
    =#
    #sort by zone
    izone_sort=sortperm(zone_num)
    mindif=zone_num[izone_sort[2]]-zone_num[izone_sort[1]]
    
    #scale the second variable
    pr_min=minimum(obj_priority)
    pr_max=maximum(obj_priority)
    
    #new quantity for sorting
    new_q=(obj_priority .- pr_min)./(pr_max-pr_min)
    new_q=new_q+zone_num
    izone_sort=sortperm(new_q,rev=:true)
    return izone_sort
end

function approx_convert_ra_dec_to_gal_b(ra,dec)
   #=Convert ra,dec in degree to galactic latitude only b
    This is only approximate and can have error of about 0.01 degree
    =#
    deg2rad=pi/180
    alpha_NGP=192.85*deg2rad #ra in deg
    delta_NGP=27.14*deg2rad  #dec in deg
    l_NCP=6.28308866829052 # in radian
    
    sinb=sin(delta_NGP)*sin(dec*deg2rad)+ cos(delta_NGP)*cos(dec*deg2rad)*cos(ra*deg2rad-alpha_NGP)
    #cosb=sqrt(1-sinb*sinb)
    #cos_sin_lb= cos(dec)*sin(ra*deg2rad-alpha_NGP)
    
    #println(sinb,' ',cosb,' ',cos_sin_lb,' ',cos_sin_lb/cosb)
    
    #get l longitude and b latitude
    b=asin(sinb)/deg2rad
    #l=l_NCP-asin(cos_sin_lb/cosb)
    return b #in degree
end


function split_Galactic_cap(ra,dec,rettype)
    #return the indices of objects in NGC or SGC
    b= @. approx_convert_ra_dec_to_gal_b(ra,dec)
    
    if(rettype=="boolean")
        ind_NGC=b.>0
        ind_SGC=b.<=0
    else
        ind_NGC=findall(b.>0)
        ind_SGC=findall(b.<=0)
    end
    #ind_SGC=deleteat!([1:size(ra);], ind_NGC)
    
    
    return ind_NGC,ind_SGC
end


function write_JLD2_file(outfile,zone_num,quant_dic)
    #=Write the structured data in jld2 format with subgroup feature
     For any tracer write following structure
      -> NGC/SGC -> groupd in ra-dec zones
    =#
    us_zone=sort(unique(zone_num))
    #append!(size(zone_num,1),fel_zone)
    #@show(fel_zone,size(fel_zone),size(us_zone))
    
    quant_write=["RA","DEC","PRIORITY","FITS_index"]

    ind_NS=Dict()
    ind_NS["NGC"],ind_NS["SGC"]= split_Galactic_cap(quant_dic["RA"],quant_dic["DEC"],"boolean")
    
    groups_order=Dict()
    groups_order["NGC"]=Dict()
    groups_order["SGC"]=Dict()
    
    
    jldopen(outfile, "a+") do file
        for zz in 1:size(us_zone,1)
            zone=us_zone[zz]
            for sky in ["NGC","SGC"]
                isel=findall((zone_num .==zone) .&  (ind_NS[sky]))
                if(size(isel,1)==0)
                    continue
                end
                
                tgroup="$(sky)/zone_$(zone)"
                if(sky in keys(file))
                    if("zone_$(zone)" in keys(file[sky]))
                        gtmp=file[sky]["zone_$(zone)"]
                    else
                        gtmp=JLD2.Group(file, tgroup)
                    end
                else
                    gtmp=JLD2.Group(file, tgroup)
                end

                if(!("zone_order" in keys(groups_order[sky])) )
                        groups_order[sky]["zone_order"]=[]
                        groups_order[sky]["nobj"]=[]
                end
                
                append!(groups_order[sky]["zone_order"],["$(tgroup)/"])
                append!(groups_order[sky]["nobj"],size(isel,1))
                
                for quant in quant_write                               
                    gtmp["$(quant)"] = quant_dic[quant][isel] 
                end #ends the quant loop
            end#ends the sky loop
        end #ends the zone loop
        
        #Now write the zone order
        groups_order["NS"]=Dict()
        groups_order["NS"]["zone_order"]=[]
        groups_order["NS"]["nobj"]=[]
        groups_order["NS"]["beg_index"]=[]
        
        count_NS=1
        for sky in ["NGC","SGC"]
            groups_order[sky]["beg_index"]=[]
            count=1
            #println("$(sky) $(size(groups_order[sky]["zone_order"])) $(size(groups_order[sky]["nobj"])))")
            for tt in 1:size(groups_order[sky]["zone_order"],1)
                append!(groups_order[sky]["beg_index"],count)
                count=count+groups_order[sky]["nobj"][tt]
                
                #Also append for count_NS
                append!(groups_order["NS"]["zone_order"],[groups_order[sky]["zone_order"][tt]])
                append!(groups_order["NS"]["nobj"],groups_order[sky]["nobj"][tt])
                append!(groups_order["NS"]["beg_index"],count_NS)
                count_NS=count_NS+groups_order[sky]["nobj"][tt]
            end
            
            file[sky]["order"]=groups_order[sky]
        end
        file["order"]=groups_order["NS"]
            
    end #closes the JLD2 file
    
    
end


function embed_sphere(ra, dec)
    #= Adapted from desimode.io """Embed `ra`, `dec` to a uniform sphere in three dimensions.
    Note that intentionallly x,y,z are kept at three arrays for efficiency as Julia uses column-major 
    and three arrays will be easier to manipulate compared to one
    """
    =# 
    
    deg2rad=pi/180.0
    
    phi = deg2rad.*ra
    theta = deg2rad.*(90.0 .- dec)
    r = @. sin(theta)
    nobj=size(ra,1)
    x_unit =@.  r * cos(phi)
    y_unit =@.  r * sin(phi)
    z_unit =@. cos(theta)
    
    return x_unit,y_unit,z_unit
end


function embed_sphere!(ra,dec,x_unit,y_unit,z_unit)
    #= Adapted from desimode.io """Embed `ra`, `dec` to a uniform sphere in three dimensions.
    Note that intentionallly x,y,z are kept at three arrays for efficiency as Julia uses column-major 
    and three arrays will be easier to manipulate compared to one
    This version uses pre-allocated arrrays to store data and hence modifies its argument
    """
    =# 
    
    deg2rad=pi/180.0
    
    phi = deg2rad.*ra
    theta = deg2rad.*(90.0 .- dec)
    r = @. sin(theta)
    nobj=size(ra,1)
    @. x_unit = r * cos(phi)
    @. y_unit = r * sin(phi)
    @. z_unit = cos(theta)
end

function embed_sphere_xyzmatrix(ra,dec)
    #= calles embed sphere but returns a matrix which is sometimes convenient=#
    nra=size(ra,1)
    txyz=zeros(nra,3)
    embed_sphere!(ra, dec,view(txyz,1:nra,1),view(txyz,1:nra,2),view(txyz,:,3))
    
    return txyz
end

function load_DESI_tiles(tile_file,program,pass,sky)
    #=load the tiles from a tile file and select for program and pass if given=#
    
    #open fits file
    fin=FITS(tile_file)
    println("Loading tiles from: $(tile_file)")
    
    nrow=size(read(fin[2],"RA"),1)
    
    tiles_dic=Dict()
    quant_need=["TILEID","RA","DEC","PASS"]
    
    sel_dic=Dict()
    if(program!=nothing)
        sel_dic["PROGRAM"]=program
    end
    if(pass!=nothing)
        sel_dic["PASS"]=pass
    end
    
    if(program==nothing && pass==nothing)
        indsel=1:nrow
    else
        indbool=trues(nrow)
        for tkey in keys(sel_dic)
            tval=read(fin[2],tkey)
            indbool .= indbool .& (tval .==sel_dic[tkey])
        end
        indsel=findall(indbool)
    end
    
    for qq in quant_need
        tiles_dic[qq]=keepat!(read(fin[2],qq),indsel)
    end
        
    #Now apply sky selection if needed
    if(sky!=nothing)
        ind_NS=Dict()
        ind_NS["NGC"],ind_NS["SGC"]=split_Galactic_cap(tiles_dic["RA"],tiles_dic["DEC"],"index")
        for qq in quant_need
            tiles_dic[qq]=tiles_dic[qq][ind_NS[sky]]
        end
    end
   
    close(fin)
    
    return tiles_dic
end

function deg2rad(rdeg)
    return pi*rdeg/180
end

function Assign_tile_id(tiles_dic,tile_radius_indeg,ra_tracer, dec_tracer)
    #= assign the object to the tile id, assumes tiles are non-overlapping=
    1) get xyz_unit for tile centers and ra, dec_tracer
    2) define the scales , assign zones, select tiles
    3) scan the object determine dot product and assign tiles 
    =#
    
    
    nobj=size(ra_tracer,1)
    ntile=size(tiles_dic["RA"],1)
    
    #min_RA=minimum(ra_tracer)-5.0*tile_radius_indeg
    #determine the radius in dot product unit
    radius_in_costheta=cos(deg2rad(tile_radius_indeg))
    #@show radius_in_costheta
    
    #determine the tiles centers in unit sphere
    tiles_xyz= embed_sphere_xyzmatrix(tiles_dic["RA"],tiles_dic["DEC"])
    
    #determine the tracers on unit sphere
    tracer_xyz=embed_sphere_xyzmatrix(ra_tracer,dec_tracer)
    
    
    count_assign=zeros(Int32,nobj)
    tileid_assign=zeros(Int32,nobj)
    
    #scan the two
    #rdist=zeros(Float64,nobj)
    for tt in 1:ntile
        rdist=costheta_uv(tiles_xyz[tt,:],tracer_xyz)
        indsel=findall(rdist .> radius_in_costheta)
        #@show(tiles_xyz[tt,:],size(indsel),nobj,minimum(rdist),maximum(rdist))
        if(size(indsel,1)==0)
            continue
        end
        tileid_assign[indsel] .= tiles_dic["TILEID"][tt]
        count_assign[indsel] .= count_assign[indsel] .+ 1
    end
    @assert maximum(count_assign)<=1
    return tileid_assign
end


function traverse_group(jldfile,group)
    #= Traverse a JLD2 group for all the file within by loading order
    =#
    #jldopen("mocks/JLD2_data/LRG.jld2", "r") do file
    #    tgroup=file["$(group)/order"]
    #end
    
    return load(jldfile,"$(group)/order")
end

function check_variable_exists(jldfile,var_path)
    #= checks if a given variable path exists in the file=#
    var_arr=split(var_path,"/")
    file=jldopen(jldfile,"r")
    res=true
    tfile=file
    for tcol in var_arr
        if(!(tcol in keys(tfile)))
            res=false
            break
        end
        tfile=tfile[tcol]
    end
        
    close(file)
    return res
end

function load_groupzone_column(jldfile,group_zone,column)
    #=loads a column for a particular group and zone
    column should be a single column name=#
    
    this_col=load(jldfile,"$(group_zone)/$(column)")
    
    return this_col
end

function Append_groupzone_xyz_unit(jldfile,group_zone)
    #=Add xyz in unit sphere for a particular group and zone
    First check for if the column name exists=#
    columns=["X_UNIT","Y_UNIT","Z_UNIT"]
   
    for col_name in columns 
        if(check_variable_exists(jldfile,"$(group_zone)$(col_name)"))
            println("exists: $(group_zone)$(col_name)")
            return
        end
    end
    
    #loads the ra,dec
    ra=load_groupzone_column(jldfile,group_zone,"RA")
    dec=load_groupzone_column(jldfile,group_zone,"DEC")
    #ra,dec=load(jldfile,"$(group)/RA","$(group)/DEC")
    
    #embeds the sphere
    unit_xyz=embed_sphere_xyzmatrix(ra, dec)
    
    #Now write this to the file
    jldopen(jldfile, "a+") do file
        for cc in 1:size(columns,1) 
            col_name=columns[cc]
            write(file,"$(group_zone)/$(col_name)",unit_xyz[:,cc])
        end
    end
end


function load_group_column(jldfile,group,columns)
    #=loads a column for a particular group scanning all the zones
    columns should be Array of available properties=#
    
    
    #load the zones in group
    jld_order=traverse_group(jldfile,group)
    
    #initialize the ditionary
    col_dic=Dict()
    for col in columns
        col_dic[col]=[]
    end
    
    for iz in 1:size(jld_order["zone_order"],1)
        tg=jld_order["zone_order"][iz]
        
        
        for col in columns
            append!(col_dic[col],load_groupzone_column(jldfile,tg,col))
        end
            
    end
    return col_dic
end


function Indexing_column(jldfile,group,columns,mask)
    #=Generate index file for groups
    writes a ditionary for each minimum and maximum values=#
    
    #load the zones in group
    jld_order=traverse_group(jldfile,group)
    
    #initialize the ditionary
    minmax_dic=Dict()
    for col in columns
        minmax_dic[col]=zeros(size(jld_order["zone_order"],1),2)
    end
    
    for iz in 1:size(jld_order["zone_order"],1)
        tg=jld_order["zone_order"][iz]
        
        for col in columns
            this_col=load_groupzone_column(jldfile,tg,col)
            if(mask!=nothing)
                ind_sub=falses(size(this_col,1))
                for mm in mask
                    ind_sub .= ind_sub .|| (this_col .== mm)
                end
                indsel=findall(@. !ind_sub)
                keepat!(this_col,indsel)
            end
            if(size(this_col,1)==0)#This will happen only when mask is available
                minmax_dic[col][iz,:]=[mask[1],mask[1]]
            else
                minmax_dic[col][iz,:]=[minimum(this_col),maximum(this_col)]
            end
        end    
    end
    
    #Now check if indexing exists
    update_dictionary_jld2file(jldfile,group,minmax_dic,"indexing")
    
    
    return minmax_dic
    
end

function LISTING_column(jldfile,group,outname,columns_use)
    #=For each zone generate the unique list of values in a particular colums
    This should be only used for columns which has small number of unique values
    for example TILEID=#
    
    #load the zones in group
    jld_order=traverse_group(jldfile,group)
    
    #initialize the ditionary
    unique_dic=Dict()
    unique_dic[outname]=Dict()
    #for col in columns
    #    minmax_dic[col]=zeros(size(jld_order["zone_order"],1),2)
    #end
    
    for iz in 1:size(jld_order["zone_order"],1)
        tg=jld_order["zone_order"][iz]
        
        comb_col=[]
        for col in columns_use
            this_col=load_groupzone_column(jldfile,tg,col)
            append!(comb_col,unique(this_col))
        end
        unique_dic[outname][iz]=unique(comb_col)
    end
    
    #Now check if indexing exists
    update_dictionary_jld2file(jldfile,group,unique_dic,"listing")
    
    
    return unique_dic
    
end



function update_dictionary_jld2file(jldfile,group,new_index,key)
    #= check if dictionary exists then load and merge with new dictionary and write the new one
    first delete the old one before writing=#
    
    dic_exists=check_variable_exists(jldfile,"$(group)/$(key)")
    old_dic=Dict()
    
    if( dic_exists)
        old_dic=load(jldfile,"$(group)/$(key)")
    end
    
    
    #Include old indexing if not updated in new_index     
    for tcol in keys(old_dic)
        if(tcol in keys(new_index))
            continue
        end
        new_index[tcol]=old_dic[tcol]
    end 
   
    jldopen(jldfile, "a+") do file
        #delete if needed or just write
        if(dic_exists)
            Base.delete!(file,"$(group)/$(key)")
        end
        write(file,"$(group)/$(key)",new_index)
    end
end

function Append_column(jldfile,group,columns)
    #=Add the given propety if doesn't already exists
    A reference group needs to be given which can empty string or NGC/SGC=#
    
    #load the zones in group
    jld_order=traverse_group(jldfile,group)
    
    for iz in 1:size(jld_order["zone_order"],1)
        tg=jld_order["zone_order"][iz]
        
        for col in columns
            if(col=="XYZ_UNIT")
                Append_groupzone_xyz_unit(jldfile,tg)
            else
                println("Not Implemented Append for $(col)")
            end
        end    
    end
end


function Append_TileID(jldfile,tile_file,group,program,npass,tile_radius_indeg,overwrite)
    #= Assign tile id for each pass to a given tracer
    =#
    
    if("SGC" in split(group,"/"))
        sky="SGC"
    elseif("NGC" in split(group,"/"))
        sky="NGC"
    else
        sky=nothing
    end
       
    
    #load tiles for all the passes
    tiles_dic=Dict()
    for pass in 0:npass-1
        tiles_dic[pass]=load_DESI_tiles(tile_file,program,pass,sky)
    end
            
    
    #load the zones in group
    jld_order=traverse_group(jldfile,group)
    
    for iz in 1:size(jld_order["zone_order"],1)
        tg=jld_order["zone_order"][iz]
        
        ra_tracer=load_groupzone_column(jldfile,tg,"RA")
        dec_tracer=load_groupzone_column(jldfile,tg,"DEC")
        
        for pass in 0:npass-1
            tile_pass_exists=check_variable_exists(jldfile,"$(tg)TILEID_PASS$(pass)")
            
            if(tile_pass_exists & !overwrite)
                println("$(tg)/TILEID_PASS$(pass) exists with no-overwrite")
                continue
            end
            #@show pass,tg,tile_pass_exists & !overwrite, tile_pass_exists
            
            tileid_assign=Assign_tile_id(tiles_dic[pass],tile_radius_indeg,ra_tracer, dec_tracer)
            jldopen(jldfile, "a+") do file
                #delete if needed or just write
                if(tile_pass_exists)
                    Base.delete!(file,"$(tg)TILEID_PASS$(pass)")
                end
                write(file,"$(tg)TILEID_PASS$(pass)",tileid_assign)
            end
        end    
    end
end


""" load all objects in a particular tile of given tracer
    also returns indices with reference to the group which can be used to load 
    or write back additional properties later
    if the tile_pass>0 then also look for object assignemnt in lower pass tiles and return a collision
    array
"""
function Load_tracers_intile(config,tracer,tile_id,tile_pass,group,columns;numobs_need=1,num_fba=1)
    #= load all objects in a particular tile
    also returns indices with reference to the group which can be used to load 
    or write back additional properties later
    =#
    
    jldfile=config["target"][tracer]["JLDfile" ]


    #load the zones in group
    jld_order=traverse_group(jldfile,group)

    #load the listing to identify the zones which has tile
    listing=load(jldfile,"$(group)/listing")


    indsel=[]
    #Dictionary to store the tracer properties
    tracer_tile=Dict()
    tracer_tile["ref_group"]=group
    tracer_tile["TILEID"]=tile_id
    tracer_tile["PASS"]=tile_pass
    tracer_tile["index"]=[]
    for col in columns
        tracer_tile[col]=[]
    end

    #list of tile which will have over lap in this region
    collided_tiles=zeros(Int32,0)

    for iz in 1:size(jld_order["zone_order"],1)
        tg=jld_order["zone_order"][iz]

        #if tileid not in listing then continue
        if( !(tile_id in listing["TILEID"][iz]))
            continue
        end


        #Now load the tile_id is in this zone
        this_id_arr=load(jldfile,"$(tg)TILEID_PASS$(tile_pass)")
        #It is important to make sure the reference groups are used consistently for indices to work
        indsel_this=findall(this_id_arr .== tile_id)

        append!(tracer_tile["index"], indsel_this.+ jld_order["beg_index"][iz] )

        for col in columns
            tcol_val=keepat!(load(jldfile,"$(tg)$(col)"),indsel_this)
            append!(tracer_tile[col],tcol_val)
        end

        #find all the tiles colliding with lower passes passes starts with 0
        if(tile_pass>0)
            for pre_pass in 0:tile_pass-1
                pre_id_arr=keepat!(load(jldfile,"$(tg)TILEID_PASS$(pre_pass)"),indsel_this)
                utiles=unique(pre_id_arr)
                for ut in utiles
                    if((!(ut in collided_tiles)) & (ut>0))
                        append!(collided_tiles,ut)
                    end
                end
                #println(pre_pass,' ',utiles)
            end
            #println("colide: ",collided_tiles)
        end
    end #end of iz loop
    
    #Now load the information whether object is already assigned
    if(size(collided_tiles,1)>0)
        count_assigned=count_assignment_tracer(config,tracer,tracer_tile,collided_tiles,num_fba) 
        #Now populate the collided boolean array:
        #one set to true will not be considered for assignements
        tracer_tile["collided"] = count_assigned .>= numobs_need
    else
        tracer_tile["collided"] = zeros(Bool,size(tracer_tile["index"],1),num_fba)
    end

    return tracer_tile

end

"""Counts the number of time any object is assigned in the given list of tiles
Only counts for the objects passed in the dictionary
fba_jldfile : JLD2 file storing assignmnet for the individual tiles
tracer_tile: A dictionary loaded for a tile along with its index
tile_list: list of tileID to count the assignmnet in
num_fba: number of fba realizations was run
"""
function count_assignment_tracer(config,tracer,tracer_tile,tile_list,num_fba)
    #total number od objects 
    nobj_tile=size(tracer_tile["index"],1)
    #index to sort the objects in tile by index field
    tile_isort=sortperm(tracer_tile["index"])

    #allocate memory
    count_assigned=zeros(nobj_tile,num_fba)

    #scan the colliding tile and count number of assignment
    for col_tile in tile_list
        #FBA jldfile for the tile
	fba_jldfile=TileFBA_FileName(config,col_tile)

	index_pre_tile=load(fba_jldfile,"$(tracer)_index")
	bool_ass=load(fba_jldfile,"$(tracer)_Bool_ass")
        #make sure the tile was assigned with consistent num_fba
	@assert num_fba==size(bool_ass,2) "inconsistent num_fba \nfor tile: $(col_tile) in $(fba_jldfile)"

        #find the sorting permutation
        isort_perm=sortperm(index_pre_tile)

        #rearrange with sorted index for this tile
        index_pre_tile=index_pre_tile[isort_perm]
        bool_ass=bool_ass[isort_perm,:]

        #append the counts #follow the logic of sorted indices
        ctile=1
        for ti in 1:size(index_pre_tile,1)
            while(index_pre_tile[ti]>tracer_tile["index"][tile_isort[ctile]])
                ctile +=1
                if(ctile>nobj_tile)
                    break
                end
            end

            if(ctile>nobj_tile)
                break
            end

            if(index_pre_tile[ti]==tracer_tile["index"][tile_isort[ctile]])
                count_assigned[tile_isort[ctile],:] .+= bool_ass[ti,:]
            end
        end 
    end #end of the col_tile loop
    return count_assigned
end# end of function



function pre_process_FITS2JLD_tracer(fname,priority,priority_frac,outfile)
    #= Pre-process a given tracer in following step
    priority is array containing list of priority
    priority_frac is simply the fraction for each priority should sum to 1
      1) read the fits file
      2) select y5 (2^1) foot-print and n(z) (2^0)
      3)  create zone ins ra-dec 20binx20bin
      3) split in NGC/SGC and next steps for each zone
        a) In each zone sort them according to priority (if multiple priority exists)
        b) Create a group and write each zone in seperate group
    =#
    
    fin=FITS(fname)
    println("Loading file: $(fname)")
    
    #apply footprint and nz selection
    status=read(fin[2],"STATUS")
    indsel=findall((status .& 2^1) .* (status .& 2^0).>1)
    #Only for testing to reduce size
    #indsel=indsel[1:1000:size(indsel,1)]
    nobj=size(indsel,1)
    println("nobj selection (footprint and nz): $(nobj)")
    
    #read ra,dec of selection
    ra=keepat!(read(fin[2],"RA"),indsel)
    dec=keepat!(read(fin[2],"DEC"),indsel)
    
    println("Assigning zones and priority and sorting")
    #Assign zone number for each object
    zone_num=assign_sky_zones(ra,dec,20,20)
    #asign_priority
    obj_priority=assign_priority(nobj,priority,priority_frac)
    
    #sort by zone and then priority within zone in reverse order (higher priority first)
    indsel_pr_zone=indsel_pr_zone=sort_by_zone_priority(indsel,zone_num,obj_priority)  
    #for each priority further sorting can be done here if needed
    
    
    #re-ordering everything
    zone_num=zone_num[indsel_pr_zone]
    
    quant_dic=Dict()
    quant_dic["RA"]=ra[indsel_pr_zone]
    quant_dic["DEC"]=dec[indsel_pr_zone]
    quant_dic["PRIORITY"]=obj_priority[indsel_pr_zone]
    quant_dic["FITS_index"]=indsel[indsel_pr_zone]
    

    println("Writing to JLD2 file: $(outfile)")
    write_JLD2_file(outfile,zone_num,quant_dic)


    
    close(fin)
    
end


"""
find the true index same as inbuilt findall but does this inplace avoiding new memory allocation
returns the number of elements in the output index_arr to be used
"""
function find_index!(bool_array,index_arr)
    nel=0
    @inbounds for ii in 1:size(bool_array,1)
        if(bool_array[ii])
            nel +=1
            index_arr[nel]=ii
        end
    end
    @assert nel<=size(index_arr,1) "The input index_array $(size(index_arr)) is smaller than needed for output $(nel)"
    return nel
end

"""
find the true index same as inbuilt findall but does this inplace avoiding new memory allocation
returns the number of elements in the output index_arr to be used
It performs an AND operation between the two inputs
"""
function find_index_AND!(bool_array1,bool_array2,index_arr)
    nel=0
    @inbounds for ii in 1:size(bool_array1,1)
        if(bool_array1[ii] & bool_array2[ii])
            nel +=1
            index_arr[nel]=ii
        end
    end
    @assert nel<=size(index_arr,1) "The input index_array $(size(index_arr)) is smaller than needed for output $(nel)"
    return nel
end


"""
find the true index same as inbuilt findall but does this inplace avoiding new memory allocation
returns the number of elements in the output index_arr to be used
It performs OR operation between the two inputs
"""
function find_index_OR!(bool_array1,bool_array2,index_arr)
    nel=0
    @inbounds for ii in 1:size(bool_array1,1)
        if(bool_array1[ii] || bool_array2[ii])
            nel +=1
            index_arr[nel]=ii
        end
    end
    @assert nel<=size(index_arr,1) "The input index_array $(size(index_arr)) is smaller than needed for output $(nel)"
    return nel
end

using Dates

"""
Converts a dictionary of events to the msg string desribing the time use
msg: any string to be prepended
events_dic: dictionar with array of events and times
depends in the Dates Package
example:
edic=Dict()
edic["events"]=["full","prepare","Assing"]
edic["times"] = [now(UTC)]
sleep(0.5)
append!(edic["times"],[now(UTC)])
sleep(0.2)
append!(edic["times"],[now(UTC)])

convert_events_to_msg(edic,"hello ")
wil give output:
"hello  full:0.709 sec(prepare:0.504 sec,Assing:0.205 sec,)"
"""
function convert_events_to_msg(events_dic,msg)
    
    nevent=length(events_dic["events"])
    if(nevent==0)
        return ""
    end
    
    #tdiff=string(Dates.canonicalize(Dates.CompoundPeriod(events_dic["times"][end] - events_dic["times"][1])))
    tdiff=Dates.DateTime(events_dic["times"][end]) - Dates.DateTime(events_dic["times"][1])
    tdiff=string(tdiff/Millisecond(1000), " sec")
    msg=string(msg," ",events_dic["events"][1],":",tdiff)
    
    msg_local=""
    for ii in 2:nevent
        #tdiff=string(Dates.canonicalize(Dates.CompoundPeriod(events_dic["times"][ii] - events_dic["times"][ii-1])))
        tdiff=Dates.DateTime(events_dic["times"][ii]) - Dates.DateTime(events_dic["times"][ii-1])
        tdiff=string(tdiff/Millisecond(1000), " sec,")
        msg_local=string(msg_local,events_dic["events"][ii],":",tdiff)
    end
    
    if(msg_local=="")
        return msg
    else
        msg_local=string(msg,"(",msg_local,")")
        return msg_local
    end
end
