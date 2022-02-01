#Tile focused function prepare, and sort a tile etc.
"""
The grid_sale is defined as the twice the maximum of sum of theta and phi arm
The factor of 2 is used so that we can also obtain look at all potential collision
otherwise this factor of 2 is not needed for only to determine all possible assignments
"""
function get_grid_scale(fp_dic)
    grid_scale=2.0*maximum(fp_dic["LENGTH_R1"].+fp_dic["LENGTH_R2"])
    return grid_scale
end

"""sorts the targets by zone and then priority
"""
function sort_targets_ZonePriority!(targets_dic)
    #get the zone
    zone_xy=zeros(eltype(targets_dic["zone_x"]),targets_dic["ntarget"])
    get_zonexy!(targets_dic["zone_x"],targets_dic["zone_y"],targets_dic["fp_bound"],zone_xy)
    
    #Actual sorting, remapping and zone_map consistinting of indices range
    sort_ZoneXp_Apply!(zone_xy,targets_dic["PRIORITY"],targets_dic)
    
end

function sort_hardware_ZonePriority!(fp_dic,fp_bound)
    #get the zone
    nfp=size(fp_dic["zone_x"])
    zone_xy=zeros(eltype(fp_dic["zone_x"]),nfp)
    get_zonexy!(fp_dic["zone_x"],fp_dic["zone_y"],fp_bound,zone_xy)
    
    #Actual sorting, remapping and zone_map consistinting of indices range
    device_code=zeros(Int8,nfp)
    #POS is 0 and everything else is one as we only want to keep all positionr together 
    #at the end and since it does reverse sort on xp
    device_code[findall(fp_dic["DEVICE_TYPE"] .!= "POS")] .= 1
    
    sort_ZoneXp_Apply!(zone_xy,device_code,fp_dic)
end


"""
sort the values based on combination of zone and XP
such that zones are together and XP is reverse ordered for each zone
This assumes zones are integers
It then apply the new sorting to everything in tdic which is array and of same size as Xp
Finally generate a map dictionary which when zone is specified gives the minimum and maximum index
One needs to be careful that zone_x and zone_y are not modified inside the function
as this function changes its arguments, one needs to make sure none of the Array is modified/sorted twice
even if they are refered within tdic and also provided directy
"""
function sort_ZoneXp_Apply!(zone_xy,Xp,tdic)
    
    #count number of entries
    nxp=size(Xp,1)
    
    
    #get the range of Xp
    xp_min=minimum(Xp)
    xp_max=maximum(Xp)

    #new quantity for sorting
    #The factor of 0.5 is to make sure the maximum value of priority is lower than next zone
    #It can be any value below 1.0, 0.5 is just a conveninet choice
    #This sort of decide the resolution of values in XP
    new_q=(Xp .- xp_min).*(0.5/(xp_max-xp_min))
    new_q=new_q+zone_xy
    izone_sort=sortperm(new_q,rev=:true)
    
    
    modified_ptr=[] #keep track of al the arrays sorted already to not apply sorting twice
    #Now apply this to all of the target arrays
    for tkey in keys(tdic)
        if(tdic[tkey] isa Union{BitArray,Array})
            #this_ptr=pointer_from_objref(tdic[tkey])
            this_ptr=repr(UInt64(pointer_from_objref(tdic[tkey])))
            if((size(tdic[tkey],1)==nxp) & (!(this_ptr in modified_ptr)) )
	        if(tdic[tkey] isa Union{BitMatrix,Matrix} )
	            tdic[tkey] .= tdic[tkey][izone_sort,:]
		else
	            tdic[tkey] .= tdic[tkey][izone_sort]
		end
		append!(modified_ptr,[this_ptr])
                #println("sorting $(tkey) $(this_ptr)")
	    elseif(!(tkey in ["fp_bound","file_map","ntarget"]))
                println("Warning: not sorting? ",tkey,size(tdic[tkey],1), this_ptr)
            end
        end
    end
    
    #Now generate the map from of indices for each zone
    #this_ptr=pointer_from_objref(zone_xy)
    this_ptr=repr(UInt64(pointer_from_objref(zone_xy)))
    if(!(this_ptr in modified_ptr) )
        zone_xy .= zone_xy[izone_sort]
        append!(modified_ptr,[this_ptr])
	#println("sorting zone_xy $(this_ptr)")
    end
        
    uszone=sort(unique(zone_xy),rev=false)
    first_zone=indexin(uszone,zone_xy)
    last_zone=(nxp+1) .- indexin(uszone,reverse(zone_xy))

    #generate the zone index map
    nunique_zone=size(uszone,1)

    tdic["zone_map"]=Dict()
    uu=1
    tdic["zone_map"][uszone[uu]]=[first_zone[uu],last_zone[uu]]
    max_per_zone=tdic["zone_map"][uszone[uu]][2]-tdic["zone_map"][uszone[uu]][1]
    min_per_zone=max_per_zone
    for uu in 2:nunique_zone
        tdic["zone_map"][uszone[uu]]=[first_zone[uu],last_zone[uu]]
        this_zone=last_zone[uu]-first_zone[uu]+1
        if(this_zone<min_per_zone)
            min_per_zone=this_zone
        elseif(this_zone>max_per_zone)
            max_per_zone=this_zone
        end
    end
    
    tdic["zone_map_keys"]=keys(tdic["zone_map"])
    tdic["count_maxper_zone"]=max_per_zone
    tdic["count_minper_zone"]=min_per_zone
end

function get_zonexy!(zone_x,zone_y,fp_bound,zone_xy)
    @. zone_xy=zone_x+ (zone_y*fp_bound.nx)
end

function get_zonexy(zone_x,zone_y,fp_bound)
    zone_xy=zone_x+ (zone_y*fp_bound.nx)
    return zone_xy
end

function split_zonexy(zone_xy,fp_bound)
    zone_x=zone_xy%fp_bound.nx
    zone_y=Int((zone_xy-zone_x)/fp_bound.nx)
    return zone_x,zone_y
end



"""
Loads all the tracer in a tile and attach a target_identification column to it
"""
function load_targets_intile(config,tile_id,tile_pass,group;
      columns=["X_UNIT","Y_UNIT","Z_UNIT","PRIORITY"],num_fba=1)
    targets_dic=Dict()
    targets_dic["file_map"]=Dict()
    
    col_accumulate=["index","filekey","collided"]
    for col in columns
        if(!(col in col_accumulate))
            append!(col_accumulate,[col])
        end
    end
    
    targets_dic["filekey"]=zeros(Int8,0)
    
    
  
    fno=1
    for tracer in keys(config["target"])
        #fname=config["target"][tracer]["JLDfile" ]
	#fba_jldfile=config["target"][tracer]["FBA_JLDfile" ]
	#fba_jldfile=TileFBA_FileName(config,tile_id)
        targets_dic["file_map"][fno]=[tracer,config["target"][tracer]["JLDfile" ]]
        
        tracer_tile=Load_tracers_intile(config,tracer,tile_id,tile_pass,group,columns;
			numobs_need=config["target"][tracer]["Num_obs"],num_fba=num_fba)
        
        nobj=size(tracer_tile["index"],1)
        for col in col_accumulate
            if(col=="filekey")
                targets_dic[col]=append!(targets_dic[col],(zeros(Int8,nobj) .+ fno))
            elseif(col in keys(targets_dic))
                #targets_dic[col]=append!(targets_dic[col],tracer_tile[col])
		targets_dic[col]=[targets_dic[col];tracer_tile[col]]
            else
                targets_dic[col]=tracer_tile[col]
            end
        end
        #println(tracer,' ',unique(targets_dic["PRIORITY"]),nobj,fno)
        fno +=1
	#println("load_targets_intile: ",tile_id,' ',num_fba,size(targets_dic["collided"]),
	#	' ',size(tracer_tile["collided"]))
    end
    targets_dic["ntarget"]=size(targets_dic["index"],1)
    return targets_dic
end


""" This prepare a tile for assignment given the config dictionary
This function execute the following steps
1) load targets in a given tile tile
"""
function Prepare_tile(config,tile_index,tiles_dic;group="",tile_date="2019-09-16T00:00:00",
        plate_scale=nothing,fp_dic=nothing,exc_dic=nothing,num_fba=1)
    #get the tile centre 
    telra=tiles_dic["RA"][tile_index]
    teldec=tiles_dic["DEC"][tile_index]
    tile_id=tiles_dic["TILEID"][tile_index]
    tile_pass=tiles_dic["PASS"][tile_index]


    #load the plate_scale
    if(plate_scale==nothing)
        plate_scale=load_platescale(config["focal_plane"]["focalplane_dir"])
    end
    
    #Load the focal plane if needed
    if((fp_dic==nothing) || (exc_dic==nothing))
        fp_dic,exc_dic=load_hw_full_FocalPlane!(config["focal_plane"];date=tile_date)
    elseif((fp_dic["tile_date"]!= tile_date) || (exc_dic["tile_date"]!= tile_date) )
        fp_dic,exc_dic=load_hw_full_FocalPlane!(config["focal_plane"];date=tile_date)
    end
    
    #Load target
    targets_dic=load_targets_intile(config,tile_id,tile_pass,group;
                    columns=["X_UNIT","Y_UNIT","Z_UNIT","PRIORITY"],num_fba=num_fba)

    #convert to focal plane
    targets_dic["X_FP"], targets_dic["Y_FP"]= xyz_unit_2_xy_focalplane(telra, teldec, 
    targets_dic["X_UNIT"],targets_dic["Y_UNIT"],targets_dic["Z_UNIT"],plate_scale)

    #define the grid scale
    grid_scale=get_grid_scale(fp_dic)

    #get the focal plane bound
    targets_dic["fp_bound"]=get_bound_focal_plane(targets_dic["X_FP"], targets_dic["Y_FP"],grid_scale;pad=1.0)

    #Assign zones
    targets_dic["zone_x"],targets_dic["zone_y"]=Assign_xy_zone(targets_dic["X_FP"], targets_dic["Y_FP"],targets_dic["fp_bound"])

    #Assifn zones for hardware in fp
    fp_dic["zone_x"],fp_dic["zone_y"]=Assign_xy_zone(fp_dic["OFFSET_X"], fp_dic["OFFSET_Y"],targets_dic["fp_bound"])

    #sort targets
    sort_targets_ZonePriority!(targets_dic)
    #sort devices
    sort_hardware_ZonePriority!(fp_dic,targets_dic["fp_bound"])
    
    return targets_dic,fp_dic,exc_dic
end


"""Runs the assignment for a given tile
num_fba: number of fba realization to generate
max_vertices maximum number of vertices any of the positioner polygon can have
max_ngb_pos: maximum number of negibohuring positioner any positioner can have
"""
function Assign_Tile(config,fp_dic,exc_dic,targets_dic;num_fba=3,max_vertices=30,max_ngb_pos=9)
    #Make sure the num_fba is propagated correctly
    @assert num_fba==size(targets_dic["collided"],2) "num_fba is not propagated correctl $(num_fba), $(size(targets_dic["collided"],2))"

    #generate two robot positioner to allocate memory to be used for efficiency
    pos_cur= POSRobotBody(max_vertices) #Assuming we will never have more than 30 vertices
    
    max_pos_zone=Integer(9)*fp_dic["count_maxper_zone"]
    pos_ngb_dic=Dict("index_central_positioner"=>0,"index"=>zeros(Int32,1max_pos_zone),
        "pos"=>[],"pos_index"=>[])
    for ii in 1:max_ngb_pos
        append!(pos_ngb_dic["pos"],[POSRobotBody(max_vertices)])
        append!(pos_ngb_dic["pos_index"],Int32(0))
    end
    
    #declare the strucutre to load for this particular zone with size of maximum possible
    #THis is so that only one time memory allocation is needed for efficiency
    #For any zone there can be 8 other neigbouring zone and hence a factor of 9
    max_target=Integer(9)*targets_dic["count_maxper_zone"]
    max_npriority=size(unique(targets_dic["PRIORITY"]),1)
    zone_target=Zone_Target(max_target,max_npriority)
    #to store the index values
    pot_ass_index=zeros(Int32,max_target)
    #pot_coll_index=zeros(Int32,max_target)
    #TO store either assigned or collided indx
    pot_ass_coll_index=zeros(Int32,max_target)
    # To store in_poly and in_circle behaviours
    in_poly=zeros(Bool,max_target)
    in_circle=zeros(Bool,max_target)
    
    #boolean for collision
    #collided=zeros(Bool,targets_dic["ntarget"],num_fba)
    #If something is already assigned needed number then it will be set to true
    #so that we do not assign this, for pass=0 this will be false everywhere in the begining
    collided = targets_dic["collided"] 

    #integer for potential collision
    potential_coll=zeros(Bool,targets_dic["ntarget"],num_fba)
    
    #dictionary for collision locations
    #integer for Assignment
    position_assigned=zeros(Int32,size(fp_dic["OFFSET_X"],1),num_fba)
    if(config["OUTPUT"]["ThetaPhi_pos"])
        theta_assigned=zeros(Float64,size(fp_dic["OFFSET_X"],1),num_fba)
        phi_assigned=zeros(Float64,size(fp_dic["OFFSET_X"],1),num_fba)
    end
    
    local_ass=zeros(Int64,num_fba) #To keep the assignent in local index
    ind_ri=zeros(Int32,num_fba)
    
    
    #Pre-allocating some memory to avoid doing this millions of time for moving positioners
    move_diff_v=zeros(Float64,2)
    move_cos_sin=zeros(Float64,2)
    move_cos_sin_sum=zeros(Float64,2)
    move_cstheta=zeros(Float64,2)
    move_csphi=zeros(Float64,2)
    
    #for reproducibility we set the random seed
    Random.seed!(config["RandomSeed"]);
    
    for zone_xy in keys(fp_dic["zone_map"])
        #get the zone
        zone_x,zone_y=split_zonexy(zone_xy,targets_dic["fp_bound"])
    
        #load the targets in zones around
        load_target_zone(zone_target,zone_x,zone_y,targets_dic["fp_bound"],targets_dic)
        #unique sorted priority for this zone
        #us_priority=sort(unique(zone_target.PRIORITY[1:zone_target.ntarget[1]]),rev=true)
        
        #@show zone_target.ntarget,zone_target.usp
        
        tmp=0
        #loop over positioner in this zone
        for pp in fp_dic["zone_map"][zone_xy][1]:fp_dic["zone_map"][zone_xy][2]
            #load the info for this positioner
            #if((fp_dic["DEVICE_TYPE"][pp]!="POS") || (fp_dic["LENGTH_R1"][pp]==0.0) )
            if(!valid_positionter(fp_dic,pp;pos_type="assignment"))
                continue
            else #load positioner
                load_positioner!(pos_cur,pp,fp_dic,exc_dic)
                #load the negihbours of thie positioner
                load_positioner_neighbours!(pos_cur,pp,targets_dic["fp_bound"],fp_dic,exc_dic,pos_ngb_dic)
                #pos_x=fp_dic["OFFSET_X"][pp]
                #pos_y=fp_dic["OFFSET_Y"][pp]
            end
            
            #maximum distance to which positioner reach
            #max_r2=(fp_dic["LENGTH_R1"][pp]+fp_dic["LENGTH_R2"][pp])^2 
            #minimum distance to which positioner can reach
            #min_r2=(fp_dic["LENGTH_R1"][pp]-fp_dic["LENGTH_R2"][pp])^2 
            
            coll_r2=(sqrt(pos_cur.arms2_lim[2])+pos_cur.body_minmax[2])^2 
            
            
            
            #scan for potential assignment and potential target
            for ii in 1:zone_target.ntarget[1]
                rad2=(zone_target.X_FP[ii]-pos_cur.xy[1])^2
                rad2 +=(zone_target.Y_FP[ii]-pos_cur.xy[2])^2
                if(rad2>coll_r2)#no collision
                    zone_target.pot_ass[ii]=false
                    zone_target.pot_coll[ii]=false
                elseif(rad2>pos_cur.arms2_lim[2]) #maximum arms length
                    zone_target.pot_ass[ii]=false
                    zone_target.pot_coll[ii]=true
                elseif(rad2>=pos_cur.arms2_lim[1]) #minimum arms length
                    zone_target.pot_ass[ii]=true
                    zone_target.pot_coll[ii]=false
                    #get theta_phi
                    zone_target.theta_pos[ii],zone_target.phi_pos[ii]=Hardware_xy_to_thetaphi(pos_cur,
                        zone_target.X_FP[ii],zone_target.Y_FP[ii],rad2)
                    #check if within theta phi range for the positioner
                    if(outside_theta_phi_range(zone_target.theta_pos[ii],zone_target.phi_pos[ii],
                                pos_cur.body_theta_limit,pos_cur.head_theta_limit) )
                        zone_target.pot_ass[ii]=false
                        zone_target.pot_coll[ii]=true
                    end
                else
                    zone_target.pot_ass[ii]=false
                    zone_target.pot_coll[ii]=false
                end
            end
                
            #a = @allocated begin
            #end; if a > 0 println("mem alloc: ",' ',a/1024/1024) end

            #priority_sel=Dict()
            #for up in us_priority
            #     priority_sel[up]=findall((zone_target.PRIORITY[1:zone_target.ntarget[1]].==up) .& 
            #        zone_target.pot_ass[1:zone_target.ntarget[1]])
            #end
            #@show us_priority,priority_sel
            
            #possible_ass=zeros(Int64,100) #To store the actual possiblity
            
            
            #pot_ass_index=findall(zone_target.pot_ass[1:zone_target.ntarget[1]])
            #pot_coll_index=findall(zone_target.pot_coll[1:zone_target.ntarget[1]])
            npot_ass=find_index!(view(zone_target.pot_ass,1:zone_target.ntarget[1]),pot_ass_index)
            #npot_coll=find_index!(zone_target.pot_coll[1:zone_target.ntarget[1]],pot_coll_index)
            npot_ass_coll=find_index_OR!(view(zone_target.pot_ass,1:zone_target.ntarget[1]), 
                view(zone_target.pot_coll,1:zone_target.ntarget[1]),pot_ass_coll_index)
            
           
                
            coll_map_dic=Dict()
            #loop over realizations
            for ri in 1:num_fba
                count=0
                for uu in 1:zone_target.ntarget[2] #do this in priority order high to low  
                    #this_index=findall( zone_target.usp_bool_sel[1:zone_target.ntarget[1],uu] .& 
                    #    zone_target.pot_ass[1:zone_target.ntarget[1]])
                    #println(' ',uu,' ',size(this_index),' ',
                    #    size(findall(zone_target.usp_bool_sel[1:zone_target.ntarget[1],uu])) )
                    for izone in pot_ass_index[1:npot_ass] #scan all the object in potential assignment
                        #@show pp,uu,ri,izone,zone_target.usp_bool_sel[izone,uu]
                        if(!zone_target.usp_bool_sel[izone,uu])
                            continue
                        end
                        
                        #check if this object is collided
                        if(collided[zone_target.index[izone],ri])
                            continue
                        end
                        
                        
                        #@show ri,uu,izone
                        #check if this object is in potential collision
                        if(potential_coll[zone_target.index[izone],ri])
                            #Asuuming no-collision and should change if found collided
                            actual_collision=false
                            
                            #Now check if this is actual collision with another positioner
                            #scan over neighbour positioner assigned
                            for ingb in 1:pos_ngb_dic["count_ngb"]
                                #check if positioner is assigned to a target
                                if(position_assigned[pos_ngb_dic["pos_index"][ingb],ri]==0)
                                    continue
                                end
                                
                                
                                #Target to which positioner is assigned
                                ingb_index=position_assigned[pos_ngb_dic["pos_index"][ingb],ri]
                                #check if there is an entry for this collision test
                                if( (izone,ingb,ingb_index) in keys(coll_map_dic) )
                                    actual_collision=coll_map_dic[(izone,ingb,ingb_index)]
                                else #now since it is not tested already check for this collision
                                    #move the neighbour to the assigned location
                                    if(config["OUTPUT"]["ThetaPhi_pos"])#use pre-computed values
                                        th_ngb=theta_assigned[pos_ngb_dic["pos_index"][ingb],ri]
                                        phi_ngb=phi_assigned[pos_ngb_dic["pos_index"][ingb],ri]
                                    else #compute the assigned theta_phi
                                        itr=position_assigned[pos_ngb_dic["pos_index"][ingb],ri]
                                        th_ngb,phi_ngb=Hardware_xy_to_thetaphi(pos_ngb_dic["pos"][ingb],
                                            targets_dic["X_FP"][itr],targets_dic["Y_FP"][itr])
                                    end
                                    #Now move the positioner
                                    move_positioner_to_theta_phi!(pos_ngb_dic["pos"][ingb],th_ngb,phi_ngb,
                                        diff_v=move_diff_v,cos_sin=move_cos_sin,cos_sin_sum=move_cos_sin_sum,
                                        cstheta=move_cstheta,csphi=move_csphi)
                                    
                                    
                                    #move current positioner to this target
                                    move_positioner_to_theta_phi!(pos_cur,
                                        zone_target.theta_pos[izone] ,zone_target.phi_pos[izone],
                                        diff_v=move_diff_v,cos_sin=move_cos_sin,cos_sin_sum=move_cos_sin_sum,
                                        cstheta=move_cstheta,csphi=move_csphi)
                                    
                                    
                                    #check for the collision
                                    actual_collision=positioner_collided(pos_cur,pos_ngb_dic["pos"][ingb])
                                    #store this in local map
                                    coll_map_dic[(izone,ingb,ingb_index)]=actual_collision
                                    
                                end
                                
                                
                                if(actual_collision)
                                    break
                                end
                            end #end of loop for scanning the neighbouring positioner
                            
                            #check of the two positioners are collided
                            if(actual_collision)
                                continue
                            end
                        end
                        
                        #This object can be assigned
                        count=count+1
                        #possible_ass[count]=izone
                        zone_target.possible_ass[count]=izone
                    end #end of izone loop
                    
                    if(count>0)
                        #@show count,uu
                        izone_sel=zone_target.possible_ass[rand(1:count)]
                        local_ass[ri]=izone_sel
                        position_assigned[pp,ri]=zone_target.index[izone_sel]
                        if(config["OUTPUT"]["ThetaPhi_pos"])
                            theta_assigned[pp,ri]= zone_target.theta_pos[izone_sel] 
                            phi_assigned[pp,ri]= zone_target.phi_pos[izone_sel] 
                        end
                        break #to break the priority loop because we got assignment at this priority
                    else
                        local_ass[ri]=0
                        position_assigned[pp,ri]=0
                        #@show "nothing can be assigned",pp,ri,count
                    end
                end #end of uu loop
                #if(local_ass[ri]==0)
                #    @show "nothing can be assigned",pp,ri,count
                #end
            end #end of ri loop
            
            #break
            #udate the collision and potential collision
            uass=sort(unique(local_ass[local_ass .> 0]))
            
            #pot_ass_coll_index=sort([pot_ass_index[1:npot_ass];pot_coll_index[1:npot_coll]])
            #in_poly=zeros(Bool,size(pot_ass_coll_index,1))
            #in_circle=zeros(Bool,size(pot_ass_coll_index,1))
            for izone in uass
                #move positioner to this target
                move_positioner_to_theta_phi!(pos_cur,zone_target.theta_pos[izone],zone_target.phi_pos[izone],
                    diff_v=move_diff_v,cos_sin=move_cos_sin,cos_sin_sum=move_cos_sin_sum,
                    cstheta=move_cstheta,csphi=move_csphi)
                
                fill!(in_poly,false)
                fill!(in_circle,false)
                
                #determine the collisions and potential collision with body
                PointsInPoly_and_PaddedCircle!(view(zone_target.X_FP,pot_ass_coll_index[1:npot_ass_coll]),
                   view(zone_target.Y_FP,pot_ass_coll_index[1:npot_ass_coll]),
                    view(pos_cur.body_on_target,1:pos_cur.nseg[1],:),
                    pos_cur.head_minmax[2];max_rad=pos_cur.body_minmax[2],
                    in_poly=in_poly,in_circle=in_circle)
                
                #determine the collisions and potential collision with head
                PointsInPoly_and_PaddedCircle!(view(zone_target.X_FP,pot_ass_coll_index[1:npot_ass_coll]),
                   view(zone_target.Y_FP,pot_ass_coll_index[1:npot_ass_coll]),
                    view(pos_cur.head_on_target,1:pos_cur.nseg[2],:),
                    pos_cur.body_minmax[2];max_rad=pos_cur.head_minmax[2],
                    in_poly=in_poly,in_circle=in_circle)


                #ind_ri=findall(local_ass .== izone)
                nri=find_index!(local_ass .== izone,ind_ri)
                
                collided[zone_target.index[pot_ass_coll_index[in_poly]],ind_ri[1:nri]] .= true
                potential_coll[zone_target.index[pot_ass_coll_index[in_circle]],ind_ri[1:nri]] .= true
            end
            #@show position_assigned[pp,:],targets_dic["PRIORITY"][position_assigned[pp,:]]
            #tmp=tmp+1
            if(false)  
                #plot_moved_positioner(pos_cur)
                #plot_zones(pos_cur,zone_target,fp_dic)
                ri=2
                nc=zone_target.ntarget[1]
                plot(zone_target.X_FP[1:nc],zone_target.Y_FP[1:nc],seriestype=:scatter,color=:black)
                
                plot!([targets_dic["X_FP"][position_assigned[pp,ri]]],
                    [targets_dic["Y_FP"][position_assigned[pp,ri]]],seriestype=:scatter,color=:red)
                
                plot_circle(pos_cur.xy[1],pos_cur.xy[2],sqrt(pos_cur.arms2_lim[2]),color=:red)
                plot_circle(pos_cur.xy[1],pos_cur.xy[2],sqrt(coll_r2),color=:blue) 
                
                #plot the collided
                ipot_col=findall(potential_coll[:,ri])
                icol=findall(collided[:,ri])
                @show icol,targets_dic["X_FP"][icol]
                plot!(targets_dic["X_FP"][icol],targets_dic["Y_FP"][icol],seriestype=:scatter,color=:blue)
                plot!(targets_dic["X_FP"][ipot_col],targets_dic["Y_FP"][ipot_col],seriestype=:scatter,color=:green)
                #|> display
                
                move_positioner_to_theta_phi!(pos_cur,theta_assigned[pp,ri] ,phi_assigned[pp,ri],
                    diff_v=move_diff_v,cos_sin=move_cos_sin,cos_sin_sum=move_cos_sin_sum,
                    cstheta=move_cstheta,csphi=move_csphi)
                                    
                plot_moved_positioner(pos_cur;linecolor=:black,bodycolor=:blue,headcolor=:red,alpha=0.2,hold=true) |> display
            end
            
            
        end
        #break
    end
    
    if(config["OUTPUT"]["ThetaPhi_pos"])
        return position_assigned,theta_assigned,phi_assigned
    else
        return position_assigned,zeros(Float64,2),zeros(Float64,2)
    end
end
    
   
"""
takes the assignment, splits in different tracers and write a group in the JLD2 file for each tile
This is expected to be post-processed after each pass

first scans the file_map dictionary in targets_dic
for each key in the dictionary
    select all the object with non-zero assignments and write them to file
"""
function TileAssignment_To_JLD2File(tile_id,config,num_fba,fp_dic,targets_dic,position_assigned;
        sum_ass_global=zeros(Int16,2),ithis=nothing,target_assignment_global=nothing,
    location_assignment_global=nothing)
    
    ntarget=size(targets_dic["X_FP"],1)
    npos=size(position_assigned,1)
    
    #declare an array to transfer the assignment
    if(size(sum_ass_global,1)<ntarget)
        sum_ass=zeros(Int16,ntarget)
    else
        sum_ass=view(sum_ass_global,1:ntarget)
        fill!(sum_ass,0)
    end
    
    if(ithis==nothing)
        ithis=zeros(Int32,ntarget)
    end
    
    if(target_assignment_global==nothing)
        target_assignment=zeros(Bool,ntarget,num_fba)
        location_assignment=zeros(Int16,ntarget,num_fba) .- 1
    else
        target_assignment=view(target_assignment_global,1:ntarget,:)
        location_assignment=view(location_assignment_global,1:ntarget,:)
        fill!(target_assignment,false)
        fill!(location_assignment,-1)
    end
    
    
    
    for ri in 1:num_fba
        #println(ri,' ',maximum(position_assigned[:,ri]),' ',ntarget)
        @inbounds for ii in 1:npos
            if(position_assigned[ii,ri]>0)
                target_assignment[position_assigned[ii,ri],ri] = true
                location_assignment[position_assigned[ii,ri],ri]=fp_dic["LOCATION"][ii]
                sum_ass[position_assigned[ii,ri]] += 1
            end
        end
    end
    
    #histogram(sum_ass[sum_ass .>0],alpha=0.1)
    msg_out=""
    for tkey in keys(targets_dic["file_map"])
        tracer=targets_dic["file_map"][tkey][1]
        nthis=find_index_AND!((targets_dic["filekey"].==tkey),(sum_ass .>0),ithis)
        #println(tkey,targets_dic["file_map"][tkey],size(ithis,1))
        
        #histogram!(sum_ass[ithis],label=targets_dic["file_map"][tkey][1],alpha=0.3)
        
        #need to collect the assignment array and index array
        #index_arr=targets_dic["index"][ithis]
        #assignment_bool=target_assignment[ithis,:]
        
        #sort by index_arr
        isort=sortperm(targets_dic["index"][ithis[1:nthis]])
        #index_arr=index_arr[isort]
        #assignment_bool=assignment_bool[isort,:]
        
        #Now plot for testing
        #plot!(targets_dic["X_FP"][ithis][isort],targets_dic["Y_FP"][ithis][isort],
        #    seriestype=:scatter,label=targets_dic["file_map"][tkey][1],color=:auto,alpha=0.5) #|> display
        
        #get the filename to write
        fba_jldfile=TileFBA_FileName(config,tile_id)
	#fba_jldfile=config["target"][tracer]["FBA_JLDfile"]

        jldopen(fba_jldfile, "a+") do file
	   write(file,"$(tracer)_Bool_ass",target_assignment[ithis[isort],:])
	   write(file,"$(tracer)_location_ass",location_assignment[ithis[isort],:])
	   write(file,"$(tracer)_index",targets_dic["index"][ithis[isort]])
        end
       
	#number of assignment for the first realization
	msg_out=string(msg_out,tracer,":",sum(target_assignment[ithis[isort],1])," ")
        #println("written: ",fba_jldfile)
    end
    #println("sum: ",sum(sum_ass .==0),' ',minimum(sum_ass),' ',maximum(sum_ass))
    return msg_out
end

"""
Generate the file_name for a tile_output
"""
function TileFBA_FileName(config,tile_id)
    outdir=string(config["OUTPUT"]["JLD2_dir"],config["OUTPUT"]["FBA-tag"],"/")
    outfile=string(outdir,"TILEID_",tile_id,".jld2")
    return outfile
end

"""
file to write the fits output of each zone including the information about FBA
"""
function ZoneFITS_FileName(config,zone,tracer,group)
    outdir=string(config["OUTPUT"]["FITS_dir"],config["OUTPUT"]["FBA-tag"],"/tmp/")
    outfile=string(outdir,group,"_",tracer,"_zone_",zone,".fits.gz")
    return outfile
end

"""Runs the full assignmnet from disk-to-disk for a single tile
tile_index : index if the tile in tiles_dic to be assigned
verbose sets the verbosity of returned message
verbose=0 empty msg
verbos=1 overall time of the function
verbose>1 time for all steps
"""
function Run_Single_Tile(config,tile_index,tiles_dic;group="",tile_date="2019-09-16T00:00:00",
    plate_scale=nothing,fp_dic=nothing,exc_dic=nothing,verbose=0)
    
    events_dic=Dict()
    events_dic["events"]=[]
    events_dic["times"]=[]
    
    
    tile_id=tiles_dic["TILEID"][tile_index]
    if(verbose>0)
        append!(events_dic["times"],[now(UTC)])
        append!(events_dic["events"],["Time"])
    end
    
    #load data for the tile and focal plane
    targets_dic,fp_dic,exc_dic=Prepare_tile(config,tile_index,tiles_dic;group=group,tile_date=tile_date,
        plate_scale=plate_scale,fp_dic=fp_dic,exc_dic=exc_dic,num_fba=config["NumFBARealization"])
    if(verbose>1)
        append!(events_dic["times"],[now(UTC)])
        append!(events_dic["events"],["Prep"])
	#count the number of objects need assignment
	need_ass=targets_dic["ntarget"]-sum(targets_dic["collided"][:,1])
    end
   
    #println("Run Single: ",tile_id,' ',config["NumFBARealization"],size(targets_dic["collided"]))

    #run the assignment
    position_assigned,theta_assigned,phi_assigned=Assign_Tile(config,fp_dic,exc_dic,targets_dic,
        num_fba=config["NumFBARealization"],max_vertices=30,max_ngb_pos=9)
    
    if(verbose>1)
        append!(events_dic["times"],[now(UTC)])
        append!(events_dic["events"],["Ass"])
    end
    
    #write the data to file
    nass_str=TileAssignment_To_JLD2File(tile_id,config,config["NumFBARealization"],fp_dic,targets_dic,position_assigned)
    
    if(verbose>0)
        append!(events_dic["times"],[now(UTC)])
        if(verbose>1)
            append!(events_dic["events"],["Out"])
        end
    end

    msg_out=string(tile_id)
    #count the number of targets in tile vs number need assignment for the first realization
    if(verbose>1)
       msg_out=string(msg_out," targ:",need_ass,"/",targets_dic["ntarget"]," ",nass_str)
    end

    return convert_events_to_msg(events_dic,msg_out)
end

"""Runs the full assignmnet from disk-to-disk for a array of tiles
tile_index_arr : Array of indices if the tile in tiles_dic to be assigned
"""
function Run_Many_Tile(config,tile_index_arr,tiles_dic;group="",tile_date="2019-09-16T00:00:00",
    plate_scale=nothing,fp_dic=nothing,exc_dic=nothing,verbose=0,pre_msg="")
    
    #load plate_scale if needed
    if(plate_scale==nothing)
        plate_scale=load_platescale(config["focal_plane"]["focalplane_dir"])
    end
    
    #load the data for first tile to use things multiple times
    #load data for the tile and focal plane
    tile_index=tile_index_arr[1]
    targets_dic,fp_dic,exc_dic=Prepare_tile(config,tile_index,tiles_dic;group=group,tile_date=tile_date,
        plate_scale=plate_scale,fp_dic=fp_dic,exc_dic=exc_dic,num_fba=config["NumFBARealization"])
    
    #Now iterate over rest of the tile
    for tile_index in tile_index_arr[1:end]
        #output file name for this tile
        fba_jldfile=TileFBA_FileName(config,tiles_dic["TILEID"][tile_index])
	#check if this file exists then continue
	if(isfile(fba_jldfile))
	   println(pre_msg," TILEID:",tiles_dic["TILEID"][tile_index], " exists")
	   continue
	end

        #println(now(UTC)," Begin TILEID:",tiles_dic["TILEID"][tile_index])
        ret_msg=Run_Single_Tile(config,tile_index,tiles_dic;group=group,tile_date=tile_date,
        plate_scale=plate_scale,fp_dic=fp_dic,exc_dic=exc_dic,verbose=verbose)
        if(verbose>0)
            #println(now(UTC),pre_msg," TILEID:",ret_msg)
            println(pre_msg," TILEID:",ret_msg)
        end
    end        
end
 
