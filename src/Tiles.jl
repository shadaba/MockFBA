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
        if(isa(tdic[tkey],Array))
            #this_ptr=pointer_from_objref(tdic[tkey])
            this_ptr=repr(UInt64(pointer_from_objref(tdic[tkey])))
            if((size(tdic[tkey],1)==nxp) & (!(this_ptr in modified_ptr)) )
                tdic[tkey]=tdic[tkey][izone_sort]
                append!(modified_ptr,[this_ptr])
                #println("sorting $(tkey) $(this_ptr)")
            else
                println("not sorting ",tkey,size(tdic[tkey],1), this_ptr)
            end
        end
    end
    
    #Now generate the map from of indices for each zone
    #this_ptr=pointer_from_objref(zone_xy)
    this_ptr=repr(UInt64(pointer_from_objref(zone_xy)))
    if(!(this_ptr in modified_ptr) )
        zone_xy=zone_xy[izone_sort]
        append!(modified_ptr,[this_ptr])
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
function load_targets_intile(config_targ,tile_id,tile_pass,group;columns=["X_UNIT","Y_UNIT","Z_UNIT","PRIORITY"])
    targets_dic=Dict()
    targets_dic["file_map"]=Dict()
    
    col_accumulate=["index","filekey"]
    for col in columns
        if(!(col in col_accumulate))
            append!(col_accumulate,[col])
        end
    end
    
    targets_dic["filekey"]=zeros(Int8,0)
    
    
  
    fno=1
    for tracer in keys(config_targ)
        fname=config_targ[tracer]["JLDfile" ]
        targets_dic["file_map"][fno]=fname
        
        tracer_tile=MockFBA.Load_tracers_intile(fname,tile_id,tile_pass,group,columns)
        
        nobj=size(tracer_tile["index"],1)
        for col in col_accumulate
            if(col=="filekey")
                targets_dic[col]=append!(targets_dic[col],(zeros(Int8,nobj) .+ fno))
            elseif(col in keys(targets_dic))
                targets_dic[col]=append!(targets_dic[col],tracer_tile[col])
            else
                targets_dic[col]=tracer_tile[col]
            end
        end
        println(tracer,' ',unique(targets_dic["PRIORITY"]),nobj,fno)
        fno +=1
    end
    targets_dic["ntarget"]=size(targets_dic["index"],1)
    return targets_dic
end


"""Runs the assignment for a given tile
num_fba: number of fba realization to generate
max_vertices maximum number of vertices any of the positioner polygon can have
max_ngb_pos: maximum number of negibohuring positioner any positioner can have
"""
function Assign_Tile(fp_dic,targets_dic;num_fba=3,max_vertices=30,max_ngb_pos=9)
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
    
    #boolean for collision
    collided=zeros(Bool,targets_dic["ntarget"],num_fba)
    #integer for potential collision
    potential_coll=zeros(Bool,targets_dic["ntarget"],num_fba)
    
    #dictionary for collision locations
    #integer for Assignment
    position_assigned=zeros(Int32,size(fp_dic["OFFSET_X"],1),num_fba)
    theta_assigned=zeros(Float64,size(fp_dic["OFFSET_X"],1),num_fba)
    phi_assigned=zeros(Float64,size(fp_dic["OFFSET_X"],1),num_fba)
    
    local_ass=zeros(Int64,num_fba) #To keep the assignent in local index
    
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

            #priority_sel=Dict()
            #for up in us_priority
            #     priority_sel[up]=findall((zone_target.PRIORITY[1:zone_target.ntarget[1]].==up) .& 
            #        zone_target.pot_ass[1:zone_target.ntarget[1]])
            #end
            #@show us_priority,priority_sel
            
            #possible_ass=zeros(Int64,100) #To store the actual possiblity
            
            pot_ass_index=findall(zone_target.pot_ass[1:zone_target.ntarget[1]])
            pot_coll_index=findall(zone_target.pot_coll[1:zone_target.ntarget[1]])
            
            coll_map_dic=Dict()
            #loop over realizations
            for ri in 1:num_fba
                count=0
                for uu in 1:zone_target.ntarget[2] #do this in priority order high to low  
                    this_index=findall( zone_target.usp_bool_sel[1:zone_target.ntarget[1],uu] .& 
                        zone_target.pot_ass[1:zone_target.ntarget[1]])
                    #println(' ',uu,' ',size(this_index),' ',
                    #    size(findall(zone_target.usp_bool_sel[1:zone_target.ntarget[1],uu])) )
                    for izone in pot_ass_index #scan all the object in potential assignment
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
                                    move_positioner_to_theta_phi!(pos_ngb_dic["pos"][ingb],
                                        theta_assigned[pos_ngb_dic["pos_index"][ingb],ri] ,
                                        theta_assigned[pos_ngb_dic["pos_index"][ingb],ri])
                                    
                                    #move current positioner to this target
                                    move_positioner_to_theta_phi!(pos_cur,
                                        zone_target.theta_pos[izone] ,zone_target.phi_pos[izone]) 
                                    
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
                        theta_assigned[pp,ri]= zone_target.theta_pos[izone_sel] 
                        phi_assigned[pp,ri]= zone_target.phi_pos[izone_sel] 
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
            
            pot_ass_coll_index=sort([pot_ass_index;pot_coll_index])
            in_poly=zeros(Bool,size(pot_ass_coll_index,1))
            in_circle=zeros(Bool,size(pot_ass_coll_index,1))
            for izone in uass
                #move positioner to this target
                move_positioner_to_theta_phi!(pos_cur,zone_target.theta_pos[izone],zone_target.phi_pos[izone])
                
                fill!(in_poly,false)
                fill!(in_circle,false)
                
                #determine the collisions and potential collision with body
                PointsInPoly_and_PaddedCircle!(zone_target.X_FP[pot_ass_coll_index],
                   zone_target.Y_FP[pot_ass_coll_index],pos_cur.body_on_target[1:pos_cur.nseg[1],:],
                    pos_cur.head_minmax[2];max_rad=pos_cur.body_minmax[2],
                    in_poly=in_poly,in_circle=in_circle)
                
                #determine the collisions and potential collision with head
                PointsInPoly_and_PaddedCircle!(zone_target.X_FP[pot_ass_coll_index],
                   zone_target.Y_FP[pot_ass_coll_index],pos_cur.head_on_target[1:pos_cur.nseg[2],:],
                    pos_cur.body_minmax[2];max_rad=pos_cur.head_minmax[2],
                    in_poly=in_poly,in_circle=in_circle)


                ind_ri=findall(local_ass .== izone)
                
                collided[zone_target.index[pot_ass_coll_index[in_poly]],ind_ri] .= true
                potential_coll[zone_target.index[pot_ass_coll_index[in_circle]],ind_ri] .= true
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
                
                move_positioner_to_theta_phi!(pos_cur,
                                        theta_assigned[pp,ri] ,phi_assigned[pp,ri]) 
                                    
                plot_moved_positioner(pos_cur;linecolor=:black,bodycolor=:blue,headcolor=:red,alpha=0.2,hold=true) |> display
            end
            
            
        end
        #break
    end
    return position_assigned,theta_assigned,phi_assigned
end
    
    

#num_fba=100
#position_assigned,theta_assigned,phi_assigned=Assign_Tile(fp_dic,targets_dic;num_fba=num_fba,max_vertices=30,max_ngb_pos=9)

#print("Finished")
