#Contain the Zone_target structure to load hold and manipulate the targets in a zone
"""loads the targets in the zone and keep track of whether a newzone in needed to be loaded
"""
struct Zone_Target
    #Loads the targets relvant for a particular zone
    #Any thing stored in two element array will always asume body first and head second
    max_target::Integer #maximum number of target
    max_usp::Integer #total maximum number of unique priority for the full tile
    
    zones::Array{Integer,1}
    #first entry for actual number of target and second entry for number of unique priority
    ntarget::Array{Integer,1} 
    
    
    index::Array{Integer,1}
    X_FP::Array{Float64,1} #To store the x-position on focal plane of target
    Y_FP::Array{Float64,1} #To store y-position on focal plane of target
    PRIORITY::Array{Integer,1}
    
    #boolean array
    pot_ass::Array{Bool,1} #To set for the potential assignment
    pot_coll::Array{Bool,1} #To set for potential collision 
    
    #To store theta,phi for a particular positioner
    theta_pos::Array{Float64,1}
    phi_pos::Array{Float64,1}
    
    #unique sorted priority for this tile
    usp::Array{Integer,1}
    usp_bool_sel::Array{Integer,2} #To store the boolean selection of each priority
    
    #To store the possible assignments
    possible_ass::Array{Int32,1}
    
    function  Zone_Target(maximum_target::Integer,nunique_priority::Integer)
        new(maximum_target,nunique_priority,zeros(3),zeros(2),zeros(maximum_target),zeros(maximum_target),zeros(maximum_target),
        zeros(maximum_target),zeros(Bool,maximum_target),zeros(Bool,maximum_target),
        zeros(maximum_target),zeros(maximum_target),
        zeros(nunique_priority),zeros(Bool,maximum_target,nunique_priority),
        zeros(Int32,maximum_target))
    end
end


function load_target_zone(zone_target::Zone_Target,zone_x,zone_y,fp_bound,targets_dic)
    #This zone is already loaded so nothing to be done
    if((zone_target.zones[1]==zone_x) & (zone_target.zones[2]==zone_y))
        return
    else
        zone_target.zones[1]=zone_x
        zone_target.zones[2]=zone_y
        zone_target.zones[3]=get_zonexy(zone_x,zone_y,fp_bound)
    end
    
    #indsel=zeros(Int32,1000)
    
    count_targ=1
    for zx in maximum([zone_x-1,1]):minimum([zone_x+1,fp_bound.nx])
        for zy in maximum([zone_y-1,1]):minimum([zone_y+1,fp_bound.ny])
            zone_xy=get_zonexy(zx,zy,fp_bound)
            
            if(!(zone_xy in targets_dic["zone_map_keys"]))
                continue
            end
            imin=targets_dic["zone_map"][zone_xy][1]
            imax=targets_dic["zone_map"][zone_xy][2]
            #@show zx,zy,zone_xy,targets_dic["zone_map"][zone_xy]
            #@show targets_dic["zone_x"][imin:imax]
            count_this=imax-imin
            zone_target.index[count_targ:count_targ+count_this] .= collect(imin:imax)
            count_targ=count_targ+count_this+1
        end
    end
    
    #set everything else to 0
    if(count_targ<=zone_target.max_target)
        zone_target.index[count_targ:end] .= 0
    end
    #Due to Julia indexing convetion 1 need to be subtracte
    count_targ=count_targ-1
    
    tmp_ind=view(zone_target.index,1:count_targ)
    #sort the indices
    sort!(tmp_ind)
    
    zone_target.ntarget[1]=count_targ
    
    zone_target.X_FP[1:count_targ]=targets_dic["X_FP"][tmp_ind]
    zone_target.Y_FP[1:count_targ]=targets_dic["Y_FP"][tmp_ind]
    zone_target.PRIORITY[1:count_targ]=targets_dic["PRIORITY"][tmp_ind]
    
    us_priority=sort(unique(zone_target.PRIORITY[1:zone_target.ntarget[1]]),rev=true)
    nup=size(us_priority,1)
    zone_target.ntarget[2]=nup
    zone_target.usp[1:nup] .= us_priority
    

    imin=1;imax=zone_target.ntarget[1]
    for uu in 1:zone_target.ntarget[2]
        zone_target.usp_bool_sel[imin:imax,uu] .= (zone_target.PRIORITY[imin:imax].== us_priority[uu])
    end

end
