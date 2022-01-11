#Function to deal with positioners
struct POSRobotBody
    #positioner body is theta arm and head is for phi arm
    #Any thing stored in two element array will always asume body first and head second
    max_seg::Integer #maximum number of segements this object can store
    
    xy::Array{Float64,1} #To store the xy position of the positioner
    arms::Array{Float64,1} #To store the theta.phi arm of the positioner
    arms2_lim::Array{Float64,1} #To store the square of minimum and maximum limit of the robot arm
    
    nseg::Array{Integer,1} #Stores the number of segment in body and head
    restricted_angels::Array{Bool,1} #Store if any angular restriction is needed
    
    body_minmax::Array{Float64,1} #Stores the minimum and maximum length of the verted for body
    body_orig::Array{Float64, 2} #Stores the orginal polygon of body
    head_minmax::Array{Float64,1} #Stores the minimum and maximum length of the verted for head
    head_orig::Array{Float64, 2} #Stores the original polygon of head
    
    axis_orig::Array{Float64,1} #To store the orignal axis
    
    #To get the theta limit and phi limit: Angles are stored in radian
    body_theta_limit::Array{Float64,1} #Store the limit of angles for body MIN_T,MAX_T,OFFSET_T
    head_theta_limit::Array{Float64,1} #Store the limit of angles for head MIN_P,MAX_P,OFFSET_P
    
    body_on_target::Array{Float64, 2} #Stores the body polygon on target
    head_on_target::Array{Float64, 2} #Stores the head polygon on target
    axis_on_target::Array{Float64,1} #This goes throgh the same transformation as head to keep tract of its axis
    

    function  POSRobotBody(maximum_seg::Integer)
        new(maximum_seg,zeros(2),zeros(2),zeros(2),zeros(2),
            zeros(Bool,2),zeros(2),zeros(maximum_seg,2),zeros(2),
            zeros(maximum_seg,2),zeros(2),zeros(3),zeros(3),
            zeros(maximum_seg,2),zeros(maximum_seg,2),zeros(2))
    end
end   


"""
Initialize the pos:POSRobotBody with the position at index in fp_dic
"""
function load_positioner!(pos::POSRobotBody,index,fp_dic,exc_dic)
    
    #get the exclusion names
    excl_name=fp_dic["EXCLUSION"][index]
    excl_body_name="$(excl_name)-theta-segments"
    excl_head_name="$(excl_name)-phi-segments"
    
    #Assign number of segments
    pos.nseg[1]=size(exc_dic[excl_body_name],1)
    pos.nseg[2]=size(exc_dic[excl_head_name],1)
    
    #Assign polygon for body
    Assign_polygon_with_maxmin_vertex!(pos.max_seg,pos.nseg[1],
        exc_dic[excl_body_name],pos.body_orig,pos.body_minmax)
    
    #Assign polygon for head
    Assign_polygon_with_maxmin_vertex!(pos.max_seg,pos.nseg[2],
        exc_dic[excl_head_name],pos.head_orig,pos.head_minmax)
    
    #get the angular limit of the robots body
    pos.body_theta_limit[1]=deg2rad(fp_dic["MIN_T"][index])
    pos.body_theta_limit[2]=deg2rad(fp_dic["MAX_T"][index])
    pos.body_theta_limit[3]=deg2rad(fp_dic["OFFSET_T"][index])
    
    #get the angular limit of the robots head
    pos.head_theta_limit[1]=deg2rad(fp_dic["MIN_P"][index])
    pos.head_theta_limit[2]=deg2rad(fp_dic["MAX_P"][index])
    pos.head_theta_limit[3]=deg2rad(fp_dic["OFFSET_P"][index])
    
    #Allocate xy
    pos.xy[1]=fp_dic["OFFSET_X"][index]
    pos.xy[2]=fp_dic["OFFSET_Y"][index]
    
    #Assign arms
    pos.arms[1]=fp_dic["LENGTH_R1"][index]
    pos.arms[2]=fp_dic["LENGTH_R2"][index]
    
    #arms limits
    pos.arms2_lim[1]=(pos.arms[1]-pos.arms[2])^2 #minimum
    pos.arms2_lim[2]=(pos.arms[1]+pos.arms[2])^2 #maximum
    
end


function Assign_polygon_with_maxmin_vertex!(max_seg,nseg,polygon,pos_poly,minmax)
    cent=[mean(polygon[:,1]),mean(polygon[:,2])]
    min_vert2=100000
    max_vert2=0
    
    for ii in 1:nseg
        pos_poly[ii,1]=polygon[ii,1]
        pos_poly[ii,2]=polygon[ii,2]
        rad2=(polygon[ii,1]-cent[1])^2 + (polygon[ii,2]-cent[2])^2
        if(max_vert2<rad2)
            max_vert2=rad2
        end
        if(min_vert2>rad2)
            min_vert2=rad2
        end
    end
    
    for ii in nseg+1:max_seg
        pos_poly[ii,1]= 0
        pos_poly[ii,2]= 0
    end
    
    minmax[1]=sqrt(min_vert2)
    minmax[2]=sqrt(max_vert2) 
end

"""
Applies the tanslation on poly_in by x,y and return result in poly_out
"""
function poly_translate!(poly_in,poly_out,x,y)
    for ii in 1:size(poly_in,1)
        poly_out[ii,1]=poly_in[ii,1]+x
        poly_out[ii,2]=poly_in[ii,2]+y
    end 
end

function point_translate!(pt_in,pt_out,x,y)
    pt_out[1]= pt_in[1] + x
    pt_out[2]= pt_in[2] + y
end

"""Rotates the points in poly around the axis
by cos(theta),sin(theta)
"""
function poly_rotate_axis!(poly_in,poly_out,cstheta,axis)
    for ii in 1:size(poly_in,1)
        poly_out[ii,:] .=point_rotate_axis(poly_in[ii,:],cstheta,axis)
    end  
end


function point_rotate_axis(point_in,cstheta,axis)
    diff_v=point_in .- axis
    if((diff_v[1]==0) & (diff_v[2]==0))
        return point_in[1],point_in[2]
    end

    rad=sqrt(sum(diff_v.^2))

    #cos_sin of difference vector
    cos_sin=diff_v ./rad

    #sum of angles cos(A+B) and sin (A+B) rules
    cos_sin_sum=[cos_sin[1]*cstheta[1] - cos_sin[2]*cstheta[2],
                 cos_sin[1]*cstheta[2] + cos_sin[2]*cstheta[1]  ]

    #rotation and translation
    return [axis[1]+(rad*cos_sin_sum[1]), axis[2]+(rad*cos_sin_sum[2])]
end



"""moves the positioner to a theta phi position
Adapted from desihub/fiberassign
"""
function move_positioner_to_theta_phi!(pos::POSRobotBody,theta,phi)
    cstheta = [cos(theta),sin(theta)]
    csphi = [cos(phi),sin(phi)]
    

    #// Move the phi polygon into the fully extended position along the X axis.
    #moving to LENGHT_R1,0
    nT=pos.nseg[1]
    nP=pos.nseg[2]
    
    poly_translate!(view(pos.head_orig,1:nP,:),view(pos.head_on_target,1:nP,:),pos.arms[1],0) 
    #translate the axis
    point_translate!(pos.axis_orig,pos.axis_on_target,pos.arms[1],0)
    
    #@show pos.head_on_target[1:pos.nseg[2],2] .== (pos.head_orig[1:pos.nseg[2],2] .+ 0)
    #@show "1: ",pos.axis_on_target,pos.arms[1]

    #// Rotate fully extended positioner an angle of theta about the origin.
    poly_rotate_axis!(view(pos.body_orig,1:nT,:),view(pos.body_on_target,1:nT,:),cstheta,[0.0,0.0])
    poly_rotate_axis!(view(pos.head_on_target,1:nP,:),view(pos.head_on_target,1:nP,:),cstheta,[0.0,0.0])
    #rotate the axis
    pos.axis_on_target .= point_rotate_axis(pos.axis_on_target,cstheta,[0.0,0.0])
    
    #@show "2: ",pos.axis_on_target,cstheta

    #// Rotate just the phi arm an angle phi about the theta arm center.
    #shpphi.rotation(csphi);
    poly_rotate_axis!(view(pos.head_on_target,1:nP,:),view(pos.head_on_target,1:nP,:),csphi,pos.axis_on_target)

  
    #// Translate the whole positioner to the center.
    poly_translate!(view(pos.head_on_target,1:nP,:),view(pos.head_on_target,1:nP,:),pos.xy[1],pos.xy[2])
    poly_translate!(view(pos.body_on_target,1:nT,:),view(pos.body_on_target,1:nT,:),pos.xy[1],pos.xy[2]) 
    
    #Not needed just for plotting
    #point_translate!(pos.axis_on_target,pos.axis_on_target,pos.xy[1],pos.xy[2])
end

"""Determines if a positioner is valid of not
fp_dic : dictionary contaning the status of focal plane
pp: index of positioner considered
pos_type: purpose of validity testing
           assignment: whether this positioner can be assigned to a target
           neighbour_collision: whether this positioner should be considered for neighbour collision
"""
function valid_positionter(fp_dic,pp;pos_type="assignment")
    if(pos_type=="assignment")
        return ((fp_dic["DEVICE_TYPE"][pp]=="POS") & (fp_dic["LENGTH_R1"][pp]!=0.0))
    elseif(pos_type=="neighbour_collision")
        return ((fp_dic["DEVICE_TYPE"][pp]=="POS") & (fp_dic["LENGTH_R1"][pp]!=0.0))
    end
end

"""
loads the positioner negibours given pindex and zones of the central positioner
pos_ngb_dic: a dictionary of structure to hold the negibours
"""
function load_positioner_neighbours!(pos_cur,pindex,fp_bound,fp_dic,exc_dic,pos_ngb_dic)
    #This zone is already loaded so nothing to be done
    if(pos_ngb_dic["index_central_positioner"]==pindex)
        #already this positioner neighbours are loaded
        return
    else
        zone_x=fp_dic["zone_x"][pindex]
        zone_y=fp_dic["zone_y"][pindex]
        zone_xy=get_zonexy(zone_x,zone_y,fp_bound)
        pos_ngb_dic["index_central_positioner"]=pindex
    end
    
    #@show zone_x,zone_y,zone_xy
    
    index
    count_pos=1
    for zx in maximum([zone_x-1,1]):minimum([zone_x+1,fp_bound.nx])
        for zy in maximum([zone_y-1,1]):minimum([zone_y+1,fp_bound.ny])
            zone_xy=get_zonexy(zx,zy,fp_bound)
            
            #@show zx,zy,zone_xy
            
            if(!(zone_xy in fp_dic["zone_map_keys"]))
                continue
            end
            imin=fp_dic["zone_map"][zone_xy][1]
            imax=fp_dic["zone_map"][zone_xy][2]
            #@show zx,zy,zone_xy,targets_dic["zone_map"][zone_xy]
            #@show targets_dic["zone_x"][imin:imax]
            count_this=imax-imin
            pos_ngb_dic["index"][count_pos:count_pos+count_this] .= collect(imin:imax)
            count_pos=count_pos+count_this+1
        end
    end
    
    
    
    #set everything else to 0
    if(count_pos<=size(pos_ngb_dic["index"],1)) 
        pos_ngb_dic["index"][count_pos:end] .= 0
    end
    #Due to Julia indexing convetion 1 need to be subtracte
    count_pos=count_pos-1
    
    
    #Determine the maximum distance to look for within positioner
    #This is defined as the maximum arams length + the maximum vertex of head
    #The above is conservative and specific choice will not matter
    lookup_rad2=(pos_cur.arms[1]+pos_cur.arms[2]+pos_cur.head_minmax[1])^2
    lookup_rad2=4*lookup_rad2 #factor of 4 is for twice quare the distance
    #@show lookup_rad2,pos_cur.arms[1],pos_cur.arms[2],pos_cur.head_minmax[2]
    
    #scan for positioners within look_up radius
    count_ngb=0
    for ii in pos_ngb_dic["index"][1:count_pos]
        if((ii==pindex) || (!valid_positionter(fp_dic,ii;pos_type="neighbour_collision")))#To not count self positioner and invalid positioner
            continue
        end
        rad2=(fp_dic["OFFSET_X"][ii]-fp_dic["OFFSET_X"][pindex])^2
        rad2 += (fp_dic["OFFSET_Y"][ii]-fp_dic["OFFSET_Y"][pindex])^2
        #@show ii,rad2,lookup_rad2,rad2<lookup_rad2
        
        if(rad2<lookup_rad2)
            count_ngb+=1
            #load this positioner info
            load_positioner!(pos_ngb_dic["pos"][count_ngb],ii,fp_dic,exc_dic)
            pos_ngb_dic["pos_index"][count_ngb]=ii
        end
    end
    pos_ngb_dic["count_ngb"]=count_ngb
    #@show count_ngb
    
    
    if(false)#To make plots
        tind=pos_ngb_dic["index"][1:count_pos]
        plot([fp_dic["OFFSET_X"][pindex]],[fp_dic["OFFSET_Y"][pindex]],
            seriestype=:scatter,color=:green,markersize=10) 
        for ii in 1:count_ngb
            xy=pos_ngb_dic["pos"][ii].xy
            plot!([xy[1]],[xy[2]],seriestype=:scatter,color=:red,markersize=10)
        end

        plot!(fp_dic["OFFSET_X"][tind],fp_dic["OFFSET_Y"][tind],seriestype=:scatter,color=:black) |> display
    end
    
end



"""Compare the moved positioner positioner to check if either the body or head-boy is collided
body-body collision is not considered as it is not possible in the focal plane
pos1: first positioner of type POSRobotBody
pos2: second positioner of type POSRobotBody
"""
function positioner_collided(pos1,pos2)
    
    head1=view(pos1.head_on_target,1:pos1.nseg[2],:)
    head2=view(pos2.head_on_target,1:pos2.nseg[2],:)
    
    #check for the head1-head2 collision
    if(PolyIntersect_boundingbox(head1,head2))
        if(PolygonsIntersect(head1,head2))
            return true
        end
    end
    
    body1=view(pos1.body_on_target,1:pos1.nseg[1],:)
    body2=view(pos2.body_on_target,1:pos2.nseg[1],:)
    
    #head1-body2 collision
    if(PolyIntersect_boundingbox(head1,body2))
        if(PolygonsIntersect(head1,body2))
            return true
        end
    end
    
    #head2-body1 collision
    if(PolyIntersect_boundingbox(head2,body1))
        if(PolygonsIntersect(head2,body1))
            return true
        end
    end
    
    return false
end

