#contains geomtry related functions
"""point in polygon
"""
function pointxy_inpoly_boundingbox(x,y,poly)
    poly_min=minimum(seg,dims=1)
    poly_max=maximum(seg,dims=1)
    ind=(x .> poly_max[1]) .|| (x .< poly_min[1]) .|| (y.>poly_max[2]) .|| (y.<poly_min[2])
end

"""compute the distance of xy from center and return boolean array for each radius limit
    radius is an array of radius which must be sorted

    return
"""
function point_in_circle!(x,y;x_cent::Number=0,y_cent::Number=0,radius1=0,radius2=0,radius3=0,
       dist2=nothing,sel_rad1=nothing,sel_rad2=nothing,sel_rad3=nothing)
    
    nobj=size(x,1)
    if(dist2==nothing)
        dist2=zeros(nobj)
    end
    
    @inbounds for ii in 1:nobj
        dist2[ii] = (x[ii]-x_cent)^2 + (y[ii]-y_cent)^2
    end
    
    
    if(radius1==0)
        return
    else
        @. sel_rad1 = dist2 < (radius1*radius1)
    end
    
    if(radius2==0)
        return
    else
        @. sel_rad2 = (!sel_rad1) & (dist2 < (radius2*radius2))
    end
    
    if(radius3==0)
        return
    else
        @. sel_rad3 = (!(sel_rad2 || sel_rad1)) & (dist2 < (radius3*radius3))
    end
    
end

"""
Estimates the winding number for point with respect to a given polygon
"""
function winding_count(pin,poly;dist=nothing)
    #Number of segments
    nseg=size(poly,1)-1
    
    #check if the point intersect y-axis
    ipos_y=@. poly[:,2]>pin[2]
    ineg_y=@. poly[:,2]<pin[2]
    
    is_up = @. ineg_y[begin:end-1] & ipos_y[begin+1:end]
    is_down = @. ipos_y[begin:end-1] & ineg_y[begin+1:end]
    
    #get is left
    a2=poly[begin:end-1,2] .- poly[begin+1:end,2]
    b2=poly[begin+1:end,1] .- poly[begin:end-1,1]
    c2=(poly[begin:end-1,1] .* poly[begin+1:end,2]) .- (poly[begin+1:end,1] .*poly[begin:end-1,2])
    
    if(dist==nothing)
        dist=zeros(nseg)
    end
    
    dist .=(a2 * pin[1]) + (b2 * pin[2]) + c2;
    is_left= @. dist<0
    is_right=@. dist>0
    
    wn = size(findall(is_up .& is_left),1) - size(findall(is_down .& is_right),1)
    return wn
    #@show wn
    
    #@show is_up,is_down
    
    if(false) 
        plot(poly[:,1],poly[:,2],ls=:solid,color=:red)
        for ii in 1:nseg
            if(is_left[ii])
                plot!(poly[ii:ii+1,1],poly[ii:ii+1,2],ls=:solid,color=:blue)
            elseif(is_right[ii])
                plot!(poly[ii:ii+1,1],poly[ii:ii+1,2],ls=:solid,color=:green)
            end
        end
        plot!([pin[1]],[pin[2]],seriestype=:scatter,color=:black)
        plot!(poly[:,1],poly[:,2],seriestype=:scatter,color=:black) |> display
    end
    
end

"""distance of point p2 from line segment defined by p1,p2"""
function line_distance(p1,p2,p3)
    dist=(p1[2]-p2[2])*p3[1] 
    dist +=  (p2[1]-p1[1])* p3[2]
    dist += ((p1[1]*p2[2]) - (p2[1]*p1[2]))
    return dist
end


""" Determines if a point is inside the polygon or not
Argument: pin: point in x,y plane of type array and length 2
 poly is array of polygon vertices with first and last vertex being same to draw all the edges

"""
function InPolyWinding_loop(pin,poly)
    #Number of segments
    nseg=size(poly,1)-1
    
    wn=0
    for ii in 1:nseg
        #println(ii)
        if((poly[ii,2]<pin[2]) & (poly[ii+1,2]>pin[2])) #is up
            dist=line_distance(view(poly,ii,:),view(poly,ii+1,:),pin)
            if(dist<0) # is left?
                wn +=1
            end
        elseif((poly[ii,2]>pin[2]) & (poly[ii+1,2]<pin[2])) # is down
            dist=line_distance(view(poly,ii,:),view(poly,ii+1,:),pin)
            if(dist>0) # is right?
                wn -=1
            end
        end
    end
    
    return wn!=0
end
        
    
"""
return the points in side the polygon and point inside the padded circle
The padded circle is defined with center at the mean of polygon vertices and 
radius as the sum of maximum vertex length and a given padding (pad_circle)
The function must only set things to true and not set things to false
as it will be using arrays preallocated which might have been true for other reasons
Except when something is inpoly we want to et in circle to false
"""
function PointsInPoly_and_PaddedCircle!(x,y,poly,pad_circle;max_rad=nothing,in_poly=nothing,in_circle=nothing)
    nvert=size(poly,1)
    nobj=size(x,1)
    
    #memory allocation if not pre-allocated
    if(in_poly==nothing)
        in_poly=zeros(Bool,nobj)
    end
    if(in_circle==nothing)
        in_circle=zeros(Bool,nobj)
    end

    
    circ_cent=[mean(poly[:,1]),mean(poly[:,2])]
    
    #determine the maximum distance of polygon vertex from the center
    if(max_rad==nothing)
        max_rad2=maximum_radius2_polygon(poly,circ_cent)
        max_rad=sqrt(max_rad2)
    else
        max_rad2=max_rad^2
    end
   
    
    pad_circ_rad2= (max_rad+pad_circle)^2
    
    
    #calculate the distance of each point from the circle center
    @inbounds for ii in 1:nobj
        if(in_poly[ii])
            continue
        end
        #calculate the distance
        dist2 = (x[ii]-circ_cent[1])^2 + (y[ii]-circ_cent[2])^2
        
        if(dist2>pad_circ_rad2) #outside the padded circle
            continue
        elseif(dist2>max_rad2)#Determine if in padded circle
            in_circle[ii] = true
            continue
        elseif(InPolyWinding_loop([x[ii],y[ii]],poly)) #checks if truly inside the polygon
            in_poly[ii] = true 
            #If something is inside the polygon then we must set in circle to false
            in_circle[ii] = false
        else
            in_circle[ii]=true
        end
    end
    
    #To make plots the result for testing and debugging
    if(false)
        plot(x[in_poly],y[in_poly],seriestype=:scatter,color=:red)
        plot!(x[in_circle],y[in_circle],seriestype=:scatter,color=:green)
        out_all= .! (in_circle .|| in_poly)
        #@show sum(out_all),sum(in_circle),sum(in_poly),size(in_poly)
        plot!(x[out_all],y[out_all],seriestype=:scatter,color=:grey)
        
        plot!([circ_cent[1]-sqrt(pad_circ_rad2),circ_cent[1]+sqrt(pad_circ_rad2)],
        [circ_cent[2],circ_cent[2]],s=:solid,color=:red,lw=2)
        
        plot!(poly[:,1],poly[:,2],ls=:solid,color=:black,lw=2) |> display
        
    end
end

"""Estimates the maximum radius of the polygon
"""
function maximum_radius2_polygon(poly,cent)
    max_rad2=0
    @inbounds for ii in 1:size(poly,1)-1
        this_rad2= (poly[ii,1]-cent[1])^2 +  (poly[ii,2]-cent[2])^2 
        if(this_rad2>max_rad2)
            max_rad2=this_rad2
        end
    end
    return max_rad2
end
 

""" checks if two polygon intersects
First check for all the vertices in poly1 and check if any in poly2
Then rpeat to all vertices is in poly2 for poly1
If none of them are inside then polygon do not intersects
This function will not work for axis aligned polygon where they fully or
partially overlaps and the vetices of one is one the edges of another.
"""
function PolygonsIntersect(poly1,poly2)
    for ii in 1:size(poly1,1)
        if(InPolyWinding_loop(poly1[ii,:],poly2))
            return true
        end
    end
    
    for ii in 1:size(poly2,1)
        if(InPolyWinding_loop(poly2[ii,:],poly1))
            return true
        end
    end
    return false
end

"""Check if the bounding box of two input polygon is non-intersecting
just touhing is considered intersecting to be conservative

example:
poly1=zeros(5,2)
poly2=zeros(5,2)

poly1[:,1] .= [0,0,1,1,0]
poly1[:,2] .= [0,1,1,0,0]

dx=0
dy=0.8
poly2[:,1] .= [0,0,1,1,0] .+ dx
poly2[:,2] .= [0,1,1,0,0] .+ dy

PolygonsIntersect(poly1,poly2)
PolyIntersect_boundingbox(poly1,poly2)

plot(Shape(poly1[:,1],poly1[:,2]),color=:green)
plot!(Shape(poly2[:,1],poly2[:,2]),color=:red) |> display

"""
function PolyIntersect_boundingbox(poly1,poly2)
    poly1_min=minimum(poly1,dims=1)
    poly1_max=maximum(poly1,dims=1)
    
    poly2_min=minimum(poly2,dims=1)
    poly2_max=maximum(poly2,dims=1)
    
    if( ((poly1_max[1]<poly2_min[1]) || (poly2_max[1]<poly1_min[1])) ||
        ((poly1_max[2]<poly2_min[2]) || (poly2_max[2]<poly1_min[2]))
        )
        return false
    else
        return true
    end
end

