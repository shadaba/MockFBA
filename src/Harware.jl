#To deal with DESI hardware sepcific functionality
"""Generate the bounding grid of focal plane
"""
struct focal_plane_Bound
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    dx_zone::Float64
    dy_zone::Float64
    nx::Int32
    ny::Int32
end


""" Adapted from desimodel.io 
https://github.com/desihub/desimodel/blob/d5f7f871873c547892a40ece30890b8610f53768/py/desimodel/io.py#L336
Loads platescale.txt.
    Returns
    -------
    :class:`~numpy.recarray`\n
        The data table read from the file.
    Notes
    -----
    The returned object has these columns:
    radius
        Radius from center of focal plane [mm].
    theta
        Radial angle that has a centroid at this radius [deg].
    radial_platescale
        Meridional (radial) plate scale [um/arcsec].
    az_platescale:
        Sagittal (azimuthal) plate scale [um/arcsec].
    arclength:
        Unknown description.

"""
function load_platescale(focalplane_dir)

    plate_scale_file="$(focalplane_dir)platescale.txt"
    plate_scale_columns=Dict("radius"=>1 ,"theta"=>2 ,"radial_platescale"=>7,
                   "az_platescale"=>8)
    
    rzs_file="$(focalplane_dir)rzsn.txt"
    rzs_column=Dict("R"=>1 ,"Z"=>2, "S"=>3, "N"=>4)

    #load and setup the plate scale disctionary
    platescale_dic=Dict()
    tab1=readdlm(plate_scale_file, Float64,comments=true)
    for col in keys(plate_scale_columns)
        platescale_dic[col]=tab1[:,plate_scale_columns[col]]
    end
    
    # One of the header column is not commented hence skipstart
    tab1=readdlm(rzs_file, Float64,comments=true,skipstart=8)
    
    #This function take yas first input and then x as second!!
    #https://htmlpreview.github.io/?https://github.com/PumasAI/DataInterpolations.jl/blob/v2.0.0/example/DataInterpolations.html
    rz_itp = QuadraticInterpolation(tab1[:,rzs_column["Z"]],tab1[:,rzs_column["R"]])
    
    platescale_dic["arclength"]=rz_itp.(platescale_dic["radius"])

    return platescale_dic
end

"""Adapted from desimodel.io
Returns an array of radii in mm given an array of radii in degrees using the platescale data
    relative to the center of the focal plane as (0,0). Supports scalar and vector inputs.\n
    Parameters
    ----------
    theta : :class:`float` or array-like
        An array that represents the angle from the center of the focal plane.
    Returns
    -------
    :class:`float` or array-like
        Radii in mm.
"""
function get_radius_mm(platescale::Dict{Any,Any},theta::Vector)
    
    #platescale = load_platescale(focalplane_dir)
    
    # Uses a quadratic one-dimensional interpolation to approximate the radius in degrees versus radius in mm
    #This function take yas first input and then x as second!!
    #https://htmlpreview.github.io/?https://github.com/PumasAI/DataInterpolations.jl/blob/v2.0.0/example/DataInterpolations.html
    rz_itp = QuadraticInterpolation(platescale["radius"],platescale["theta"])
    
    radius = rz_itp.(theta)
    
    #This to plot the interpolation for testing
    #scatter(platescale["theta"], platescale["radius"], label="input data")
    #plot!(rz_itp) 
    #plot!(theta,radius) |> display
    
    return radius
end

function get_radius_mm(focalplane_dir::String,theta)
    platescale = load_platescale(focalplane_dir)
    radius=get_radius_mm(platescale,theta)
    return radius
end


    

""" Adapted from desimodel.focalplane.geometry \n
Returns arrays of the x, y positions of given celestial objects
    on the focal plane given an arbitrary telescope pointing in RA and Dec and
    arrays of the `ra` and `dec` of celestial objects in the sky.\n
    Parameters
    ----------
    telra : :class:`float`
        The telescope's RA pointing in degrees.
    teldec : :class:`float`
        The telescope's Dec pointing in degrees.
    ra : array-like
        An array of RA values for locations in the sky.
    dec : array-like
        An array of Dec values for locations in the sky.
    Returns
    -------
    tuple
        The x, y positions corrsponding to `ra`, `dec`.
    Notes
    -----
    Implements the Haversine formula.
"""
function xyz_unit_2_xy_focalplane(telra, teldec, x_unit,y_unit,z_unit,platescale::Dict{Any,Any})
    # Inclination is 90 degrees minus the declination in degrees
    nobj=size(x_unit,1)
    
    telra_rad=deg2rad(telra)
    teldec_rad=deg2rad(teldec)

    # Clockwise rotation around the z-axis by the right ascension of the tile center
    #Note that this uses column major order which mean
    #the first three values are in first column and not first row
    rarotate = SMatrix{3,3}(cos(telra_rad),-sin(telra_rad),0,
                            sin(telra_rad),cos(telra_rad),0,
                            0,0,1)
    

    # Counter-Clockwise rotation around y axis by declination of the tile center
    decrotate = SMatrix{3,3}(cos(teldec_rad),0,-sin(teldec_rad),
                            0,1,0,
                            sin(teldec_rad),0,cos(teldec_rad))
        

    dec_ra_rotate=decrotate*rarotate
    

    x_new = dec_ra_rotate[1,1].*x_unit .+ dec_ra_rotate[1,2].*y_unit .+ dec_ra_rotate[1,3].*z_unit
    y_new = dec_ra_rotate[2,1].*x_unit .+ dec_ra_rotate[2,2].*y_unit .+ dec_ra_rotate[2,3].*z_unit
    z_new = dec_ra_rotate[3,1].*x_unit .+ dec_ra_rotate[3,2].*y_unit .+ dec_ra_rotate[3,3].*z_unit

    newteldec = 0
    newtelra = 0
    
    ra_rad = atan.(y_new, x_new)
    
    dec_rad = (pi / 2) .- acos.(z_new) #assuming x,y,z on unit sphere so norm is 1
    
    radius_rad= (sin.(0.5 .* dec_rad).^2) .+ ((sin.(0.5 .* ra_rad).^2).*cos.(dec_rad))
    @. radius_rad = 2.0 .* asin.(sqrt.(radius_rad))
    
    #Finally convert this in degree but keep the variable name same for memory
    radius_rad .= radius_rad .* (180/pi)

    q_rad = atan.(z_new, -y_new)

    radius_mm = get_radius_mm(platescale,radius_rad)

    x_focalplane =@. radius_mm * cos(q_rad)
    y_focalplane =@. radius_mm * sin(q_rad)

    return x_focalplane, y_focalplane
end

function radec2xy_focalplane(telra, teldec, ra, dec,platescale::Dict{Any,Any})
    
    #Convert dec to inclination angle in radian
    x_unit,y_unit,z_unit = embed_sphere(ra,dec)
    
    x_focalplane, y_focalplane= xyz_unit_2_xy_focalplane(telra, teldec, x_unit,y_unit,z_unit,platescale)
    return x_focalplane, y_focalplane
end

"""
Estimates the focal plane boundaries for the efficient neighbour finding zone definition
"""
function get_bound_focal_plane(xfocal,yfocal,grid_scale;pad=0.2)
    xlim=[minimum(xfocal)-pad,maximum(xfocal)+pad]
    ylim=[minimum(yfocal)-pad,maximum(yfocal)+pad]
    
    
    nx=floor(Int,(xlim[2]-xlim[1])/grid_scale)
    ny=floor(Int,(ylim[2]-ylim[1])/grid_scale)
    
    dx_zone=(xlim[2]-xlim[1])/(nx-1)
    dy_zone=(ylim[2]-ylim[1])/(ny-1)
    
    fp_bound=focal_plane_Bound(xlim[1],xlim[2],ylim[1],ylim[2],dx_zone,dy_zone,nx,ny)
    
    return fp_bound
end

function Assign_xy_zone(xfocal,yfocal,fp_bound)
    nobj=size(xfocal,1)
    
    zone_x=zeros(Int16,nobj)
    zone_y=zeros(Int16,nobj)
    
    @. zone_x = floor(Int16,(xfocal-fp_bound.xmin)/fp_bound.dx_zone)
    @. zone_y = floor(Int16,(yfocal-fp_bound.ymin)/fp_bound.dy_zone)
    
    
    if(false)#To make some plots for visualization
        plot([0],[0],seriestype=:scatter,color=:green)
        for ii in 0:fp_bound.nx-1
            for jj in 0:fp_bound.ny-1
                indsel= findall((zone_x .== ii) .& (zone_y .== jj))
                if(size(indsel,1)>0)
                    println(ii,' ',jj,' ',size(indsel,1))
                    plot!(xfocal[indsel],yfocal[indsel],seriestype=:scatter,color=:auto)
                end
            end
        end
                    
        plot!([0],[0],seriestype=:scatter,color=:blue) |> display
    end
    
    return zone_x,zone_y
end

"""
Loads the status of focal plane based on given date
look for following file with appropriate dates
    desi-exclusion_2019-09-16T00:00:00.yaml  
    desi-focalplane_2019-09-16T00:00:00.ecsv 
    desi-state_2019-09-16T00:00:00.ecsv
"""
function load_focal_plane_hardware(config;date="2019-09-16T00:00:00")
    
    harware_fp=Dict()
    
    file_fp="$(focalplane_dir)desi-focalplane_$(date).ecsv"
    df=CSV.read(file_fp, DataFrame;delim=" ",comment="#",header=1)
    
    
    return df

end

"""Load the focal plane file from jld2 format
if jld2 format doesn't exists then load the ascii version
and convert it to jld2 seperating into two group
one called base consisting of column we should need for fba
and other called extra contain everything else
It only loads the base part on any future call if jld2 file exists
"""
function load_hw_focalplane_file(config;date="2019-09-16T00:00:00")

    base_cols=["OFFSET_X","OFFSET_Y","LOCATION","MIN_T","MAX_T","MIN_P","MAX_P",
                    "OFFSET_T","OFFSET_P","LENGTH_R1","LENGTH_R2","DEVICE_TYPE"]
    
    fp_jld2="$(config["focalplane_dir_jld2"])desi-focalplane_$(date).jld2"
    

    if(isfile(fp_jld2))
        #fp_base=load(fp_jld2,"base")
        fp_base=read_columns_jld2(fp_jld2,base_cols)
    else
        file_fp="$(config["focalplane_dir"])desi-focalplane_$(date).ecsv"
        if(! isfile(file_fp))
            warn("File no Found $(file_fp)")
        else#Load the csv file convert to jld2 file
            df=CSV.read(file_fp, DataFrame;delim=" ",comment="#",header=1)
            
            DataFrame_TO_JLD2(df,fp_jld2,base_cols)
            fp_base=read_columns_jld2(fp_jld2,base_cols)
        end
    end
    return fp_base
end

function load_hw_state_file(config;date="2019-09-16T00:00:00")

    base_cols=["LOCATION", "STATE", "EXCLUSION"]
    st_jld2="$(config["focalplane_dir_jld2"])desi-state_$(date).jld2"
    

    if(isfile(st_jld2))
        st_base=read_columns_jld2(st_jld2,base_cols)
    else
        file_st="$(config["focalplane_dir"])desi-state_$(date).ecsv"
        if(! isfile(file_st))
            warn("File no Found $(file_st)")
        else#Load the csv file convert to jld2 file
            df=CSV.read(file_st, DataFrame;delim=" ",comment="#",header=1)
            #Write the JLD2 file
            DataFrame_TO_JLD2(df,st_jld2,base_cols)
            #st_base=load(st_jld2,"base")
            st_base=read_columns_jld2(st_jld2,base_cols)
        end
    end
    return st_base
end


"""transforms the nested ditionary loaded from yaml to flat dictionary
"""
function transform_dictionary_flatJLD(jldfile,dicin)
    jldopen(jldfile,"a+",compress=LZ4FrameCompressor()) do file
        for tkey in keys(dicin)
            for pkey in keys(dicin[tkey])
                for lkey in keys(dicin[tkey][pkey])
                    nel=length(desi_exc[tkey][pkey][lkey])
                    if(nel==0)
                        continue
                    elseif(nel>1)
                        @warn(" Found more than one element $(tkey) > $(pkey) > $(lkey) : $(nel)\n
                            Ignoring elements after the first one please chekc the raw file for errors")
                    end
                    
                    tarr=desi_exc[tkey][pkey][lkey][1]
                    nseg=length(tarr)
                    if(lkey=="circles")
                        out_arr=[tarr[1][1],tarr[1][2],tarr[2][1]]
                    else
                        out_arr=zeros(nseg,2)
                        for ii in 1:nseg
                            out_arr[ii,:]=tarr[ii]
                        end
                    end
                    write(file,"$(tkey)-$(pkey)-$(lkey)",out_arr,compress=true)
                end
            end
        end
    end
end

"""load the flat dictionary of array and converts it into a StaticArray
"""
function load_exclusion(jldfile)
    tdic=load(jldfile)

    for tkey in keys(tdic)
        nsize=size(tdic[tkey])
        if(length(nsize)==2)
            tdic[tkey]=SMatrix{nsize[1],nsize[2],Float32}(tdic[tkey])
        else
            tdic[tkey]=SVector{nsize[1],Float32}(tdic[tkey])
        end
        #println(tkey,nsize[1])
    end
    return tdic
end


function load_hw_exclusion_file(config;date="2019-09-16T00:00:00")

    exc_jld2="$(config["focalplane_dir_jld2"])desi-exclusion_$(date).jld2"
    
    
    if(isfile(exc_jld2))
        exc=load_exclusion(exc_jld2)
    else
        file_exc="$(config["focalplane_dir"])desi-exclusion_$(date).yaml"
        #loads the yaml file
        desi_exc=YAML.load_file(file_exc)
        #flattens the ditionary, convert to array and compress and save
        transform_dictionary_flatJLD(exc_jld2,desi_exc)
        
        #finally load the new format and covert them to static arrays for speed
        exc=load_exclusion(exc_jld2)
    end
    return exc
end



function DataFrame_TO_JLD2(df,jldfile,base_cols;compress=false)
    #save the base property
    #tdic_base=Dict()
    #tdic_extra=Dict()
    
    jldopen(jldfile, "a+",compress=LZ4FrameCompressor()) do file
        for col in names(df)
            if(col in base_cols)
                #tdic_base[col]=df[!,col]
                write(file,"$(col)",df[!,col],compress=true)
            else
                #tdic_extra[col]=df[!,col]
                write(file,"$(col)",df[!,col],compress=true)
            end
        end

    
        #delete if needed or just write
        #write(file,"extra",tdic_extra,compress=compress)
    end
            
    println("saved base and extra properties in $(jldfile)")
end


function read_columns_jld2(jldfile,columns)
    tdic=Dict()
    jldopen(jldfile,"r") do file
        for col in columns
            tdic[col]=read(file,col)
        end
    end
    return tdic
end

"""
To load the full status of focal plane on certain date
first loads the focalplane, state and exclusion library
It then combines the focal plane and state together in one coherent structure
example: 
config_fp=Dict("focalplane_dir" => "focalplane/",
            "focalplane_dir_jld2" => "focalplane_jld2/",
            "date" => "2019-09-16T00:00:00")

tdic=load_hw_focalplane_file(config_fp;date="2019-09-16T00:00:00")
"""
function load_hw_full_FocalPlane!(config;date="2019-09-16T00:00:00",
      default_prop=Dict("STATE"=> 0 ,"EXCLUSION"=>"default") )
    #load the focal plane
    fp_dic=load_hw_focalplane_file(config;date=date)
    #load the state
    state_dic=load_hw_state_file(config;date=date)
    #load the exclusion polygons
    exclusion_dic=load_hw_exclusion_file(config;date=date)
    
    #Now append a state column to the original fp_dic
    nobj=size(fp_dic["OFFSET_X"],1)
    fp_dic["STATE"]=zeros(Int8,nobj)
    fp_dic["STATE"] .= default_prop["STATE"]
    fp_dic["EXCLUSION"]=[default_prop["EXCLUSION"] for ii=1:nobj]#zeros(String,nobj)
    #In the future we might want to update other positioner info like R1,R2 or theta,phi limits
    
    
    isort_state=sortperm(state_dic["LOCATION"])
    isort_fp=sortperm(fp_dic["LOCATION"])
    
    nst=size(isort_state,1)
    ist=1
    for ii in 1:nobj
        if(ist>nst)#fill default
            break
        elseif(fp_dic["LOCATION"][isort_fp[ii]]==state_dic["LOCATION"][isort_state[ist]])
            for prop in ["STATE","EXCLUSION"]
               fp_dic[prop][isort_fp[ii]]= state_dic[prop][isort_state[ist]]
            end
        #elseif(fp_dic["LOCATION"][isort_fp[ii]]<state_dic["LOCATION"][isort_state[ist]])
        #    #fill the default
        #    for prop in ["STATE","EXCLUSION"]
        #       fp_dic[prop][isort_fp[ii]]= default_state[prop]
        #    end
        end
        ist +=1
    end
  
    fp_dic["tile_date"] =date
    exclusion_dic["tile_date"]=date 
    return fp_dic,exclusion_dic
end




"""convert positioner x,y to theta phi for the positioner
"""
function Hardware_xy_to_thetaphi(pos,x,y,rad2)
    dx=x-pos.xy[1]
    dy=y-pos.xy[2]
    
    if(isapprox(rad2,pos.arms2_lim[2]))
        # We are at the maximum arm extension.  Force phi angle to zero
        # and compute theta.
        phi = 0.0;
        theta = atan(dy, dx)
    elseif (isapprox(rad2,pos.arms2_lim[1]))
        # We are at the limit of the arm folded inwards.  Force phi angle
        # to PI and compute theta.
        phi = pi;
        theta = atan(dy,dx)
    else #We are on neither limit.
        theta_r2=pos.arms[1]*pos.arms[1]
        phi_r2=pos.arms[2]*pos.arms[2]
        
        # Use law of cosines to compute "opening" angle at the "elbow".
        opening = acos((theta_r2 + phi_r2 - rad2)/ (2.0 * pos.arms[1]*pos.arms[2]))

        # The PHI angle is just the supplement of this.
        phi = pi - opening

        # Compute the theta angle.
        # Use law of cosines to compute angle from theta arm to the line from
        # the origin to the X/Y position.
        nrm_offset = sqrt(rad2)
        txy = acos((theta_r2 + rad2 - phi_r2)/ (2 * pos.arms[1] * nrm_offset))
        theta = atan(dy, dx) - txy
    end
    
    return theta,phi
end

function Hardware_xy_to_thetaphi(pos,x,y)
    rad2=(x-pos.xy[1])^2 + (y-pos.xy[2])^2
    return Hardware_xy_to_thetaphi(pos,x,y,rad2)
end

"""
Takes the difference of two angle
"""
function angle_diff(hi, low) 
    twopi = 2.0 * pi
    # range reduction to [-Pi, Pi)
    if (hi >= pi) 
        hi -= twopi
    elseif(hi < -pi) 
        hi += twopi
    end
    
    if (low >= pi) 
        low -= twopi
    elseif(low < -pi) 
        low += twopi
    end
    
    diff = hi - low
    if (diff > 0.0) 
        return diff
    else
        return twopi + diff
    end
end
                        
"""
theta_lim consists of MIN_T,MAX_T,OFFSET_T
phi_lim consists of MIN_P, MAX_P, OFFSET_P
"""
function outside_theta_phi_range(theta,phi,theta_lim,phi_lim)
    #= Adapted from dsihub/fibreassign
    // The definitions of the theta / phi angles and ranges can be found in
    // DESI-0899.

    // Theta angle.  The theta offset loaded from desimodel is already in global
    // focalplane coordinates (not petal-local).  The theta min / max values are
    // relative to this offset (not the coordinate system).  Compute the angle
    // difference rotating both directions from theta_zero and see if we can reach
    // it at least one of those ways.  Note that theta_min is negative to indicate
    // a clockwise rotation from theta_zero.
    =#
    diff_hi = angle_diff(theta, theta_lim[3])
    diff_lo = angle_diff(theta_lim[3], theta)

    #@show "theta",diff_hi,diff_lo,theta_lim
    #=
    // Since calling code may add theta_offset (aka theta_zero here),
    // which we then subtract off, we want to put a little margin on
    // the min/max checks.  (In most cases where we are passing
    // theta/phi values through from the desimodel inputs and we need
    // to add offsets, we also ignore this range check, so this is
    // maybe overly cautious.)
    =#
    eps = 1e-12;

    if ((   ( diff_hi > (theta_lim[2] + eps)) || ( diff_hi < (theta_lim[1] - eps)))
        & ((-diff_lo > (theta_lim[2] + eps)) || (-diff_lo < (theta_lim[1] - eps))) 
        ) 
        return true
    end

    #=
    // Phi angle.  The phi offset (phi_zero) is relative to the coordinate system
    // defined by the theta arm along the X axis.  The phi min / max values are
    // relative to this offset (not the theta-arm coordinate system).  Note that
    // a negative phi_min indicates a clockwise rotation from phi_zero.
    =#
    diff_hi = angle_diff(phi, phi_lim[3])
    diff_lo = angle_diff(phi_lim[3], phi)
    #@show "phi",diff_hi,diff_lo,phi_lim
    if ( (   ( diff_hi > (phi_lim[2] + eps)) || ( diff_hi < (phi_lim[1] - eps)) )
        & ((-diff_lo > (phi_lim[2] + eps)) || (-diff_lo < (phi_lim[1] - eps)))
        ) 
        return true
    end

    #// not outside range
    return false
end



