#Any plotting functionality should be here

"""
Function to make plot of positioners
Takes input a POSRobotBody structure and plots the on_target location of positioner
"""
function plot_moved_positioner(pos;linecolor=:black,bodycolor=:blue,headcolor=:red,alpha=0.2,hold=false)
    nT=pos.nseg[1]
    nP=pos.nseg[2]
    
    #plots center
    if(hold)
        plot!([pos.xy[1]],[pos.xy[2]],seriestype=:scatter,color=:black,markersize = 5) 
    else
        plot([pos.xy[1]],[pos.xy[2]],seriestype=:scatter,color=:black,markersize = 5) 
    end
    #plot([pos.xy[1],pos.axis_on_target[1]],[pos.xy[2],pos.axis_on_target[2]],seriestype=:scatter,color=:black,markersize = 5) 
    
    #plot(pos.body_on_target[1:nT,1],pos.body_on_target[1:nT,2],linestyle=:solid)
    #plots body
    plot!(pos.body_on_target[1:nT,1],pos.body_on_target[1:nT,2],
        seriestype=[:shape,],lw=0.5,c=bodycolor,linecolor=linecolor,
        legend=false,fillalpha=alpha)
    
    plot!(pos.head_on_target[1:nP,1],pos.head_on_target[1:nP,2],
        seriestype=[:shape,],lw=0.5,c=headcolor,linecolor=linecolor,
        legend=false,fillalpha=alpha)
    
    #plot!(pos.head_orig[1:nP,1],pos.head_orig[1:nP,2],
    #    seriestype=[:shape,],lw=0.5,c=:green,linecolor=linecolor,
    #    legend=false,fillalpha=alpha)
    
    plot_circle(pos.xy[1],pos.xy[2],pos.arms[1],color=:black,alpha=0.1)
    plot_circle(pos.xy[1],pos.xy[2],pos.arms[1]+pos.arms[2],color=:green,alpha=0.05)
end


"""plots the assignment in a tile
Currently assumes targets with three different priority
Should be generaziled later if this tool is used
"""
function plot_tile_assignment(position_assigned,targets_dic,fp_dic;ri=1,outfile="test.png")
    
    iass=position_assigned[:,ri] .>0
    iass=sort(position_assigned[iass,ri])
    
    bool_fib=zeros(Bool,size(targets_dic["X_FP"],1))
    bool_fib[iass] .=true
    
    #plot all target and assigned objects
    p1=plot(targets_dic["X_FP"][bool_fib],targets_dic["Y_FP"][bool_fib],seriestype=:scatter,color=:black,
        ma=0.5,size=(1500,1500),axis=nothing,title="Assigned targets")
    #plot un-assigned object
    p2=plot(targets_dic["X_FP"][.!bool_fib],targets_dic["Y_FP"][.!bool_fib],seriestype=:scatter,color=:black,
        ma=0.1,size=(1500,1500),axis=nothing,title="Un-Assigned targets")
    #plot un-assigned positioner
    pos_un_ass=position_assigned[:,ri] .== 0
    p3=plot(fp_dic["OFFSET_X"][pos_un_ass],fp_dic["OFFSET_Y"][pos_un_ass],seriestype=:scatter,color=:red,
        ma=0.5,size=(1500,1500),axis=nothing,title="Un-Assigned Positioners")
    
    
    #split by priority for top three priority
    usp=sort(unique(targets_dic["PRIORITY"]),rev=:true)
    bp=findall(bool_fib .& (targets_dic["PRIORITY"] .== usp[1]))
    ncount=size(bp,1)
    #@show size(bp),sum(targets_dic["PRIORITY"] .== usp[1]) 
    p4=plot(targets_dic["X_FP"][bp],targets_dic["Y_FP"][bp],seriestype=:scatter,color=:red,
        ma=1.0,size=(1500,1500),axis=nothing,title="PRIORITY=$(usp[1]) ,count=$(ncount)")
    
    bp=findall(bool_fib .& (targets_dic["PRIORITY"] .== usp[2]))
    ncount=size(bp,1)
    p5=plot(targets_dic["X_FP"][bp],targets_dic["Y_FP"][bp],seriestype=:scatter,color=:blue,
        ma=1.0,size=(1500,1500),axis=nothing,title="PRIORITY=$(usp[2]) ,count=$(ncount)")
    #@show size(bp),sum(targets_dic["PRIORITY"] .== usp[2])
    
    bp=findall(bool_fib .& (targets_dic["PRIORITY"] .== usp[3]))
    ncount=size(bp,1)
    p6=plot(targets_dic["X_FP"][bp],targets_dic["Y_FP"][bp],seriestype=:scatter,color=:green,
        ma=0.5,size=(1500,1500),axis=nothing,title="PRIORITY=$(usp[3]) ,count=$(ncount)")
    #@show size(bp),sum(targets_dic["PRIORITY"] .== usp[3])
    
    plot(p1, p2,p3,p4,p5,p6, layout = (2, 3), #title=["Assigned","UnAssigned"],
        legend = false,aspect_ratio=1)
    savefig(outfile)
    println("Generated: $(outfile)")
    
end

"""plots a circle with given properties
"""
function plot_circle(x,y,rad;color=:blue,linecolor=:black,alpha=0.2)
    theta=LinRange(0,2*pi,500)
    plot!(x .+ rad.*sin.(theta), y.+ rad.*cos.(theta),
        seriestype=[:shape,],lw=0.5,c=color,linecolor=linecolor,
        legend=false,fillalpha=alpha,aspect_ratio=1)
end


"""
plots the target in the zone and positioner
as as any auxiliary informatio
"""
function plot_zones(pos,zone_target,fp_dic)
    nc=zone_target.ntarget[1]
    plot(zone_target.X_FP[1:nc],zone_target.Y_FP[1:nc],seriestype=:scatter,color=:black)
    ifp=fp_dic["zone_map"][zone_target.zones[3]]
    #plot!(fp_dic["OFFSET_X"][ifp[1]:ifp[2]],fp_dic["OFFSET_Y"][ifp[1]:ifp[2]],
    #    seriestype=:scatter,color=:green,markersize=10)
    
    iass=findall(zone_target.pot_ass[1:nc])
    plot!(zone_target.X_FP[iass],zone_target.Y_FP[iass],seriestype=:scatter,color=:red)
    icoll=findall(zone_target.pot_coll[1:nc])
    plot!(zone_target.X_FP[icoll],zone_target.Y_FP[icoll],seriestype=:scatter,color=:blue) 
    
    
end
