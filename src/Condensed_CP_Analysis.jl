module Condensed_CP_Analysis

import Interpolations: scale, interpolate, BSpline, Linear, Quadratic, Cubic, Flat, Periodic, OnGrid, gradient, hessian, Gridded
import LinearAlgebra: dot, cross, normalize, norm
import LoggingExtras, LoggingFormats
import NLsolve
import Parsers
import Interact
import Roots: find_zeros, Newton
using StaticArrays
using Meshes, MeshViz
import GLMakie as Mke
import SplitApplyCombine: invert, flatten
import NaNStatistics: movmean
using Plots
using LaTeXStrings
using Assignment
import CSV: CSV
import DataFrames: DataFrame
using ImageFiltering

function set_str_to_int_array(set_str)
    out = []
    for i in split(replace(set_str, "[" => "", "]" => ""), ",")
        if occursin("-", i)
            start, stop = [parse(Int, j) for j in split(i, "-")]
            for j in start:stop
                push!(out, j)
            end
        else
            push!(out, parse(Int,i))
        end
    end
    return out
end

function import_dat(fname)
    @info "Importing data from $fname"
    f = open(fname)

    out = Dict()

    # header
    out["title"] = strip(replace(join(split(readline(f), "=")[2:end]), "\"" => ""))
    out["variables"] = [strip(replace(join(split(readline(f), "=")[2:end]), "\"" => ""))]
    line = ""
    while !eof(f)
        line = readline(f)
        if !startswith(line, "\"")
            break
        end
        push!(out["variables"], strip(replace(line, "\"" => "")))
    end

    # dataset aux data
    while !startswith(line, "ZONE") && !eof(f)
        sl = split(line)
        aux1 = sl[1]
        sl2 = split(join(sl[2:end]), "=")
        aux2 = sl2[1]
        auxval = replace(sl2[end], "\"" => "")
        if aux1 in keys(out)
            out[aux1][aux2] = auxval
        else
            out[aux1] = Dict(aux2 => auxval)
        end
        line = readline(f)
    end


    # now get all the zones
    out["zones"] = Dict()

    while !eof(f)
        if !startswith(line, "ZONE")
            line = readline(f)
            continue
        end
        # get zone header
        line = replace(line, "ZONE " => "")
        zone = Dict()
        while occursin("=", line) && !eof(f)
            props = split(line, ", ")
            for prop in props
                if occursin("AUXDATA", prop)
                    sl = split(prop)
                    aux1 = strip(sl[1])
                    sl2 = split(join(sl[2:end]), "=")
                    aux2 = strip(sl2[1])
                    auxval = replace(sl2[end], "\"" => "")
                    if aux1 in keys(zone)
                        zone[aux1][aux2] = auxval
                    else
                        zone[aux1] = Dict(aux2 => auxval)
                    end
                else
                    ps = split(prop, "=")
                    zone[strip(ps[1])] = strip(replace(join(ps[2:end],"="), "\"" => ""))
                end
            end
            line = readline(f)
        end

        # prepare some header data for processing
        num_elem_lines, rem_elem_col = 0,0
        val_per_line = 5
        if zone["ZONETYPE"] == "Ordered"
            num_points = 1
            for i in ["I","J","K"]
                zone[i] = parse(Int, zone[i])
                num_points *= zone[i]
            end
            zone["NUMPOINTS"] = num_points
        elseif zone["ZONETYPE"] == "FETriangle"
            for i in ["Nodes","Elements"]
                zone[i] = parse(Int, zone[i])
            end
            zone["NUMPOINTS"] = zone["Nodes"]

            num_elem_lines, rem_elem_col = divrem(zone["Elements"], val_per_line)
            if rem_elem_col > 0
                num_elem_lines += 1
            end
        else
            break
        end

        cell_centered_var_nums = []
        if "VARLOCATION" in keys(zone) && occursin("CELLCENTERED", zone["VARLOCATION"])
            m = match(r"\[(.*)\]",zone["VARLOCATION"])
            cell_centered_var_nums = set_str_to_int_array(m[1])
        end

        zone["cell_centered_var_nums"] = cell_centered_var_nums

        num_lines, rem_col = divrem(zone["NUMPOINTS"], val_per_line)
        if rem_col > 0
            num_lines += 1
        end
        if zone["ZONETYPE"] != "FETriangle"
            num_elem_lines, rem_elem_col = num_lines, rem_col
        end

        # get the zone data at points/Nodes
        zone["variables"] = Dict()
        zone["DT"] = split(strip(replace(zone["DT"], "(" => "", ")" => "")))
        for (vi,dt) in enumerate(zone["DT"])
            if dt == "BIT"
                tmp_num_lines = (vi in cell_centered_var_nums ? num_elem_lines : num_lines)
                for li in 1:tmp_num_lines
                    line = readline(f)
                end
            else
                vals = Vector{Float32}(undef,(vi in cell_centered_var_nums ? zone["Elements"] : zone["NUMPOINTS"]))
                i = 1
                tmp_num_lines = (vi in cell_centered_var_nums ? num_elem_lines : num_lines)
                for li in 1:tmp_num_lines
                    for val in split(line)
                        vals[i] = parse(Float32, val)
                        i += 1
                    end
                    line = readline(f)
                end
                zone["variables"][out["variables"][vi]] = vals
            end
        end

        # save matrix of xyz values for convenience
        zone["nodes_xyz"] = invert([zone["variables"][dir] for dir in ["X","Y","Z"]])

        # get Elements
        if zone["ZONETYPE"] == "FETriangle"
            elements = Array{Int32}(undef, zone["Elements"], 3)
            for li in 1:zone["Elements"]
                for (j,e) in enumerate(split(line))
                    elements[li, j] = parse(Int32, e)
                end
                line = readline(f)
            end
            zone["element_list"] = elements

            # make mesh
            points = Point.(zone["variables"]["X"], zone["variables"]["Y"], zone["variables"]["Z"])
            tris = connect.([Tuple(elements[i,:]) for i in 1:size(elements,1)], Triangle)
            mesh = SimpleMesh(points, tris)
            zone["mesh"] = mesh
        end
        
        out["zones"][zone["T"]] = zone
    end

    close(f)
    return out
end

# angle difference of vectors in plane given plan normal n
# (from https://stackoverflow.com/a/33920320/2620767)
angle(a, b, n) = atan((a×b)⋅n, a⋅b)

function pitick(start, stop, denom; mode=:text)
    a = Int(cld(start, π/denom))
    b = Int(fld(stop, π/denom))
    tick = range(a*π/denom, b*π/denom; step=π/denom)
    ticklabel = piticklabel.((a:b) .// denom, Val(mode))
    tick, ticklabel
end

function piticklabel(x::Rational, ::Val{:text})
    iszero(x) && return "0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return S * N * "π"
    S * N * "π/" * repr(d)
end

function piticklabel(x::Rational, ::Val{:latex})
    iszero(x) && return L"0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return L"%$S%$N\pi"
    L"%$S\frac{%$N\pi}{%$d}"
end

function sphere_slice_zone_point_theta(zone)
    origin = [(minimum(zone["variables"][dir]) + maximum(zone["variables"][dir])) / 2 for dir in ["X","Y","Z"]]
    points = zone["nodes_xyz"]
    θ0 = normalize(points[1] - origin) # use first point to define the zero theta
    norm_vec = normalize(cross(points[1] - origin, points[div(end,3)] - origin))
    return [angle(normalize(points[i] - origin), θ0, norm_vec) for i in eachindex(points)], origin, norm_vec
end

function sphere_slice_analysis(file_path; smoothing_factor=0, var_check_str="INS: ", spacing=π/100)
    sys = import_dat(file_path)
    @info "Processing file $(sys["title"])"
    smoothing_str = (smoothing_factor > 0 ? " smoothing n=$smoothing_factor" : "")
    for (zk, zv) in sys["zones"]
        if zv["ZONETYPE"] ≠ "Ordered" || zv["I"] == 1 || zv["J"] > 1 || zv["K"] > 1
            continue
        end
        @info "Processing zone $zk"

        θ, origin, norm_vec = sphere_slice_zone_point_theta(zv)
        ignore_ind = Set{Int}([i for i in eachindex(θ) if θ[i] == θ[mod1(i+1,end)]])
        @debug "Ignoring duplicate point indices" sort(collect(ignore_ind)) θ

        ind = [i for i in eachindex(θ) if i ∉ ignore_ind]
        θ = θ[ind]
        order = sortperm(θ)
        θ = θ[order]
        θreg = -π : spacing : π
        θreg_padded = -2π : spacing : 2π
        plot_spacing = π/500
        θreg_plot = θ[1]+9plot_spacing : π/500 : θ[end]-9plot_spacing

        # need to get XYZ values of CPs found in the θ space
        xyz_of_theta = []
        for vk in ["X","Y","Z"]
            vv = zv["variables"][vk]
            θtmp = []
            vtmp = []
            for offset in [-2π, 0, 2π]
                θtmp = vcat(θtmp, θ .+ offset)
                vtmp = vcat(vtmp, vv[order])
            end
            itp_gridded = interpolate((θtmp,), vtmp, Gridded(Linear()))
            v = [itp_gridded[th] for th in θreg_padded]
            itp_cubic = scale(interpolate(v, BSpline(Cubic(Periodic(OnGrid())))), θreg_padded)
            push!(xyz_of_theta, itp_cubic)
        end
        var_cp_data = Dict()
        var_itp = Dict()
        var_itp_k = Dict()
        var_itp_g = Dict()
        var_itp_h = Dict()

        for (vk, vv) in zv["variables"]
            if !occursin(var_check_str, vk)
                continue
            end

            @info "Processing variable $vk"

            vk_short = replace(vk, var_check_str => "")

            var_cp_data[vk_short] = Dict()
            
            θtmp = Vector{Float32}()
            vtmp = Vector{Float32}()
            for offset in [-2π, 0, 2π]
                θtmp = vcat(θtmp, θ .+ offset)
                vtmp = vcat(vtmp, vv[order])
            end
            # first interpolate to regular grid, then use higher order interpolation for derivatives
            itp_gridded = interpolate((θtmp,), vtmp, Gridded(Linear()))
            # itp_k_gridded = interpolate((θ,), movmean(vv[order], max(1, smoothing_factor)), Gridded(Linear()))
            itp_k_gridded = interpolate((θtmp,), imfilter(vtmp, ImageFiltering.Kernel.gaussian((max(0,smoothing_factor),))), Gridded(Linear()))
            # function values on regular grid
            v = [itp_gridded[th] for th in θreg_padded]
            v_k = [itp_k_gridded[th] for th in θreg_padded]
            # quadratic and cubic interpolators using the interpolated grid values
            itp_quadratic = scale(interpolate(v, BSpline(Quadratic(Periodic(OnGrid())))), θreg_padded)
            itp_cubic = scale(interpolate(v, BSpline(Cubic(Periodic(OnGrid())))), θreg_padded)
            itp_k_quadratic = scale(interpolate(v_k, BSpline(Quadratic(Periodic(OnGrid())))), θreg_padded)
            itp_k_cubic = scale(interpolate(v_k, BSpline(Cubic(Periodic(OnGrid())))), θreg_padded)
            # print some values and compare the interpolation error
            for i in 5:length(θ)÷3:length(θ)-5
                θmid = (θ[i] + θ[mod1(i+1, end)]) / 2
                gridded = itp_gridded(θ[i])
                quad = itp_quadratic(θ[i])
                quad_error = quad - vv[order][i]
                cub = itp_cubic(θ[i])
                cub_error = cub - vv[order][i]
                @debug "Test interp for f($(θ[i])) = $(vv[order][i])" gridded quad quad_error cub cub_error itp_gridded(θmid) itp_quadratic(θmid) itp_cubic(θmid)
            end

            itp = itp_quadratic
            itp_k = itp_k_quadratic
            var_itp[vk] = itp
            var_itp_k[vk] = itp_k
            # now find roots using quadratic or cubic interpolation
            f(x) = gradient(itp, x)[1]
            g(x) = hessian(itp, x)[1]
            f_k(x) = gradient(itp_k, x)[1]
            g_k(x) = hessian(itp_k, x)[1]
            var_itp_g[vk] = f_k
            var_itp_h[vk] = g_k
            @debug "test grad hess" f(0) f_k(0) g(0) g_k(0)
            cps = nothing
            try
                cps = find_zeros(f_k, -π, π)
            catch e
                @error "Failed to run variable $vk for zone $zk. Skipping!"
                continue
            end
            @debug "Zeros" cps

            # save CP data and derivative values
            cp_info = []
            for (ri, r) in enumerate(cps)
                cp = Dict([
                    "#" => ri,
                    "variable" => vk,
                    "plane" => zk,
                    "theta" => r,
                    "x" => xyz_of_theta[1](r),
                    "y" => xyz_of_theta[2](r),
                    "z" => xyz_of_theta[3](r),
                    "f(xyz)" => itp(r),
                    "df/dtheta" => f_k(r),
                    "d2f/dtheta2" => g_k(r)   
                ])
                push!(cp_info, cp)
            end
            var_cp_data[vk_short]["cp_info"] = cp_info

            # make plot of function and CPs
            if smoothing_factor > 0
                p = plot(θreg_plot, [itp.(θreg_plot), itp_k.(θreg_plot)], title="$(sys["title"]) $zk", label=["raw" "gaussian$smoothing_str"], xtick=pitick(-π,π,4; mode=:latex))
                plot!(legend=:outerbottom, legendcolumns=3)
            else
                p = plot(θreg_plot, itp.(θreg_plot), title="$(sys["title"]) $zk", xtick=pitick(-π,π,4; mode=:latex))
                plot!(legend=:outerbottom, legendcolumns=2)
            end
            scatter!(p, cps, itp_k.(cps),label="CPs", mc=:cyan, ms=5, ma=0.4)
            for (ri,r) in enumerate(cps)
                annotate!(r, itp_k.(r), text(ri, :blue, 4))
            end
            xlabel!(p, "θ")
            ylabel!(p, vk)

            var_cp_data[vk_short]["plot_f"] = p

            # and grad and hess
            p = plot(θreg_plot, f_k.(θreg_plot), title="$(sys["title"]) $zk ∂/∂θ", xtick=pitick(-π,π,4; mode=:latex))
            scatter!(p, cps, f_k.(cps),label="CPs", mc=:cyan, ms=5, ma=0.4)
            for (ri,r) in enumerate(cps)
                annotate!(r, f_k.(r), text(ri, :blue, 4))
            end
            plot!(legend=:outerbottom, legendcolumns=2)
            xlabel!(p, "θ")
            ylabel!(p, "∂$vk/∂θ")
            var_cp_data[vk_short]["plot_g"] = p

            p = plot(θreg, g_k.(θreg), title="$(sys["title"]) $zk ∂²/∂θ²", xtick=pitick(-π,π,4; mode=:latex))
            scatter!(p, cps, g_k.(cps),label="CPs", mc=:cyan, ms=5, ma=0.4)
            for (ri,r) in enumerate(cps)
                annotate!(r, g_k.(r), text(ri, :blue, 4))
            end
            plot!(legend=:outerbottom, legendcolumns=2)
            xlabel!(p, "θ")
            ylabel!(p, "∂²$vk/∂θ²")
            var_cp_data[vk_short]["plot_h"] = p
        end
        sys["zones"][zk]["condensed_cp_info"] = var_cp_data
        sys["zones"][zk]["itp"] = var_itp
        sys["zones"][zk]["itp_k"] = var_itp_k
        sys["zones"][zk]["g_itp"] = var_itp_g
        sys["zones"][zk]["h_itp"] = var_itp_h
    end

    # Now we have all the CPs for all the INS variables of all the zones.
    # Next to match CPs across the zones by solving the linear sum assignment problem across zones pairwise
    var_cp_info = Dict()
    zones = [zk for (zk, zv) in sys["zones"] if "condensed_cp_info" ∈ keys(zv)]
    for vk in keys(sys["zones"][zones[1]]["condensed_cp_info"])
        @debug "Matching CPs for variable: $vk"
        cp_info = []
        for (i,zk1) in enumerate(zones[2:end])
            r1 = [[cp[v] for v in ["x","y","z"]] for cp in sys["zones"][zk1]["condensed_cp_info"][vk]["cp_info"]]
            for j in 1:i
                @debug "Matching CPs between zones $zk1 and $(zones[j])"
                r2 = [[cp[v] for v in ["x","y","z"]] for cp in sys["zones"][zones[j]]["condensed_cp_info"][vk]["cp_info"]]
                costs = [norm(ri - rj) for ri in r1, rj in r2]
                sol = find_best_assignment(costs)
                if length(sol.row4col) > length(sol.col4row)
                    costs = [costs[sol.col4row[i], i] for i in eachindex(sol.col4row)]
                    cost_order = partialsortperm(costs, 1:2)
                    @debug "closest two CP distances between zones" costs[cost_order]
                    for cj in cost_order
                        push!(cp_info, Dict(
                            zk1 => sys["zones"][zk1]["condensed_cp_info"][vk]["cp_info"][sol.col4row[cj]],
                            zones[j] => sys["zones"][zones[j]]["condensed_cp_info"][vk]["cp_info"][cj],
                            "distance" => costs[cj]
                        ))
                    end
                else
                    costs = [costs[i, sol.row4col[i]] for i in eachindex(sol.row4col)]
                    cost_order = partialsortperm(costs, 1:2)
                    @debug "closest two CP distances between zones" costs[cost_order]
                    for cj in cost_order
                        push!(cp_info, Dict(
                            zk1 => sys["zones"][zk1]["condensed_cp_info"][vk]["cp_info"][cj],
                            zones[j] => sys["zones"][zones[j]]["condensed_cp_info"][vk]["cp_info"][sol.row4col[cj]],
                            "distance" => costs[cj]
                        ))
                    end
                end
            end 
        end
        var_cp_info[vk] = cp_info
    end
    sys["condensed_cp_info"] = var_cp_info

    # Now CPs are matched between zones. Just need to print out all the data and plots.
    base_dir = dirname(file_path)
    out_dir = "$base_dir/out"
    mkpath(out_dir)
    cd(out_dir)

    # First CPs matched between different zones. 
    cp_info = Dict(["# matched" => [], "distance" => []])
    for zone_pair_cps in values(sys["condensed_cp_info"])
        cp_num = 1
        for zone_pair in zone_pair_cps
            for (zk,cp) in zone_pair
                if zk == "distance"
                    continue
                end
                push!(cp_info["distance"], zone_pair["distance"])
                push!(cp_info["# matched"], cp_num)
                for (var_name,value) in cp
                    if var_name in keys(cp_info)
                        push!(cp_info[var_name], value)
                    else
                        cp_info[var_name] = [value]
                    end
                end
            end
            cp_num += 1
        end
    end
    df = DataFrame(cp_info)
    sys["condensed_cp_dataframe_matched"] = df
    CSV.write("$(sys["title"]) matched CPs$smoothing_str.csv",df)

    # now all cps for each zone, and plots
    cp_info = Dict()
    for (zk,zv) in sys["zones"]
        if "condensed_cp_info" in keys(zv)
            for (vk, vv) in zv["condensed_cp_info"]
                plot_file_suffix = "$smoothing_str.pdf"
                plots_keys = [
                    ("plot_f", "itp", "$(sys["title"]) condensed CPs zone $zk var $vk"),
                    ("plot_g", "g_itp", "$(sys["title"]) condensed CPs zone $zk var $vk first deriv"),
                    ("plot_h", "h_itp", "$(sys["title"]) condensed CPs zone $zk var $vk second deriv"),
                ]
                for (pk, itp, ppath) in plots_keys
                    if pk ∉ keys(vv)
                        @info "Skipping plot $ppath, variable $vk because plot wasn't found"
                        continue
                    end
                    savefig(vv[pk], "$ppath$plot_file_suffix")
                    # then zoomed versions at cps
                    old_xlims = xlims(vv[pk])
                    old_ylims = ylims(vv[pk])
                    theta_last = -4π/8
                    for (cpi,cp) in enumerate(vv["cp_info"])
                        for (var_name,value) in cp
                            if var_name in keys(cp_info)
                                push!(cp_info[var_name], value)
                            else
                                cp_info[var_name] = [value]
                            end
                        end
                        if cp["theta"] < 4π/8 && cp["theta"] - theta_last >= π/8
                            try
                                new_xlims = (max(old_xlims[1],cp["theta"] - π/6), min(old_xlims[2], cp["theta"] + π/6))
                                test_vals = [zv[itp][var_check_str * vk](i) for i in new_xlims[1]:(new_xlims[2]-new_xlims[1])/200:new_xlims[2]]
                                new_ylims = (minimum(test_vals), maximum(test_vals))
                                buffer = (new_ylims[2] - new_ylims[1]) / 20
                                new_ylims = (new_ylims[1] - buffer, new_ylims[2] + buffer)
                                xlims!(vv[pk], new_xlims)
                                ylims!(vv[pk], new_ylims)
                                xticks!(vv[pk], pitick(new_xlims[1],new_xlims[2],16; mode=:latex))
                                savefig(vv[pk], "$(replace(ppath, "condensed CPs " => "")) condensed CP $cpi$plot_file_suffix")
                                theta_last = cp["theta"]
                            catch e
                                @error "error printing subplot"
                            end
                        end
                    end
                    xlims!(vv[pk], old_xlims)
                    ylims!(vv[pk], old_ylims)
                    xticks!(vv[pk], pitick(-π,π,4; mode=:latex))
                end
            end
        end
    end
    df = DataFrame(cp_info)
    sys["condensed_cp_dataframe_all"] = df
    CSV.write("$(sys["title"]) CPs$smoothing_str.csv",df)

    return sys
end

end
