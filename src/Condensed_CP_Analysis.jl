module Condensed_CP_Analysis

import Interpolations: scale, interpolate, BSpline, Linear, Quadratic, Cubic, Flat, Periodic, OnGrid, gradient, hessian, Gridded
import LinearAlgebra: dot, cross, normalize, norm
import LoggingExtras, LoggingFormats
import NLsolve
import Parsers
import Interact
import Plots: plot, plot!, scatter, scatter!
import Roots: find_zeros, Newton
using StaticArrays
using Meshes, MeshViz
import GLMakie as Mke
import SplitApplyCombine: invert
import NaNStatistics: movmean

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

function sphere_slice_zone_point_theta(zone, zero_theta_vec)
    θ0 = zero_theta_vec
    origin = [(minimum(zone["variables"][dir]) + maximum(zone["variables"][dir])) / 2 for dir in ["X","Y","Z"]]
    points = zone["nodes_xyz"]
    norm_vec = normalize(cross(points[1] - origin, points[div(end,3)] - origin))
    return [angle(normalize(points[i] - origin), θ0, norm_vec) for i in eachindex(points)], origin, norm_vec
end

function sphere_slice_analysis(file_path, zero_theta_vec; smoothing_factor=1)
    sys = import_dat(file_path)
    @info "Processing file $(sys["title"])"
    for (zk, zv) in sys["zones"]
        if zv["ZONETYPE"] ≠ "Ordered" || zv["I"] == 1 || zv["J"] > 1 || zv["K"] > 1
            continue
        end
        @info "Processing zone $zk"

        θ, origin, norm_vec = sphere_slice_zone_point_theta(zv, zero_theta_vec)
        ignore_ind = Set{Int}([i for i in eachindex(θ) if θ[i] == θ[mod1(i+1,end)]])
        @debug "Ignoring duplicate point indices" sort(collect(ignore_ind)) θ

        ind = [i for i in eachindex(θ) if i ∉ ignore_ind]
        θ = θ[ind]
        order = sortperm(θ)
        θ = θ[order]
        θreg = θ[1] : 2π/length(θ) : θ[end]

        # need to get XYZ values of CPs found in the θ space
        xyz_of_theta = []
        for vk in ["X","Y","Z"]
            vv = zv["variables"][vk]
            itp_gridded = interpolate((θ,), vv[order], Gridded(Linear()))
            v = [itp_gridded[th] for th in θreg]
            itp_cubic = scale(interpolate(v, BSpline(Cubic(Flat(OnGrid())))), θreg)
            push!(xyz_of_theta, itp_cubic)
        end

        for (vk, vv) in zv["variables"]
            if !occursin("INS: Electron", vk)
                continue
            end

            @info "Processing variable $vk"

            # first interpolate to regular grid, then use higher order interpolation for derivatives
            itp_gridded = interpolate((θ,), vv[order], Gridded(Linear()))
            # function values on regular grid
            v = [itp_gridded[th] for th in θreg]
            # quadratic and cubic interpolators using the interpolated grid values
            itp_quadratic = scale(interpolate(v, BSpline(Quadratic(Flat(OnGrid())))), θreg)
            itp_cubic = scale(interpolate(v, BSpline(Cubic(Flat(OnGrid())))), θreg)
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
            # now find roots using quadratic or cubic interpolation
            f(x) = gradient(itp, x)[1]
            g(x) = hessian(itp, x)[1]
            @info "test grad hess" f(0) g(0)
            cub_zeros = find_zeros(f, -3, 3)
            @info "Zeros" cub_zeros

            p = plot(θreg, [itp.(θreg), movmean(itp.(θreg), max(1,smoothing_factor))], title=zk, label=[vk "moving average n=$smoothing_factor"])
            return p
            break
        end
        break
    end
end

end
