module Condensed_CP_Analysis

import Interpolations
import LinearAlgebra
import LoggingExtras, LoggingFormats
import NLsolve
import Parsers
import Interact
import Plots

using Interact

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

    while startswith(line, "ZONE") && !eof(f)
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
        end

        # save matrix of xyz values for convenience
        nodes = Array{Float32}(undef, zone["NUMPOINTS"], 3)
        for (idir,dir) in enumerate(["X", "Y", "Z"])
            vals = zone["variables"][dir]
            for i in 1:zone["NUMPOINTS"]
                nodes[i, idir] = vals[i]
            end
        end

        zone["nodes_xyz"] = nodes
        
        out["zones"][zone["T"]] = zone
    end

    close(f)
    return out
end

function data_explorer(sys)
    @manipulate for zone = collect(keys(sys["zones"]))
        zone
    end
end

end
