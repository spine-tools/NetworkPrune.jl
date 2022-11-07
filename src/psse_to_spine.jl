#############################################################################
# Copyright (C) 2017 - 2022  Spine Project
#
# This file is part of NetworkPrune.
#
# NetworkPrune is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# NetworkPrune is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#############################################################################


struct FilteredIO <: IO
    inner::IO
    flt
end

Base.readlines(io::FilteredIO) = filter!(io.flt, readlines(io.inner))

"""
    psse_to_spine(psse_path, db_url)

Parse the psse raw file (`psse_path`) using PowerModels.jl
and create a SpineOpt model at `db_url` using `nodes`, `units` and `connections`.
"""
function psse_to_spine(psse_path, db_url::String; skip=(), bus_codes=Dict())
    pm_data = open(psse_path) do io
        filtered_io = FilteredIO(io, line -> !startswith(line, "@!") && !isempty(strip(line)))
        PowerModels.parse_psse(filtered_io)
    end
    psse_to_spine(pm_data, db_url; skip=skip, bus_codes=bus_codes)
end
function psse_to_spine(ps_system::Dict, db_url::String; skip=(), bus_codes=Dict())

    objects = []
    object_groups = []
    relationships = []
    object_parameter_values = []
    relationship_parameter_values = []

    object_parameters = [
        ("node", "minimum_voltage"),
        ("node", "psse_bus_name"),
        ("node", "bus_code"),
        ("node", "voltage")
    ]

    commodity_name = "elec"

    baseMVA = ps_system["baseMVA"]

    areas = []
    zones = []

    push!(objects, ["commodity", commodity_name])

    push!(object_parameter_values, ("commodity", commodity_name, "commodity_lodf_tolerance", 0.1))
    push!(object_parameter_values, ("commodity", commodity_name, "commodity_physics", "commodity_physics_ptdf"))
    push!(object_parameter_values, ("commodity", commodity_name, "commodity_ptdf_flow_tolerance", 0.1))
    push!(object_parameter_values, ("commodity", commodity_name, "commodity_ptdf_threshold", 0.0001))
    push!(object_parameter_values, ("commodity", commodity_name, "commodity_slack_penalty", 1000))

    node_demand = Dict{Int64,Float64}()
    node_name = Dict{Int64,String}()

    for b in ps_system["bus"]
        data = b[2]
        i = data["bus_i"]
        node_demand[i] = 0
        psse_bus_name = strip(data["name"])
        bus_code_raw = isempty(bus_codes) ? psse_bus_name[1:3] : get(bus_codes, psse_bus_name, nothing)
        if ismissing(bus_code_raw) || isnothing(bus_code_raw)
            @warn "couldn't find bus code for $psse_bus_name"
            bus_code_raw = "MISSING_BUS_CODE"
        end
        bus_code = replace(bus_code_raw, "-" => "_")
        voltage_level = Int(round(data["base_kv"], digits=0))
        node_name[i] = name = join([i, bus_code, voltage_level], "_")
        push!(objects, ("node", name))
        push!(object_parameter_values, ("node", name, "voltage", data["base_kv"]))
        push!(object_parameter_values, ("node", name, "psse_bus_name", psse_bus_name))        
        push!(object_parameter_values, ("node", name, "bus_code", bus_code))        

        if data["bus_type"] == 3
            push!(object_parameter_values, ("node", name, "node_opf_type", "node_opf_type_reference"))
        end

        area_name = string("area_", data["area"])
        zone_name = string("zone_", data["zone"])
        if !(area_name in areas)
            push!(areas, area_name)
            push!(objects, ("node", area_name))
            push!(object_parameter_values, ("node", area_name, "minimum_voltage", 110.0))
        end        
        push!(object_groups, ["node", area_name, name])

        if !(zone_name in zones)
            push!(zones, zone_name)
            push!(objects, ("node", zone_name))
        end        
        push!(object_groups, ["node", zone_name, name])

        obj_list = []
        push!(obj_list, name)
        push!(obj_list, "elec")
        push!(relationships,("node__commodity", obj_list))
    end

    for b in ("branch" in skip ? () : ps_system["branch"])
        data = b[2]
        if data["br_status"] in (0, 1)
            from_bus_name = node_name[data["f_bus"]]
            to_bus_name =  node_name[data["t_bus"]]
            ckt = rstrip(string(data["source_id"][4]))
            name = string(from_bus_name, "__", to_bus_name, "_", ckt)
            connection = ("connection", name)
            push!(objects, connection)
            br_r = data["br_r"]
            push!(object_parameter_values, ("connection", name, "connection_resistance", isempty(br_r) ? 0 : br_r))
            push!(object_parameter_values, ("connection", name, "connection_reactance", data["br_x"]))
            push!(object_parameter_values, ("connection", name, "connection_monitored", 1))
            push!(object_parameter_values, ("connection", name, "connection_contingency", 1))
            push!(object_parameter_values, ("connection", name, "connection_availability_factor", 1.0))
            obj_list = []
            push!(obj_list, name)
            push!(obj_list, to_bus_name)
            push!(relationships, ("connection__to_node", obj_list))
            rate_a = round(get(data, "rate_a", 0) * baseMVA, digits=2)
            if haskey(data, "rate_b")
                rate_b = round(data["rate_b"] * baseMVA, digits=2)
            else
                rate_b = rate_a
            end

            push!(relationship_parameter_values, ("connection__to_node", obj_list, "connection_capacity", rate_a))
            push!(
                relationship_parameter_values,
                ("connection__to_node", obj_list, "connection_emergency_capacity", rate_b)
            )

            obj_list = []
            push!(obj_list, name)
            push!(obj_list, from_bus_name)
            push!(relationships, ("connection__from_node", obj_list))
            push!(relationship_parameter_values, ("connection__from_node", obj_list, "connection_capacity", rate_a))
            push!(
                relationship_parameter_values,
                ("connection__from_node", obj_list, "connection_emergency_capacity", rate_b)
            )
        end
    end

    for dc in ("dcline" in skip ? () : ps_system["dcline"])
        data = dc[2]
        from_bus_name = node_name[data["source_id"][2]]
        to_bus_name = node_name[data["source_id"][3]]
        ckt = 1
        name = string(from_bus_name, "__", to_bus_name, "_", ckt)
        connection = ("connection", name)
        pmaxt = round(data["pmaxt"]  * baseMVA, digits=2)
        pmaxf = round(data["pmaxf"]  * baseMVA, digits=2)

        push!(objects, connection)
        push!(object_parameter_values, ("connection", name, "connection_resistance", 0))
        push!(object_parameter_values, ("connection", name, "connection_reactance", 0.0001))
        push!(object_parameter_values, ("connection", name, "connection_monitored", 0))
        push!(object_parameter_values, ("connection", name, "connection_contingency", 0))
        push!(object_parameter_values, ("connection", name, "connection_availability_factor", 1.0))

        obj_list = []
        push!(obj_list, name)
        push!(obj_list, to_bus_name)
        push!(relationships,("connection__to_node",obj_list))
        push!(relationship_parameter_values, ("connection__to_node", obj_list, "connection_capacity", pmaxt))
        push!(relationship_parameter_values, ("connection__to_node", obj_list, "connection_emergency_capacity", pmaxt))

        obj_list = []
        push!(obj_list, name)
        push!(obj_list, from_bus_name)
        push!(relationships,("connection__from_node",obj_list))
        push!(relationship_parameter_values, ("connection__from_node", obj_list, "connection_capacity", pmaxf))
        push!(relationship_parameter_values, ("connection__from_node", obj_list, "connection_emergency_capacity", pmaxf))
    end

    gen_ids = []
    for g in ("gen" in skip ? () : ps_system["gen"])
        data = g[2]
        bus_name = node_name[data["gen_bus"]]
        name = string(ps_system["bus"][string(data["gen_bus"])]["name"][1:3], "_", rstrip(data["source_id"][3]))
        if name in gen_ids
            sub_index = 1
            new_name = string(name, "_", sub_index)
            while new_name in gen_ids
                sub_index = sub_index + 1
                new_name = string(name, "_", sub_index)
            end
            name = new_name
        end
        push!(gen_ids, name)

        unit = ("unit", name)
        push!(objects, unit)
        push!(object_parameter_values, ("unit", name, "number_of_units", 1))
        push!(object_parameter_values, ("unit", name, "online_variable_type", "unit_online_variable_type_binary"))
        push!(object_parameter_values, ("unit", name, "unit_availability_factor", 1))

        obj_list = []
        push!(obj_list, name)
        push!(obj_list, bus_name)
        push!(relationships, ("unit__to_node", obj_list))

        push!(
            relationship_parameter_values,
            ("unit__to_node", obj_list, "unit_capacity", round(data["pmax"] * baseMVA, digits=2))
        )
        pmin = round(data["pmin"] / data["pmax"], digits=4)
        if pmin isa Number
            push!(relationship_parameter_values, ("unit__to_node", obj_list, "minimum_operating_point", pmin))
        end
    end

    for l in ps_system["load"]
        data = l[2]
        if data["status"] == 1
            node_demand[data["load_bus"]] = node_demand[data["load_bus"]] + data["pd"]
        end
    end

    for b in ps_system["bus"]
        data = b[2]
        name = node_name[data["bus_i"]]
        load_contender = round(baseMVA * node_demand[data["bus_i"]], digits=4)
        if load_contender > 0
            push!(object_parameter_values, ("node", name, "demand", load_contender))
        end
    end
    @info "writing PSSE data to $(db_url)"

    added, err_log = import_data(db_url, SpineOpt.template(), "Load SpineOpt template")
    if !isempty(err_log)
        @error join(err_log, "\n")
    end

    @info "importing data to $(db_url)"

    comment = "Import powersystems to Spine"
    added, err_log = import_data(
        db_url,
        comment;
        objects=objects,
        object_groups=object_groups,
        object_parameters=object_parameters,
        relationships=relationships,
        object_parameter_values=object_parameter_values,
        relationship_parameter_values=relationship_parameter_values
    )
    if !isempty(err_log)
        @error join(err_log, "\n")
    end

    @info "data imported to $(db_url)"

end
