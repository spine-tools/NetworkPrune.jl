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


function prune_network(db_url::String, prunned_db_url::String; alternative="Base")
    con__mon = []  # set of monitored and contingent line tuples that must be considered as defined by the
                   # connection_monitored and connection_contingency parameters
    monitored_lines = []
    ptdf_conn_n = Dict()  # ptdfs returned by PowerSystems.jl
    lodf_con_mon = Dict()  # lodfs calcuated based on ptdfs returned by PowerSystems.jl
    net_inj_nodes = []  # set of nodes with demand or generation

    objects = []
    object_parameter_values = []
    relationships = []
    relationship_parameter_values = []

    using_spinedb(db_url, SpineOpt)

    @info "pruning system at $(db_url)"    

    for c in commodity()        
        if commodity_physics(commodity=c) in (:commodity_physics_ptdf, :commodity_physics_lodf)

            @info "Processing network for commodity $(c) with network_physics $(commodity_physics(commodity=c))"
            n_islands, island_node = islands()
            @info "Your network consists of $(n_islands) islands"
            if n_islands > 1
                @warn "Your network consists of multiple islands, this may end badly."
                print(island_node)
            end

            net_inj_nodes = get_net_inj_nodes()  # returns list of nodes that have demand and/or generation

            @info "calculating ptdfs"
            ptdf_conn_n = calculate_ptdfs()

            if commodity_physics(commodity = c) == :commodity_physics_lodf
                con__mon = Tuple{Object,Object}[]
                @info "calculating lodfs" lodf_con_mon = calculate_lodfs(ptdf_conn_n, con__mon)
            end
        end
    end

    inj_nodes = []
    comm_nodes = []

    traversed = Dict{Object,Bool}()
    node__new_nodes = Dict{Object,Array}()
    min_v = 0
    for c in commodity()
        if commodity_physics(commodity=c) in (:commodity_physics_lodf, :commodity_physics_ptdf)
            for n in node__commodity(commodity=c)
                push!(comm_nodes, n)
                if !isempty(unit__to_node(node=n)) || demand(node=n) > 0
                    for ng in groups(n)
                        if minimum_voltage(node=ng) != nothing
                            min_v = minimum_voltage(node=ng)
                            break
                        end
                    end
                    if voltage(node=n) < min_v
                        push!(inj_nodes, n)
                    end
                end
            end
        end
    end

    # @info "Writing ptdf diagnostic file"
    # write_ptdfs(ptdf_conn_n, comm_nodes)

    @info "Traversing nodes"
    for n in inj_nodes
        node__new_nodes[n] = []
    end

    for n in inj_nodes
        for n2 in comm_nodes
            traversed[n2] = false
        end
        min_v = 0
        for ng in groups(n)
            if ! (minimum_voltage(node=ng) == nothing)
                min_v = minimum_voltage(node=ng)
                break
            end
        end
        traverse(n, n, traversed, node__new_nodes, min_v, ptdf_conn_n, 1)
    end

    to_prune_object_keys = []
    nodes_pruned = 0
    connections_pruned=0
    units_moved = 0
    units_distributed = 0
    demands_moved = 0
    demands_distributed = 0    
    for n in comm_nodes
        for ng in groups(n)            
            if ! (minimum_voltage(node=ng) == nothing)
                min_v = minimum_voltage(node=ng)
                if voltage(node=n) < min_v
                    push!(to_prune_object_keys, (n.class_name, n.name))                    
                    nodes_pruned += 1
                    for conn in connection__to_node(node=n)
                        if !((conn.class_name, conn.name) in to_prune_object_keys)
                            push!(to_prune_object_keys, (conn.class_name, conn.name))
                            connections_pruned += 1
                        end
                    end
                    for conn in connection__from_node(node=n)
                        if !((conn.class_name, conn.name) in to_prune_object_keys)
                            push!(to_prune_object_keys, (conn.class_name, conn.name))
                            connections_pruned += 1
                        end
                    end
                end
                break
            end            
        end
    end

    @info "Writing output"
    write_node__new_nodes(node__new_nodes)

    new_demand_dict = Dict{Object,Float64}()
    new_gen_dict = Dict{Object,Float64}()
    gens_to_move = Dict{Object,Object}()

    for (n, data) in node__new_nodes
        if demand(node=n) == nothing
            demand_to_shift = 0
        else
            demand_to_shift = demand(node=n)
        end
        if size(data, 1) == 1  # only one connected higher voltage node, move all the demand here
            if demand_to_shift > 0
                if haskey(new_demand_dict, data[1])
                    new_demand_dict[data[1][1]] = new_demand_dict[data[1][1]] + demand_to_shift
                else
                    new_demand_dict[data[1][1]] = demand_to_shift
                end
                demands_moved += 1
            end
            for u in unit__to_node(node=n)
                gens_to_move[u] = data[1][1]
                units_moved += 1
            end
        else
            gen_to_shift = 0
            for u in unit__to_node(node=n)
                gen_to_shift = gen_to_shift + unit_capacity(unit=u, node=n)
                push!(to_prune_object_keys, (u.class_name, u.name))                    
                units_distributed += 1
            end
            for (n2, ptdf) in data
                ptdf = abs(ptdf)
                if demand_to_shift > 0
                    if haskey(new_demand_dict, n2)
                        new_demand_dict[n2] = new_demand_dict[n2] + demand_to_shift * ptdf
                    else
                        new_demand_dict[n2] = demand_to_shift * ptdf
                    end
                    demands_distributed += 1
                end
                if haskey(new_gen_dict, n2)
                    new_gen_dict[n2] = new_gen_dict[n2] + gen_to_shift * ptdf
                else
                    new_gen_dict[n2] = gen_to_shift * ptdf
                end
            end
        end
    end

    # Update demand parameter of higher voltage nodes to add demand of pruned nodes

    for (n, new_demand) in new_demand_dict
        if demand(node=n) == nothing
            updated_demand = new_demand
        else
            updated_demand = demand(node=n) + new_demand
        end
        push!(object_parameter_values, ("node", string(n), "demand", new_demand))
    end

    for (u, new_node) in gens_to_move
        obj_list=[]
        push!(obj_list,string(u))
        push!(obj_list,string(new_node))
        push!(relationships,("unit__to_node",obj_list))
        for old_node in unit__to_node(unit=u)
            push!(
                relationship_parameter_values,
                ("unit__to_node", obj_list, "unit_capacity", unit_capacity(unit=u, node=old_node))
            )
            if !(minimum_operating_point(unit=u, node=old_node) == nothing)
                push!(
                    relationship_parameter_values,
                    ("unit__to_node", obj_list, "minimum_operating_point", minimum_operating_point(unit=u, node=old_node))
                )
            end
        end
    end

    for (n, new_gen) in new_gen_dict
        if new_gen > 0
            unit_name = string("U_DIST_", n)
            push!(objects, ("unit", unit_name))
            push!(object_parameter_values, ("unit", unit_name, "number_of_units", 1))
            push!(object_parameter_values, ("unit", unit_name, "online_variable_type", "unit_online_variable_type_binary"))
            push!(object_parameter_values, ("unit", unit_name, "unit_availability_factor", 1))

            obj_list=[]
            push!(obj_list,unit_name)
            push!(obj_list,string(n))
            push!(relationships,("unit__to_node",obj_list))
            push!(relationship_parameter_values,("unit__to_node", obj_list, "unit_capacity", new_gen))
        end
    end

    all_data = run_request(db_url, "export_data")
    run_request(prunned_db_url, "import_data", (all_data, ""))
    to_prune_object_ids = [
        x["id"]
        for x in run_request(prunned_db_url, "query", ("ext_object_sq",))["ext_object_sq"]
        if (Symbol(x["class_name"]), Symbol(x["name"])) in to_prune_object_keys
    ]
    run_request(prunned_db_url, "call_method", ("cascade_remove_items",), Dict(:object => to_prune_object_ids))

    @info "committing removed items to the void..."
    comment = "Network pruning: low voltage nodes removed"
    run_request(prunned_db_url, "call_method", ("commit_session", comment))

    object_parameter_values = [(opv..., alternative) for opv in object_parameter_values]
    relationship_parameter_values = [(opv..., alternative) for opv in relationship_parameter_values]

    added, err_log = import_data(
        prunned_db_url,
        "";
        objects=objects,
        relationships=relationships,
        object_parameter_values=object_parameter_values,
        relationship_parameter_values=relationship_parameter_values,
        alternatives=[alternative]
    )
    @info "added $(added) items"
    for err in err_log
        @info "import error: " err.msg
    end

    comment = "Network pruning: demand and generation shifts"
    run_request(prunned_db_url, "call_method", ("commit_session", comment))

    @info "network pruned successfully" nodes_pruned connections_pruned units_moved units_distributed demands_moved demands_distributed
end

function traverse(n_t, n, traversed, node__new_nodes, min_v, ptdf_conn_n, ptdf)
    # @info "traversing $(n) for $(n_t) with voltage $(voltage(node=n)) and min voltage $(min_v)"
    traversed[n] = true
    if voltage(node=n) >= min_v
        push!(node__new_nodes[n_t], (n, ptdf) )
    else
        for conn in connection__from_node(node=n)
            for n2 in connection__to_node(connection=conn)
                if !traversed[n2]
                    traverse(n_t, n2, traversed, node__new_nodes, min_v, ptdf_conn_n, ptdf_conn_n[(conn,n_t)] )
                end
            end
        end
        for conn in connection__to_node(node=n)
            for n2 in connection__from_node(connection=conn)
                if !traversed[n2]
                    traverse(n_t, n2, traversed, node__new_nodes, min_v, ptdf_conn_n, ptdf_conn_n[(conn,n_t)])
                end
            end
        end
    end
end

function write_node__new_nodes(node__new_nodes)
    io = open("node_mapping.csv", "w")
    print(io, "node,V,total_gens,total_generation,total_demand,node1,V1,ptdf1,node2,V2,ptdf2,node3,V3,ptdf3,node4,V4,ptdf4,node5,V5,ptdf5\n")
    for n in node__new_nodes
        total_generators=0
        total_generation=0
        total_demand=0
        for u in unit__to_node(node=n[1])
            total_generators = total_generators + 1
            total_generation = total_generation + unit_capacity(unit=u, node=n[1])
        end
        if demand(node=n[1]) == nothing
            total_demand = 0
        else
            total_demand = demand(node=n[1])
        end
        print(io, string(n[1]), ",", voltage(node=n[1]), ",", total_generators, ",", total_generation, ",", total_demand)
        if size(n[2],1) == 1
            print(io, ",", string(n[2][1][1]), ",", voltage(node=n[2][1][1]), ",", 1)
        else
            for n2 in n[2]
                print(io, ",", string(n2[1]), ",", voltage(node=n2[1]), ",", string(abs(n2[2])))
            end
        end
        print(io, "\n")
    end
    close(io)
end


"""
    islands()

Determines the number of islands in a commodity network - used for diagnostic purposes

"""
function islands()
    visited_d = Dict{Object,Bool}()
    island_node = Dict{Int64,Array}()
    island = 0

    for c in commodity()
        if commodity_physics(commodity=c) in(:commodity_physics_lodf, :commodity_physics_ptdf)
            for n in node__commodity(commodity=c)
                visited_d[n] = false
            end
        end
    end

    for c in commodity()
        if commodity_physics(commodity=c) in(:commodity_physics_lodf, :commodity_physics_ptdf)
            for n in node__commodity(commodity=c)
                if !visited_d[n]
                    island = island + 1
                    island_node[island] = Object[]
                    visit(n, island, visited_d, island_node)
                end
            end
        end
    end
    return island, island_node
end

"""
    visit()

Function called recursively to visit nodes in the network to determine number of islands

"""
function visit(n, island, visited_d, island_node)
    visited_d[n] = true
    push!(island_node[island], n)
    for conn in connection__from_node(node=n)
        for n2 in connection__to_node(connection=conn)
            if !visited_d[n2]
                visit(n2, island, visited_d, island_node)
            end
        end
    end
    for conn in connection__to_node(node=n)
        for n2 in connection__from_node(connection=conn)
            if !visited_d[n2]
                visit(n2, island, visited_d, island_node)
            end
        end
    end
end

"""
    check_x()

Check for low reactance values

"""
function check_x()
    @info "Checking reactances"
    for conn in connection()
        if conn_reactance(connection=conn) < 0.0001
            @info "Low reactance may cause problems for line " conn
        end
    end
end


"""
    calculate_ptdfs()

Returns a dict indexed on tuples of (connection, node) containing the ptdfs of the system currently in memory.

"""
function calculate_ptdfs()
    ps_busses = Bus[]
    ps_lines = Line[]

    node_ps_bus = Dict{Object,Bus}()
    i = 1
    for c in commodity()
        if commodity_physics(commodity=c) in(:commodity_physics_lodf, :commodity_physics_ptdf)
            for n in node__commodity(commodity=c)
                if node_opf_type(node=n) == :node_opf_type_reference
                    bustype = BusTypes.REF
                else
                    bustype = BusTypes.PV
                end
                ps_bus = Bus(
                    number = i,
                    name = string(n),
                    bustype = bustype,
                    angle = 0.0,
                    #voltage = 0.0,
                    magnitude = 0.0,
                    voltage_limits = (min = 0.0, max = 0.0),
                    base_voltage = nothing,
                    area = nothing,
                    load_zone = LoadZone(nothing),
                    ext = Dict{String, Any}()
                )

                push!(ps_busses,ps_bus)
                node_ps_bus[n] = ps_bus
                i = i + 1
            end

            # InfrastructureSystems.buscheck(ps_busses)
            # InfrastructureSystems.slack_bus_check(ps_busses)

            for conn in connection()
                for n_from in connection__from_node(connection=conn)
                    for n_to in connection__to_node(connection=conn)
                        for c in node__commodity(node=n_from)
                            if commodity_physics(commodity=c) in(:commodity_physics_lodf, :commodity_physics_ptdf)
                                ps_arc = Arc(node_ps_bus[n_from], node_ps_bus[n_to])
                                new_line = Line(;
                                    name = string(conn),
                                    available = true,
                                    active_power_flow = 0.0,
                                    reactive_power_flow = 0.0,
                                    arc = ps_arc,
                                    r = connection_resistance(connection=conn),
                                    x = max(connection_reactance(connection=conn), 0.00001),
                                    b = (from=0.0, to=0.0),
                                    rate = 0.0,
                                    angle_limits = (min = 0.0, max = 0.0)
                                )
                                push!(ps_lines,new_line)
                            end   # in case there are somehow multiple commodities
                            break
                        end

                    end
                end
            end
        end
    end

    ps_ptdf=PowerSystems.PTDF(ps_lines,ps_busses)

    ptdf=Dict{Tuple{Object,Object},Float64}()

    for c in commodity()
        if commodity_physics(commodity=c) in(:commodity_physics_lodf, :commodity_physics_ptdf)
            for n in node__commodity(commodity=c)
                for conn in connection()
                    ptdf[(conn,n)] = ps_ptdf[string(conn), node_ps_bus[n].number]
                end
            end
        end
    end

    return ptdf

    # buildlodf needs to be updated to account for cases
    # lodfs=PowerSystems.buildlodf(ps_lines,ps_busses)
end


"""
    calculate_lodfs(ptdf_b_n, b_con__b_mon)

Returns lodfs for the system specified by ptdf_b_n ,b_con__b_mon as a dict of tuples: contingent_line, monitored_line

"""
# This function takes a long time.
# PowerSystems has a function that does it faster using linear algebra but doesn't handle the case of tails like
# I would like.

function calculate_lodfs(ptdf_conn_n, con__mon)
    lodf_con_mon = Dict{Tuple{Object,Object},Float64}()
    considered_contingencies = 0
    skipped = 0
    tolerance = 0
    for conn_con in connection()
        if connection_contingency(connection = conn_con) == 1
            for n_from in connection__from_node(connection=conn_con)
                for n_to in connection__to_node(connection=conn_con)
                    demoninator = 1 - (ptdf_conn_n[(conn_con, n_from)]-ptdf_conn_n[(conn_con, n_to)])
                    if abs(demoninator) < 0.001
                        demoninator = -1
                    end
                    for conn_mon in connection()
                        if connection_monitored(connection=conn_mon) == 1 && conn_con != conn_mon
                            if demoninator == -1
                                lodf_trial = (ptdf_conn_n[(conn_mon, n_from)] - ptdf_conn_n[(conn_mon, n_to)]) / demoninator
                            else
                                lodf_trial = -ptdf_conn_n[(conn_mon, n_from)]
                            end
                            c = first(indices(commodity_lodf_tolerance))
                            tolerance = commodity_lodf_tolerance(commodity=c)
                            if abs(lodf_trial) > tolerance
                                considered_contingencies = considered_contingencies + 1
                                push!(con__mon, (conn_con, conn_mon))
                                lodf_con_mon[(conn_con, conn_mon)] = lodf_trial
                            else
                                skipped = skipped + 1
                            end
                        end
                    end
                end
            end
        end
    end
    # @info "Contingencies summary " considered_contingencies skipped tolerance
    return lodf_con_mon
end

function get_net_inj_nodes()
    net_inj_nodes = []
    for c in commodity()
        if commodity_physics(commodity=c) in(:commodity_physics_lodf, :commodity_physics_ptdf)
            for n in node__commodity(commodity=c)
                for u in unit__to_node(node=n)
                    if !(n in net_inj_nodes)
                        push!(net_inj_nodes, n)
                    end
                end
                for u in unit__from_node(node=n)
                    if !(n in net_inj_nodes)
                        push!(net_inj_nodes, n)
                    end
                end
                for ng in groups(n)
                    if fractional_demand(node=n) > 0 || demand(node=n) > 0
                        if !(n in net_inj_nodes)
                            push!(net_inj_nodes, n)
                        end
                    end
                end
            end
        end
    end
    return net_inj_nodes
end


function write_ptdfs(ptdfs, net_inj_nodes)
    io = open("ptdfs.csv", "w")
    print(io, "connection,")
    for n in net_inj_nodes
        print(io, string(n), ",")
    end
    print(io, "\n")
    for conn in connection()
        print(io, string(conn), ",")
        for n in net_inj_nodes
            if haskey(ptdfs, (conn,n))
                print(io, ptdfs[(conn,n)], ",")
            else
                print(io, "_nan", ",")
            end
        end
        print(io, "\n")
    end
    close(io)
end
