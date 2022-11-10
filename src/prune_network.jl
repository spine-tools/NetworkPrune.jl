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

function prune_network(
    db_url::String, prunned_db_url::String; alternative="Base", node_mapping_file_name="node_mapping.csv"
)
    I = Module()
    @eval I using SpineInterface
    using_spinedb(db_url, I)
    comm = first(
        c for c in I.commodity()
        if I.commodity_physics(commodity=c) in (:commodity_physics_ptdf, :commodity_physics_lodf)
    )
    @info "Processing network for commodity $(comm) with network_physics $(I.commodity_physics(commodity=comm))"
    n_islands, island_node = islands(I, comm)
    @info "Your network consists of $(n_islands) islands"
    if n_islands > 1
        @warn "Your network consists of multiple islands, this may end badly."
        print(island_node)
    end
    @info "Calculating ptdfs"
    ptdf_conn_n = calculate_ptdfs(I, comm)
    lodf_con_mon, con__mon = if I.commodity_physics(commodity=comm) == :commodity_physics_lodf
        @info "Calculating lodfs" calculate_lodfs(I, ptdf_conn_n)
    else
        Dict(), []
    end
    @info "Pruning system at $(db_url)"    
    objects = []
    object_parameter_values = []
    relationships = []
    relationship_parameter_values = []
    inj_nodes = []
    comm_nodes = []
    min_v = 0
    for n in I.node__commodity(commodity=comm)
        push!(comm_nodes, n)
        # isempty(I.unit__to_node(node=n)) && (isnothing(I.demand(node=n)) || iszero(I.demand(node=n))) && continue
        # for some reason, I was only mapping nodes with load or generation, but we want to map all nodes, right?
        for ng in groups(n)
            if I.minimum_voltage(node=ng) != nothing
                min_v = I.minimum_voltage(node=ng)
                break
            end
        end
        if I.voltage(node=n) < min_v
            push!(inj_nodes, n)
        end
    end
    # @info "Writing ptdf diagnostic file"
    # write_ptdfs(I, ptdf_conn_n, comm_nodes)
    @info "Traversing nodes"
    traversed = Dict{Object,Bool}()
    node__new_nodes = Dict(n => [] for n in inj_nodes)
    for n in inj_nodes
        for n2 in comm_nodes
            traversed[n2] = false
        end
        min_v = 0
        for ng in groups(n)
            if ! (I.minimum_voltage(node=ng) == nothing)
                min_v = I.minimum_voltage(node=ng)
                break
            end
        end
        traverse(I, n, n, traversed, node__new_nodes, min_v, ptdf_conn_n, 1)
    end
    to_prune_object_keys = []
    nodes_pruned = 0
    connections_pruned = 0
    units_moved = 0
    units_distributed = 0
    demands_moved = 0
    demands_distributed = 0    
    for n in comm_nodes
        for ng in groups(n)            
            if I.minimum_voltage(node=ng) != nothing
                min_v = I.minimum_voltage(node=ng)
                if I.voltage(node=n) < min_v
                    push!(to_prune_object_keys, (n.class_name, n.name))                    
                    nodes_pruned += 1
                    for conn in I.connection__to_node(node=n)
                        if !((conn.class_name, conn.name) in to_prune_object_keys)
                            push!(to_prune_object_keys, (conn.class_name, conn.name))
                            connections_pruned += 1
                        end
                    end
                    for conn in I.connection__from_node(node=n)
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
    new_demand_dict = Dict{Object,Float64}()
    new_gen_dict = Dict{Object,Float64}()
    gens_to_move = Dict{Object,Object}()
    for (n, new_nodes) in node__new_nodes
        if I.demand(node=n) == nothing
            demand_to_shift = 0
        else
            demand_to_shift = I.demand(node=n)
        end
        if size(new_nodes, 1) == 1  # only one connected higher voltage node, move all the demand here
            n2, _ptdf = new_nodes[1]
            if demand_to_shift > 0
                if haskey(new_demand_dict, n2)
                    new_demand_dict[n2] += demand_to_shift
                else
                    new_demand_dict[n2] = demand_to_shift
                end
                demands_moved += 1
            end
            for u in I.unit__to_node(node=n)
                gens_to_move[u] = n2
                units_moved += 1
            end
        else
            gen_to_shift = 0
            for u in I.unit__to_node(node=n)
                gen_to_shift = gen_to_shift + I.unit_capacity(unit=u, node=n)
                push!(to_prune_object_keys, (u.class_name, u.name))                    
                units_distributed += 1
            end
            for (n2, ptdf) in new_nodes
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
        if I.demand(node=n) == nothing
            updated_demand = new_demand
        else
            updated_demand = I.demand(node=n) + new_demand
        end
        push!(object_parameter_values, ("node", string(n), "demand", new_demand))
    end
    for (u, new_node) in gens_to_move
        rel = [string(u), string(new_node)]
        push!(relationships, ("unit__to_node", rel))
        for old_node in I.unit__to_node(unit=u)
            push!(
                relationship_parameter_values,
                ("unit__to_node", rel, "unit_capacity", I.unit_capacity(unit=u, node=old_node))
            )
            if I.minimum_operating_point(unit=u, node=old_node) != nothing
                push!(
                    relationship_parameter_values,
                    ("unit__to_node", rel, "minimum_operating_point", I.minimum_operating_point(unit=u, node=old_node))
                )
            end
        end
    end
    for (n, new_gen) in new_gen_dict
        if new_gen > 0
            unit_name = string("U_DIST_", n)
            push!(objects, ("unit", unit_name))
            push!(object_parameter_values, ("unit", unit_name, "number_of_units", 1))
            push!(
                object_parameter_values, ("unit", unit_name, "online_variable_type", "unit_online_variable_type_binary")
            )
            push!(object_parameter_values, ("unit", unit_name, "unit_availability_factor", 1))
            rel = [unit_name, string(n)]
            push!(relationships,("unit__to_node", rel))
            push!(relationship_parameter_values,("unit__to_node", rel, "unit_capacity", new_gen))
        end
    end
    all_data = run_request(db_url, "export_data")
    run_request(prunned_db_url, "import_data", (all_data, ""))
    object_parameter_values = [(opv..., alternative) for opv in object_parameter_values]
    relationship_parameter_values = [(opv..., alternative) for opv in relationship_parameter_values]
    data_to_import = Dict(
        :objects => objects,
        :relationships => relationships,
        :object_parameter_values => object_parameter_values,
        :relationship_parameter_values => relationship_parameter_values,
        :alternatives => [alternative]
    )
    comment = "Network pruning: demand and generation shifts"
    _prune_and_import(prunned_db_url, to_prune_object_keys, data_to_import, comment)
    @info "Network pruned successfully" nodes_pruned connections_pruned units_moved units_distributed demands_moved demands_distributed
    k = 1
    while true
        @info "Trimming tails - pass $k"
        trim_tails(prunned_db_url, node__new_nodes; alternative=alternative) || break
        k += 1
    end
    @info "Writing node mapping"    
    path = joinpath(pwd(), node_mapping_file_name)
    write_node__new_nodes(I, node__new_nodes, path)
    @info "Node mapping written at $path"
end

function traverse(I, n_t, n, traversed, node__new_nodes, min_v, ptdf_conn_n, ptdf)
    # @info "traversing $(n) for $(n_t) with voltage $(voltage(node=n)) and min voltage $(min_v)"
    traversed[n] = true
    if I.voltage(node=n) >= min_v
        push!(node__new_nodes[n_t], (n, ptdf))
    else
        for conn in I.connection__from_node(node=n)
            for n2 in I.connection__to_node(connection=conn)
                if !traversed[n2]
                    traverse(I, n_t, n2, traversed, node__new_nodes, min_v, ptdf_conn_n, ptdf_conn_n[(conn, n_t)])
                end
            end
        end
        for conn in I.connection__to_node(node=n)
            for n2 in I.connection__from_node(connection=conn)
                if !traversed[n2]
                    traverse(I, n_t, n2, traversed, node__new_nodes, min_v, ptdf_conn_n, ptdf_conn_n[(conn, n_t)])
                end
            end
        end
    end
end

"""
If a node is connected to only one other node:
- If there is no demand or generation, delete then node and connection regardless of reactance
- If there is only a demand, move the demand and delete the node and connection regardless of reactance
- If there is generation:
    - If the reactance is lower or equal than 0.0001, then move the generation and delete the node and connection.
    - If the reactance is greater than 0.0001, then do nothing
"""
function trim_tails(prunned_db_url::String, node__new_nodes; alternative="Base")
    P = Module()
    @eval P using SpineInterface
    using_spinedb(prunned_db_url, P)
    to_remove = Dict()
    comm = first(
        c for c in P.commodity()
        if P.commodity_physics(commodity=c) in (:commodity_physics_ptdf, :commodity_physics_lodf)
    )
    comm_nodes = P.node__commodity(commodity=comm)
    tail_conn_next_tuples = []
    for conn in P.connection()
        nodes_from = [n for n in P.connection__from_node(connection=conn)]
        nodes_to = [n for n in P.connection__to_node(connection=conn)]
        if isempty(nodes_from) || isempty(nodes_to)
            push!(get!(to_remove, :connection, []), conn)
            continue
        end
        if length(nodes_from) == length(nodes_to) == 1
            n_from = nodes_from[1]
            n_to = nodes_to[1]
            n_from_other_conns = [
                other_conn
                for connection__node in (P.connection__to_node, P.connection__from_node)
                for other_conn in connection__node(node=n_from)
                if other_conn != conn
            ]
            n_to_other_conns = [
                other_conn
                for connection__node in (P.connection__to_node, P.connection__from_node)
                for other_conn in connection__node(node=n_to)
                if other_conn != conn
            ]
            if isempty(n_from_other_conns)
                push!(tail_conn_next_tuples, (n_from, conn, n_to))
            elseif isempty(n_to_other_conns)
                push!(tail_conn_next_tuples, (n_to, conn, n_from))
            end
        end
    end
    filter!(x -> x[1] in comm_nodes && x[3] in comm_nodes, tail_conn_next_tuples)
    rels = []
    opvs = []
    new_demands = Dict()
    for (tail_node, conn, next_node) in tail_conn_next_tuples
        tail_demand = P.demand(node=tail_node)
        units = P.unit__to_node(node=tail_node)
        react = P.connection_reactance(connection=conn)
        if isempty(units) || react <= 0.0001
            # remove tail and conn
            push!(get!(to_remove, :node, []), tail_node)            
            push!(get!(to_remove, :connection, []), conn)
            # save mapping
            push!(get!(node__new_nodes, tail_node, []), (next_node, 1))
            # move any tail demand to next
            if !isnothing(tail_demand) && !iszero(tail_demand)
                push!(get!(new_demands, next_node, []), tail_demand)
            end
            # move any units to next
            if !isempty(units)
                append!(rels, [("unit__to_node", (u.name, next_node.name)) for u in units])
            end
        end
    end
    for (next_node, tail_demands) in new_demands
        next_demand = P.demand(node=next_node)
        new_demand = isnothing(next_demand) ? 0 : next_demand
        new_demand += sum(tail_demands)
        push!(opvs, ("node", next_node.name, "demand", unparse_db_value(new_demand), alternative))
    end
    to_prune_object_keys = [(class_name, x.name) for (class_name, objects) in to_remove for x in objects]
    if isempty(to_prune_object_keys)
        @info "No tails left to trim"
        return false
    end
    data_to_import = Dict(:relationships => rels, :object_parameter_values => opvs, :alternatives => [alternative])
    comment = "Network pruning: trim tails"
    _prune_and_import(prunned_db_url, to_prune_object_keys, data_to_import, comment)
    nodes_pruned = length(get(to_remove, :node, ()))
    connections_pruned = length(get(to_remove, :connection, ()))
    units_moved = length(rels)
    demands_moved = length(opvs)
    @info "Network tails trimmed successfully" nodes_pruned connections_pruned units_moved demands_moved
    true
end

function _prune_and_import(prunned_db_url, to_prune_object_keys, data_to_import, comment)
    to_prune_object_ids = [
        x["id"]
        for x in run_request(prunned_db_url, "query", ("ext_object_sq",))["ext_object_sq"]
        if (Symbol(x["class_name"]), Symbol(x["name"])) in to_prune_object_keys
    ]
    run_request(prunned_db_url, "call_method", ("cascade_remove_items",), Dict(:object => to_prune_object_ids))
    added, err_log = import_data(prunned_db_url, ""; data_to_import...)
    @info "Added $(added) items"
    for err in err_log
        @info "Import error: " err
    end
    run_request(prunned_db_url, "call_method", ("commit_session", comment))
end

function write_node__new_nodes(I, node__new_nodes, path)
    io = open(path, "w")
    print(io, "node,V,total_gens,total_generation,total_demand,")
    print(io, "node1,V1,ptdf1,node2,V2,ptdf2,node3,V3,ptdf3,node4,V4,ptdf4,node5,V5,ptdf5\n")
    for (n, new_nodes) in node__new_nodes
        total_generators = 0
        total_generation = 0
        total_demand = 0
        for u in I.unit__to_node(node=n)
            total_generators = total_generators + 1
            total_generation = total_generation + I.unit_capacity(unit=u, node=n)
        end
        if I.demand(node=n) == nothing
            total_demand = 0
        else
            total_demand = I.demand(node=n)
        end
        print(
            io, string(n), ",", I.voltage(node=n), ",", total_generators, ",", total_generation, ",", total_demand
        )
        if size(new_nodes, 1) == 1
            n2, _ptdf = new_nodes[1]
            print(io, ",", string(n2), ",", I.voltage(node=n2), ",", 1)
        else
            for (n2, ptdf) in new_nodes
                print(io, ",", string(n2), ",", I.voltage(node=n2), ",", string(abs(ptdf)))
            end
        end
        print(io, "\n")
    end
    close(io)
end

"""
    islands(I, comm)

Determines the number of islands in a commodity network - used for diagnostic purposes
"""
function islands(I, comm)
    visited_d = Dict{Object,Bool}()
    island_node = Dict{Int64,Array}()
    island = 0
    for n in I.node__commodity(commodity=comm)
        visited_d[n] = false
    end
    for n in I.node__commodity(commodity=comm)
        if !visited_d[n]
            island = island + 1
            island_node[island] = Object[]
            visit(I, n, island, visited_d, island_node)
        end
    end
    return island, island_node
end

"""
    visit()

Function called recursively to visit nodes in the network to determine number of islands
"""
function visit(I, n, island, visited_d, island_node)
    visited_d[n] = true
    push!(island_node[island], n)
    for conn in I.connection__from_node(node=n)
        for n2 in I.connection__to_node(connection=conn)
            if !visited_d[n2]
                visit(I, n2, island, visited_d, island_node)
            end
        end
    end
    for conn in I.connection__to_node(node=n)
        for n2 in I.connection__from_node(connection=conn)
            if !visited_d[n2]
                visit(I, n2, island, visited_d, island_node)
            end
        end
    end
end

"""
    check_x()

Check for low reactance values
"""
function check_x(I)
    @info "Checking reactances"
    for conn in I.connection()
        if I.conn_reactance(connection=conn) < 0.0001
            @info "Low reactance may cause problems for line " conn
        end
    end
end

"""
    calculate_ptdfs()

Returns a dict indexed on tuples of (connection, node) containing the ptdfs of the system currently in memory.
"""
function calculate_ptdfs(I, comm)
    ps_busses = Bus[]
    ps_lines = Line[]
    node_ps_bus = Dict{Object,Bus}()
    i = 1
    for n in I.node__commodity(commodity=comm)
        if I.node_opf_type(node=n) == :node_opf_type_reference
            bustype = BusTypes.REF
        else
            bustype = BusTypes.PV
        end
        ps_bus = Bus(
            number = i,
            name = string(n),
            bustype = bustype,
            angle = 0.0,
            # voltage = 0.0,
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
    for conn in I.connection()
        for n_from in I.connection__from_node(connection=conn)
            for n_to in I.connection__to_node(connection=conn)
                if comm in vcat(I.node__commodity(node=n_from), I.node__commodity(node=n_to))
                    ps_arc = Arc(node_ps_bus[n_from], node_ps_bus[n_to])
                    new_line = Line(;
                        name = string(conn),
                        available = true,
                        active_power_flow = 0.0,
                        reactive_power_flow = 0.0,
                        arc = ps_arc,
                        r = I.connection_resistance(connection=conn),
                        x = max(I.connection_reactance(connection=conn), 0.00001),
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
    ps_ptdf = PowerSystems.PTDF(ps_lines, ps_busses)
    ptdf = Dict{Tuple{Object,Object},Float64}()
    for n in I.node__commodity(commodity=comm)
        for conn in I.connection()
            ptdf[conn, n] = ps_ptdf[string(conn), node_ps_bus[n].number]
        end
    end
    # buildlodf needs to be updated to account for cases
    # lodfs = PowerSystems.buildlodf(ps_lines,ps_busses)
    return ptdf
end

"""
    calculate_lodfs(I, ptdf_b_n)

Returns lodfs for the system specified by ptdf_b_n ,b_con__b_mon as a dict of tuples: contingent_line, monitored_line

"""
# This function takes a long time.
# PowerSystems has a function that does it faster using linear algebra but doesn't handle the case of tails like
# I would like.
function calculate_lodfs(I, ptdf_conn_n)
    lodf_con_mon = Dict{Tuple{Object,Object},Float64}()
    con__mon = Tuple{Object,Object}[]
    considered_contingencies = 0
    skipped = 0
    tolerance = 0
    for conn_con in I.connection()
        if I.connection_contingency(connection=conn_con) == 1
            for n_from in I.connection__from_node(connection=conn_con)
                for n_to in I.connection__to_node(connection=conn_con)
                    denom = 1 - (ptdf_conn_n[(conn_con, n_from)] - ptdf_conn_n[(conn_con, n_to)])
                    if abs(denom) < 0.001
                        denom = -1
                    end
                    for conn_mon in I.connection()
                        if I.connection_monitored(connection=conn_mon) == 1 && conn_con != conn_mon
                            if denom == -1
                                lodf_trial = (ptdf_conn_n[(conn_mon, n_from)] - ptdf_conn_n[(conn_mon, n_to)]) / denom
                            else
                                lodf_trial = -ptdf_conn_n[(conn_mon, n_from)]
                            end
                            c = first(indices(commodity_lodf_tolerance))
                            tolerance = I.commodity_lodf_tolerance(commodity=c)
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
    return lodf_con_mon, con__mon
end

function write_ptdfs(I, ptdfs, net_inj_nodes)
    io = open("ptdfs.csv", "w")
    print(io, "connection,")
    for n in net_inj_nodes
        print(io, string(n), ",")
    end
    print(io, "\n")
    for conn in I.connection()
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
