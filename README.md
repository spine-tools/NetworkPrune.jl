# NetworkPrune

Production cost modelling and load flow analysis are common and often separate tasks in power systems analysis. Often it is desirable to conduct both sets of analyses from the same network dataset. However, network models used in load flow and contingency analysis often have a far higher level of detail than is required by production cost and energy systems models. For example, network models for load flow analysis often model down to the level of the low voltage generator or distribution transformer. The motivation behind this repository is to provide a set of tools which can aid with this process and, from the same loadflow dataset, generate a model suitable for production cost modelling without loss of detail or accuracy.

This packages provides two main functions:
 - `prune_network`, which reduces a network model by voltage level.
 - `psse_to_spine`, which parses PSSE raw files using PowerSystems.jl and creates a minimal SpineOpt model.
  
## `prune_network`
The goal of this function is to reduce the dimensionality of the network model by removing lower voltage nodes and connections that are not necessary for production cost and energy system modelling studies. The algorithm calculates the power transfer distribution factors (PTDFs) of the network and uses tree traversal to find the nearest higher-voltage level node above the minimum voltage level specified. When only a single higher-level node is identified, generation and loads are simply moved to the nearest higher level node. When multiple paths exist to a higher voltage level, the network PTDFs are used to distribute the load and generation in the correct proportion. This approach will yield an identical DC power flow solution to the original network but will contain fewer nodes and connections - the trade-off is that the flows on the lower level network are unknown.

### Usage
To run the function on a Spine database located at `db_url` and write results on `prunned_db_url`, simply call `prune_network(db_url, prunned_db_url)`.

For example:

```julia
using NetworkPrune

db_url = "sqlite:///powersystem.sqlite"
prunned_db_url = "sqlite:///powersystem_prunned.sqlite"
prune_network(db_url, prunned_db_url)
```

Requirements for the input Spine database
 - The network is defined by nodes and connections. It is assumed each connection has a single `connection__from_node` and a single `connection__to_node relationship`
 - To be included in the network, each node must be related to a single `commodity` with `commodity_physics(commodity=c) in(:commodity_physics_ptdf, :commodity_physics_lodf)`
 - To define the minimum voltage level, a node group can be created and the parameter `minimum_voltage` specified which will set the minimum voltage level for all member nodes. 
 - A sample database can be found in the `examples` folder.

## `psse_to_spine`
This function takes as input a PSSE raw file and uses Powersystems.jl to parse it. It then creates an equivalent minimal SpineOpt system which can be used directly with network_prune.jl as described above. The function automatically creates node groups based on the zones and areas in the PSSE raw file. A minimum voltage level of 110 is currently pre-specified. This can easily be changed by updated the minimum_voltage parameter in the resulting Spine database that is created.

### Usage
To convert the PSSE raw file at `psse_path` to a Spine Database at `db_url` one simply calls `psse_to_spine(psse_path, db_url)`

For example:

```julia
using NetworkPrune

db_url = "sqlite:///powersystem.sqlite"
psse_path = "WP2019.raw"
psse_to_spine(psse_path, db_url)
```
