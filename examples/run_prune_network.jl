using NetworkPrune

db_url = "sqlite:///powersystem.sqlite"
prunned_db_url = "sqlite:///powersystem_prunned.sqlite"
prune_network(db_url)
