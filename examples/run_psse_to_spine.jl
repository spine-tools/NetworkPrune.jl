using NetworkPrune

psse_path = "WP2019.raw"
db_url = "sqlite:///powersystem.sqlite"
psse_to_spine(psse_path, db_url)