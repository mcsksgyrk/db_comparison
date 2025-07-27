import pandas as pd
from database.db_api import DuckdbAPI
from config import *
import duckdb
from apicalls.api_oop import *
import pickle

source_db = OUTPUT_DIR / "test2.duckdb"
output_db = OUTPUT_DIR / "w_tissues.duckdb"
sql_seed = PROJECT_ROOT / "database" / "duckdb_seed.sql"

source_db = duckdb.connect(source_db)

query_nodes = """
              SELECT id, name FROM node
              """
nodes_df = source_db.execute(query_nodes).fetchdf()
source_db.close()
uniprot = UniProtClient()
gene_names, fails = uniprot.batch_convert_from_uniprot_id("Ensembl", nodes_df.name.to_list(), batch_size=10)
nodes_df['ensembl_id'] = nodes_df['name'].map(gene_names)

nodes_df.to_csv("break.csv")

bgee = BGEEClient()

target_localisations = {
    'dendritic_cell': 'CL:0000451',
    'macrophage': 'CL:0000235',
    'fibroblast': 'CL:0000057',
    't_cell': 'CL:0000084',
    'epithelial_cell': 'CL:0000066'
}

bgee = BGEEClient()
results = {}
failures = []
ids = nodes_df[~nodes_df.ensembl_id.isna()].ensembl_id.to_list()
for i, id in enumerate(ids):
    try:
        res = bgee.get_expression_anat_entity(id)
        results[id] = res
    except Exception as e:
        failures.append(id)
        print(f"{id} failed: {e}")

    if (i + 1) % 100 == 0:
        with open(f'./backup/backup_{i+1}.pkl', 'wb') as f:
            pickle.dump(results, f)
        print(f"Saved at {i+1}")

flattened = {k: next(iter(v.values())) for k, v in results.items() if v}
nodes_df['uberon'] = nodes_df.ensembl_id.map(flattened)
nodes_df.to_csv("nodes_w_uberon_annot.csv")

uberons = set()
for uberon_list in nodes_df.uberon.dropna():
    uberons.update(uberon_list)

test = nodes_df.uberon.to_list()
