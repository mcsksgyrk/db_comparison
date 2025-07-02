from config import OUTPUT_DIR, SOURCE_DIR, PROJECT_ROOT
import sqlite3
from database.db_api import DuckdbAPI
import pandas as pd


source_db_path = PROJECT_ROOT / 'test_omnipath.db'
conn = sqlite3.connect(source_db_path)
query = """
    SELECT name, display_name, tax_id, type
    FROM node
    WHERE type = 'protein' AND primary_id_type = 'uniprot_id'
    """
df_ferr = pd.read_sql_query(query, conn)
df_ferr['source_database'] = 'ferr'

query = """
        SELECT interactor_a_node_name, interactor_b_node_name, layer
        FROM edge
        WHERE layer IN ('0', '1', '2', 0, 1, 2)
        """
df_ferr_edges = pd.read_sql_query(query, conn)
df_ferr_edges['source_database'] = 'ferr'

conn.close()

pr_names_set = set(df_ferr.name)
protein_edges = df_ferr_edges[
    df_ferr_edges['interactor_a_node_name'].isin(pr_names_set) &
    df_ferr_edges['interactor_b_node_name'].isin(pr_names_set)
].to_dict('records')

all_node_dict = df_ferr.to_dict('records')

sql_seed = PROJECT_ROOT / 'database' / 'duckdb_seed.sql'
output_db = OUTPUT_DIR / 'testing.duckdb'
db_end = DuckdbAPI(sql_seed, output_db, False)
for node in all_node_dict:
    print(node['name'])
    db_end.insert_or_update_node(node)
db_end.close()
