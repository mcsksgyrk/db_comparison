from config import OUTPUT_DIR, SOURCE_DIR, PROJECT_ROOT
import sqlite3
from database.db_api import DuckdbAPI
import pandas as pd

#### ferrdb ####
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
# have to check later why these values exist
df_ferr.loc[df_ferr.tax_id != 9606, 'tax_id'] = 9606
df_ferr_edges['source_db'] = 'ferr'
conn.close()
pr_names_set = set(df_ferr.name)
protein_edges = df_ferr_edges[
    df_ferr_edges['interactor_a_node_name'].isin(pr_names_set) &
    df_ferr_edges['interactor_b_node_name'].isin(pr_names_set)
]
all_node_dict = df_ferr.reset_index(drop=True).to_dict('records')
##################

#### ARN ####
sql_seed = PROJECT_ROOT / 'database' / 'duckdb_seed.sql'

output_db = OUTPUT_DIR / 'combined_db.duckdb'
db_arn_data = DuckdbAPI(sql_seed, output_db, False)

query = """
    SELECT * FROM node
    """
df_arn = db_arn_data.db.execute(query).fetchdf()
db_arn_data.close()

max_id = df_arn['id'].max()
##################
common_nodes = df_arn[
    df_arn['name'].isin(df_ferr['name'])
]

mask = df_arn['name'].isin(df_ferr['name'])

df_arn.loc[mask, 'source_database'] += '|ferr'

ferr_only = df_ferr[
    ~df_ferr['name'].isin(df_arn['name'])
].copy()
ferr_only['id'] = range(max_id + 1, max_id + 1 + len(ferr_only))

df_combined = pd.concat([df_arn, ferr_only], ignore_index=True)
########## edges #################


def edge_exists(row, existing_edges):
    return ((existing_edges['interactor_a_node_id'] == row['interactor_a_node_id']) &
        (existing_edges['interactor_b_node_id'] == row['interactor_b_node_id'])).any()


name_to_id = df_combined.set_index('name')['id'].to_dict()
protein_edges_extended = protein_edges.copy()
protein_edges_extended['interactor_a_node_id'] = protein_edges_extended['interactor_a_node_name'].map(name_to_id)
protein_edges_extended['interactor_b_node_id'] = protein_edges_extended['interactor_b_node_name'].map(name_to_id)

query = """
    SELECT * FROM edge
    """
db_arn_data = DuckdbAPI(sql_seed, output_db, False)
df_arn_edge = db_arn_data.db.execute(query).fetchdf()
db_arn_data.close()
new_edges = protein_edges_extended[
    ~protein_edges_extended.apply(lambda row: edge_exists(row, df_arn_edge), axis=1)
]
max_edge_id = df_arn_edge['id'].max()
new_edges['id'] = range(max_edge_id + 1, max_edge_id + 1 + len(new_edges))

df_final_edges = pd.concat([df_arn_edge, new_edges], ignore_index=True)

db_end = DuckdbAPI(sql_seed, OUTPUT_DIR / 'test2.duckdb')
db_end.db.execute("INSERT INTO node SELECT * FROM df_combined")
db_end.db.execute("INSERT INTO edge SELECT * FROM df_final_edges")
db_end.close()
