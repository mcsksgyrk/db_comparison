import duckdb
from pathlib import Path
import pandas as pd
from database.sqlite_db_api3 import PsimiSQL
from config import PROJECT_ROOT, OUTPUT_DIR

source_dict = Path("actual_db_arn2")
nodes_json = source_dict / "nodes_RC2_2023_06_06.json"
edges_json = source_dict / "edges_RC2_2023_06_06.json"

conn = duckdb.connect()

nodes_query = f"""
    SELECT * FROM read_json('{nodes_json}')
    WHERE taxon.id = 9606
    AND moleculeType[1].value LIKE '%protein%'
"""

pr_nodes = conn.execute(nodes_query).fetchdf()
pr_hash = set(pr_nodes.name)

edges_query = f"""
    SELECT * FROM read_json('{edges_json}')
    WHERE isDirect = true
    AND layer[1].value < 3
"""

edges = conn.execute(edges_query).fetchdf()
pr_edges = edges[
    edges['source'].isin(pr_hash) &
    edges['target'].isin(pr_hash)
]

nodes_clean = pd.DataFrame({
   'name': pr_nodes['name'],
   'primary_id_type': 'uniprot_id',
   'display_name': pr_nodes['displayedName'],
   'tax_id': pr_nodes['taxon'].apply(lambda x: x['id'] if x else None),
   'type': 'protein',
   'pathways': pr_nodes['pathways'],
   'role_in_ferroptosis': '',
   'function': '',
   'source_db': 'ARN'
})

edges_clean = pd.DataFrame({
    'interactor_a_node_name': pr_edges['source'],
    'interactor_b_node_name': pr_edges['target'],
    'layer': pr_edges['layer'].apply(lambda x: x[0]['value'] if len(x) > 0 else None),
    'interaction_types': pr_edges.apply(
        lambda row: str(row['isDirect']) + "|" + str(row['isDirected']) + "|" +
        (row['interactionType'][0]['value'] if len(row['interactionType']) > 0 else ''),
        axis=1
    ),
    'effect_on_ferroptosis': '',
    'source_db': 'ARN'
})

sql_seed = PROJECT_ROOT / "database" / "network_db_seed3.sql"
db_api = PsimiSQL(sql_seed)
for _, node_row in nodes_clean.iterrows():
    node_dict = node_row.to_dict()
    db_api.insert_node(node_dict)

for _, edge_row in edges_clean.iterrows():
    edge_dict = edge_row.to_dict()
    interactor_a = db_api.get_node_by_name(edge_dict['interactor_a_node_name'])
    interactor_b = db_api.get_node_by_name(edge_dict['interactor_b_node_name'])

    if interactor_a and interactor_b:
        db_api.insert_edge(interactor_a, interactor_b, edge_dict)
    else:
        print(f"Missing nodes for edge: {edge_dict['interactor_a_node_name']} -> {edge_dict['interactor_b_node_name']}")

db_api.save_db_to_file(str(OUTPUT_DIR / 'arn.db'))
