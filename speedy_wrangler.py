from config import OUTPUT_DIR, SOURCE_DIR, PROJECT_ROOT
import importlib
import sys
from database.db_api import DuckdbAPI
from database.arn_api import ArnAPI
importlib.reload(sys.modules['database.arn_api'])


source_db_path = OUTPUT_DIR / 'ARN.duckdb'
db_source = ArnAPI(source_db_path)
# tissues col in nodes for tissues, but nested struct,list of dictionaries,
# dictionaries values: value(str), db(str), url(str, but empty), searcheable(bool)
proteins = db_source.get_all_protein_name()
edges = db_source.get_edges_by_layer(2)


protein_edges = [edge for edge in edges
                 if edge['source'] in proteins and edge['target'] in proteins]

sources = {item['source'] for item in protein_edges}
targets = {item['target'] for item in protein_edges}
all_nodes = sources | targets

sql_seed = PROJECT_ROOT / 'database' / 'duckdb_seed.sql'
output_db = OUTPUT_DIR / 'combined_db.duckdb'
db_end = DuckdbAPI(sql_seed, output_db)

for node in list(all_nodes):
    try:
        node_dict = db_source.get_node_by_id(node)
        node_dict['source_database'] = 'ARN'
        db_end.inser_node(node_dict)
    except Exception as e:
        print(f"failed to add {node}: {e}")

for edge in protein_edges:
    try:
        interactor_a_dict = {}
        interactor_b_dict = {}
        edge_dict = {}

        source_res = db_end.get_node_by_id(edge['source'])
        target_res = db_end.get_node_by_id(edge['target'])

        interactor_a_dict['id'] = source_res['id']
        interactor_b_dict['id'] = target_res['id']
        interactor_a_dict['name'] = source_res['name']
        interactor_b_dict['name'] = source_res['name']

        edge_dict['layer'] = edge['layer']
        edge_dict['interaction_types'] = ''
        edge_dict['source_db'] = 'ARN'

        db_end.insert_edge(interactor_a_dict, interactor_b_dict, edge_dict)
    except Exception as e:
        print(f"failed to add edge: {e}")
db_end.close()
db_source.close()
