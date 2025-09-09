from config import OUTPUT_DIR, PROJECT_ROOT, SOURCE_DIR
import duckdb
from database.sqlite_db_api3 import PsimiSQL

# Create temporary DuckDB connection
temp_conn = duckdb.connect(':memory:')

# Load JSON files
print("Loading ARN JSON files...")
nodes_json = SOURCE_DIR / 'nodes_RC2_2023_06_06.json'
edges_json = SOURCE_DIR / 'edges_RC2_2023_06_06.json'

temp_conn.execute(f"""
    CREATE TABLE raw_nodes AS
    SELECT * FROM read_json_auto('{nodes_json}')
""")

temp_conn.execute(f"""
    CREATE TABLE raw_edges AS
    SELECT * FROM read_json_auto('{edges_json}')
""")

# Extract nodes with simple filtering
print("Extracting nodes...")
nodes_query = """
    SELECT
        name,
        displayedName as display_name,
        taxon,
        pathways,
        topologicalFeatures,
        moleculeType
    FROM raw_nodes
    WHERE taxon.id = '9606'
"""

df_nodes_raw = temp_conn.execute(nodes_query).fetchdf()

# Extract edges
print("Extracting edges...")
edges_query = """
    SELECT
        source,
        target,
        layer,
        isDirected,
        interactionType
    FROM raw_edges
"""

df_edges_raw = temp_conn.execute(edges_query).fetchdf()
temp_conn.close()

# Process nodes in Python to avoid DuckDB JSON issues
print("Processing nodes...")
processed_nodes = []
for idx, row in df_nodes_raw.iterrows():
    # Check if protein
    is_protein = False
    try:
        mol_types = row['moleculeType']
        if isinstance(mol_types, list):
            for mt in mol_types:
                if isinstance(mt, dict) and 'value' in mt:
                    if 'protein' in mt['value'].lower():
                        is_protein = True
                        break
    except:
        continue

    if not is_protein:
        continue

    # Extract pathways
    pathways_str = ""
    try:
        pathways = row['pathways']
        if isinstance(pathways, list):
            pathway_values = []
            for p in pathways:
                if isinstance(p, dict) and 'value' in p:
                    pathway_values.append(str(p['value']))
            pathways_str = "|".join(pathway_values)
    except:
        pathways_str = ""

    # Extract function
    function_str = ""
    try:
        features = row['topologicalFeatures']
        if isinstance(features, list):
            feature_values = []
            for f in features:
                if isinstance(f, dict) and 'value' in f:
                    feature_values.append(str(f['value']))
            function_str = "|".join(feature_values)
    except:
        function_str = ""

    node_dict = {
        'name': row['name'],
        'primary_id_type': 'uniprot_id',
        'display_name': row['display_name'],
        'tax_id': 9606,
        'type': 'protein',
        'pathways': pathways_str,
        'role_in_ferroptosis': '',
        'function': function_str,
        'source_db': 'ARN'
    }
    processed_nodes.append(node_dict)

print(f"Debug info:")
print(f"Total nodes processed: {total_nodes}")
print(f"Proteins found: {protein_found}")
print(f"Sample molecule types: {mol_type_debug[:10]}")
print(f"Final processed nodes: {len(processed_nodes)}")

# Process edges
print("Processing edges...")
processed_edges = []
for idx, row in df_edges_raw.iterrows():
    # Get layer value
    layer_val = "0"
    try:
        layer = row['layer']
        if isinstance(layer, list) and len(layer) > 0:
            if isinstance(layer[0], dict) and 'value' in layer[0]:
                layer_val = str(layer[0]['value'])
                if int(layer_val) >= 3:
                    continue
    except:
        continue

    # Get interaction type
    interaction_type = "unknown"
    try:
        int_type = row['interactionType']
        if isinstance(int_type, list) and len(int_type) > 0:
            if isinstance(int_type[0], dict) and 'value' in int_type[0]:
                interaction_type = str(int_type[0]['value'])
    except:
        interaction_type = "unknown"

    # Build interaction types string
    is_directed = str(row.get('isDirected', False)).lower()
    interaction_types = f"is_directed:{is_directed}|{interaction_type}"

    edge_dict = {
        'interactor_a_node_name': row['source'],
        'interactor_b_node_name': row['target'],
        'layer': layer_val,
        'interaction_types': interaction_types,
        'effect_on_ferroptosis': '',
        'source_db': 'ARN'
    }
    processed_edges.append(edge_dict)

# Create SQLite database
print("Creating SQLite database...")
sql_seed_path = PROJECT_ROOT / 'database' / 'network_db_seed3.sql'
output_path = OUTPUT_DIR / 'arn_converted.db'
db_api = PsimiSQL(sql_seed_path)

# Insert nodes
print(f"Inserting {len(processed_nodes)} nodes...")
for node_dict in processed_nodes:
    db_api.insert_node(node_dict)
    node_id = node_dict['id']
    db_api.insert_node_identifier(node_id, 'uniprot_id', node_dict['name'], True)

# Insert edges
print(f"Processing {len(processed_edges)} edges...")
valid_edges = 0
for edge_dict in processed_edges:
    source_node = db_api.get_node_by_any_identifier(edge_dict['interactor_a_node_name'])
    target_node = db_api.get_node_by_any_identifier(edge_dict['interactor_b_node_name'])

    if source_node and target_node:
        db_api.insert_edge(source_node, target_node, edge_dict)
        valid_edges += 1

# Save database
db_api.save_db_to_file(str(output_path))
print(f"ARN database converted to SQLite: {output_path}")
print(f"Final: {len(processed_nodes)} nodes, {valid_edges} valid edges")
