from config import OUTPUT_DIR, PROJECT_ROOT
import duckdb
from database.sqlite_db_api3 import PsimiSQL

# Connect to ARN DuckDB and use SQL to flatten everything
arn_db_path = OUTPUT_DIR / 'ARN.duckdb'
arn_conn = duckdb.connect(str(arn_db_path))

# Extract nodes with all flattening done in SQL
print("Extracting and flattening ARN nodes...")
nodes_query = """
    SELECT
        name,
        displayedName as display_name,
        CAST(taxon.id AS INTEGER) as tax_id,
        'protein' as type,
        'uniprot_id' as primary_id_type,
        COALESCE(
            list_transform(pathways, x -> COALESCE(x.value, CAST(x AS VARCHAR))),
            []
        ) as pathways_array,
        COALESCE(
            list_transform(topologicalFeatures, x -> COALESCE(x.value, CAST(x AS VARCHAR))),
            []
        ) as function_array,
        '' as role_in_ferroptosis,
        'ARN' as source_db
    FROM nodes
    WHERE EXISTS (
        SELECT 1 FROM unnest(moleculeType) AS mt(elem)
        WHERE elem.value LIKE '%protein%'
    )
    AND taxon.id = '9606'
"""

# Convert arrays to pipe-separated strings in SQL
nodes_with_strings_query = f"""
    WITH flattened AS ({nodes_query})
    SELECT
        name,
        display_name,
        tax_id,
        type,
        primary_id_type,
        array_to_string(pathways_array, '|') as pathways,
        role_in_ferroptosis,
        array_to_string(function_array, '|') as function,
        source_db
    FROM flattened
"""

df_nodes = arn_conn.execute(nodes_with_strings_query).fetchdf()

# Extract edges with flattening done in SQL
print("Extracting and flattening ARN edges...")
edges_query = """
    SELECT
        source as interactor_a_node_name,
        target as interactor_b_node_name,
        CAST(layer[1].value AS VARCHAR) as layer,
        array_to_string([
            CONCAT('is_directed:', CAST(isDirected AS VARCHAR)),
            COALESCE(interactionType[1].value, 'unknown')
        ], '|') as interaction_types,
        '' as effect_on_ferroptosis,
        'ARN' as source_db
    FROM edges
    WHERE layer[1].value < 3
"""

df_edges = arn_conn.execute(edges_query).fetchdf()
arn_conn.close()

# Create SQLite database
print("Creating SQLite database...")
sql_seed_path = PROJECT_ROOT / 'database' / 'network_db_seed3.sql'
output_path = OUTPUT_DIR / 'arn_converted_from_duck.db'
db_api = PsimiSQL(sql_seed_path)

# Insert nodes
print(f"Inserting {len(df_nodes)} nodes...")
for idx, row in df_nodes.iterrows():
    node_dict = row.to_dict()
    db_api.insert_node(node_dict)

    node_id = node_dict['id']
    db_api.insert_node_identifier(node_id, 'uniprot_id', row['name'], True)

# Insert edges
print(f"Inserting {len(df_edges)} edges...")
valid_edges = 0
for idx, row in df_edges.iterrows():
    source_node = db_api.get_node_by_any_identifier(row['interactor_a_node_name'])
    target_node = db_api.get_node_by_any_identifier(row['interactor_b_node_name'])

    if source_node and target_node:
        edge_dict = {
            'layer': row['layer'],
            'interaction_types': row['interaction_types'],
            'effect_on_ferroptosis': row['effect_on_ferroptosis'],
            'source_db': row['source_db']
        }
        db_api.insert_edge(source_node, target_node, edge_dict)
        valid_edges += 1

# Save database
db_api.save_db_to_file(str(output_path))
print(f"ARN database converted to SQLite: {output_path}")
print(f"Nodes: {len(df_nodes)}, Valid edges: {valid_edges}")
