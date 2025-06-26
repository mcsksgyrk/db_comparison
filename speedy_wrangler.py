from config import OUTPUT_DIR, SOURCE_DIR, PROJECT_ROOT
import importlib
import sys
from database.db_api import DuckdbAPI
from database.arn_api import ArnAPI
importlib.reload(sys.modules['database.arn_api'])


source_db = OUTPUT_DIR / 'ARN.duckdb'
db_source = ArnAPI(source_db)
# tissues col in nodes for tissues, but nested struct,list of dictionaries,
# dictionaries values: value(str), db(str), url(str, but empty), searcheable(bool)
# gets the protein nodes from ARN db
proteins = db_source.db.execute("""
    SELECT DISTINCT name
    FROM nodes
    WHERE EXISTS (
        SELECT 1 FROM unnest(moleculeType) AS mt(elem)
        WHERE elem.value LIKE '%protein%'
    )
    AND taxon.id = '9606';
""").fetchall()
protein_names = list({row[0] for row in proteins})
bro = []
for protein in protein_names:
    bro.append(db_source.get_node_by_id(protein))


sql_seed = PROJECT_ROOT / 'database' / 'duckdb_seed.sql'
db_end = DuckdbAPI(sql_seed)
for node in bro:
    db_end.inser_node(node)
db_end.inser_node(bro[0])

db_end.db.execute("ATTACH './outputs/nodes_inserted.duckdb' AS file_db")
db_end.db.execute("COPY FROM DATABASE memory TO file_db")
db_end.db.execute("DETACH file_db")

edges = db_source.execute("""
    SELECT
        source,
        target,
        layer[1].value as layer,
        isDirected,
        isDirect,
        interactionType,
        sourceDatabases
    FROM edges
    WHERE layer[1].value < 3;
""").fetchall()

protein_edges = [edge for edge in edges
                 if edge[0] in protein_names and edge[1] in protein_names]

# gets all core autophagy proteins
aut_core_proteins = db_source.execute("""
    SELECT name, tissues FROM nodes
    WHERE len(autophagyPhase) > 0;
""").fetchall()

db_source.close()

#sql_seed = PROJECT_ROOT / 'database' / 'duckdb_seed.sql'
#db_end = DuckdbAPI(sql_seed)
