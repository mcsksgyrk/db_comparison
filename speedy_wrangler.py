import duckdb
from config import OUTPUT_DIR, SOURCE_DIR, PROJECT_ROOT
import os
import json
from typing import List
from database.db_api import DuckdbAPI


source_db = OUTPUT_DIR / 'ARN.duckdb'
sql_seed = PROJECT_ROOT / 'database' / 'duckdb_seed.sql'

db_end = DuckdbAPI(sql_seed)

db_source = duckdb.connect(source_db)
# tissues col in nodes for tissues, but nested struct,list of dictionaries,
# dictionaries values: value(str), db(str), url(str, but empty), searcheable(bool)

proteins = db_source.execute("""
    SELECT DISTINCT name
    FROM nodes
    WHERE EXISTS (
        SELECT 1 FROM unnest(moleculeType) AS mt(elem)
        WHERE elem.value LIKE '%protein%'
    )
    AND taxon.id = '9606';
""").fetchall()

protein_names = {row[0] for row in proteins}
test_name = list(protein_names)[0]
taxon_test = db_source.execute(f"SELECT taxon FROM nodes WHERE name='{test_name}'").fetchall()
taxon_test

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
