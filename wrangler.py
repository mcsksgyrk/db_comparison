import duckdb
from config import OUTPUT_DIR, SOURCE_DIR
import os
import json
from typing import List


def query_executor(conn):
    try:
        conn.execute(f"""
            SELECT moleculeType FROM nodes WHERE moleculeType
        """)
    except Exception as e:
        print(f"Error: {e}")


def tissue_extractor(tissues_struct: List[dict]) -> List[str]:
    tissues = []
    for tissue in tissues_struct:
        tissues.append(tissue['value'])
    return tissues


source_db = OUTPUT_DIR / 'rawdata.duckdb'
conn = duckdb.connect(source_db)
query_executor(conn)
conn.fetchall()
conn.execute("DESCRIBE edges;").fetchall()


# select every edge where source and target is a protein,
# and the edge is of layer-0 -1 or -2.
conn.execute("""
    SELECT COUNT(*)
    FROM edges
    INNER JOIN nodes AS source_node ON edges.source = source_node.name
    INNER JOIN nodes AS target_node ON edges.target = target_node.name
    WHERE EXISTS (
        SELECT 1 FROM unnest(source_node.moleculeType) AS mt(elem)
        WHERE elem.value LIKE '%protein%'
    )
    AND EXISTS (
        SELECT 1 FROM unnest(target_node.moleculeType) AS mt(elem)
        WHERE elem.value LIKE '%protein%'
    )
    AND edges.layer[1].value < 3;
""")

# gets all core autophagy proteins
conn.execute("""
    SELECT name, tissues FROM nodes
    WHERE len(autophagyPhase) > 0;
""")
aut_coore = conn.fetchall()
test = aut_coore[0]
print(tissue_extractor(test[1]))
