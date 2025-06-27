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
proteins = db_source.get_all_protein_name()
edges = db_source.get_edges_by_layer(2)
edges

sql_seed = PROJECT_ROOT / 'database' / 'duckdb_seed.sql'
output_db = OUTPUT_DIR / 'combined_db.duckdb'
db_end = DuckdbAPI(sql_seed)

protein_edges = [edge for edge in edges
                 if edge['source'] in proteins and edge['target'] in proteins]
protein_edges
db_source.close()
