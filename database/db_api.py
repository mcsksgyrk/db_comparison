import duckdb
from pathlib import Path
from typing import Optional, List, Any, Dict


class DuckdbAPI:
    def __init__(self, sql_seed: Path, db_path: Optional[Path] = None):
        self.sql_seed = open(sql_seed).read()

        if db_path:
            self.db = duckdb.connect(str(db_path))
        else:
            self.db = duckdb.connect()
        self.create_schema()

    def create_schema(self) -> None:
        if self.sql_seed:
            self.db.execute(self.sql_seed)

    def read_existing_db(self, existing_db_path: Path) -> None:
        try:
            self.db.execute(f"ATTACH '{existing_db_path}' AS source_db")
            self.db.execute("ATTACH INTO edges SELECT * FROM source_db.edges")
            self.db.execute("ATTACH INTO nodes SELECT * FROM source_db.nodes")
            try:
                self.db.execute("ATTACH INTO tissue SELECT * FROM source_db.tissue")
            except Exception:
                pass
            try:
                self.db.execute("ATTACH INTO node_tissue SELECT * FROM source_db.node_tissue")
            except Exception:
                pass
        except Exception as e:
            print(f"failed to import db {e}")
            raise

    def inser_node(self, node_dict: Dict) -> None:
        query = """
                INSERT INTO node
                (name, display_name, tax_id, type, pathways, source, function, source_database)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """
        self.db.execute(query, (
                node_dict['name'],
                node_dict.get('display_name', ''),
                node_dict.get('tax_id', '9606'),
                node_dict.get('type', 'protein'),
                node_dict.get('pathways', []),
                node_dict.get('source', ''),
                node_dict.get('function', []),
                node_dict.get('source_database', '')
                ))
