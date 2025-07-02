import duckdb
from pathlib import Path
from typing import Optional, List, Any, Dict, Set


class DuckdbAPI:
    def __init__(self,
                 sql_seed: Path,
                 db_path: Optional[Path] = None,
                 create_new: bool = True):
        self.sql_seed = open(sql_seed).read()

        if db_path:
            self.db = duckdb.connect(str(db_path))
        else:
            self.db = duckdb.connect()
        if create_new:
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

    def _execute_as_dict(self, query, params=None):
        return self.db.execute(query, params or []).fetchdf().to_dict('records')

    def _convert_to_set(self, res) -> Set:
        return {row[0] for row in res}

    def get_node_by_id(self, name: str) -> Dict:
        query = "SELECT * FROM node WHERE name = ?"
        res = self._execute_as_dict(query, (name,))
        return res[0] if res else None

    def inser_node(self, node_dict: Dict) -> None:
        res = self.get_node_by_id(node_dict['name'])
        if res:
            print(f"Node {node_dict['name']} already exists")
            return None
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

    def insert_edge(self, interactor_a_dict, interactor_b_dict, edge_dict):
        query = """
                INSERT INTO edge (
                interactor_a_node_id,
                interactor_b_node_id,
                interactor_a_node_name,
                interactor_b_node_name,
                layer,
                source_db,
                interaction_types
                )
                VALUES ( ?, ?, ?, ?, ?, ?, ?)
                """

        tup = (
            interactor_a_dict['id'],
            interactor_b_dict['id'],
            interactor_a_dict['name'],
            interactor_b_dict['name'],
            edge_dict['layer'],
            edge_dict['source_db'],
            edge_dict.get('interaction_types', [])
        )

        self.db.execute(query, tup)

    def update_existing_node(self, node_id: int, node_dict: Dict) -> None:
        existing_node = self.get_node_by_id(node_dict['name'])

        updates = {}
        for field, new_value in node_dict.items():
            if field in ['id', 'name', 'tax_id']:
                continue

            old_value = existing_node.get(field, '')
            if isinstance(new_value, list) and isinstance(old_value, list):
                merged = list(set(old_value + new_value))
                if merged != old_value:
                    updates[field] = merged
            elif isinstance(new_value, str) and isinstance(old_value, str):
                if old_value and new_value and new_value not in old_value:
                    merged = f"{old_value}|{new_value}"
                    updates[field] = merged
                elif not old_value and new_value:
                    updates[field] = new_value

            if not updates:
                return None
        try:

            set_clauses = [f"{field} = ?" for field in updates.keys()]
            values = list(updates.values()) + [node_id]
            query = f"""
                    UPDATE node
                    SET {', '.join(set_clauses)}
                    WHERE id = ?
                    """
            self.db.execute(query, values)
        except Exception as e:
            print(f"failed to update {existing_node} with {node_dict}: {e}")

    def insert_or_update_node(self, node_dict: Dict) -> int:
        existing_node = self.get_node_by_id(node_dict['name'])
        if existing_node:
            self.update_existing_node(existing_node['id'], node_dict)
            return existing_node['id']
        else:
            self.inser_node(node_dict)

    def close(self) -> None:
        self.db.close()
