import duckdb
from pathlib import Path
from typing import List, Any, Dict


class ArnAPI:
    def __init__(self, db_path: Path):
        try:
            self.db = duckdb.connect(str(db_path))
        except Exception as e:
            print(f"Couldn't open databse: {e}")
            raise

    def _execute_as_dict(self, query, params=None):
        return self.db.execute(query, params or []).fetchdf().to_dict('records')

    def get_node_by_id(self, id: str, all=False) -> Dict:
        if all:
            query = """
                    SELECT * FROM nodes WHERE name = ?
                    """
        else:
            query = """
                    SELECT
                          displayedName,
                          name,
                          list_transform(topologicalFeatures, x -> x.value) as function
                    FROM nodes
                    WHERE name = ?
                    """
        res = self._execute_as_dict(query, ([id]))
        return res[0] if res else None
