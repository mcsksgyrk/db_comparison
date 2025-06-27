import duckdb
from pathlib import Path
from typing import List, Dict, Set, Optional


class ArnAPI:
    def __init__(self, db_path: Path):
        try:
            self.db = duckdb.connect(str(db_path))
        except Exception as e:
            print(f"Couldn't open databse: {e}")
            raise

    def _execute_as_dict(self, query, params=None):
        return self.db.execute(query, params or []).fetchdf().to_dict('records')

    def _convert_to_set(self, res) -> Set:
        return {row[0] for row in res}

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
        res = self._execute_as_dict(query, (id,))
        return res[0] if res else None

    def get_all_protein_name(self) -> Set:
        proteins = self.db.execute("""
                SELECT DISTINCT name
                FROM nodes
                WHERE EXISTS (
                    SELECT 1 FROM unnest(moleculeType) AS mt(elem)
                    WHERE elem.value LIKE '%protein%'
                )
                AND taxon.id = '9606';
        """).fetchall()
        return self._convert_to_set(proteins)

    def get_core_proteins(self) -> Set:
        aut_core_proteins = self.db.execute("""
            SELECT name, tissues FROM nodes
            WHERE len(autophagyPhase) > 0;
        """).fetchall()
        return self._convert_to_set(aut_core_proteins)

    def get_edges_by_layer(self, layer: Optional[int] = 2) -> List[Dict]:
        query = """
            SELECT
                source,
                target,
                layer[1].value as layer,
                isDirected,
                isDirect,
                interactionType,
                sourceDatabases
            FROM edges
            WHERE layer[1].value < ?;
            """
        res = self._execute_as_dict(query, (layer,))
        return res
