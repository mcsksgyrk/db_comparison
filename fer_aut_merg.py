from config import OUTPUT_DIR, SOURCE_DIR, PROJECT_ROOT
import sqlite3
from database.db_api import DuckdbAPI
import pandas as pd
from pathlib import Path
from typing import Tuple, Duct
import logging

logging.basicConfig(lebel=logging.INFO)
logger = logging.getLogger(__name__)


class DatabaseMerger:
    def __init__(self, ferr_db_path: Path, arn_path: Path,
                 sql_seed_path: Path, output_path: Path):
        self.ferr_db_path = ferr_db_path
        self.sql_seed_path = sql_seed_path
        self.arn_path = arn_path
        self.output_path = output_path

    def load_ferroptsis_data(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        logger.info("Loading ferroptosis data")
        with sqlite3.connect(self.ferr_db_path) as conn:
            node_query = """
                SELECT name, display_name, tax_id, type
                FROM node
                WHERE type = 'protein' AND primary_id_type = 'uniprot_id'
                """
            df_nodes = pd.read_sql_query(node_query, conn)
            df_nodes['source_database'] = 'ferr'

            # have to check later why these values exist
            df_nodes.loc[df_nodes.tax_id != 9606, 'tax_id'] = 9606

            edges_query = """
                    SELECT interactor_a_node_name, interactor_b_node_name, layer
                    FROM edge
                    WHERE layer IN ('0', '1', '2', 0, 1, 2)
                    """
            df_edges = pd.read_sql_query(edges_query, conn)
            df_edges['source_db'] = 'ferr'

        pr_names_set = set(df_nodes.name)
        df_edges = df_edges[
            df_edges['interactor_a_node_name'].isin(pr_names_set) &
            df_edges['interactor_b_node_name'].isin(pr_names_set)
        ]
        logger.info(f"Loaded {len(df_nodes)} ferroptosis nodes, {len(df_edges)} edges")
        return df_nodes, df_edges
        # all_node_dict = df_nodes.reset_index(drop=True).to_dict('records')

    def load_arn_data(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        logger.info("Loading ARN data")
        db_arn = DuckdbAPI(self.sql_seed_path, self.arn_path, create_new=False)

        try:
            df_nodes = db_arn.db.execute("SELECT * FROM node").fetchdf()
            df_edges = db_arn.db.execute("SELECT * FROM edge").fetchdf()
        finally:
            db_arn.close()

        logger.info(f"Loaded {len(df_nodes)} autophagy nodes, {len(df_edges)} edges")
        return df_nodes, df_edges

    def merge_nodes(self, df_arn: pd.DataFrame, df_ferr: pd.DataFrame) -> pd.DataFrame:
        logger.info("Merging nodes")

        overlap_mask = df_arn['name'].isin(df_ferr['name'])
        df_arn.loc[overlap_mask, 'source_database'] += '|ferr'

        max_id = df_arn['id'].max()
        ferr_only = df_ferr[
            ~df_ferr['name'].isin(df_arn['name'])
        ].copy()
        ferr_only['id'] = range(max_id + 1, max_id + 1 + len(ferr_only))

        df_combined = pd.concat([df_arn, ferr_only], ignore_index=True)
        return df_combined

    def merge_edges(self, df_arn_edges: pd.DataFrame, df_ferr_edges: pd.DataFrame,
                    name_to_id_map: dict[str, id]) -> pd.DataFrame:
        logger.info("Merging edges")

        def edge_exists(row, existing_edges):
            return ((existing_edges['interactor_a_node_id'] == row['interactor_a_node_id']) &
                (existing_edges['interactor_b_node_id'] == row['interactor_b_node_id'])).any()

        protein_edges_extended = df_ferr_edges.copy()
        protein_edges_extended['interactor_a_node_id'] = protein_edges_extended['interactor_a_node_name'].map(name_to_id_map)
        protein_edges_extended['interactor_b_node_id'] = protein_edges_extended['interactor_b_node_name'].map(name_to_id_map)

        new_edges = protein_edges_extended[
            ~protein_edges_extended.apply(lambda row: edge_exists(row, df_arn_edges), axis=1)
        ]
        if len(new_edges) > 0:
            max_edge_id = df_arn_edges['id'].max()
            new_edges = new_edges.copy()
            new_edges['id'] = range(max_edge_id + 1, max_edge_id + 1 + len(new_edges))

        df_combined_edges = pd.concat([df_arn_edges, new_edges], ignore_index=True)

        logger.info(f"Added {len(new_edges)} new edges")
        return df_combined_edges

    def save_combined_database(self, df_nodes: pd.DataFrame, df_edges: pd.DataFrame):
        logger.info(f"Saving combined database to {self.output_path}")
        db_output = DuckdbAPI(self.sql_seed_path, self.output_path, create_new=True)

        try:
            db_output.db.execute("INSERT INTO node SELECT * FROM df_nodes")
            db_output.db.execute("INSERT INTO edge SELECT * FROM df_edges")
            logger.info("Database saved successfully")
        finally:
            db_output.close()

    def merge_databases(self):
        logger.info("Merging")

        df_ferr_nodes, df_ferr_edges = self.load_ferroptsis_data()
        df_arn_nodes, df_arn_edges = self.load_arn_data()

        df_combined_nodes = self.merge_nodes(df_arn_nodes, df_ferr_nodes)
        name_to_id_map = df_combined_nodes.set_index('name')['id'].to_dict()
        df_combined_edges = self.merge_edges(df_arn_edges, df_ferr_edges, name_to_id_map)

        self.save_combined_database(df_combined_nodes, df_combined_edges)

        logger.info("Database saved")


def main():
    ferr_db_path = PROJECT_ROOT / 'test_omnipath.db'
    sql_seed_path = PROJECT_ROOT / 'database' / 'duckdb_seed.sql'
    output_path = OUTPUT_DIR / 'test2.duckdb'

    merger = DatabaseMerger(ferr_db_path, sql_seed_path, output_path)
    df_nodes, df_edges = merger.merge_databases()


if __name__ == "__main__":
    main()
