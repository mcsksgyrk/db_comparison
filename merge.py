import duckdb
import sqlite3
import pandas as pd
from config import SOURCE_DIR


def migrate_sqlite_to_duckdb(sqlite_path, duckdb_path):
    """
    Migrate SQLite database to DuckDB with new schema
    """

    sqlite_conn = sqlite3.connect(sqlite_path)
    duck_conn = duckdb.connect(duckdb_path)

    try:
        sql_seed = open(str(SOURCE_DIR / "databse/sql_seed.sql")).read()
        duck_conn.execute(sql_seed)

        print("Creating DuckDB schema...")
        print("Migrating nodes...")
        nodes_df = pd.read_sql_query("SELECT * FROM node", sqlite_conn)

        # Add pathway classification
        nodes_df['source_database'] = 'ferroptosis_db'
        nodes_df.loc[nodes_df['source'].str.contains('autophagy', case=False, na=False), 'source_database'] = 'autophagy'
        nodes_df.loc[nodes_df['source'].str.contains('ferroptosis', case=False, na=False), 'source_database'] = 'ferroptosis'

        # Insert into DuckDB
        duck_conn.register('nodes_temp', nodes_df)
        duck_conn.execute("INSERT INTO node SELECT * FROM nodes_temp")

        print("Migrating node identifiers...")
        try:
            identifiers_df = pd.read_sql_query("SELECT * FROM node_identifier", sqlite_conn)
            duck_conn.register('identifiers_temp', identifiers_df)
            duck_conn.execute("INSERT INTO node_identifier SELECT * FROM identifiers_temp")
        except pd.errors.DatabaseError:
            print("No node_identifier table in source database")

        print("Migrating edges...")
        edges_df = pd.read_sql_query("SELECT * FROM edge", sqlite_conn)

        # Add new columns
        edges_df['pathway_type'] = 'unknown'
        edges_df['confidence_score'] = None
        edges_df['is_directed'] = False
        edges_df['publications'] = None

        # Classify edges
        edges_df.loc[edges_df['source_db'].str.contains('autophagy', case=False, na=False), 'pathway_type'] = 'autophagy'
        edges_df.loc[edges_df['source_db'].str.contains('ferroptosis', case=False, na=False), 'pathway_type'] = 'ferroptosis'

        duck_conn.register('edges_temp', edges_df)
        duck_conn.execute("INSERT INTO edge SELECT * FROM edges_temp")

        print("Creating indexes...")
        duck_conn.execute("CREATE INDEX idx_node_name ON node(name)")
        duck_conn.execute("CREATE INDEX idx_node_source ON node(source_database)")
        duck_conn.execute("CREATE INDEX idx_edge_pathway ON edge(pathway_type)")
        duck_conn.execute("CREATE INDEX idx_edge_interactors ON edge(interactor_a_node_id, interactor_b_node_id)")

        print("Migration verification...")

        # Verification queries
        result = duck_conn.execute("""
            SELECT 'Nodes' as table_name, COUNT(*) as count FROM node
            UNION ALL
            SELECT 'Edges', COUNT(*) FROM edge
            UNION ALL
            SELECT 'Node IDs', COUNT(*) FROM node_identifier
        """).fetchall()

        print("\nMigration Results:")
        for row in result:
            print(f"{row[0]}: {row[1]}")

        # Pathway distribution
        pathway_dist = duck_conn.execute("""
            SELECT source_database, COUNT(*) as count
            FROM node
            GROUP BY source_database
            ORDER BY count DESC
        """).fetchall()

        print("\nNode Pathway Distribution:")
        for row in pathway_dist:
            print(f"{row[0]}: {row[1]}")

    finally:
        sqlite_conn.close()
        duck_conn.close()

    print(f"\nMigration completed! DuckDB database: {duckdb_path}")


if __name__ == "__main__":
    migrate_sqlite_to_duckdb(
        sqlite_path="old_database.db",
        duckdb_path="new_analysis.duckdb"
    )
