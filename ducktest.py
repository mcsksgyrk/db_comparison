import duckdb
import os
import json

files = ['./actual_db_arn2/nodes_RC2_2023_06_06.json', './actual_db_arn2/edges_RC2_2023_06_06.json', './actual_db_arn2/attributes_RC2_2023_06_06.json']

for file in files:
    if os.path.exists(file):
        size_mb = os.path.getsize(file) / (1024**2)
        print(f"{file}: {size_mb:.1f} MB")
    else:
        print(f"{file}: NOT FOUND")

conn = duckdb.connect('test.db')
try:
    conn.execute("""
        CREATE TABLE attributes AS
        SELECT * FROM read_json_auto('./actual_db_arn2/attributes_RC2_2023_06_06.json')
    """)

    count = conn.execute("SELECT COUNT(*) FROM attributes").fetchone()[0]
    print(f"Attributes loaded: {count} records")

    columns = conn.execute("DESCRIBE attributes").fetchall()
    print("Columns:", [col[0] for col in columns])

    sample = conn.execute("SELECT * FROM attributes LIMIT 3").fetchall()
    print("Sample data:", sample)
except Exception as e:
    print(f"Error: {e}")
    print("Let's try alternative format...")

try:
    print("Loading nodes (this may take a few minutes)...")
    conn.execute("""
        CREATE TABLE nodes AS
        SELECT * FROM read_json_auto('./actual_db_arn2/nodes_RC2_2023_06_06.json')
    """)

    count = conn.execute("SELECT COUNT(*) FROM nodes").fetchone()[0]
    print(f" Nodes loaded: {count:,} records")

    columns = conn.execute("DESCRIBE nodes").fetchall()
    print(f"Node columns ({len(columns)}): {[col[0] for col in columns[:5]]}")

except Exception as e:
    print(f" Nodes error: {e}")

try:
    print("Loading edges...")
    conn.execute("""
        CREATE TABLE edges AS
        SELECT * FROM read_json_auto('./actual_db_arn2/edges_RC2_2023_06_06.json')
    """)

    count = conn.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    print(f" Edges loaded: {count:,} records")

    columns = conn.execute("DESCRIBE edges").fetchall()
    print(f"Edge columns ({len(columns)}): {[col[0] for col in columns[:5]]}")

except Exception as e:
    print(f" Edges error: {e}")

conn.close()

conn = duckdb.connect('test.db')

tables = conn.execute("SHOW TABLES").fetchall()

conn.execute("ATTACH 'graph_database_sqlite.db' AS sqlite_db (TYPE sqlite)")

for table in tables:
    table_name = table[0]
    print(f"Exporting {table_name}...")
    conn.execute(f"CREATE TABLE sqlite_db.{table_name} AS SELECT * FROM {table_name}")

conn.close()
print(" Export complete!")
