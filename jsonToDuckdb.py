import duckdb
from config import OUTPUT_DIR, SOURCE_DIR
import os
import json

files = ['nodes_RC2_2023_06_06.json', 'edges_RC2_2023_06_06.json', 'attributes_RC2_2023_06_06.json']

for file in files:
    path = SOURCE_DIR / file
    if os.path.exists(path):
        size_mb = os.path.getsize(path) / (1024**2)
        print(f"{path}: {size_mb:.1f} MB")
    else:
        print(f"{path}: NOT FOUND")

output_db = OUTPUT_DIR / 'rawdata.duckdb'
conn = duckdb.connect(output_db)
try:
    print("Loading attributes...")
    attributes_json = SOURCE_DIR / 'attributes_RC2_2023_06_06.json'
    conn.execute(f"""
        CREATE TABLE attributes AS
        SELECT * FROM read_json_auto('{attributes_json}')
    """)

except Exception as e:
    print(f"Error: {e}")

try:
    print("Loading nodes...")
    nodes_json = SOURCE_DIR / 'nodes_RC2_2023_06_06.json'
    conn.execute(f"""
        CREATE TABLE nodes AS
        SELECT * FROM read_json_auto('{nodes_json}')
    """)

except Exception as e:
    print(f" Nodes error: {e}")

try:
    print("Loading edges...")
    edges_json = SOURCE_DIR / 'edges_RC2_2023_06_06.json'
    conn.execute(f"""
        CREATE TABLE edges AS
        SELECT * FROM read_json_auto('{edges_json}')
    """)

except Exception as e:
    print(f" Edges error: {e}")

tables = conn.execute("SHOW TABLES").fetchall()

conn.close()
print(" Export complete!")
