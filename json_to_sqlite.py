import duckdb

conn = duckdb.connect('test.db')

tables = conn.execute("SHOW TABLES").fetchall()
for table in tables:
    table_name = table[0]
    print(f"\n=== {table_name} ===")
    columns = conn.execute(f"DESCRIBE {table_name}").fetchall()
    print("Columns:", [col[0] for col in columns])

    # Show sample data to see the format
    sample = conn.execute(f"SELECT * FROM {table_name} LIMIT 2").fetchall()
    print("Sample data:", sample[0] if sample else "No data")

conn.close()
