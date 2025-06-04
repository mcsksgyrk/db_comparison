import duckdb
import sqlite3
import json

def extract_simple_values(complex_value):
    """Extract useful values from complex nested structures"""
    if complex_value is None:
        return ''

    if isinstance(complex_value, (str, int, float, bool)):
        return str(complex_value)

    if isinstance(complex_value, list):
        if len(complex_value) == 0:
            return ''

        # Handle list of dictionaries (like publications, databases)
        if isinstance(complex_value[0], dict):
            values = []
            for item in complex_value:
                if 'value' in item:
                    values.append(str(item['value']))
                elif 'name' in item:
                    values.append(str(item['name']))
            return '; '.join(values) if values else ''

        # Handle simple lists
        else:
            return '; '.join(str(item) for item in complex_value)

    if isinstance(complex_value, dict):
        # Handle single objects like taxon
        if 'name' in complex_value:
            return str(complex_value['name'])
        elif 'value' in complex_value:
            return str(complex_value['value'])
        elif 'id' in complex_value:
            return str(complex_value['id'])
        else:
            # Fallback: return all values joined
            return '; '.join(str(v) for v in complex_value.values() if v)

    return str(complex_value)

def clean_biological_data():
    """Clean your biological network database"""
    print("🧬 Cleaning biological network database...")

    # Connect to DuckDB
    duckdb_conn = duckdb.connect('test.db')

    # Create clean SQLite database
    sqlite_conn = sqlite3.connect('biological_network_clean.db')

    # Clean each table with specific logic
    tables_config = {
        'nodes': {
            'simple_columns': ['displayedName', 'fullName', 'name'],
            'extract_columns': {
                'taxon_name': 'taxon',  # Extract species name
                'pathways_list': 'pathways',
                'tissues_list': 'tissues',
                'cellular_location_major': 'majorCellularLocalization',
                'cellular_location_minor': 'minorCellularLocalization',
                'autophagy_phase': 'autophagyPhase',
                'go_annotation': 'geneOntologyAnnotation',
                'phenotype': 'phenotypeAnnotation',
                'molecule_type': 'moleculeType',
                'external_refs': 'externalReferences'
            }
        },
        'edges': {
            'simple_columns': ['sourceDisplayedName', 'targetDisplayedName', 'source', 'target', 'isDirected', 'isDirect'],
            'extract_columns': {
                'interaction_methods': 'interactionDetectionMethods',
                'confidence_score': 'confidenceScore',
                'publications': 'publications',
                'source_databases': 'sourceDatabases',
                'interaction_type': 'interactionType'
            }
        },
        'attributes': {
            'simple_columns': ['dataType', 'displayedName', 'isNode', 'key', 'searchable'],
            'extract_columns': {
                'possible_values': 'possibleValues'
            }
        }
    }

    for table_name, config in tables_config.items():
        print(f"\n📊 Processing {table_name}...")

        try:
            # Get original data
            data = duckdb_conn.execute(f"SELECT * FROM {table_name}").fetchall()
            original_columns = [desc[0] for desc in duckdb_conn.execute(f"DESCRIBE {table_name}").fetchall()]

            # Create new column structure
            new_columns = config['simple_columns'] + list(config['extract_columns'].keys())

            # Create SQLite table
            col_defs = [f'"{col}" TEXT' for col in new_columns]
            create_sql = f'CREATE TABLE "{table_name}" ({", ".join(col_defs)})'
            sqlite_conn.execute(create_sql)

            # Process each row
            cleaned_rows = []
            for row in data:
                row_dict = dict(zip(original_columns, row))
                new_row = []

                # Add simple columns as-is
                for col in config['simple_columns']:
                    value = row_dict.get(col, '')
                    new_row.append(str(value) if value is not None else '')

                # Extract complex columns
                for new_col, original_col in config['extract_columns'].items():
                    complex_value = row_dict.get(original_col, '')
                    extracted = extract_simple_values(complex_value)
                    new_row.append(extracted)

                cleaned_rows.append(new_row)

            # Insert cleaned data
            placeholders = ', '.join(['?' for _ in new_columns])
            insert_sql = f'INSERT INTO "{table_name}" VALUES ({placeholders})'
            sqlite_conn.executemany(insert_sql, cleaned_rows)

            print(f"   ✅ {len(cleaned_rows):,} records cleaned")
            print(f"   📋 Columns: {', '.join(new_columns[:4])}...")

        except Exception as e:
            print(f"   ❌ Error processing {table_name}: {e}")

    # Create summary
    print(f"\n📈 Creating analysis views...")

    # Gene-focused view
    sqlite_conn.execute("""
        CREATE VIEW genes AS
        SELECT
            name as gene_id,
            displayedName as gene_symbol,
            taxon_name as species,
            pathways_list as pathways,
            go_annotation as gene_ontology,
            cellular_location_major as location,
            molecule_type
        FROM nodes
        WHERE molecule_type LIKE '%protein%' OR molecule_type LIKE '%gene%'
    """)

    # Interaction summary
    sqlite_conn.execute("""
        CREATE VIEW interactions_summary AS
        SELECT
            sourceDisplayedName as gene1,
            targetDisplayedName as gene2,
            interaction_type,
            confidence_score,
            publications,
            source_databases
        FROM edges
        WHERE sourceDisplayedName IS NOT NULL
        AND targetDisplayedName IS NOT NULL
    """)

    sqlite_conn.commit()
    sqlite_conn.close()
    duckdb_conn.close()

    print(f"\n🎉 SUCCESS! Clean database: biological_network_clean.db")
    print(f"📋 Ready for ferroptosis analysis!")

    return 'biological_network_clean.db'

# Run the cleaning
clean_db = clean_biological_data()

# Quick verification
print(f"\n🔍 VERIFICATION:")
conn = sqlite3.connect(clean_db)
cursor = conn.cursor()

tables = cursor.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()
for table in tables:
    table_name = table[0]
    count = cursor.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()[0]
    print(f"   {table_name}: {count:,} records")

# Show sample from nodes
print(f"\n🧬 Sample genes:")
sample = cursor.execute("SELECT gene_symbol, species, pathways FROM genes LIMIT 5").fetchall()
for row in sample:
    print(f"   {row[0]}: {row[1]} | {row[2][:50]}...")

conn.close()
