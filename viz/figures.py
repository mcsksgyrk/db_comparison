import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import duckdb
import sqlite3
import pandas as pd
from config import OUTPUT_DIR, PROJECT_ROOT
from typing import List
from pathlib import Path


def get_core_arn_proteins(arn_path: Path) -> List[str]:
    arn_db = duckdb.connect(arn_path)
    query = """
            SELECT name
            FROM nodes
            WHERE len(autophagyPhase) > 0
            """
    res = arn_db.execute(query).fetchall()
    arn_db.close()
    return [val[0] for val in res]


def get_core_fer_proteins(ferr_path: Path) -> List[str]:
    conn = sqlite3.connect(ferr_path)
    curs = conn.cursor()
    query = """
            SELECT name
            FROM node
            WHERE source='KEGG'
            AND type='protein'
            """
    curs.execute(query)
    res = curs.fetchall()
    conn.close()
    return [val[0] for val in res]


db_path = OUTPUT_DIR / 'test2.duckdb'
arn_path = OUTPUT_DIR / 'ARN.duckdb'
fer_path = PROJECT_ROOT / 'test_omnipath.db'

arn_core = set(get_core_arn_proteins(arn_path))
fer_core = set(get_core_fer_proteins(fer_path))
common_core = arn_core & fer_core
db = duckdb.connect(db_path)
sources = ['ARN', 'ferr', 'ARN|ferr']
group_members = []
for source in sources:
    query_nodes = f"""
                  SELECT * FROM node WHERE source_database = '{source}'
                  """
    group_members.append(len(db.execute(query_nodes).fetchdf()))

plt.figure(figsize=(8, 6), dpi=150)
venn2(subsets=group_members,
      set_labels=('Autophagy', 'Ferroptosis'))
plt.title('Node Overlap Between Pathways')
plt.show()

edges_df = db.execute("SELECT * FROM edge").fetchdf()

""" ARN_core → FER_core  or  FER_core → ARN_core """
crosstalk = edges_df[
    (edges_df.interactor_a_node_name.isin(arn_core)) &
    (edges_df.interactor_b_node_name.isin(fer_core)) |
    (edges_df.interactor_a_node_name.isin(fer_core)) &
    (edges_df.interactor_b_node_name.isin(arn_core))
]
""" ARN_core → X ← FER_core  (X = shared target) """
arn_targets = set(edges_df[
    edges_df.interactor_a_node_name.isin(arn_core)
    ].interactor_b_node_name)

fer_targets = set(edges_df[
    edges_df.interactor_a_node_name.isin(fer_core)
    ].interactor_b_node_name)

shared_targets = arn_targets & fer_targets

shared_target_edges = edges_df[
    (edges_df.interactor_a_node_name.isin(arn_core | fer_core)) &
    (edges_df.interactor_b_node_name.isin(shared_targets))
]
""" X → ARN_core  and  X → FER_core  (X = shared regulator) """
arn_regulators = set(edges_df[
    edges_df.interactor_b_node_name.isin(arn_core)
].interactor_a_node_name)

fer_regulators = set(edges_df[
    edges_df.interactor_b_node_name.isin(fer_core)
].interactor_a_node_name)

shared_regulators_set = arn_regulators & fer_regulators

shared_regulator_edges = edges_df[
    (edges_df.interactor_a_node_name.isin(shared_regulators_set)) &
    (edges_df.interactor_b_node_name.isin(arn_core | fer_core))
]

db.close()

from apicalls.api_oop import BGEEClient, UniProtClient

target_localisations = {
    'dendritic_cell': 'CL:0000451',
    'macrophage': 'CL:0000235',
    'fibroblast': 'CL:0000057',
    't_cell': 'CL:0000084',
    'epithelial_cell': 'CL:0000066'
}
uniprot = UniProtClient()
bgee = BGEEClient()
anat_loc_dicts = dict()
for pr in common_core:
    ensembl_id, _ = uniprot.batch_convert_from_uniprot_id("Ensembl", [pr])
    res = bgee.get_expression_anat_entity(ensembl_id[pr])
    anat_loc_dicts[pr] = {
        'ensembl_id': ensembl_id[pr],
        'anatEntity': res[ensembl_id[pr]]
    }

target_uberon_ids = set(target_localisations.values())
test = {'UBERON:0000059', 'CL:0000039', 'CL:0000039'}
for k, v in anat_loc_dicts.items():
    print(set(v['anatEntity']) & test)

anat_loc_dicts['A6NCE7']
