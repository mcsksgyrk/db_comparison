import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import duckdb
import sqlite3
import pandas as pd
from config import OUTPUT_DIR, PROJECT_ROOT
from typing import List, Dict, Set
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


def find_node_edges(node, edges_df, cond_set=None):
    if cond_set is not None:
        res = edges_df[
            (edges_df.interactor_a_node_name == node) &
            (edges_df.interactor_b_node_name.isin(cond_set))
        ]
    else:
        res = edges_df[
            (edges_df.interactor_a_node_name == node)
        ]
    return res


def analyize_downstream_connectivity_nc(node_set: Set, edges_df: pd.DataFrame,
                                        arn_core: Set, fer_core: Set) -> Dict:
    shared_core = arn_core & fer_core
    arn_only_core = arn_core - shared_core
    fer_only_core = fer_core - shared_core
    res = {}
    for node in list(node_set):
        arn_conn = len(find_node_edges(node, edges_df, arn_only_core))
        fer_conn = len(find_node_edges(node, edges_df, fer_only_core))
        shared_conn = len(find_node_edges(node, edges_df, shared_core))
        all_conn = len(find_node_edges(node, edges_df))
        score = arn_conn + fer_conn + shared_conn
        crosstalk_score = min(arn_conn+shared_conn, fer_conn+shared_conn)
        res[node] = {
            'fer_conn': fer_conn,
            'arn_conn': arn_conn,
            'shared_conn': shared_conn,
            'all_conn': all_conn,
            'score': score,
            'crosstalk_score': crosstalk_score
        }
    return res


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

#from apicalls.api_oop import BGEEClient, UniProtClient
#
#target_localisations = {
#    'dendritic_cell': 'CL:0000451',
#    'macrophage': 'CL:0000235',
#    'fibroblast': 'CL:0000057',
#    't_cell': 'CL:0000084',
#    'epithelial_cell': 'CL:0000066'
#}
#uniprot = UniProtClient()
#bgee = BGEEClient()
#anat_loc_dicts = dict()
#for pr in common_core:
#    ensembl_id, _ = uniprot.batch_convert_from_uniprot_id("Ensembl", [pr])
#    res = bgee.get_expression_anat_entity(ensembl_id[pr])
#    anat_loc_dicts[pr] = {
#        'ensembl_id': ensembl_id[pr],
#        'anatEntity': res[ensembl_id[pr]]
#    }
#
#target_uberon_ids = set(target_localisations.values())
#test = {'UBERON:0000059', 'CL:0000039', 'CL:0000039'}
#for k, v in anat_loc_dicts.items():
#    print(set(v['anatEntity']) & test)

res = analyize_downstream_connectivity_nc(shared_regulators_set, edges_df,
                                          arn_core, fer_core)

df_upstream_reg = pd.DataFrame.from_dict(res, orient='index')
nodes_to_analyze = df_upstream_reg[
    (df_upstream_reg.crosstalk_score > 2) &
    (df_upstream_reg.fer_conn > 0)
].index.tolist()

shared_core = arn_core & fer_core
arn_only_core = arn_core - shared_core
fer_only_core = fer_core - shared_core

regulators = {}
for node in nodes_to_analyze:
    edges = find_node_edges(node, edges_df, fer_core | arn_core)

    if node not in regulators:
        regulators[node] = {'fer': [],
                            'arn': [],
                            'shared': []}

    for idx, edge in edges.iterrows():
        if edge.interactor_b_node_name in shared_core:
            regulators[node]['shared'].append(edge.interactor_b_node_name)
        elif edge.interactor_b_node_name in arn_only_core:
            regulators[node]['arn'].append(edge.interactor_b_node_name)
        elif edge.interactor_b_node_name in fer_only_core:
            print(f"{edge.interactor_a_node_name} effects {edge.interactor_b_node_name}")
            regulators[node]['fer'].append(edge.interactor_b_node_name)

for node in regulators:
    regulators[node]['fer'] = list(set(regulators[node]['fer']))
    regulators[node]['arn'] = list(set(regulators[node]['arn']))
    regulators[node]['shared'] = list(set(regulators[node]['shared']))

import json
target_pr_path = OUTPUT_DIR / "upstream_regulators.json"
with open(target_pr_path, "w") as f:
    json.dump(regulators, f)

