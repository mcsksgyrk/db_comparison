import gseapy as gp
from config import OUTPUT_DIR, PROJECT_ROOT
from apicalls.api_oop import UniProtClient
from typing import Set, List
import duckdb
import matplotlib.pyplot as plt
import json


def check_term_overlap(pathway1: Set, pathway2: Set) -> List[str]:
    overlap = pathway1 & pathway2
    return list(overlap)


db_path = OUTPUT_DIR / "test2.duckdb"
uniprot_client = UniProtClient()

db = duckdb.connect(db_path)

query_nodes = """
              SELECT * FROM node
              """
db_df = db.execute(query_nodes).fetchdf()

libs = gp.get_library_name(organism='human')

kegg_libs = [lib for lib in libs if "KEGG" in lib]
go_libs = [lib for lib in libs if 'GO_Biological_Process' in lib]

libraries = [
    'GO_Biological_Process_2025',
    'KEGG_2021_Human',
    'Reactome_Pathways_2024',
    'DisGeNET',
    'OMIM_Disease',
    'Jensen_DISEASES_Curated_2025',
    'MSigDB_Hallmark_2020',
    'Disease_Perturbations_from_GEO_up'
]
sources = ['ARN', 'ferr', 'ARN|ferr']
enrich_results = {}
for source in sources:
    pr_group = list(db_df[db_df.source_database == source].name)
    gene_names_dict, failed_ids = uniprot_client.batch_convert_from_uniprot_id("GeneCards",
                                                                               pr_group,
                                                                               batch_size=100)

    if failed_ids:
        print(f"failed conevrsion of ids {failed_ids}")
    enrich = gp.enrichr(gene_list=list(gene_names_dict.values()),
                        gene_sets=libraries,
                        organism='human',
                        background=None)
    enrich_results[source] = enrich

import pickle
with open('enrichment_results.pkl', 'wb') as f:
    pickle.dump(enrich_results, f)
