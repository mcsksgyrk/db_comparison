import gseapy as gp
from database.db_api import DuckdbAPI
from config import OUTPUT_DIR, PROJECT_ROOT
from apicalls.api_oop import UniProtClient
import matplotlib.pyplot as plt
from typing import Set, List
import json


def check_term_overlap(pathway1: Set, pathway2: Set) -> List[str]:
    overlap = pathway1 & pathway2
    return list(overlap)


db_path = OUTPUT_DIR / "test2.duckdb"
sql_seed = PROJECT_ROOT / 'database' / 'duckdb_seed.sql'
db_api = DuckdbAPI(sql_seed, db_path, create_new=False)
uniprot_client = UniProtClient()

df_prs = db_api.db.execute("SELECT * FROM node").fetch_df()
db_api.close()

pathway_groups = df_prs.source_database.unique()
gene_sets = dict()

for pathway in pathway_groups:
    uniprot_ids = df_prs[df_prs.source_database == pathway].name.to_list()
    gene_names_dict, failed_ids = uniprot_client.batch_convert_from_uniprot_id("GeneCards", uniprot_ids, batch_size=100)

    gene_sets[pathway] = {
        'uniprot_ids': list(gene_names_dict.keys()) if gene_names_dict else [],
        'gene_names': list(gene_names_dict.values()) if gene_names_dict else [],
        'failed_ids': failed_ids,
        'enrichment': None,
        'significant_terms': None
    }

libs = gp.get_library_name(organism='human')
go_libs = [lib for lib in libs if 'GO_Biological_Process' in lib]

for pathway, data in gene_sets.items():
    enrich = gp.enrichr(gene_list=list(data['gene_names']),
                        gene_sets='Reactome_Pathways_2024',
                        organism='human',
                        background=None)

    gene_sets[pathway]['enrichment'] = enrich
    gene_sets[pathway]['significant_terms'] = enrich.results[enrich.results['Adjusted P-value'] < 0.05]

intersections = dict()
for i in range(len(gene_sets.keys())):
    for j in range(i+1, len(gene_sets.keys())):
        pathway1 = pathway_groups[i]
        pathway2 = pathway_groups[j]

        key = pathway1+"&"+pathway2

        terms1 = set(gene_sets[pathway1]['significant_terms']['Term'])
        terms2 = set(gene_sets[pathway2]['significant_terms']['Term'])

        intersections[key] = check_term_overlap(terms1, terms2)

# ferr and ARN groups no same term
# common and arn 45
# common and ferr 220
for k, v in intersections.items():
    print(f"{k} have {len(v)} many common functions")

intersect_of_intersects = check_term_overlap(set(intersections['ARN&ARN|ferr']),
                                             set(intersections['ARN|ferr&ferr']))
intersect_of_intersects
set(intersections['ARN&ferr']) & set(intersect_of_intersects)
