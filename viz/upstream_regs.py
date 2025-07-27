import gseapy as gp
from config import OUTPUT_DIR, PROJECT_ROOT
from apicalls.api_oop import UniProtClient
from typing import Set, List
import json
import matplotlib.pyplot as plt


def check_term_overlap(pathway1: Set, pathway2: Set) -> List[str]:
    overlap = pathway1 & pathway2
    return list(overlap)


db_path = OUTPUT_DIR / "test2.duckdb"
upstream_regs_path = OUTPUT_DIR / "upstream_regulators.json"
uniprot_client = UniProtClient()

with open(upstream_regs_path, "r") as f:
    upstream_prs = json.load(f)

pr_ids = list(upstream_prs.keys())
gene_names_dict, failed_ids = uniprot_client.batch_convert_from_uniprot_id("GeneCards", pr_ids, batch_size=100)
gene_names_dict
libs = gp.get_library_name(organism='human')

kegg_libs = [lib for lib in libs if "KEGG" in lib]
go_libs = [lib for lib in libs if 'GO_Biological_Process' in lib]

libraries = [
    'GO_Biological_Process_2025',
    'KEGG_2021_Human',
    'Reactome_Pathways_2024'
]
disease_libs = [
    'DisGeNET',
    'OMIM_Disease',
    'Jensen_DISEASES_Curated_2025'
]
specialized_libs = [
    'MSigDB_Hallmark_2020',
    'Disease_Perturbations_from_GEO_up'
]

enrich = gp.enrichr(gene_list=list(gene_names_dict.values()),
                    gene_sets=disease_libs,
                    organism='human',
                    background=None)

ax = gp.barplot(enrich.res2d,title='KEGG_2021_Human', figsize=(4, 5), color='darkred')
plt.show()
