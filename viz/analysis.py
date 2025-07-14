import gseapy as gp
from gseapy import barplot, dotplot
from database.db_api import DuckdbAPI
from config import OUTPUT_DIR, PROJECT_ROOT
from apicalls.api_oop import UniProtClient
import pandas as pd


db_path = OUTPUT_DIR / "test2.duckdb"
sql_seed = PROJECT_ROOT / 'database' / 'duckdb_seed.sql'
db_api = DuckdbAPI(sql_seed, db_path, create_new=False)

df_prs = db_api.db.execute("SELECT * FROM node").fetch_df()
db_api.close()

arn_prs = df_prs[df_prs.source_database == "ARN"].name.to_list()
fer_prs = df_prs[df_prs.source_database == "ferr"].name.to_list()
com_prs = df_prs[df_prs.source_database == "ARN|ferr"].name.to_list()

uniprot_client = UniProtClient()
arn_gene = uniprot_client.batch_convert_from_uniprot_id("GeneCards", arn_prs, batch_size=100)
arn_gene[0].values()

libs = gp.get_library_name(organism='human')
go_libs = [lib for lib in libs if 'GO_Biological_Process' in lib]

# it needs gene symbols
arn_go = gp.enrichr(gene_list=list(arn_gene[0].values()),
                    gene_sets=go_libs,
                    organism='human',
                    background=None)
arn_go.results
ax = barplot(arn_go.res2d,title='Test', figsize=(4, 5), color='darkred')
ax = dotplot(arn_go.res2d, title='KEGG_2021_Human',cmap='viridis_r', size=10, figsize=(3,5))
import matplotlib.pyplot as plt
plt.show()
