import pickle
import gseapy as gp
import matplotlib.pyplot as plt

with open('enrichment_results.pkl', 'rb') as f:
    enrichr_result = pickle.load(f)

intersection_results = enrichr_result['ARN|ferr'].results
sig_diseases = intersection_results[
    (intersection_results['Adjusted P-value'] < 0.5) &
    (intersection_results['Gene_set'].isin(['DisGeNET', 'OMIM_Disease', 'Jensen_DISEASES_Curated_2025']))
]

print("Shared disease associations:")
for _, row in sig_diseases.nsmallest(10, 'Adjusted P-value').iterrows():
    print(f"{row['Term']} (p={row['Adjusted P-value']:.2e})")

ibd_terms = intersection_results[
    intersection_results['Term'].str.contains('IBD|inflammatory bowel|Crohn|colitis', case=False, na=False)
]
sig_ibd_terms = ibd_terms[ibd_terms['Adjusted P-value'] < 0.05]
interesting = sig_ibd_terms[sig_ibd_terms['Overlap'].str.split('/').str[0].astype(int) >= 3]
interesting.Term
sig_ibd_terms.Term

gp.dotplot(sig_ibd_terms,
               title='IBD Terms - Autophagy-Ferroptosis Intersection',
               figsize=(8, 6),
               size_factor=200)

plt.show()
# Built-in bar plot
gp.barplot(sig_ibd_terms,
               title='IBD Enrichment',
               figsize=(10, 6))
plt.show()
