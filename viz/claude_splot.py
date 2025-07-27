import pickle
import gseapy as gp
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Load enrichment results
with open('enrichment_results.pkl', 'rb') as f:
    enrichr_result = pickle.load(f)

# Define consistent IBD search terms for all groups
IBD_SEARCH_PATTERN = 'IBD|inflammatory bowel|Crohn|colitis|ulcerative colitis|intestinal inflammation'

print("IBD ANALYSIS: AUTOPHAGY-FERROPTOSIS INTERSECTION")
print("=" * 60)

# Analyze intersection results
intersection_results = enrichr_result['ARN|ferr'].results

# Get significant disease associations
sig_diseases = intersection_results[
    (intersection_results['Adjusted P-value'] < 0.05) &
    (intersection_results['Gene_set'].isin(['DisGeNET', 'OMIM_Disease', 'Jensen_DISEASES_Curated_2025']))
]

print(f"Significant disease associations (p < 0.05): {len(sig_diseases)}")
print("Top 10 most significant:")
for _, row in sig_diseases.nsmallest(10, 'Adjusted P-value').iterrows():
    print(f"  {row['Term']} (p={row['Adjusted P-value']:.2e}, overlap={row['Overlap']})")

# Search for IBD terms
ibd_terms = intersection_results[
    intersection_results['Term'].str.contains(IBD_SEARCH_PATTERN, case=False, na=False)
]

print(f"\nAll IBD-related terms found: {len(ibd_terms)}")
for _, row in ibd_terms.iterrows():
    print(f"  {row['Term']} (p={row['Adjusted P-value']:.2e}, overlap={row['Overlap']})")

# Filter for significant IBD terms
sig_ibd_terms = ibd_terms[ibd_terms['Adjusted P-value'] < 0.05]
print(f"\nSignificant IBD terms (p < 0.05): {len(sig_ibd_terms)}")

if len(sig_ibd_terms) > 0:
    print("Significant IBD terms:")
    for _, row in sig_ibd_terms.iterrows():
        overlap_nums = row['Overlap'].split('/')
        genes_in_term = int(overlap_nums[0])
        total_genes = int(overlap_nums[1])
        print(f"  {row['Term']}")
        print(f"    P-value: {row['Adjusted P-value']:.2e}")
        print(f"    Overlap: {genes_in_term}/{total_genes} genes")
        print(f"    Gene set: {row['Gene_set']}")
        genes_display = row['Genes'][:100] + "..." if len(row['Genes']) > 100 else row['Genes']
        print(f"    Genes: {genes_display}")

    # Filter for high-confidence terms (>=3 genes)
    interesting = sig_ibd_terms[sig_ibd_terms['Overlap'].str.split('/').str[0].astype(int) >= 3]
    print(f"\nHigh-confidence IBD terms (>=3 genes): {len(interesting)}")
    for term in interesting['Term']:
        print(f"  {term}")

    # Create visualizations
    print("\nCreating visualizations...")

    # Custom dot plot
    plt.figure(figsize=(12, 8))
    y_pos = np.arange(len(sig_ibd_terms))
    gene_counts = [int(overlap.split('/')[0]) for overlap in sig_ibd_terms['Overlap']]

    scatter = plt.scatter(
        -np.log10(sig_ibd_terms['Adjusted P-value']),
        y_pos,
        s=[count * 50 for count in gene_counts],
        c=-np.log10(sig_ibd_terms['Adjusted P-value']),
        cmap='viridis',
        alpha=0.7,
        edgecolors='black'
    )

    plt.yticks(y_pos, sig_ibd_terms['Term'])
    plt.xlabel('-log10(Adjusted P-value)')
    plt.title('IBD Terms - Autophagy-Ferroptosis Intersection\n(Dot size = number of genes)')
    plt.colorbar(scatter, label='-log10(Adjusted P-value)')
    plt.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='p=0.05')
    plt.legend()
    plt.tight_layout()
    plt.savefig('filename.png', dpi=300, bbox_inches='tight')
    plt.show()

    # Bar plot
    plt.figure(figsize=(12, 8))
    sorted_terms = sig_ibd_terms.sort_values('Adjusted P-value')
    y_pos = np.arange(len(sorted_terms))
    bars = plt.barh(y_pos, -np.log10(sorted_terms['Adjusted P-value']),
                   color='steelblue', alpha=0.7)

    plt.yticks(y_pos, sorted_terms['Term'])
    plt.xlabel('-log10(Adjusted P-value)')
    plt.title('IBD Enrichment - Autophagy-Ferroptosis Intersection')
    plt.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='p=0.05')

    for i, (bar, pval) in enumerate(zip(bars, sorted_terms['Adjusted P-value'])):
        plt.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
                f'{pval:.2e}', va='center', fontsize=10)

    plt.legend()
    plt.tight_layout()
    plt.show()

else:
    print("No significant IBD terms found")

# Compare with individual groups
print("\nComparison with individual groups:")

# Collect IBD results for all groups
group_ibd_results = {}
for group_name in ['ARN', 'ferr', 'ARN|ferr']:
    if group_name in enrichr_result:
        if group_name == 'ARN|ferr':
            group_results = intersection_results
        else:
            group_results = enrichr_result[group_name].results

        group_ibd = group_results[
            group_results['Term'].str.contains(IBD_SEARCH_PATTERN, case=False, na=False)
        ]
        group_sig_ibd = group_ibd[group_ibd['Adjusted P-value'] < 0.05]

        group_ibd_results[group_name] = {
            'total': len(group_ibd),
            'significant': len(group_sig_ibd),
            'terms': group_sig_ibd['Term'].tolist()
        }

        print(f"\n{group_name} group:")
        print(f"  Total IBD terms: {len(group_ibd)}")
        print(f"  Significant IBD terms: {len(group_sig_ibd)}")

        if len(group_sig_ibd) > 0:
            print("  Significant terms:")
            for _, row in group_sig_ibd.iterrows():
                print(f"    {row['Term']} (p={row['Adjusted P-value']:.2e})")

# Create comparison figure
print("\nCreating group comparison figure...")
plt.figure(figsize=(10, 6))

groups = list(group_ibd_results.keys())
total_counts = [group_ibd_results[group]['total'] for group in groups]
sig_counts = [group_ibd_results[group]['significant'] for group in groups]
non_sig_counts = [total - sig for total, sig in zip(total_counts, sig_counts)]

x = np.arange(len(groups))
width = 0.6

# Stack non-significant and significant IBD terms
bars1 = plt.bar(x, non_sig_counts, width, label='Non-significant IBD terms', alpha=0.7, color='lightgray')
bars2 = plt.bar(x, sig_counts, width, bottom=non_sig_counts, label='Significant IBD terms', alpha=0.7, color='darkblue')

plt.xlabel('Groups')
plt.ylabel('Number of IBD terms')
plt.title('IBD Terms: Significant vs Non-significant')
plt.xticks(x, groups)
plt.legend()

# Add value labels
for i, (total, sig, non_sig) in enumerate(zip(total_counts, sig_counts, non_sig_counts)):
    if sig > 0:
        plt.text(i, non_sig + sig/2, str(sig), ha='center', va='center', fontweight='bold', color='white')
    if non_sig > 0:
        plt.text(i, non_sig/2, str(non_sig), ha='center', va='center', fontweight='bold', color='black')
    plt.text(i, total + 0.2, f'Total: {total}', ha='center', va='bottom', fontsize=10)

plt.tight_layout()
plt.show()

print(f"\nSummary: Found {len(sig_ibd_terms)} significant IBD terms in ARN|ferr intersection")
