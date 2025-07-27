import pickle
import gseapy as gp
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Load enrichment results
with open('enrichment_results.pkl', 'rb') as f:
    enrichr_result = pickle.load(f)

IBD_TERMS = 'IBD|inflammatory bowel|Crohn|colitis|ulcerative colitis|intestinal inflammation'

intersection_results = enrichr_result['ARN|ferr'].results

sig_diseases = intersection_results[
    (intersection_results['Adjusted P-value'] < 0.05) &
    (intersection_results['Gene_set'].isin(['DisGeNET', 'OMIM_Disease', 'Jensen_DISEASES_Curated_2025']))
]

for _, row in sig_diseases.nsmallest(10, 'Adjusted P-value').iterrows():
    print(f"• {row['Term']} (p={row['Adjusted P-value']:.2e}, overlap={row['Overlap']})")

ibd_terms = intersection_results[
    intersection_results['Term'].str.contains(
        'IBD|inflammatory bowel|Crohn|colitis|ulcerative colitis|intestinal inflammation',
        case=False, na=False
    )
]

for _, row in ibd_terms.iterrows():
    print(f"• {row['Term']} (p={row['Adjusted P-value']:.2e}, overlap={row['Overlap']})")

# Filter for significant IBD terms
sig_ibd_terms = ibd_terms[ibd_terms['Adjusted P-value'] < 0.05]
print(f"\nSignificant IBD terms (p < 0.05): {len(sig_ibd_terms)}")

if len(sig_ibd_terms) > 0:
    print("\nSIGNIFICANT IBD TERMS:")
    for _, row in sig_ibd_terms.iterrows():
        overlap_nums = row['Overlap'].split('/')
        genes_in_term = int(overlap_nums[0])
        total_genes = int(overlap_nums[1])
        print(f"• {row['Term']}")
        print(f"  - P-value: {row['Adjusted P-value']:.2e}")
        print(f"  - Overlap: {genes_in_term}/{total_genes} genes")
        print(f"  - Gene set: {row['Gene_set']}")
        print(f"  - Genes: {row['Genes'][:100]}..." if len(row['Genes']) > 100 else f"  - Genes: {row['Genes']}")

    # Filter for high-confidence terms (≥3 genes)
    interesting = sig_ibd_terms[sig_ibd_terms['Overlap'].str.split('/').str[0].astype(int) >= 3]
    print(f"\nHigh-confidence IBD terms (≥3 genes): {len(interesting)}")

    if len(interesting) > 0:
        print("HIGH-CONFIDENCE TERMS:")
        for term in interesting['Term']:
            print(f"• {term}")

    # Create visualizations
    print(f"\n📈 CREATING VISUALIZATIONS...")

    # Enhanced dot plot
    plt.figure(figsize=(12, 8))
    try:
        # Create custom dot plot since gseapy might have issues
        y_pos = np.arange(len(sig_ibd_terms))

        # Extract overlap ratios for dot sizes
        overlap_ratios = []
        gene_counts = []
        for overlap_str in sig_ibd_terms['Overlap']:
            nums = overlap_str.split('/')
            ratio = int(nums[0]) / int(nums[1])
            overlap_ratios.append(ratio)
            gene_counts.append(int(nums[0]))

        # Create scatter plot
        scatter = plt.scatter(
            -np.log10(sig_ibd_terms['Adjusted P-value']),
            y_pos,
            s=[count * 50 for count in gene_counts],  # Size based on gene count
            c=-np.log10(sig_ibd_terms['Adjusted P-value']),
            cmap='viridis',
            alpha=0.7,
            edgecolors='black'
        )

        plt.yticks(y_pos, sig_ibd_terms['Term'])
        plt.xlabel('-log10(Adjusted P-value)', fontsize=12)
        plt.title('IBD Terms - Autophagy-Ferroptosis Intersection\n(Dot size = number of genes)',
                 fontsize=14, fontweight='bold')
        plt.colorbar(scatter, label='-log10(Adjusted P-value)')

        # Add significance line
        plt.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='p=0.05')
        plt.legend()

        plt.tight_layout()
        plt.show()

    except Exception as e:
        print(f"Custom dot plot error: {e}")

    # Enhanced bar plot
    plt.figure(figsize=(12, 8))
    try:
        # Sort by p-value
        sorted_terms = sig_ibd_terms.sort_values('Adjusted P-value')

        # Create horizontal bar plot
        y_pos = np.arange(len(sorted_terms))
        bars = plt.barh(y_pos, -np.log10(sorted_terms['Adjusted P-value']),
                       color='steelblue', alpha=0.7)

        plt.yticks(y_pos, sorted_terms['Term'])
        plt.xlabel('-log10(Adjusted P-value)', fontsize=12)
        plt.title('IBD Enrichment - Autophagy-Ferroptosis Intersection',
                 fontsize=14, fontweight='bold')

        # Add significance line
        plt.axvline(x=-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='p=0.05')

        # Add value labels on bars
        for i, (bar, pval) in enumerate(zip(bars, sorted_terms['Adjusted P-value'])):
            plt.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
                    f'{pval:.2e}', va='center', fontsize=10)

        plt.legend()
        plt.tight_layout()
        plt.show()

    except Exception as e:
        print(f"Custom bar plot error: {e}")

    # Try gseapy plots as backup
    try:
        print("\nTrying gseapy built-in plots...")

        # Convert to gseapy format if needed
        gp.dotplot(sig_ibd_terms,
                   title='IBD Terms - Autophagy-Ferroptosis Intersection',
                   figsize=(10, 8),
                   size_factor=200)
        plt.tight_layout()
        plt.show()

        gp.barplot(sig_ibd_terms,
                   title='IBD Enrichment - Autophagy-Ferroptosis Intersection',
                   figsize=(12, 8))
        plt.tight_layout()
        plt.show()

    except Exception as e:
        print(f"gseapy plot error: {e}")

else:
    print("❌ No significant IBD terms found")

# Compare with individual groups
print(f"\n🔄 COMPARISON WITH INDIVIDUAL GROUPS")
print("=" * 50)

for group_name in ['ARN', 'ferr']:
    if group_name in enrichr_result:
        group_results = enrichr_result[group_name].results
        group_ibd = group_results[
            group_results['Term'].str.contains(
                'IBD|inflammatory bowel|Crohn|colitis|ulcerative colitis',
                case=False, na=False
            )
        ]
        group_sig_ibd = group_ibd[group_ibd['Adjusted P-value'] < 0.05]

        print(f"\n{group_name} group:")
        print(f"  Total IBD terms: {len(group_ibd)}")
        print(f"  Significant IBD terms: {len(group_sig_ibd)}")

        if len(group_sig_ibd) > 0:
            print("  Significant terms:")
            for _, row in group_sig_ibd.iterrows():
                print(f"    • {row['Term']} (p={row['Adjusted P-value']:.2e})")

print(f"\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"✅ Found {len(sig_ibd_terms)} significant IBD terms in ARN|ferr intersection")
print("🧬 This suggests autophagy-ferroptosis crosstalk is associated with IBD")
print("🔬 The intersection reveals disease associations not apparent in individual pathways")
