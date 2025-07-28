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

# Search for IBD terms
ibd_terms = intersection_results[
    intersection_results['Term'].str.contains(IBD_SEARCH_PATTERN, case=False, na=False)
]

# Filter for significant IBD terms
sig_ibd_terms = ibd_terms[ibd_terms['Adjusted P-value'] < 0.05].copy()

if len(sig_ibd_terms) > 0:
    print(f"Found {len(sig_ibd_terms)} significant IBD terms")

    # Extract additional metrics for better visualization
    gene_counts = []
    gene_ratios = []
    fold_enrichments = []

    for _, row in sig_ibd_terms.iterrows():
        # Parse overlap (e.g., "12/150" -> 12 genes found, 150 total in pathway)
        overlap_parts = row['Overlap'].split('/')
        genes_found = int(overlap_parts[0])
        genes_in_pathway = int(overlap_parts[1])

        gene_counts.append(genes_found)
        gene_ratios.append(genes_found / genes_in_pathway)

        # Calculate fold enrichment (if available in your data)
        # If not available, we can approximate it
        if 'Odds Ratio' in row:
            fold_enrichments.append(row['Odds Ratio'])
        elif 'Combined Score' in row:
            # Combined score is often related to fold enrichment
            fold_enrichments.append(row['Combined Score'] / 100)  # Normalize
        else:
            # Approximate fold enrichment from p-value and overlap
            # This is a rough approximation - replace with actual values if available
            fold_enrichment = genes_found * (-np.log10(row['Adjusted P-value'])) / 10
            fold_enrichments.append(max(1.1, fold_enrichment))

    # Add calculated metrics to dataframe
    sig_ibd_terms['Gene_Count'] = gene_counts
    sig_ibd_terms['Gene_Ratio'] = gene_ratios
    sig_ibd_terms['Fold_Enrichment'] = fold_enrichments

    # Sort by p-value for better display
    sig_ibd_terms = sig_ibd_terms.sort_values('Adjusted P-value')

    # ENHANCED DOTPLOT - VERSION 1: Use Fold Enrichment for Color
    plt.figure(figsize=(14, 10))

    y_pos = np.arange(len(sig_ibd_terms))

    # NOW USING DIFFERENT METRICS FOR EACH VISUAL ELEMENT!
    scatter = plt.scatter(
        -np.log10(sig_ibd_terms['Adjusted P-value']),  # X-axis: Statistical significance
        y_pos,                                         # Y-axis: Terms
        s=sig_ibd_terms['Gene_Count'] * 30,           # Size: Number of genes (practical significance)
        c=sig_ibd_terms['Fold_Enrichment'],           # Color: Fold enrichment (effect size) - DIFFERENT!
        cmap='YlOrRd',                                # Yellow to Red colormap
        alpha=0.8,
        edgecolors='black',
        linewidth=0.8
    )

    # Customize plot
    plt.yticks(y_pos, sig_ibd_terms['Term'], fontsize=11)
    plt.xlabel('-log10(Adjusted P-value)', fontsize=14, fontweight='bold')
    plt.ylabel('IBD Terms', fontsize=14, fontweight='bold')
    plt.title('IBD Terms - Autophagy-Ferroptosis Intersection\n' +
              '(Color = Fold Enrichment, Size = Gene Count)',
              fontsize=16, fontweight='bold')

    # Add significance threshold line
    plt.axvline(x=-np.log10(0.05), color='gray', linestyle='--', alpha=0.6, linewidth=1)
    plt.text(-np.log10(0.05) + 0.1, len(sig_ibd_terms)-1, 'p=0.05', rotation=90,
             va='top', ha='left', fontsize=9, color='gray')

    # Add colorbar for fold enrichment (DIFFERENT from x-axis!)
    cbar = plt.colorbar(scatter)
    cbar.set_label('Fold Enrichment', fontsize=12, fontweight='bold')

    # Add size legend for gene counts
    max_genes = max(gene_counts)
    size_values = [3, 8, max_genes] if max_genes > 8 else [1, 3, max_genes]
    size_labels = [f'{val} genes' for val in size_values]

    size_legend_elements = []
    for size_val in size_values:
        size_legend_elements.append(
            plt.scatter([], [], s=size_val*30, c='gray', alpha=0.6,
                       edgecolors='black', linewidth=0.8)
        )

    size_legend = plt.legend(size_legend_elements, size_labels,
                           title='Gene Count', loc='upper right',
                           title_fontsize=11, fontsize=10, frameon=True)
    size_legend.get_frame().set_facecolor('white')
    size_legend.get_frame().set_alpha(0.8)

    # Add grid and improve layout
    plt.grid(True, alpha=0.3, axis='x')
    plt.tight_layout()
    plt.savefig('enhanced_ibd_dotplot_v1.png', dpi=300, bbox_inches='tight')
    plt.show()

    # ALTERNATIVE VERSION 2: Use Gene Ratio for Color
    plt.figure(figsize=(14, 10))

    scatter2 = plt.scatter(
        -np.log10(sig_ibd_terms['Adjusted P-value']),  # X-axis: Statistical significance
        y_pos,                                         # Y-axis: Terms
        s=sig_ibd_terms['Gene_Count'] * 30,           # Size: Number of genes
        c=sig_ibd_terms['Gene_Ratio'],                # Color: Gene ratio (pathway coverage) - DIFFERENT!
        cmap='viridis',                               # Different colormap
        alpha=0.8,
        edgecolors='black',
        linewidth=0.8
    )

    plt.yticks(y_pos, sig_ibd_terms['Term'], fontsize=11)
    plt.xlabel('-log10(Adjusted P-value)', fontsize=14, fontweight='bold')
    plt.ylabel('IBD Terms', fontsize=14, fontweight='bold')
    plt.title('IBD Terms - Autophagy-Ferroptosis Intersection\n' +
              '(Color = Gene Ratio, Size = Gene Count)',
              fontsize=16, fontweight='bold')

    # Add significance threshold
    plt.axvline(x=-np.log10(0.05), color='gray', linestyle='--', alpha=0.6, linewidth=1)
    plt.text(-np.log10(0.05) + 0.1, len(sig_ibd_terms)-1, 'p=0.05', rotation=90,
             va='top', ha='left', fontsize=9, color='gray')

    # Colorbar for gene ratio
    cbar2 = plt.colorbar(scatter2)
    cbar2.set_label('Gene Ratio\n(Found/Total in Pathway)', fontsize=12, fontweight='bold')

    # Size legend
    size_legend2 = plt.legend(size_legend_elements, size_labels,
                             title='Gene Count', loc='upper right',
                             title_fontsize=11, fontsize=10, frameon=True)

    plt.grid(True, alpha=0.3, axis='x')
    plt.tight_layout()
    plt.savefig('enhanced_ibd_dotplot_v2.png', dpi=300, bbox_inches='tight')
    plt.show()

    # PRINT SUMMARY OF IMPROVEMENTS
    print("\n" + "="*60)
    print("DOTPLOT ENHANCEMENT SUMMARY")
    print("="*60)
    print("ORIGINAL PROBLEM:")
    print("- X-axis: -log10(Adjusted P-value)")
    print("- Color:  -log10(Adjusted P-value)  <- SAME AS X-axis (redundant!)")
    print("- Size:   Gene count")
    print()
    print("ENHANCED VERSION 1:")
    print("- X-axis: -log10(Adjusted P-value) (statistical significance)")
    print("- Color:  Fold Enrichment (effect size)")
    print("- Size:   Gene Count (clearly labeled)")
    print()
    print("ENHANCED VERSION 2:")
    print("- X-axis: -log10(Adjusted P-value) (statistical significance)")
    print("- Color:  Gene Ratio (pathway coverage)")
    print("- Size:   Gene Count (clearly labeled)")
    print()
    print("Now each visual element provides DIFFERENT information!")
    print("="*60)

    # Print detailed results
    print(f"\nDETAILED RESULTS:")
    print(f"{'Term':<50} {'P-value':<12} {'Genes':<8} {'Ratio':<8} {'Fold Enrich':<12}")
    print("-" * 100)

    for _, row in sig_ibd_terms.iterrows():
        term_short = row['Term'][:47] + "..." if len(row['Term']) > 50 else row['Term']
        print(f"{term_short:<50} {row['Adjusted P-value']:<12.2e} {row['Gene_Count']:<8} {row['Gene_Ratio']:<8.3f} {row['Fold_Enrichment']:<12.2f}")

else:
    print("No significant IBD terms found")

# Original comparison code continues...
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

# Create comparison figure
print("\nCreating group comparison figure...")
plt.figure(figsize=(10, 6))

groups = list(group_ibd_results.keys())
total_counts = [group_ibd_results[group]['total'] for group in groups]
sig_counts = [group_ibd_results[group]['significant'] for group in groups]
non_sig_counts = [total - sig for total, sig in zip(total_counts, sig_counts)]

x = np.arange(len(groups))
width = 0.6

bars1 = plt.bar(x, non_sig_counts, width, label='Non-significant IBD terms', alpha=0.7, color='lightgray')
bars2 = plt.bar(x, sig_counts, width, bottom=non_sig_counts, label='Significant IBD terms', alpha=0.7, color='darkblue')

plt.xlabel('Groups')
plt.ylabel('Number of IBD terms')
plt.title('IBD Terms: Significant vs Non-significant')
plt.xticks(x, groups)
plt.legend()

for i, (total, sig, non_sig) in enumerate(zip(total_counts, sig_counts, non_sig_counts)):
    if sig > 0:
        plt.text(i, non_sig + sig/2, str(sig), ha='center', va='center', fontweight='bold', color='white')
    if non_sig > 0:
        plt.text(i, non_sig/2, str(non_sig), ha='center', va='center', fontweight='bold', color='black')
    plt.text(i, total + 0.2, f'Total: {total}', ha='center', va='bottom', fontsize=10)

plt.tight_layout()
plt.savefig('ibd_group_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"\nSummary: Found {len(sig_ibd_terms) if len(sig_ibd_terms) > 0 else 0} significant IBD terms in ARN|ferr intersection")
