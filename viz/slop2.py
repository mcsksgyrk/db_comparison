import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import json
from matplotlib.patches import FancyBboxPatch
import pandas as pd
from apicalls.api_oop import UniProtClient
from config import *


def get_gene_name(id):
    return uniprot_to_gene[id]

# A0 poster configuration
plt.rcParams.update({
    'font.size': 18,
    'axes.titlesize': 26,
    'axes.labelsize': 22,
    'xtick.labelsize': 18,
    'ytick.labelsize': 18,
    'legend.fontsize': 20,
    'figure.figsize': [20, 14]
})

json_path = OUTPUT_DIR / "upstream_regulators.json"
with open(json_path, 'rb') as f:
    upstream_data = json.load(f)

uniprot_ids = set()
for reg, pw in upstream_data.items():
    uniprot_ids.add(reg)
    for p, id in pw.items():
        uniprot_ids.update(id)

uniprot = UniProtClient()
uniprot_to_gene, fails = uniprot.batch_convert_from_uniprot_id("GeneCards", list(uniprot_ids), 1)

ferroptosis_targets = set()
for reg_data in upstream_data.values():
    ferroptosis_targets.update(reg_data['fer'])

ferroptosis_targets = sorted(list(ferroptosis_targets))
regulators = sorted(list(upstream_data.keys()))

# Create regulation matrix
regulation_matrix = []
for reg in regulators:
    row = []
    for target in ferroptosis_targets:
        if target in upstream_data[reg]['fer']:
            row.append(1)  # Regulated
        else:
            row.append(0)  # Not regulated
    regulation_matrix.append(row)

regulation_df = pd.DataFrame(regulation_matrix,
                             index=[get_gene_name(reg) for reg in regulators],
                             columns=[get_gene_name(target) for target in ferroptosis_targets])

# Create the heatmap
fig3, ax3 = plt.subplots(figsize=(16, 12))

# Custom colormap: white for no regulation, red for regulation
import matplotlib.colors as mcolors
colors = ['white', '#D32F2F']  # White -> Dark Red
n_bins = 2
cmap = mcolors.ListedColormap(colors)

# Create heatmap
im = ax3.imshow(regulation_matrix, cmap=cmap, aspect='auto', alpha=0.8)

ax3.set_xticks(range(len(ferroptosis_targets)))
ax3.set_yticks(range(len(regulators)))
ax3.set_xticklabels([get_gene_name(target) for target in ferroptosis_targets],
                   rotation=45, ha='right', fontsize=14)
ax3.set_yticklabels([get_gene_name(reg) for reg in regulators], fontsize=14)

# Add grid
ax3.set_xticks(np.arange(len(ferroptosis_targets))+0.5, minor=True)
ax3.set_yticks(np.arange(len(regulators))+0.5, minor=True)
ax3.grid(which='minor', color='gray', linestyle='-', linewidth=1, alpha=0.3)

# Add text annotations
for i in range(len(regulators)):
    for j in range(len(ferroptosis_targets)):
        if regulation_matrix[i][j] == 1:
            ax3.text(j, i, '●', ha='center', va='center',
                    fontsize=20, color='white', fontweight='bold')

ax3.set_xlabel('Ferroptosis Target Genes', fontsize=18, fontweight='bold')
ax3.set_ylabel('Upstream Regulators', fontsize=18, fontweight='bold')
ax3.set_title('Ferroptosis Regulation Matrix',
              fontsize=22, fontweight='bold')
ax3.title.set_position([0.5, 1.08])  # Position title higher

# Add summary statistics
total_regulations = sum(sum(row) for row in regulation_matrix)
ax3.text(len(ferroptosis_targets), len(regulators)-1,
         f'Total regulations: {total_regulations}\n' +
         f'Target genes: {len(ferroptosis_targets)}\n' +
         f'Regulators: {len(regulators)}',
         fontsize=12, bbox=dict(boxstyle="round,pad=0.5", facecolor='lightblue', alpha=0.8))

plt.tight_layout()

plt.savefig('ferroptosis_regulation_matrix.pdf', bbox_inches='tight',
            facecolor='white', pad_inches=0.3)
plt.show()
