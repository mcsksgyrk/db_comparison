import matplotlib.pyplot as plt


def plot_pathway_integration_summary(arn_core, fer_core, crosstalk, shared_targets, shared_regulators_set, edges_df):
    """Create comprehensive integration summary"""

    # Calculate integration metrics
    core_overlap = len(arn_core & fer_core)
    total_unique_cores = len(arn_core | fer_core)
    integration_score = (core_overlap + len(crosstalk) + len(shared_targets) + len(shared_regulators_set)) / 4

    # Create figure
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 16))

    # 1. Core overlap bar chart
    overlap_data = ['ARN Only', 'Shared', 'FER Only']
    overlap_counts = [len(arn_core) - core_overlap, core_overlap, len(fer_core) - core_overlap]
    overlap_colors = ['#4ECDC4', '#96CEB4', '#FF6B6B']

    bars1 = ax1.bar(overlap_data, overlap_counts, color=overlap_colors, alpha=0.8)
    ax1.set_title('Core Protein Distribution', fontweight='bold')
    ax1.set_ylabel('Number of Proteins')

    y_max = max(overlap_counts)
    ax1.set_ylim(0, y_max * 1.3)
    for bar, count in zip(bars1, overlap_counts):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + y_max*0.05,
                str(count), ha='center', fontweight='bold')

    # 2. Interaction types
    interaction_types = ['Direct\nCrosstalk', 'Shared\nTargets', 'Shared\nRegulators']
    interaction_counts = [len(crosstalk), len(shared_targets), len(shared_regulators_set)]
    interaction_colors = ['#FF6B6B', '#4ECDC4', '#45B7D1']

    bars2 = ax2.bar(interaction_types, interaction_counts, color=interaction_colors, alpha=0.8)
    ax2.set_title('Crosstalk Mechanisms', fontweight='bold')
    ax2.set_ylabel('Number of Proteins/Interactions')

    y_max2 = max(interaction_counts)
    ax2.set_ylim(0, y_max2 * 1.3)
    for bar, count in zip(bars2, interaction_counts):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + y_max2*0.05,
                str(count), ha='center', fontweight='bold')

    # 3. NETWORK STATS
    arn_core_set = set(arn_core)
    fer_core_set = set(fer_core)

    # Calculate network statistics
    all_nodes = set(edges_df['interactor_a_node_name']) | set(edges_df['interactor_b_node_name'])
    total_nodes = len(all_nodes)
    total_edges = len(edges_df)

    # ARN network stats
    arn_edges = edges_df[
        (edges_df['interactor_a_node_name'].isin(arn_core_set)) |
        (edges_df['interactor_b_node_name'].isin(arn_core_set))
    ]
    arn_nodes = set(arn_edges['interactor_a_node_name']) | set(arn_edges['interactor_b_node_name'])

    # FER network stats
    fer_edges = edges_df[
        (edges_df['interactor_a_node_name'].isin(fer_core_set)) |
        (edges_df['interactor_b_node_name'].isin(fer_core_set))
    ]
    fer_nodes = set(fer_edges['interactor_a_node_name']) | set(fer_edges['interactor_b_node_name'])

    # Average degree calculation
    node_degrees = {}
    for _, row in edges_df.iterrows():
        node_a, node_b = row['interactor_a_node_name'], row['interactor_b_node_name']
        node_degrees[node_a] = node_degrees.get(node_a, 0) + 1
        node_degrees[node_b] = node_degrees.get(node_b, 0) + 1

    avg_degree_total = sum(node_degrees.values()) / len(node_degrees) if node_degrees else 0

    # Create network stats visualization
    stat_categories = ['Total\nNodes', 'Total\nEdges', 'ARN\nNetwork', 'FER\nNetwork', 'Avg\nDegree']
    stat_values = [
        total_nodes,
        total_edges,
        len(arn_nodes),
        len(fer_nodes),
        int(avg_degree_total)
    ]
    stat_colors = ['#95A5A6', '#95A5A6', '#4ECDC4', '#FF6B6B', '#9B59B6']

    bars3 = ax3.bar(stat_categories, stat_values, color=stat_colors, alpha=0.8)
    ax3.set_title('Network Statistics', fontweight='bold')
    ax3.set_ylabel('Count')

    y_max3 = max(stat_values)
    ax3.set_ylim(0, y_max3 * 1.3)
    for bar, value in zip(bars3, stat_values):
        ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + y_max3*0.05,
                str(value), ha='center', fontweight='bold')

    # 4. SUMMARY TEXT - MAKE EVEN SHORTER TO FIT COMPLETELY
    ax4.axis('off')
    # VERY SHORT TEXT that definitely fits
    summary_text = f"""SUMMARY

Core Proteins:
• Autophagy: {len(arn_core)}
• Ferroptosis: {len(fer_core)}
• Shared: {core_overlap}

Network:
• Nodes: {total_nodes}
• Edges: {total_edges}

Crosstalk:
• Direct: {len(crosstalk)}
• Targets: {len(shared_targets)}"""

    # Position at TOP with more margin from bottom
    ax4.text(0.1, 0.85, summary_text, transform=ax4.transAxes, fontsize=11,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle="round,pad=0.4", facecolor="lightgray", alpha=0.8))

    # MAIN TITLE
    plt.suptitle('Autophagy-Ferroptosis Network Integration Analysis',
                 fontsize=16, fontweight='bold', y=0.94)

    # INCREASE bottom margin to prevent cutoff
    plt.subplots_adjust(left=0.1, bottom=0.15, right=0.9, top=0.85,
                       wspace=0.3, hspace=0.6)
    plt.show()
