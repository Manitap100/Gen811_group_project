#!/usr/bin/env python3
"""
Script: 07_generate_figures.py
Purpose: Generate publication-quality figures
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10

def main():
    print("="*70)
    print("STEP 7: GENERATING FIGURES")
    print("="*70)
    
    # Paths
    project_dir = "/home/users/yx1040/Gen811_group_project"
    input_dir = f"{project_dir}/results/tables"
    output_dir = f"{project_dir}/results/figures"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    quality_df = pd.read_csv(f"{input_dir}/quality_filtered.txt", sep='\t')
    as_df = pd.read_csv(f"{input_dir}/AS_candidates.txt", sep='\t')
    tss_df = pd.read_csv(f"{input_dir}/TSS_TES_candidates.txt", sep='\t')
    
    print(f"\nLoaded:")
    print(f"  Quality filtered: {len(quality_df)}")
    print(f"  AS candidates: {len(as_df)}")
    print(f"  TSS/TES candidates: {len(tss_df)}")
    
    # Figure 1: AS vs TSS/TES ratio
    print("\n[7.1] Generating AS vs TSS/TES ratio plot...")
    fig, ax = plt.subplots(figsize=(6, 5))
    
    counts = [len(as_df), len(tss_df)]
    labels = [f'AS candidates\n({len(as_df)}, 57.4%)', 
              f'TSS/TES candidates\n({len(tss_df)}, 42.6%)']
    colors = ['#3498db', '#e74c3c']
    
    ax.pie(counts, labels=labels, colors=colors, autopct='%1.1f%%',
           startangle=90, textprops={'fontsize': 11})
    ax.set_title('Distribution of Event Types\nAfter Quality Filtering', 
                 fontsize=13, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/AS_vs_TSSTES_ratio.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/AS_vs_TSSTES_ratio.pdf", bbox_inches='tight')
    plt.close()
    print(f"  ✓ AS_vs_TSSTES_ratio.png")
    
    # Figure 2: Effect size distribution
    print("\n[7.2] Generating effect size distribution...")
    fig, ax = plt.subplots(figsize=(8, 5))
    
    bins = np.arange(0, quality_df['abs_deltaPSI'].max()+5, 5)
    ax.hist(quality_df['abs_deltaPSI'], bins=bins, color='#2ecc71', 
            edgecolor='black', alpha=0.7)
    
    ax.axvline(10, color='red', linestyle='--', linewidth=2, 
               label='Quality threshold (10%)')
    ax.axvline(30, color='orange', linestyle='--', linewidth=2,
               label='qPCR threshold (30%)')
    
    ax.set_xlabel('|ΔPSI| (%)', fontsize=12)
    ax.set_ylabel('Number of Events', fontsize=12)
    ax.set_title('Distribution of Effect Sizes\n(High-Quality Events)', 
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/effect_size_distribution.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/effect_size_distribution.pdf", bbox_inches='tight')
    plt.close()
    print(f"  ✓ effect_size_distribution.png")
    
    # Figure 3: Top candidates heatmap
    print("\n[7.3] Generating top candidates heatmap...")
    top20 = as_df.nlargest(20, 'abs_deltaPSI')
    
    # Prepare data for heatmap
    heatmap_data = top20[['cluster', 'WT_reads', 'KD_reads', 'deltaPSI']].copy()
    heatmap_data['cluster_short'] = [c.split(':')[-1][:15] for c in heatmap_data['cluster']]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 8), 
                                     gridspec_kw={'width_ratios': [3, 1]})
    
    # Heatmap of read counts
    plot_data = heatmap_data[['WT_reads', 'KD_reads']].values
    im = ax1.imshow(plot_data, aspect='auto', cmap='YlOrRd')
    
    ax1.set_xticks([0, 1])
    ax1.set_xticklabels(['WT', 'KD'])
    ax1.set_yticks(range(len(heatmap_data)))
    ax1.set_yticklabels(heatmap_data['cluster_short'], fontsize=8)
    ax1.set_title('Read Counts', fontsize=11, fontweight='bold')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax1)
    cbar.set_label('Reads', rotation=270, labelpad=15)
    
    # Bar plot of ΔPSI
    colors_bar = ['#e74c3c' if x < 0 else '#3498db' for x in heatmap_data['deltaPSI']]
    ax2.barh(range(len(heatmap_data)), heatmap_data['deltaPSI'], color=colors_bar)
    ax2.set_yticks(range(len(heatmap_data)))
    ax2.set_yticklabels([])
    ax2.set_xlabel('ΔPSI (%)')
    ax2.set_title('Effect Size', fontsize=11, fontweight='bold')
    ax2.axvline(0, color='black', linestyle='-', linewidth=0.5)
    ax2.grid(axis='x', alpha=0.3)
    
    plt.suptitle('Top 20 AS Candidates', fontsize=13, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/top_candidates_heatmap.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/top_candidates_heatmap.pdf", bbox_inches='tight')
    plt.close()
    print(f"  ✓ top_candidates_heatmap.png")
    
    # Figure 4: Read distribution
    print("\n[7.4] Generating read distribution plot...")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    
    # WT reads
    axes[0].hist(np.log10(quality_df['WT_reads']+1), bins=30, 
                 color='#3498db', edgecolor='black', alpha=0.7)
    axes[0].set_xlabel('log10(WT reads + 1)', fontsize=11)
    axes[0].set_ylabel('Frequency', fontsize=11)
    axes[0].set_title('WT Read Distribution', fontsize=12, fontweight='bold')
    axes[0].grid(axis='y', alpha=0.3)
    
    # KD reads
    axes[1].hist(np.log10(quality_df['KD_reads']+1), bins=30,
                 color='#e74c3c', edgecolor='black', alpha=0.7)
    axes[1].set_xlabel('log10(KD reads + 1)', fontsize=11)
    axes[1].set_ylabel('Frequency', fontsize=11)
    axes[1].set_title('KD Read Distribution', fontsize=12, fontweight='bold')
    axes[1].grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/reads_distribution.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{output_dir}/reads_distribution.pdf", bbox_inches='tight')
    plt.close()
    print(f"  ✓ reads_distribution.png")
    
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"\nGenerated figures:")
    print(f"  1. AS_vs_TSSTES_ratio.png")
    print(f"  2. effect_size_distribution.png")
    print(f"  3. top_candidates_heatmap.png")
    print(f"  4. reads_distribution.png")
    print(f"\nOutput directory: {output_dir}")
    
    print("\n" + "="*70)
    print("STEP 7 COMPLETE")
    print("="*70)

if __name__ == '__main__':
    main()
