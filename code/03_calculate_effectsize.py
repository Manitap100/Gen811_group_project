#!/usr/bin/env python3
"""
Script: 03_calculate_effectsize.py
Purpose: Calculate PSI and ΔPSI for each splicing event
"""

import pandas as pd
import gzip
import os
import sys

def main():
    print("="*70)
    print("STEP 3: EFFECT SIZE CALCULATION")
    print("="*70)
    
    # Paths
    project_dir = "/home/users/yx1040/Gen811_group_project"
    input_dir = f"{project_dir}/data/leafcutter"
    output_dir = f"{project_dir}/results/tables"
    os.makedirs(output_dir, exist_ok=True)
    
    counts_file = f"{input_dir}/leafcutter_perind.counts.gz"
    diff_file = f"{input_dir}/diff_splicing_results.txt"
    
    # Check files
    if not os.path.exists(counts_file):
        print(f"ERROR: {counts_file} not found")
        sys.exit(1)
    if not os.path.exists(diff_file):
        print(f"ERROR: {diff_file} not found")
        sys.exit(1)
    
    print("\n[3.1] Loading differential results...")
    diff_results = pd.read_csv(diff_file, sep='\s+')
    
    # Find FDR column
    fdr_col = None
    for col in ['p.adjust', 'FDR', 'fdr']:
        if col in diff_results.columns:
            fdr_col = col
            break
    
    if not fdr_col:
        print("ERROR: Cannot find FDR column")
        print(f"Columns: {diff_results.columns.tolist()}")
        sys.exit(1)
    
    sig_clusters = diff_results[diff_results[fdr_col] < 0.05].copy()
    print(f"  Total testable: {len(diff_results)}")
    print(f"  Significant (FDR<0.05): {len(sig_clusters)}")
    
    print("\n[3.2] Loading junction counts...")
    counts = []
    with gzip.open(counts_file, 'rt') as f:
        header = f.readline().strip().split()
        samples = header[1:]
        print(f"  Samples: {samples}")
        
        for line in f:
            parts = line.strip().split()
            junction_id = parts[0]
            values = [int(x) for x in parts[1:]]
            
            try:
                junc_parts = junction_id.split(':')
                chrom = junc_parts[0]
                start = int(junc_parts[1])
                end = int(junc_parts[2])
                cluster_info = junc_parts[3]
                
                cluster_parts = cluster_info.rsplit('_', 1)
                cluster_id = cluster_parts[0]
                strand = cluster_parts[1] if len(cluster_parts) > 1 else '+'
                cluster = f"{chrom}:{cluster_id}"
                
                counts.append({
                    'junction_full': junction_id,
                    'cluster': cluster,
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    **{samples[i]: values[i] for i in range(len(samples))}
                })
            except:
                continue
    
    counts_df = pd.DataFrame(counts)
    print(f"  Junctions: {len(counts_df)}")
    print(f"  Clusters: {counts_df['cluster'].nunique()}")
    
    # Define groups
    wt_samples = [s for s in samples if 'WT' in s]
    kd_samples = [s for s in samples if 'KD' in s and s not in ['KD3', 'KD4']]
    print(f"  WT samples: {wt_samples}")
    print(f"  KD samples: {kd_samples}")
    
    print("\n[3.3] Calculating PSI and ΔPSI...")
    results = []
    
    for cluster_name in sig_clusters['cluster']:
        cluster_juncs = counts_df[counts_df['cluster'] == cluster_name]
        if len(cluster_juncs) == 0:
            continue
        
        cluster_wt_sum = cluster_juncs[wt_samples].sum().sum()
        cluster_kd_sum = cluster_juncs[kd_samples].sum().sum()
        
        if cluster_wt_sum == 0 or cluster_kd_sum == 0:
            continue
        
        for idx, junc in cluster_juncs.iterrows():
            wt_reads = junc[wt_samples].sum()
            kd_reads = junc[kd_samples].sum()
            
            psi_wt = 100 * wt_reads / cluster_wt_sum
            psi_kd = 100 * kd_reads / cluster_kd_sum
            delta_psi = psi_kd - psi_wt
            
            results.append({
                'cluster': cluster_name,
                'chrom': junc['chrom'],
                'start': junc['start'],
                'end': junc['end'],
                'strand': junc['strand'],
                'WT_reads': int(wt_reads),
                'KD_reads': int(kd_reads),
                'total_reads': int(cluster_wt_sum + cluster_kd_sum),
                'PSI_WT': psi_wt,
                'PSI_KD': psi_kd,
                'deltaPSI': delta_psi,
                'abs_deltaPSI': abs(delta_psi)
            })
    
    all_junctions = pd.DataFrame(results)
    
    print("\n[3.4] Selecting representative junction per cluster...")
    max_per_cluster = all_junctions.loc[
        all_junctions.groupby('cluster')['abs_deltaPSI'].idxmax()
    ]
    
    # Merge with p-values
    p_col = 'p' if 'p' in sig_clusters.columns else 'p.value'
    if p_col in sig_clusters.columns:
        max_per_cluster = max_per_cluster.merge(
            sig_clusters[['cluster', p_col, fdr_col]],
            on='cluster'
        )
        max_per_cluster.rename(columns={p_col: 'p_value', fdr_col: 'fdr'}, inplace=True)
    
    print(f"  Clusters: {len(max_per_cluster)}")
    
    # Save
    output_file = f"{output_dir}/effect_sizes.txt"
    max_per_cluster.to_csv(output_file, sep='\t', index=False, float_format='%.3f')
    print(f"\n✓ Saved: {output_file}")
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Clusters analyzed: {len(max_per_cluster)}")
    print(f"\n|ΔPSI| distribution:")
    print(f"  Mean:   {max_per_cluster['abs_deltaPSI'].mean():.2f}%")
    print(f"  Median: {max_per_cluster['abs_deltaPSI'].median():.2f}%")
    print(f"  Range:  {max_per_cluster['abs_deltaPSI'].min():.2f}% - {max_per_cluster['abs_deltaPSI'].max():.2f}%")
    
    print(f"\nTop 10 by |ΔPSI|:")
    print("-"*70)
    top10 = max_per_cluster.nlargest(10, 'abs_deltaPSI')
    print(top10[['cluster', 'deltaPSI', 'total_reads', 'WT_reads', 'KD_reads']].to_string(index=False))
    
    print("\n" + "="*70)
    print("STEP 3 COMPLETE")
    print("="*70)

if __name__ == '__main__':
    main()
