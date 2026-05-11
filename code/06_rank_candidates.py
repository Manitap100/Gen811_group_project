#!/usr/bin/env python3
"""
Script: 06_rank_candidates.py
Purpose: Rank AS candidates and identify validation targets
"""

import pandas as pd
import os

def main():
    print("="*70)
    print("STEP 6: RANKING AS CANDIDATES")
    print("="*70)
    
    # Paths
    project_dir = "/home/users/yx1040/Gen811_group_project"
    input_dir = f"{project_dir}/results/tables"
    output_dir = input_dir
    
    input_file = f"{input_dir}/AS_candidates.txt"
    df = pd.read_csv(input_file, sep='\t')
    
    print(f"\nTotal AS candidates: {len(df)}")
    
    # Sort by effect size
    df_sorted = df.sort_values('abs_deltaPSI', ascending=False)
    
    print("\n[6.1] Top 20 AS candidates by |ΔPSI|:")
    print("="*70)
    top20 = df_sorted.head(20)
    print(top20[['cluster', 'deltaPSI', 'total_reads', 'WT_reads', 'KD_reads']].to_string(index=False))
    
    print("\n[6.2] Candidates by effect size:")
    print("-"*70)
    for threshold in [50, 40, 30, 20, 15]:
        n = (df['abs_deltaPSI'] >= threshold).sum()
        print(f"  |ΔPSI| ≥ {threshold}%: {n:2d} candidates")
    
    print("\n[6.3] qPCR-ready candidates:")
    print("-"*70)
    qpcr_ready = df[
        (df['total_reads'] >= 200) &
        (df['WT_reads'] >= 10) &
        (df['KD_reads'] >= 10) &
        (df['abs_deltaPSI'] >= 20)
    ].sort_values('abs_deltaPSI', ascending=False)
    
    print(f"\nCriteria: total≥200, WT≥10, KD≥10, |ΔPSI|≥20%")
    print(f"Qualified: {len(qpcr_ready)} candidates")
    
    if len(qpcr_ready) > 0:
        print("\nqPCR-ready candidates:")
        print(qpcr_ready[['cluster', 'deltaPSI', 'total_reads', 'WT_reads', 'KD_reads']].to_string(index=False))
    
    # Save
    df_sorted.to_csv(f"{output_dir}/ranked_AS_candidates.txt", sep='\t', index=False, float_format='%.2f')
    if len(qpcr_ready) > 0:
        qpcr_ready.to_csv(f"{output_dir}/qPCR_candidates.txt", sep='\t', index=False, float_format='%.2f')
    
    print(f"\n✓ Saved:")
    print(f"  ranked_AS_candidates.txt")
    if len(qpcr_ready) > 0:
        print(f"  qPCR_candidates.txt ({len(qpcr_ready)} candidates)")
    
    print("\n" + "="*70)
    print("STEP 6 COMPLETE")
    print("="*70)

if __name__ == '__main__':
    main()
