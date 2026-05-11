#!/usr/bin/env python3
"""
Script: 04_quality_filter.py
Purpose: Apply quality filters to identify high-confidence events
"""

import pandas as pd
import os

def main():
    print("="*70)
    print("STEP 4: QUALITY FILTERING")
    print("="*70)
    
    # Paths
    project_dir = "/home/users/yx1040/Gen811_group_project"
    input_dir = f"{project_dir}/results/tables"
    output_dir = input_dir
    
    input_file = f"{input_dir}/effect_sizes.txt"
    df = pd.read_csv(input_file, sep='\t')
    
    print(f"\nInput: {len(df)} clusters")
    
    print("\n[4.1] Applying quality filters...")
    print("\nCriteria:")
    print("  1. |ΔPSI| ≥ 10%")
    print("  2. WT_reads ≥ 20")
    print("  3. KD_reads ≥ 20")
    print("  4. total_reads ≥ 100")
    
    MIN_DPSI = 10
    MIN_WT = 20
    MIN_KD = 20
    MIN_TOTAL = 100
    
    quality_pass = df[
        (df['abs_deltaPSI'] >= MIN_DPSI) &
        (df['WT_reads'] >= MIN_WT) &
        (df['KD_reads'] >= MIN_KD) &
        (df['total_reads'] >= MIN_TOTAL)
    ].copy()
    
    print(f"\nPassed: {len(quality_pass)} / {len(df)} ({100*len(quality_pass)/len(df):.1f}%)")
    
    # Save
    output_file = f"{output_dir}/quality_filtered.txt"
    quality_pass.to_csv(output_file, sep='\t', index=False, float_format='%.2f')
    print(f"\n✓ Saved: {output_file}")
    
    # Statistics
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    print(f"\nEffect size distribution:")
    for threshold in [50, 40, 30, 20, 15, 10]:
        n = (quality_pass['abs_deltaPSI'] >= threshold).sum()
        pct = 100 * n / len(quality_pass)
        print(f"  |ΔPSI| ≥ {threshold}%: {n:3d} ({pct:5.1f}%)")
    
    print(f"\nRead depth:")
    print(f"  Mean total:   {quality_pass['total_reads'].mean():.0f}")
    print(f"  Median total: {quality_pass['total_reads'].median():.0f}")
    print(f"  Mean WT:      {quality_pass['WT_reads'].mean():.0f}")
    print(f"  Mean KD:      {quality_pass['KD_reads'].mean():.0f}")
    
    print("\n" + "="*70)
    print("STEP 4 COMPLETE")
    print("="*70)

if __name__ == '__main__':
    main()
