#!/usr/bin/env python3
"""
Script: 05_validate_splice_sites.py
Purpose: Check for canonical splice sites and classify events
"""

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import os

def main():
    print("="*70)
    print("STEP 5: SPLICE SITE VALIDATION")
    print("="*70)
    
    # Paths
    project_dir = "/home/users/yx1040/Gen811_group_project"
    input_dir = f"{project_dir}/results/tables"
    output_dir = input_dir
    genome_file = f"{project_dir}/data/reference/ToxoDB-68_TgondiiME49_Genome.fasta"
    
    input_file = f"{input_dir}/quality_filtered.txt"
    df = pd.read_csv(input_file, sep='\t')
    
    print(f"\nInput: {len(df)} quality-filtered clusters")
    
    print("\n[5.1] Loading genome...")
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    print(f"  ✓ Loaded {len(genome)} chromosomes")
    
    print("\n[5.2] Checking splice sites...")
    has_splice = []
    no_splice = []
    
    for idx, row in df.iterrows():
        chrom = row['chrom']
        start = row['start']
        end = row['end']
        strand = row['strand']
        
        if chrom not in genome:
            no_splice.append(idx)
            continue
        
        try:
            if strand == '+':
                acceptor = str(genome[chrom].seq[start-2:start]).upper()
                donor = str(genome[chrom].seq[end:end+2]).upper()
            else:
                acceptor_raw = str(genome[chrom].seq[start-2:start])
                donor_raw = str(genome[chrom].seq[end-2:end])
                acceptor = str(Seq(acceptor_raw).reverse_complement()).upper()
                donor = str(Seq(donor_raw).reverse_complement()).upper()
            
            if (acceptor == "AG") and (donor in ["GT", "GC"]):
                has_splice.append(idx)
            else:
                no_splice.append(idx)
        except:
            no_splice.append(idx)
    
    with_splice = df.loc[has_splice].copy()
    without_splice = df.loc[no_splice].copy()
    
    print(f"\n" + "="*70)
    print("RESULTS")
    print("="*70)
    print(f"\nWith canonical splice sites: {len(with_splice)} ({100*len(with_splice)/len(df):.1f}%)")
    print(f"  → AS candidates")
    print(f"\nWithout canonical splice sites: {len(without_splice)} ({100*len(without_splice)/len(df):.1f}%)")
    print(f"  → TSS/TES candidates")
    
    # Save
    with_splice.to_csv(f"{output_dir}/AS_candidates.txt", sep='\t', index=False, float_format='%.2f')
    without_splice.to_csv(f"{output_dir}/TSS_TES_candidates.txt", sep='\t', index=False, float_format='%.2f')
    
    print(f"\n✓ Saved:")
    print(f"  AS_candidates.txt ({len(with_splice)} events)")
    print(f"  TSS_TES_candidates.txt ({len(without_splice)} events)")
    
    print(f"\nTop 10 AS candidates:")
    print("-"*70)
    top10 = with_splice.nlargest(10, 'abs_deltaPSI')
    print(top10[['cluster', 'deltaPSI', 'total_reads', 'WT_reads', 'KD_reads']].to_string(index=False))
    
    print("\n" + "="*70)
    print("STEP 5 COMPLETE")
    print("="*70)

if __name__ == '__main__':
    main()
