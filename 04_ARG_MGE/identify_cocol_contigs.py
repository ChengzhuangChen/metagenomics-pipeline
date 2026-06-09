#!/usr/bin/env python3
# ============================================================
# Identify ARG-MGE Co-localized Contigs (10kb window)
# For Nature article supplementary materials
#
# This script identifies contigs where ARGs and MGEs are located
# within 10 kb of each other, indicating potential horizontal
# gene transfer (HGT).
#
# ARGs and MGEs within 10 kb on the same contig were regarded
# as co-occurring.
# ============================================================

import argparse
import os
import glob
from collections import defaultdict


def load_orf_positions_from_gff(gff_dir):
    """
    Load ORF positions from GFF files.

    Args:
        gff_dir: Directory with Prodigal GFF files

    Returns:
        Dict: {contig_id: [(orf_id, start, end), ...]}
    """
    orf_positions = defaultdict(list)

    for gff_file in glob.glob(os.path.join(gff_dir, "*.gff")):
        sample_name = os.path.basename(gff_file).replace('.gff', '')

        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue

                # Skip header lines
                if fields[0].startswith('#'):
                    continue

                contig_id = fields[0]
                start = int(fields[3])
                end = int(fields[4])

                # Parse GFF attributes for ORF ID
                attributes = {}
                for attr in fields[8].split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attributes[key] = value

                # Get ID attribute (could be 'ID' or 'locus_tag')
                orf_id = attributes.get('ID', '')
                if not orf_id:
                    orf_id = attributes.get('locus_tag', '')

                if orf_id:
                    # Add sample prefix if not present
                    if not orf_id.startswith(sample_name):
                        orf_id = f"{sample_name}_{orf_id}"
                    if not contig_id.startswith(sample_name):
                        contig_id = f"{sample_name}_{contig_id}"

                    orf_positions[contig_id].append({
                        'orf_id': orf_id,
                        'start': start,
                        'end': end,
                        'midpoint': (start + end) // 2
                    })

    print(f"Loaded positions for {len(orf_positions)} contigs")
    return orf_positions


def load_arg_orfs_from_csv(csv_file):
    """
    Load ARG ORF IDs from DIAMOND results.

    Args:
        csv_file: Path to merged ARG DIAMOND results

    Returns:
        Set of ARG ORF IDs
    """
    arg_orfs = set()

    with open(csv_file, 'r') as f:
        header = f.readline()  # Skip header

        for line in f:
            if not line.strip():
                continue

            fields = line.strip().split('\t')
            if len(fields) > 0:
                orf_id = fields[0]
                arg_orfs.add(orf_id)

    return arg_orfs


def load_mge_orfs_from_csv(csv_file):
    """
    Load MGE ORF IDs from DIAMOND results.

    Args:
        csv_file: Path to merged MGE DIAMOND results

    Returns:
        Set of MGE ORF IDs
    """
    mge_orfs = set()

    with open(csv_file, 'r') as f:
        header = f.readline()  # Skip header

        for line in f:
            if not line.strip():
                continue

            fields = line.strip().split('\t')
            if len(fields) > 0:
                orf_id = fields[0]
                mge_orfs.add(orf_id)

    return mge_orfs


def find_cocol_contigs_with_distance(arg_csv, mge_csv, gff_dir, output_file, stats_file, distance_threshold=10000):
    """
    Find contigs where ARGs and MGEs are within specified distance.

    Args:
        arg_csv: Path to merged ARG DIAMOND results
        mge_csv: Path to merged MGE DIAMOND results
        gff_dir: Directory with GFF files for position info
        output_file: Output file for co-localized results
        stats_file: Output file for statistics
        distance_threshold: Distance in bp (default: 10,000 = 10kb)
    """
    print(f"Loading ARG results from: {arg_csv}")
    arg_orfs = load_arg_orfs_from_csv(arg_csv)
    print(f"  Found {len(arg_orfs)} ARG ORFs")

    print(f"Loading MGE results from: {mge_csv}")
    mge_orfs = load_mge_orfs_from_csv(mge_csv)
    print(f"  Found {len(mge_orfs)} MGE ORFs")

    print(f"Loading ORF positions from GFF files...")
    orf_positions = load_orf_positions_from_gff(gff_dir)

    # Build ORF lookup
    orf_to_contig = {}
    for contig, orfs in orf_positions.items():
        for orf in orfs:
            orf_to_contig[orf['orf_id']] = {
                'contig': contig,
                'start': orf['start'],
                'end': orf['end'],
                'midpoint': orf['midpoint']
            }

    # Find ARG contigs and MGE contigs
    arg_contigs = set()
    mge_contigs = set()

    for orf_id in arg_orfs:
        if orf_id in orf_to_contig:
            arg_contigs.add(orf_to_contig[orf_id]['contig'])

    for orf_id in mge_orfs:
        if orf_id in orf_to_contig:
            mge_contigs.add(orf_to_contig[orf_id]['contig'])

    print(f"\nContigs with ARGs: {len(arg_contigs)}")
    print(f"Contigs with MGEs: {len(mge_contigs)}")

    # Find co-localized (within 10kb) ARGs and MGEs
    # Strategy: For each contig with both ARGs and MGEs,
    # check if any ARG-MGE pair is within distance_threshold

    cocol_results = []  # [(contig, arg_orf, mge_orf, distance), ...]
    cocol_contigs = set()

    # Get contigs with both ARGs and MGEs
    shared_contigs = arg_contigs & mge_contigs

    print(f"\nAnalyzing {len(shared_contigs)} contigs with both ARGs and MGEs...")

    for contig in shared_contigs:
        # Get all ORFs on this contig
        contig_orfs = orf_positions.get(contig, [])

        # Separate into ARG and MGE ORFs
        arg_orfs_on_contig = [o for o in contig_orfs if o['orf_id'] in arg_orfs]
        mge_orfs_on_contig = [o for o in contig_orfs if o['orf_id'] in mge_orfs]

        # Check each ARG-MGE pair
        for arg_orf in arg_orfs_on_contig:
            for mge_orf in mge_orfs_on_contig:
                # Calculate distance between midpoints
                distance = abs(arg_orf['midpoint'] - mge_orf['midpoint'])

                if distance <= distance_threshold:
                    cocol_results.append({
                        'contig': contig,
                        'arg_orf': arg_orf['orf_id'],
                        'arg_start': arg_orf['start'],
                        'arg_end': arg_orf['end'],
                        'mge_orf': mge_orf['orf_id'],
                        'mge_start': mge_orf['start'],
                        'mge_end': mge_orf['end'],
                        'distance': distance
                    })
                    cocol_contigs.add(contig)

    # Write results
    print(f"\nWriting results to: {output_file}")
    with open(output_file, 'w') as f:
        f.write("Contig\tARG_ORF\tARG_Start\tARG_End\tMGE_ORF\tMGE_Start\tMGE_End\tDistance_bp\tDistance_kb\n")
        for r in sorted(cocol_results, key=lambda x: (x['contig'], x['distance'])):
            f.write(f"{r['contig']}\t{r['arg_orf']}\t{r['arg_start']}\t{r['arg_end']}\t")
            f.write(f"{r['mge_orf']}\t{r['mge_start']}\t{r['mge_end']}\t")
            f.write(f"{r['distance']}\t{r['distance']/1000:.2f}\n")

    # Write statistics
    total_arg = len(arg_orfs)
    total_mge = len(mge_orfs)
    total_cocol = len(cocol_contigs)
    total_pairs = len(cocol_results)

    # Calculate distances
    distances = [r['distance'] for r in cocol_results]
    avg_distance = sum(distances) / len(distances) if distances else 0
    median_distance = sorted(distances)[len(distances)//2] if distances else 0

    cocol_rate_arg = (total_cocol / len(arg_contigs) * 100) if arg_contigs else 0
    cocol_rate_mge = (total_cocol / len(mge_contigs) * 100) if mge_contigs else 0

    with open(stats_file, 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("ARG-MGE Co-localization Analysis (10 kb window)\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Distance threshold: {distance_threshold} bp ({distance_threshold/1000} kb)\n\n")

        f.write("Basic Statistics:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total ARG ORFs detected: {total_arg}\n")
        f.write(f"Total MGE ORFs detected: {total_mge}\n")
        f.write(f"Contigs with ARGs: {len(arg_contigs)}\n")
        f.write(f"Contigs with MGEs: {len(mge_contigs)}\n")
        f.write(f"Contigs with both ARGs and MGEs: {len(shared_contigs)}\n\n")

        f.write("Co-localization Results (within 10 kb):\n")
        f.write("-" * 40 + "\n")
        f.write(f"ARG-MGE pairs within {distance_threshold/1000} kb: {total_pairs}\n")
        f.write(f"Co-localized contigs: {total_cocol}\n")
        f.write(f"Co-localization rate (of ARG contigs): {cocol_rate_arg:.2f}%\n")
        f.write(f"Co-localization rate (of MGE contigs): {cocol_rate_mge:.2f}%\n\n")

        f.write("Distance Statistics:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Average distance: {avg_distance:.2f} bp ({avg_distance/1000:.2f} kb)\n")
        f.write(f"Median distance: {median_distance} bp ({median_distance/1000:.2f} kb)\n")
        f.write(f"Min distance: {min(distances) if distances else 0} bp\n")
        f.write(f"Max distance: {max(distances) if distances else 0} bp\n")

    print(f"\nResults:")
    print(f"  ARG-MGE pairs within {distance_threshold/1000} kb: {total_pairs}")
    print(f"  Co-localized contigs: {total_cocol}")
    print(f"  Co-localization rate (of ARG contigs): {cocol_rate_arg:.2f}%")
    print(f"  Average distance: {avg_distance:.2f} bp")
    print(f"\nOutput files:")
    print(f"  Results: {output_file}")
    print(f"  Statistics: {stats_file}")

    # Also output simple list of co-localized contigs
    list_file = output_file.replace('.tsv', '_contigs.txt')
    with open(list_file, 'w') as f:
        for contig in sorted(cocol_contigs):
            f.write(f"{contig}\n")
    print(f"  Contig list: {list_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Identify ARG-MGE co-localized contigs (within 10 kb)"
    )
    parser.add_argument(
        "--arg-csv", "-a",
        required=True,
        help="Merged ARG DIAMOND results CSV"
    )
    parser.add_argument(
        "--mge-csv", "-m",
        required=True,
        help="Merged MGE DIAMOND results CSV"
    )
    parser.add_argument(
        "--gff-dir", "-g",
        required=True,
        help="Directory with Prodigal GFF files"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output TSV file for co-localization results"
    )
    parser.add_argument(
        "--stats", "-s",
        required=True,
        help="Output file for statistics"
    )
    parser.add_argument(
        "--distance", "-d",
        type=int,
        default=10000,
        help="Distance threshold in bp (default: 10000 = 10 kb)"
    )

    args = parser.parse_args()

    find_cocol_contigs_with_distance(
        args.arg_csv,
        args.mge_csv,
        args.gff_dir,
        args.output,
        args.stats,
        args.distance
    )


if __name__ == "__main__":
    main()