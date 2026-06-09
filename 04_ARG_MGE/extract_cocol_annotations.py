#!/usr/bin/env python3
# ============================================================
# Extract Annotations for Co-localized Contigs
# For Nature article supplementary materials
#
# This script extracts protein sequences and gene annotations
# for contigs that contain both ARGs and MGEs.
# ============================================================

import argparse
import os
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def load_contig_names(contigs_file):
    """
    Load contig names from file.

    Args:
        contigs_file: Path to file with contig names

    Returns:
        Set of contig names
    """
    contigs = set()
    with open(contigs_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                contigs.add(line)
    return contigs


def extract_proteins(contigs_file, faa_dir, output_faa):
    """
    Extract protein sequences for co-localized contigs.

    Args:
        contigs_file: File with list of contig names
        faa_dir: Directory with Prodigal FAA files
        output_faa: Output FASTA file
    """
    print(f"Extracting proteins from: {faa_dir}")

    cocol_contigs = load_contig_names(contigs_file)
    protein_count = 0
    contig_count = 0

    with open(output_faa, 'w') as out_f:
        for faa_file in glob.glob(os.path.join(faa_dir, "*.faa")):
            sample_name = os.path.basename(faa_file).replace('.faa', '')

            for record in SeqIO.parse(faa_file, 'fasta'):
                # Extract contig ID from ORF ID
                # Format: contigID_orfNUMBER
                parts = record.id.rsplit('_', 1)
                if len(parts) == 2 and parts[1].startswith('orf'):
                    contig_id = parts[0]
                else:
                    contig_id = '_'.join(parts[:-1]) if len(parts) > 1 else parts[0]

                if contig_id in cocol_contigs:
                    # Add sample prefix if not present
                    if not record.id.startswith(sample_name):
                        record.id = f"{sample_name}_{record.id}"

                    SeqIO.write(record, out_f, 'fasta')
                    protein_count += 1

    print(f"  Extracted {protein_count} proteins")
    return protein_count


def extract_genes(contigs_file, gff_dir, output_gff):
    """
    Extract gene annotations for co-localized contigs.

    Args:
        contigs_file: File with list of contig names
        gff_dir: Directory with Prodigal GFF files
        output_gff: Output GFF file
    """
    print(f"Extracting genes from: {gff_dir}")

    cocol_contigs = load_contig_names(contigs_file)
    gene_count = 0

    with open(output_gff, 'w') as out_gff:
        for gff_file in glob.glob(os.path.join(gff_dir, "*.gff")):
            sample_name = os.path.basename(gff_file).replace('.gff', '')

            with open(gff_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    fields = line.strip().split('\t')
                    if len(fields) < 9:
                        continue

                    contig_id = fields[0]

                    if contig_id in cocol_contigs:
                        # Add sample prefix if not present
                        if not contig_id.startswith(sample_name):
                            fields[0] = f"{sample_name}_{contig_id}"

                        out_gff.write('\t'.join(fields) + '\n')
                        gene_count += 1

    print(f"  Extracted {gene_count} gene annotations")
    return gene_count


def create_annotation_summary(contigs_file, arg_csv, mge_csv, faa_dir, output_summary):
    """
    Create comprehensive annotation summary for co-localized contigs.

    Args:
        contigs_file: File with list of contig names
        arg_csv: Merged ARG DIAMOND results
        mge_csv: Merged MGE DIAMOND results
        faa_dir: Directory with FAA files
        output_summary: Output TSV file
    """
    print("Creating annotation summary...")

    cocol_contigs = load_contig_names(contigs_file)

    # Load ARG annotations
    arg_data = {}
    with open(arg_csv, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) > 1:
                orf_id = fields[0]
                parts = orf_id.rsplit('_', 1)
                contig_id = parts[0] if len(parts) > 1 else orf_id
                if contig_id in cocol_contigs:
                    arg_data[orf_id] = {
                        'sseqid': fields[1] if len(fields) > 1 else '',
                        'identity': fields[2] if len(fields) > 2 else '',
                        'evalue': fields[10] if len(fields) > 10 else ''
                    }

    # Load MGE annotations
    mge_data = {}
    with open(mge_csv, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) > 1:
                orf_id = fields[0]
                parts = orf_id.rsplit('_', 1)
                contig_id = parts[0] if len(parts) > 1 else orf_id
                if contig_id in cocol_contigs:
                    mge_data[orf_id] = {
                        'sseqid': fields[1] if len(fields) > 1 else '',
                        'identity': fields[2] if len(fields) > 2 else '',
                        'evalue': fields[10] if len(fields) > 10 else ''
                    }

    # Write summary
    with open(output_summary, 'w') as out:
        out.write("Contig_ID\tORF_ID\tAnnotation_Type\tHit_ID\tIdentity\tEvalue\n")

        # Write ARG annotations
        for orf_id, data in arg_data.items():
            parts = orf_id.rsplit('_', 1)
            contig_id = parts[0] if len(parts) > 1 else orf_id
            out.write(f"{contig_id}\t{orf_id}\tARG\t{data['sseqid']}\t{data['identity']}\t{data['evalue']}\n")

        # Write MGE annotations
        for orf_id, data in mge_data.items():
            parts = orf_id.rsplit('_', 1)
            contig_id = parts[0] if len(parts) > 1 else orf_id
            out.write(f"{contig_id}\t{orf_id}\tMGE\t{data['sseqid']}\t{data['identity']}\t{data['evalue']}\n")

    print(f"  Summary written to: {output_summary}")
    return len(arg_data), len(mge_data)


def main():
    parser = argparse.ArgumentParser(
        description="Extract annotations for ARG-MGE co-localized contigs"
    )
    parser.add_argument(
        "--contigs", "-c",
        required=True,
        help="File with list of co-localized contig names"
    )
    parser.add_argument(
        "--faa-dir", "-f",
        required=True,
        help="Directory with Prodigal FAA files"
    )
    parser.add_argument(
        "--gff-dir", "-g",
        required=True,
        help="Directory with Prodigal GFF files"
    )
    parser.add_argument(
        "--arg-csv", "-a",
        help="Merged ARG DIAMOND results (for summary)"
    )
    parser.add_argument(
        "--mge-csv", "-m",
        help="Merged MGE DIAMOND results (for summary)"
    )
    parser.add_argument(
        "--output-dir", "-o",
        required=True,
        help="Output directory"
    )

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Extract proteins
    output_faa = os.path.join(args.output_dir, "cocol_proteins.faa")
    extract_proteins(args.contigs, args.faa_dir, output_faa)

    # Extract genes
    output_gff = os.path.join(args.output_dir, "cocol_genes.gff")
    extract_genes(args.contigs, args.gff_dir, output_gff)

    # Create summary if both ARG and MGE CSVs provided
    if args.arg_csv and args.mge_csv:
        output_summary = os.path.join(args.output_dir, "cocol_annotation_summary.tsv")
        create_annotation_summary(
            args.contigs, args.arg_csv, args.mge_csv,
            args.faa_dir, output_summary
        )

    print("\nAnnotation extraction complete!")


if __name__ == "__main__":
    main()