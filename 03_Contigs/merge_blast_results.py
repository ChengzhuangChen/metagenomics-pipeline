#!/usr/bin/env python3
# ============================================================
# Merge BLAST/DIAMOND Results
# For Nature article supplementary materials
#
# This script processes multiple sample BLAST/DIAMOND result files,
# adds sample prefixes to query IDs, and merges them into a single file.
# ============================================================

import glob
import os
import csv
import argparse


def process_and_merge_files(source_dir, output_dir, suffix, sample_prefix):
    """
    Process BLAST result files and merge them.

    Args:
        source_dir: Directory containing BLAST result files
        output_dir: Output directory for merged results
        suffix: File suffix to match (e.g., "_diamond_blastp_arg.txt")
        sample_prefix: Prefix to add to sample IDs in output
    """
    os.makedirs(output_dir, exist_ok=True)

    processed_files = []
    empty_files = []

    # Find all matching files
    file_pattern = os.path.join(source_dir, f"*{suffix}")
    files = glob.glob(file_pattern)

    if not files:
        print(f"No files matching '{suffix}' found in {source_dir}")
        return

    # Process each file
    for file_path in files:
        file_name = os.path.basename(file_path)

        # Extract sample ID from filename
        sample_id = file_name.replace(suffix, "")

        # Define CSV output path
        csv_file_name = f"{sample_id}_processed.csv"
        csv_file_path = os.path.join(output_dir, csv_file_name)

        # Process file
        has_data = False
        with open(file_path, "r") as infile, open(csv_file_path, "w", newline='') as outfile:
            csv_writer = csv.writer(outfile, delimiter=',')

            for line_num, line in enumerate(infile):
                stripped_line = line.strip()
                if not stripped_line:
                    continue

                columns = stripped_line.split('\t')

                if line_num == 0:
                    # Write header
                    csv_writer.writerow(columns)
                else:
                    # Mark as having data
                    has_data = True
                    if len(columns) > 0:
                        # Add sample prefix to query ID
                        columns[0] = f"{sample_prefix}{columns[0]}"
                        csv_writer.writerow(columns)

        if not has_data:
            empty_files.append(file_name)
            print(f"Note: {file_name} has only header, no data rows")
        else:
            print(f"Generated CSV: {csv_file_name}")

        processed_files.append(csv_file_path)

    # Merge all processed CSV files
    if processed_files:
        merged_file = os.path.join(output_dir, f"all_samples_merged{suffix.replace('.txt', '.csv')}")

        with open(merged_file, "w", newline='') as outfile:
            csv_writer = csv.writer(outfile, delimiter=',')
            header_written = False

            for file in processed_files:
                with open(file, "r") as infile:
                    csv_reader = csv.reader(infile, delimiter=',')
                    for row_num, row in enumerate(csv_reader):
                        if row_num == 0 and not header_written:
                            csv_writer.writerow(row)
                            header_written = True
                        elif row_num > 0:
                            csv_writer.writerow(row)

        print(f"\nMerged all files to: {merged_file}")

        if empty_files:
            print(f"Note: {len(empty_files)} files had no data rows")

    else:
        print(f"No files found in {source_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Merge BLAST/DIAMOND results from multiple samples"
    )
    parser.add_argument(
        "--source-dir", "-i",
        required=True,
        help="Directory containing BLAST result files"
    )
    parser.add_argument(
        "--output-dir", "-o",
        required=True,
        help="Output directory for merged results"
    )
    parser.add_argument(
        "--suffix", "-s",
        required=True,
        help="File suffix to match (e.g., '_diamond_blastp_arg.txt')"
    )
    parser.add_argument(
        "--sample-prefix", "-p",
        default="",
        help="Prefix to add to sample IDs (default: none)"
    )

    args = parser.parse_args()

    process_and_merge_files(
        args.source_dir,
        args.output_dir,
        args.suffix,
        args.sample_prefix
    )


if __name__ == "__main__":
    main()