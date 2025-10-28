import subprocess
import pandas as pd
import os
import shutil
import datetime

def safe_remove(filepath):
    """Remove a file if it exists."""
    try:
        if os.path.exists(filepath):
            os.remove(filepath)
    except Exception as e:
        print(f"Warning: Could not remove file {filepath}. Error: {e}")

# function: Filters and prepares the header for the SAM file to be used with ngsLCA.
def filter_header(header_full, header_new, header_subset):
    """
    Filters and prepares the header for the SAM file to be used with ngsLCA.

    Parameters:
    - header_full (str): Path to the file containing the full header.
    - header_new (str): Path to the file containing new headers to include.
    - header_subset (str): Output path for the filtered header file.
    """
    df1 = pd.read_csv(header_new, header=None)
    df1 = df1[0].unique()
    df1 = [f"SN:{value}" for value in df1]

    header = pd.read_csv(header_full, header=None, sep='\t')
    header = header[header[1].isin(df1)]

    header.to_csv(header_subset, header=False, index=False, sep='\t')


# function: Processes BAM files to prepare for ngsLCA.
def process_bam_files(bam_file, bname, thread, output_base):
    """
    Processes BAM files to prepare for ngsLCA, including extracting headers and converting to sorted BAM.

    Parameters:
    - bam_file (str): Path to the input BAM file.
    - bname (str): Base name for the output files.
    - thread (int): Number of threads to use for processing.
    - output_base (str): Base path for the output files.

    Returns:
    - str: Path to the processed BAM file, or None if processing is skipped.
    """
    header_subset_1_file = f"{output_base}.{bname}.header_subset.1.txt"
    header_subset_2_file = f"{output_base}.{bname}.header_subset.2.txt"
    header_subset_3_file = f"{output_base}.{bname}.header_subset.3.txt"
    header_new_file = f"{output_base}.{bname}.header_new.txt"
    alignment_file = f"{output_base}.{bname}.alignment.txt"
    processed_bam = f"{output_base}.{bname}.bam"
    
    # Check if BAM file contains any mapped reads
    result = subprocess.run(f'samtools view {bam_file} | head -n 1', shell=True, capture_output=True, text=True)
    if not result.stdout.strip():
        print(f'No reads mapped for file {bname}. Skipping.')
        return None
    
    # Extract header information from BAM
    subprocess.run(f"samtools view -H {bam_file} | grep @HD > {header_subset_1_file}", shell=True, check=True)
    subprocess.run(f"samtools view -H {bam_file} | grep @SQ > {header_subset_2_file}", shell=True, check=True)
    subprocess.run(f"samtools view -H {bam_file} | grep @PG > {header_subset_3_file}", shell=True, check=True)

    # Extract alignment information
    subprocess.run(f"samtools view {bam_file} > {alignment_file}", shell=True, check=True)

    # Remove duplicate @SQ lines using awk
    awk_command = f'''awk '!seen[$0]++' {header_subset_2_file} > {header_subset_2_file}.tmp'''
    subprocess.run(awk_command, shell=True, check=True)
    safe_remove(header_subset_2_file)
    os.rename(header_subset_2_file + ".tmp", header_subset_2_file)

    # Extract unique values from alignment_file
    extract_cmd = f"cat {alignment_file} | cut -f3 | sort -u | uniq > {header_new_file}"
    subprocess.run(extract_cmd, shell=True, check=True)

    # Filter header: rename header_subset_2_file for filtering, then restore
    os.rename(header_subset_2_file, header_subset_2_file + ".tmp")
    filter_header(header_subset_2_file + ".tmp", header_new_file, header_subset_2_file)
    subprocess.run(
        f'cat {header_subset_1_file} {header_subset_2_file} {header_subset_3_file} {alignment_file} | samtools view -@ {thread} -b -o {processed_bam}',
        shell=True, check=True)

    # Remove intermediate files
    safe_remove(header_subset_1_file)
    safe_remove(header_subset_2_file)
    safe_remove(header_subset_3_file)
    safe_remove(header_new_file)
    safe_remove(alignment_file)
    safe_remove(header_subset_2_file + ".tmp")
    
    return processed_bam


# function: Merges and sorts BAM files, preparing for ngsLCA.
def merge_sort(bam_path, output_base, thread):
    """
    Merges and sorts BAM files, preparing for ngsLCA analysis.

    Parameters:
    - bam_path (str): Directory containing BAM files to merge and sort.
    - output_base (str): Base name for the output sorted BAM file.
    - thread (int): Number of threads to use for processing.

    Returns:
    - str: Path to the sorted BAM file.
    """
    # Check the number of BAM files in bam_path that start with the basename of output_base
    bam_files = [f for f in os.listdir(bam_path) if f.endswith('.bam')
                 and f.startswith(os.path.basename(output_base))]
    if len(bam_files) == 1:
        # If only one BAM file, skip merging and just sort
        merged_bam = os.path.join(bam_path, bam_files[0])
    else:
        # Merge multiple BAM files
        merge_cmd = f'''samtools merge -n -@ {thread} -o {output_base}.merged.bam {' '.join([f'{bam_path}/{f}' for f in bam_files])}'''
        subprocess.run(merge_cmd, shell=True, check=True)
        merged_bam = f"{output_base}.merged.bam"

    # Sort the BAM file
    sort_cmd = f"samtools sort -n -@ {thread} -O bam -o {output_base}.sorted.bam {merged_bam}"
    subprocess.run(sort_cmd, shell=True, check=True)

    return f"{output_base}.sorted.bam"


# function: call ngsLCA to run eDNA pipeline
def run_ngsLCA(bam, names_dmp, nodes_dmp, acc2tax, minedit, maxedit, output):
    """
    Executes ngsLCA for taxonomic identification based on specified parameters.

    Parameters:
    - bam (str): Path to the input BAM file for ngsLCA.
    - names_dmp, nodes_dmp, acc2tax (str): Paths to NCBI taxonomy and accession to taxonomy mapping files.
    - minedit, maxedit (int): Minimum and maximum edit distances for ngsLCA identification.
    - output (str): Base path for the output file generated by ngsLCA.

    Returns:
    - str: Path to the LCA output file.
    """
    ngsLCA_cmd = f"ngsLCA -names {names_dmp} -nodes {nodes_dmp} -acc2tax {acc2tax} -editdistmin {minedit} -editdistmax {maxedit} -bam {bam} -outnames {output}.min{minedit}max{maxedit}"
    subprocess.run(ngsLCA_cmd, shell=True, check=True)
    return f"{output}.min{minedit}max{maxedit}.lca"

if __name__ == '__main__':
    import argparse

    usage = """ python run_ngsLCA.py -b bam_folder --thread 16 --names names.dmp --nodes nodes.dmp -a acc2tax.txt --minedit 0 --maxedit 2 -o output """

    parser = argparse.ArgumentParser(description=usage)
    parser.add_argument("-b", "--bam_folder", dest="bam_folder", action="store", nargs='?',
                        help="Input folder containing BAM files of different libraries.", metavar="FOLDER", required=True)
    parser.add_argument("--thread", dest="thread", action="store", nargs='?', default=1,
                        help="Number of thread (default: 1). ", metavar="STRING")
    parser.add_argument("--names", dest="names", action="store", nargs='?',
                        help="NCBI taxonomy names.dmp", metavar="FILE", required=True)
    parser.add_argument("--nodes", dest="nodes", action="store", nargs='?',
                        help="NCBI taxonomy nodes.dmp", metavar="FILE", required=True)
    parser.add_argument("-a", "--acc2tax", dest="acc2tax", action="store", nargs='?',
                        help="accession2taxid.txt, the match file of genome accession and taxid.", metavar="FILE",
                        required=True)
    parser.add_argument("--minedit", dest="minedit", action="store", nargs='?', default=0,
                        help="Minimum edit distance to assign a sequence in ngsLCA (default: 0).", metavar="STRING",
                        required=True)
    parser.add_argument("--maxedit", dest="maxedit", action="store", nargs='?', default=2,
                        help="Maximum edit distance to assign a sequence in ngsLCA.", metavar="STRING", required=True)
    parser.add_argument("-o", "--output", dest="output", action="store", nargs='?',
                        help="Output prefix.", metavar="STRING", required=True)

    args = parser.parse_args()
    
    # Create a temporary directory for intermediate files
    tmp_dir = f'{args.output}_ngsLCA_temp'
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)
    output_base = os.path.join(tmp_dir, os.path.basename(args.output))

    print("Processing BAM files...", datetime.datetime.now())
    for bam_file in os.listdir(args.bam_folder):
        print(f"Checking file: {bam_file}")

        if bam_file.startswith(".") or not bam_file.endswith(".bam"):
            print(f"Skipping invalid or hidden file: {bam_file}")
            continue

        full_path = os.path.join(args.bam_folder, bam_file)

        if not os.path.isfile(full_path):
            print(f"Skipping invalid path: {full_path}")
            continue

        if os.path.getsize(full_path) == 0:
            print(f"Skipping empty BAM file: {full_path}")
            continue
            
        result = process_bam_files(full_path, bam_file.replace(".bam", ""), int(args.thread), output_base)
        if result is None:
            print(f"Processing of {bam_file} failed or was skipped.")
    print("BAM files processing completed.", datetime.datetime.now())

    print("Merging and sorting BAM files...", datetime.datetime.now())
    prepared_bam = merge_sort(tmp_dir, output_base, int(args.thread))
    print("Merging and sorting completed.", datetime.datetime.now())

    print("Running taxonomic identification with ngsLCA...", datetime.datetime.now())
    lca_file = run_ngsLCA(prepared_bam, args.names, args.nodes, args.acc2tax, args.minedit, args.maxedit, args.output)
    
    # Clean up temporary directory
    shutil.rmtree(tmp_dir)