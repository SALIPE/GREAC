import argparse
import glob
import os

# Script to get CASTOR-KRFE input HIV example and separate files classes for comparison.



def gather_input_files(paths):
    """
    Given a list of file or directory paths, return all existing FASTA files.
    """
    fasta_files = []
    for p in paths:
        if os.path.isdir(p):
            fasta_files.extend(glob.glob(os.path.join(p, '*.fasta')))
            fasta_files.extend(glob.glob(os.path.join(p, '*.fa')))
        elif os.path.isfile(p):
            fasta_files.append(p)
        else:
            print(f"Warning: '{p}' not found, skipping.")
    return fasta_files

def split_fasta_by_class(input_files, output_dir):
    # Create base output directory if needed
    os.makedirs(output_dir, exist_ok=True)

    # Keep file handles open for each class
    handles = {}

    try:
        for infile in input_files:
            with open(infile) as f:
                header = None
                seq_lines = []
                for line in f:
                    line = line.rstrip('\n')
                    if line.startswith('>'):
                        if header:
                            write_record(header, seq_lines, handles, output_dir)
                        header = line
                        seq_lines = []
                    else:
                        seq_lines.append(line)
                if header:
                    write_record(header, seq_lines, handles, output_dir)
    finally:
        # Close all open handles
        for h in handles.values():
            h.close()


def extract_class_from_header(header):
    """
    Extract the class identifier from a FASTA header of the form '>ID|CLASS'.
    Returns 'unknown' if the pattern is not found.
    """
    parts = header[1:].split('|')
    return parts[1] if len(parts) >= 2 else 'unknown'


def write_record(header, seq_lines, handles, output_dir):
    cl = extract_class_from_header(header)
    # Prepare directory for this class
    class_dir = os.path.join(output_dir, cl)
    os.makedirs(class_dir, exist_ok=True)

    # Open handle for this class if not already
    if cl not in handles:
        out_path = os.path.join(class_dir, f"{cl}.fasta")
        handles[cl] = open(out_path, 'a')
    h = handles[cl]

    # Write record
    h.write(header + '\n')
    h.write('\n'.join(seq_lines) + '\n')


if __name__ == '__main__':
    home = '/home/salipe/Desktop/rrm-genomic-extractor/comparison_scripts'
    training_file = f'{home}/castor_hiv_data/training_data.fasta'
    testing_file = f'{home}/castor_hiv_data/testing_data.fasta'
    fasta_files = gather_input_files([training_file,testing_file])
    split_fasta_by_class(fasta_files, f'{home}/castor_hiv_data/variants') 
