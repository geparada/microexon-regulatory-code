#!/usr/bin/env python

from pkg_resources import resource_filename
from pangolin.model import *
import argparse
import csv
import numpy as np
import torch
from Bio import SeqIO
from Bio.Seq import Seq

Genome = {}

def Genomictabulator(fasta):
    with open(fasta) as f:
        for chrfa in SeqIO.parse(f, "fasta"):
            Genome[chrfa.id] = str(chrfa.seq.upper())

def apply_variant(seq, pos, ref, alt):
    """
    Apply a variant to a sequence using efficient string slicing.

    Parameters:
    - seq: str
        The original sequence.
    - pos: int
        The 0-based position in the sequence where the variant occurs.
    - ref: str
        The reference allele.
    - alt: str
        The alternate allele.

    Returns:
    - mutated_seq: str
        The mutated sequence.
    """
    variant_start = pos
    variant_end = pos + len(ref)

    # Verify that the reference allele matches the sequence at the position
    if seq[variant_start:variant_end].upper() != ref.upper():
        raise ValueError(
            f"Reference allele '{ref}' does not match the sequence at position {pos}: "
            f"found '{seq[variant_start:variant_end]}' instead."
        )

    # Construct the mutated sequence using string slicing and concatenation
    mutated_seq = seq[:variant_start] + alt.upper() + seq[variant_end:]

    return mutated_seq

def load_models():
    """
    Load Pangolin models once and return them.
    """
    # Define the models and tissues
    model_nums = [0, 2, 4, 6]
    model_names = ["Heart", "Liver", "Brain", "Testis"]

    models = []
    for i in model_nums:
        model_set = []
        for j in range(1, 6):
            model = Pangolin(L, W, AR)
            if torch.cuda.is_available():
                model.cuda()
                weights = torch.load(resource_filename("pangolin", f"models/final.{j}.{i}.3"))
            else:
                weights = torch.load(resource_filename("pangolin", f"models/final.{j}.{i}.3"),
                                     map_location=torch.device('cpu'))
            model.load_state_dict(weights)
            model.eval()
            model_set.append(model)
        models.append((model_names[model_nums.index(i)], model_set))
    return models

def pangolin_PSI(seq, strand, models):
    """
    Predict PSI using Pangolin models for a single sequence and strand.
    """
    model_nums = [0, 2, 4, 6]
    model_names = ["Heart", "Liver", "Brain", "Testis"]

    IN_MAP = np.asarray([[0, 0, 0, 0],
                         [1, 0, 0, 0],
                         [0, 1, 0, 0],
                         [0, 0, 1, 0],
                         [0, 0, 0, 1]])
    INDEX_MAP = {0:1, 1:2, 2:4, 3:5, 4:7, 5:8, 6:10, 7:11}

    def one_hot_encode(seq, strand):
        seq = seq.upper().replace('A', '1').replace('C', '2')
        seq = seq.replace('G', '3').replace('T', '4').replace('N', '0')
        if strand == '+':
            seq = np.asarray(list(map(int, list(seq))))
        elif strand == '-':
            seq = np.asarray(list(map(int, list(seq[::-1]))))
            seq = (5 - seq) % 5  # Reverse complement
        return IN_MAP[seq.astype('int8')]

    seq_encoded = one_hot_encode(seq, strand).T
    seq_tensor = torch.from_numpy(np.expand_dims(seq_encoded, axis=0)).float()

    if torch.cuda.is_available():
        seq_tensor = seq_tensor.to(torch.device("cuda"))

    results = []
    for tissue, model_set in models:
        model_num = model_nums[model_names.index(tissue)]
        score = []
        # Average across 5 models
        for model in model_set:
            with torch.no_grad():
                output = model(seq_tensor)
                score.append(output[0][INDEX_MAP[model_num], :].cpu().numpy())
        mean_score = np.mean(score, axis=0)[0]
        std_score = np.std(score, axis=0)[0]
        results.append([tissue, mean_score, std_score])
    return results

def main():
    parser = argparse.ArgumentParser(description='Predict delta PSI for variants using Pangolin.')
    parser.add_argument('--genome', required=True, help='Path to genome FASTA file')
    parser.add_argument('--variant_file', required=True, help='Path to variant file')
    parser.add_argument('--event', required=True, help='Event identifier')
    parser.add_argument('--mic_anno', required=True, help='Path to microexon annotation file')
    parser.add_argument('--output', required=True, help='Output file path')
    args = parser.parse_args()

    # Load genome
    Genomictabulator(args.genome)

    # Load models once
    models = load_models()

    # Load microexon annotation
    human_event_coords = dict()
    with open(args.mic_anno) as file:
        reader = csv.DictReader(file, delimiter=",")
        for row in reader:
            lengthDiff = int(row["lengthDiff"])
            chrom = row["chrom"]
            strand = row["strand"]
            C1 = int(row["upIntStart"])
            C2 = int(row["dnIntEnd"])
            ME_start = int(row['upIntEnd'])
            ME_end = int(row['dnIntStart'])
            if strand == "-":
                C1 = int(row['upIntEnd'])
                C2 = int(row['dnIntStart'])
                ME_start = int(row['dnIntEnd'])
                ME_end = int(row['upIntStart'])
            human_event_coords[row["event"]] = (C1, C2, ME_start, ME_end-1, chrom, strand, lengthDiff)

    event = args.event

    # Retrieve event coordinates
    if event not in human_event_coords:
        print(f"Event {event} not found in microexon annotation. Exiting.")
        return

    C1, C2, ME_start, ME_end, chrom_event, strand, lengthDiff = human_event_coords[event]

    chrom_seq = Genome[chrom_event]
    chrom_length = len(chrom_seq)

    # Initialize output
    with open(args.output, 'w') as out, open(args.variant_file) as variant_file:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(['event', 'var_ID', 'chrom', 'pos', 'ref', 'alt', 'strand', 'ss', 'tissue', 'WT_PSI_mean', 'WT_PSI_sd', 'MUT_PSI_mean', 'MUT_PSI_sd', 'dPSI', 'norm_dPSI'])

        variant_reader = csv.DictReader(variant_file, delimiter='\t')  # Adjust delimiter if needed

        # Precompute wild-type sequences and PSI for ss3 and ss5
        WT_res = {}

        for ss_type in ['3ss', '5ss']:
            # Determine the splice site position
            if ss_type == '3ss':
                ss_pos = ME_start if strand == '+' else ME_end
            elif ss_type == '5ss':
                ss_pos = ME_end if strand == '+' else ME_start

            # Extract wild-type sequence around the splice site
            seq_start = ss_pos - 5000
            seq_end = ss_pos + 5001

            # Adjust positions for 0-based indexing -- ##corrected
            seq_start_idx = max(seq_start, 0)
            seq_end_idx = min(seq_end, chrom_length)

            wt_seq = chrom_seq[seq_start_idx:seq_end_idx]

            # Handle boundary conditions by padding with 'N's if necessary
            if seq_start - 1 < 0:
                padding = 'N' * abs(seq_start)
                wt_seq = padding + wt_seq
            if seq_end - 1 > chrom_length:
                padding = 'N' * (seq_end - chrom_length)
                wt_seq = wt_seq + padding

            # Convert to uppercase
            wt_seq = wt_seq.upper()

            # Reverse complement if on the negative strand
            if strand == '-':
                wt_seq = str(Seq(wt_seq).reverse_complement())

            # Predict PSI for wild-type sequence
            WT_values = pangolin_PSI(wt_seq, '+', models)  # After reverse complement, treat as '+'

            for result in WT_values:
                tissue, PSI_mean, PSI_sd = result
                WT_res[(ss_type, tissue)] = (PSI_mean, PSI_sd)

        for row in variant_reader:
            var_ID = row['variant']

            # Parse variant information
            try:
                chrom_var, pos_var, ref_var, alt_var = var_ID.split(':')
                pos_var = int(pos_var)  # 1-based position
            except ValueError:
                print(f"Skipping invalid variant format: {var_ID}")
                continue

            # Ensure that the variant chromosome matches the event chromosome
            if chrom_var != chrom_event:
                print(f"Variant chromosome {chrom_var} does not match event chromosome {chrom_event} for event {event}. Skipping.")
                continue

            # Convert positions to 0-based index for chrom_seq
            pos_var0 = pos_var - 1  # 0-based index

            # Apply variant to the chromosome sequence
            try:
                mut_chrom_seq = apply_variant(chrom_seq, pos_var0, ref_var, alt_var)
            except Exception as e:
                print(f"Error applying variant {var_ID} to chromosome sequence: {e}")
                continue

            # Calculate variant length change
            variant_length_change = len(alt_var) - len(ref_var)

            for ss_type in ['3ss', '5ss']:
                # Determine the splice site position
                if ss_type == '3ss':
                    ss_pos = ME_start if strand == '+' else ME_end
                elif ss_type == '5ss':
                    ss_pos = ME_end if strand == '+' else ME_start

                # Determine the new splice site position after applying the variant
                if pos_var < ss_pos:
                    ss_pos_mut = ss_pos + variant_length_change
                else:
                    ss_pos_mut = ss_pos

                # Extract mutant sequence around the new splice site position
                seq_start_mut = ss_pos_mut - 5000
                seq_end_mut = ss_pos_mut + 5001

                # Adjust positions for 0-based indexing
                seq_start_idx_mut = max(seq_start_mut, 0)
                seq_end_idx_mut = min(seq_end_mut, len(mut_chrom_seq))

                mut_seq = mut_chrom_seq[seq_start_idx_mut:seq_end_idx_mut]

                # Handle boundary conditions by padding with 'N's if necessary
                if seq_start_mut < 0:
                    padding = 'N' * abs(seq_start_mut)
                    mut_seq = padding + mut_seq
                if seq_end_mut > len(mut_chrom_seq):
                    padding = 'N' * (seq_end_mut - len(mut_chrom_seq))
                    mut_seq = mut_seq + padding

                # Convert to uppercase
                mut_seq = mut_seq.upper()

                # Reverse complement if on the negative strand
                if strand == '-':
                    mut_seq = str(Seq(mut_seq).reverse_complement())

                # Predict PSI for mutant sequence
                try:
                    MUT_values = pangolin_PSI(mut_seq, '+', models)  # After reverse complement, treat as '+'
                except Exception as e:
                    print(f"Error predicting MUT PSI for event {event}, variant {var_ID}: {e}")
                    continue

                # Process results and calculate delta PSI
                for result in MUT_values:
                    tissue, MUT_PSI_mean, MUT_PSI_sd = result
                    if (ss_type, tissue) in WT_res:
                        WT_PSI_mean, WT_PSI_sd = WT_res[(ss_type, tissue)]
                        dPSI = MUT_PSI_mean - WT_PSI_mean
                        if WT_PSI_mean != 0:
                            norm_dPSI = dPSI / WT_PSI_mean
                        else:
                            norm_dPSI = 0
                        writer.writerow([event, var_ID, chrom_event, pos_var, ref_var, alt_var, strand, ss_type, tissue, WT_PSI_mean, WT_PSI_sd, MUT_PSI_mean, MUT_PSI_sd, dPSI, norm_dPSI])
                    else:
                        print(f"No WT PSI value for tissue {tissue}, ss {ss_type} in event {event}")

if __name__ == "__main__":
    main()