import os
import torch
import copy
import argparse
import csv
import numpy as np
from pyfasta import Fasta
from deltasplice.constant import (
    default_model_paths,
    model,
    EL,
    SeqTable,
    repdict,
    IN_MAP,
)
from Bio.Seq import Seq

# Adjust the default model paths if necessary
default_model_paths = ['/path/to/DeltaSplice/' + x for x in default_model_paths]

# Function to load models onto the CPU
def load_models():
    Models = [copy.deepcopy(model) for _ in default_model_paths]
    for m, b in zip(Models, default_model_paths):
        # Load the model state dict onto the CPU
        m.load_state_dict(torch.load(b, map_location=torch.device('cpu')))
        m.eval()  # Set the model to evaluation mode
    return Models

# Function to load the reference genome once
def load_genome(genome_name):
    reference_genome = Fasta(genome_name)
    return reference_genome

# Function to apply a variant to a sequence
def apply_variant(seq, pos, ref, alt):
    """
    Apply a variant to a sequence.

    Parameters:
    - seq: The original sequence (as a string).
    - pos: The 0-based position in the sequence where the variant occurs.
    - ref: The reference allele.
    - alt: The alternate allele.

    Returns:
    - The mutated sequence (as a string).
    """
    seq_list = list(seq)

    variant_start = pos
    variant_end = pos + len(ref)

    # Verify that the reference allele matches the sequence
    if seq[variant_start:variant_end].upper() != ref.upper():
        raise ValueError(f"Reference allele does not match the sequence at position {pos}.")

    # Replace the reference allele with the alternate allele
    seq_list[variant_start:variant_end] = list(alt.upper())

    mutated_seq = ''.join(seq_list)
    return mutated_seq

# Prediction function running on CPU
def predict_ssu(seq, Models, ss):
    # Convert sequence to uppercase
    seq = seq.upper()

    # Convert bases to indices
    try:
        seq_indices = [SeqTable[base] for base in seq]
    except KeyError as e:
        raise ValueError(f"Invalid base '{e.args[0]}' found in sequence.")

    # Convert indices to one-hot encoding
    seq_array = IN_MAP[seq_indices][:, :4]

    # Convert to tensor on CPU
    seq_tensor = torch.tensor(seq_array, dtype=torch.float32)[None]  # Shape: (1, L, 4)

    # Predict using the loaded models
    pred = 0
    for m in Models:
        with torch.no_grad():
            pred += m.predict({"X": seq_tensor}, use_ref=False)["single_pred_psi"]
    pred = (pred / len(Models))[0]
    pred = pred[pred.shape[0] // 2]

    # Return only the requested SSU
    if ss == "3ss":
        return pred[1].item()  # Acceptor SSU
    elif ss == "5ss":
        return pred[2].item()  # Donor SSU
    else:
        raise ValueError("Invalid ss parameter. Must be '3ss' or '5ss'.")

def get_window_seq(ss_pos, strand, chrom_seq):
    """
    Extracts the sequence window around the splice site.

    Parameters:
    - ss_pos: int
        Splice site position (1-based coordinate).
    - strand: str
        Strand of the gene ('+' or '-').
    - chrom_seq: str
        The sequence of the chromosome.

    Returns:
    - window_seq: str
        The extracted sequence window around the splice site.
    """
    # Assuming EL is defined globally
    chrom_length = len(chrom_seq)

    # Calculate start and end positions for the sequence window
    seq_start = ss_pos - (EL // 2)
    seq_end = ss_pos + (EL // 2) + 1  # Include the base at seq_end

    # Adjust positions for 0-based indexing (Python uses 0-based indexing)
    seq_start_idx = max(seq_start - 1, 0)
    seq_end_idx = min(seq_end - 1, chrom_length)

    # Extract the sequence window
    window_seq = chrom_seq[seq_start_idx:seq_end_idx]

    # Handle boundary conditions by padding with 'N's if necessary
    if seq_start - 1 < 0:
        padding = 'N' * abs(seq_start - 1)
        window_seq = padding + window_seq
    if seq_end - 1 > chrom_length:
        padding = 'N' * (seq_end - 1 - chrom_length)
        window_seq = window_seq + padding

    # Convert to uppercase
    window_seq = window_seq.upper()

    # Reverse complement if on the negative strand
    if strand == '-':
        window_seq = str(Seq(window_seq).reverse_complement())

    return window_seq

def main():
    parser = argparse.ArgumentParser(description='Predict delta SSU for variants using DeltaSplice.')
    parser.add_argument('--genome', required=True, help='Path to genome FASTA file')
    parser.add_argument('--variant_file', required=True, help='Path to variant file')
    parser.add_argument('--var_column', required=True, help='Column name for variant IDs (chr:pos:ref:alt)')
    parser.add_argument('--event_column', default='event', help='Column name for event IDs')
    parser.add_argument('--mic_anno', required=True, help='Path to microexon annotation file')
    parser.add_argument('--output', required=True, help='Output file path')
    args = parser.parse_args()

    # Load models and genome
    Models = load_models()
    reference_genome = load_genome(args.genome)

    # Read microexon annotation and build human_event_coords dictionary
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

            human_event_coords[row["event"]] = (C1, C2, ME_start, ME_end, chrom, strand, lengthDiff)

    # Initialize output
    with open(args.output, 'w') as out, open(args.variant_file) as variant_file:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(['event', 'var_ID', 'chrom', 'pos', 'ref', 'alt', 'strand', 'ss', 'WT_SSU', 'MUT_SSU', 'delta_SSU'])

        variant_reader = csv.DictReader(variant_file, delimiter='\t')  # Adjust delimiter if needed

        for row in variant_reader:
            print(row)
            var_ID = row[args.var_column]
            event = row[args.event_column]

            # Parse variant information
            try:
                chrom_var, pos_var, ref_var, alt_var = var_ID.split(':')
                pos_var = int(pos_var)  # 1-based position
            except ValueError:
                print(f"Skipping invalid variant format: {var_ID}")
                continue

            # Retrieve event coordinates
            if event not in human_event_coords:
                print(f"Event {event} not found in microexon annotation. Skipping.")
                continue

            C1, C2, ME_start, ME_end, chrom_event, strand, lengthDiff = human_event_coords[event]

            # Ensure that the variant chromosome matches the event chromosome
            if chrom_var != chrom_event:
                print(f"Variant chromosome {chrom_var} does not match event chromosome {chrom_event} for event {event}. Skipping.")
                continue

            # For each splice site (3ss and 5ss)
            for ss_type in ['3ss', '5ss']:
                # Determine the splice site position
                if ss_type == '3ss':
                    ss_pos = ME_start if strand == '+' else ME_end
                elif ss_type == '5ss':
                    ss_pos = ME_end if strand == '+' else ME_start

                chrom_seq = reference_genome[chrom_event]
                chrom_length = len(chrom_seq)

                # Get the wild-type sequence
                wt_seq = get_window_seq(ss_pos, strand, chrom_seq)

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

                # Determine the new splice site position after applying the variant
                if strand == '+':
                    if pos_var < ss_pos:
                        ss_pos_mut = ss_pos + variant_length_change
                    else:
                        ss_pos_mut = ss_pos
                elif strand == '-':
                    if pos_var > ss_pos:
                        ss_pos_mut = ss_pos - variant_length_change
                    else:
                        ss_pos_mut = ss_pos
                else:
                    raise ValueError(f"Invalid strand '{strand}' for event {event}")

                # Check if the variant spans the splice site
                variant_start = pos_var
                variant_end = pos_var + len(ref_var) - 1

                if strand == '+':
                    spans_ss = (variant_start <= ss_pos) and (variant_end >= ss_pos)
                else:  # Negative strand
                    spans_ss = (variant_start >= ss_pos) and (variant_end <= ss_pos)

                # Predict SSU for wild-type sequence
                try:
                    WT_SSU = predict_ssu(wt_seq, Models, ss_type)
                except Exception as e:
                    print(f"Error predicting WT SSU for event {event}, variant {var_ID}: {e}")
                    continue

                if spans_ss:
                    # Variant spans the splice site; evaluate over a range
                    MUT_SSUs = []
                    # Define the range around ss_pos_mut
                    for offset in range(-3, variant_length_change + 4):
                        if strand == '+':
                            test_ss_pos = ss_pos_mut + offset
                        else:
                            test_ss_pos = ss_pos_mut - offset
                        # Get the mutant sequence at this position
                        try:
                            mut_seq = get_window_seq(test_ss_pos, strand, mut_chrom_seq)
                            MUT_SSU = predict_ssu(mut_seq, Models, ss_type)
                            MUT_SSUs.append(MUT_SSU)
                        except Exception as e:
                            continue  # Skip positions that cause errors

                    if MUT_SSUs:
                        # Select the maximum MUT_SSU value
                        MUT_SSU = max(MUT_SSUs)
                    else:
                        print(f"No valid MUT_SSU predictions for event {event}, variant {var_ID}")
                        continue
                else:
                    # Get the mutant sequence
                    mut_seq = get_window_seq(ss_pos_mut, strand, mut_chrom_seq)

                    # Predict SSU for mutant sequence using mut_seq
                    try:
                        MUT_SSU = predict_ssu(mut_seq, Models, ss_type)
                    except Exception as e:
                        print(f"Error predicting MUT SSU for event {event}, variant {var_ID}: {e}")
                        continue

                # Calculate delta SSU
                delta_SSU = MUT_SSU - WT_SSU

                # Write results
                writer.writerow([event, var_ID, chrom_event, pos_var, ref_var, alt_var, strand, ss_type, WT_SSU, MUT_SSU, delta_SSU])

if __name__ == "__main__":
    main()