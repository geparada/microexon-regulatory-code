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
    Fapath,
    EL,
    SeqTable,
    repdict,
    IN_MAP,
)
from Bio import SeqIO
from Bio.Seq import Seq

default_model_paths = ['EXTERNAL_DATA/DeltaSplice/' + x for x in default_model_paths ]

# Function to load models once
def load_models():
    Models = [copy.deepcopy(model) for _ in default_model_paths]
    for m, b in zip(Models, default_model_paths):
        m.load_state_dict(torch.load(b))
    return Models

# Function to load the reference genome once
def load_genome(genome_name):
    reference_genome = Fasta(genome_name)
    return reference_genome

# Prediction function with an additional parameter 'ss'
#def predict_ssu(chrom, pos, strand, reference_genome, Models, ss):
def predict_ssu(seq, Models, ss):
    # seq_start = pos - EL // 2
    # seq_end = seq_start + EL + 1
    # chrom_length = len(reference_genome[chrom])

    # # Handle boundary conditions
    # seq = reference_genome[chrom][max(seq_start, 0) : min(seq_end, chrom_length)]
    # if seq_start < 0:
    #     seq = "N" * abs(seq_start) + seq
    # if seq_end > chrom_length:
    #     seq = seq + "N" * (seq_end - chrom_length)
    # seq = seq.upper()

    # # Reverse complement if negative strand
    # if strand == "-":
    #     seq = [repdict[base] for base in seq][::-1]
    # seq = IN_MAP[[SeqTable[base] for base in seq]][:, :4]

    # Predict using the loaded models
    pred = 0
    for m in Models:
        pred += m.predict({"X": torch.tensor(seq)[None]}, use_ref=False)["single_pred_psi"]
    pred = (pred / len(Models))[0]
    pred = pred[pred.shape[0] // 2]

    # Return only the requested SSU
    if ss == "3ss":
        return pred[1]  # Acceptor SSU
    elif ss == "5ss":
        return pred[2]  # Donor SSU
    else:
        raise ValueError("Invalid ss parameter. Must be '3ss' or '5ss'.")

def main():
    parser = argparse.ArgumentParser(description='Run DeltaSplice ISM analysis.')
    parser.add_argument('--genome', required=True, help='Path to genome FASTA file')
    parser.add_argument('--mic-anno', required=True, help='Path to microexon annotation file')
    parser.add_argument('--event', required=True, help='Event identifier')
    parser.add_argument('--output', required=True, help='Output file path')
    args = parser.parse_args()

    # Load models and genome
    Models = load_models()
    reference_genome = load_genome(args.genome)

    # Read microexon annotation
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

            human_event_coords[row["event"]] = (C1, C2, ME_start, ME_end - 1, chrom, strand, lengthDiff)

    # Get event coordinates
    C1, C2, ME_start, ME_end, chrom, strand, lengthDiff = human_event_coords[args.event]

    # Define window size
    w = 300  # Adjust window size as needed

    # Initialize output
    with open(args.output, 'w') as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(['ss', 'rel_pos', 'ref', 'alt', 'WT_SSU', 'MUT_SSU', 'delta_SSU'])

        # For each splice site (3ss and 5ss)
        for i, ss_type in enumerate(['3ss', '5ss']):
            # Determine the splice site position
            if ss_type == '3ss':
                if strand == '+':
                    ss_pos = ME_start
                else:
                    ss_pos = ME_end
            elif ss_type == '5ss':
                if strand == '+':
                    ss_pos = ME_end
                else:
                    ss_pos = ME_start

            # Extract sequence around the splice site
            # seq_start = ss_pos - w    # incorrect
            # seq_end = ss_pos + w + 1
            seq_start = ss_pos - EL // 2  # The perseption field window needs to be determined by the EL variable
            seq_end = seq_start + EL + 1

            seq = reference_genome[chrom][max(seq_start, 0):min(seq_end, len(reference_genome[chrom]))]
            if seq_start < 0:
                seq = 'N' * abs(seq_start) + seq
            if seq_end > len(reference_genome[chrom]):
                seq = seq + 'N' * (seq_end - len(reference_genome[chrom]))
            

            # Reverse complement if negative strand
            if strand == '-':
                seq = str(Seq(''.join(seq)).reverse_complement())

            seq = seq.upper()
            wt_seq = IN_MAP[[SeqTable[base] for base in seq]][:, :4]

            # Predict SSU for wild-type sequence
            #WT_SSU = predict_ssu(chrom, ss_pos, strand, reference_genome, Models, ss_type)
            WT_SSU = predict_ssu(wt_seq, Models, ss_type)

            if i==0:
                pos_range = list(range(-w, w+1+lengthDiff))
            elif i==1:
                pos_range = list(range(-w-lengthDiff, w+1))

            # Perform mutations within the window
            #for rel_pos in range(-w, w + 1): #incorrect
            for rel_pos in pos_range:
                abs_pos = ss_pos + rel_pos if strand == '+' else ss_pos - rel_pos
                if abs_pos < 0 or abs_pos >= len(reference_genome[chrom]):
                    continue

                ref_base = seq[(EL // 2) + rel_pos]  # Index adjusted for 0-based Python indexing
                for alt_base in set(['A', 'C', 'G', 'T']) - set(ref_base):
                    # Create mutant sequence
                    seq_list = list(seq)
                    seq_list[(EL // 2) + rel_pos] = alt_base
                    mut_seq = ''.join(seq_list)

                    # Predict SSU for mutant sequence
                    #MUT_SSU = predict_ssu(chrom, ss_pos, strand, reference_genome, Models, ss_type)
                    mut_seq = IN_MAP[[SeqTable[base] for base in mut_seq]][:, :4]

                    MUT_SSU = predict_ssu(mut_seq, Models, ss_type)

                    # Calculate delta SSU
                    delta_SSU = MUT_SSU - WT_SSU

                    # Write results
                    writer.writerow([ss_type, rel_pos, ref_base, alt_base, WT_SSU, MUT_SSU, delta_SSU])

if __name__ == "__main__":
    main()
