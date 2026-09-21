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
from Bio import SeqIO

# Adjust the default model paths if necessary
default_model_paths = ['EXTERNAL_DATA/DeltaSplice/' + x for x in default_model_paths]

# Function to load models once
def load_models():
    Models = [copy.deepcopy(model) for _ in default_model_paths]
    for m, b in zip(Models, default_model_paths):
        m.load_state_dict(torch.load(b))
    return Models

# Prediction function
def predict_ssu(seq, Models, ss):
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
    parser = argparse.ArgumentParser(description='Process sib event fastas with DeltaSplice')
    parser.add_argument('--sib_event', required=True, help='Path to genome FASTA file')
    parser.add_argument('--output', required=True, help='Output file path')
    args = parser.parse_args()

    # Load models
    Models = load_models()

    # Read the FASTA file
    with open(args.sib_event, "r") as file, open(args.output, 'w') as out:

        writer = csv.writer(out, delimiter="\t")
        writer.writerow(['sibbling', 'event', 'ss', 'ssu'])

        for record in SeqIO.parse(file, "fasta"):
            sib, event, ss = record.id.split("|")  # The ID of the sequence
            event_seq = str(record.seq)  # The sequence as a string
            event_seq = IN_MAP[[SeqTable[base] for base in event_seq]][:, :4]

            ssu = predict_ssu(event_seq, Models, ss)
            out_row = [sib, event, ss, ssu]

            writer.writerow(out_row)

if __name__ == "__main__":
    main()




