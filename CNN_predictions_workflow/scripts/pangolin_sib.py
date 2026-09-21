from pkg_resources import resource_filename
from pangolin.model import *
import argparse
import csv
import numpy as np
import torch
from Bio import SeqIO
from Bio.Seq import Seq

# Function to load models once
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

# Prediction function
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
        #std_score = np.std(score, axis=0)[0]
        results.append([tissue, mean_score])
    return results
    
def main():
    parser = argparse.ArgumentParser(description='Process sib event fastas with Pangolin')
    parser.add_argument('--sib_event', required=True, help='Path to genome FASTA file')
    parser.add_argument('--output', required=True, help='Output file path')
    args = parser.parse_args()

    # Load models once
    models = load_models()

    # Read the FASTA file
    with open(args.sib_event, "r") as file, open(args.output, 'w') as out:

        writer = csv.writer(out, delimiter="\t")
        writer.writerow(['sibbling', 'event', 'ss', 'tissue', 'ssu'])

        for record in SeqIO.parse(file, "fasta"):
            sib, event, ss = record.id.split("|")  # The ID of the sequence
            event_seq = str(record.seq)  # The sequence as a string
            sib_results = pangolin_PSI(event_seq, '+', models)

            for result in sib_results:
                tissue, ssu = result
                out_row = [sib, event, ss, tissue, ssu]

                writer.writerow(out_row)

if __name__ == "__main__":
    main()