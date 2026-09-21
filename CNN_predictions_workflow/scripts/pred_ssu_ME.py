import os
import torch
import copy
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

# Function to load models once
def load_models():
    Models = [copy.deepcopy(model) for _ in default_model_paths]
    for m, b in zip(Models, default_model_paths):
        m.load_state_dict(torch.load(b))
    return Models

# Function to load the reference genome once
def load_genome(genome_name):
    reference_genome = Fasta(os.path.join(Fapath, genome_name + ".fa"))
    return reference_genome

# Prediction function
def predict_ssu(chrom, pos, strand, reference_genome, Models, ss):
    seq_start = pos - EL // 2
    seq_end = seq_start + EL + 1
    chrom_length = len(reference_genome[chrom])

    # Handle boundary conditions
    seq = reference_genome[chrom][max(seq_start, 0) : min(seq_end, chrom_length)]
    if seq_start < 0:
        seq = "N" * abs(seq_start) + seq
    if seq_end > chrom_length:
        seq = seq + "N" * (seq_end - chrom_length)
    seq = seq.upper()

    # Reverse complement if negative strand
    if strand == "-":
        seq = [repdict[base] for base in seq][::-1]
    seq = IN_MAP[[SeqTable[base] for base in seq]][:, :4]

    # Predict using the loaded models
    pred = 0
    for m in Models:
        pred += m.predict({"X": torch.tensor(seq)[None]}, use_ref=False)["single_pred_psi"]
    pred = (pred / len(Models))[0]
    pred = pred[pred.shape[0] // 2]

    # Return the predictions as a list
    #return [pred[1], pred[2]]

    # Return only the requested SSU
    if ss == "3ss":
        return pred[1]  # Acceptor SSU
    elif ss == "5ss":
        return pred[2]  # Donor SSU
    else:
        raise ValueError("Invalid ss parameter. Must be '3ss' or '5ss'.")

default_model_paths = ['EXTERNAL_DATA/DeltaSplice/' + x for x in default_model_paths ]

# Load models and genome once
Models = load_models()
reference_genome = load_genome('/scratch/ssd004/datasets/brainemo/Genome/hg38/hg38')  # Replace 'hg19' with your genome version

# Predict for a given chrom, position, and strand
chrom = 'chr10'
position = 102862159
strand = '+'

prediction_3ss = predict_ssu(chrom, position, strand, reference_genome, Models, '3ss')
prediction_5ss = predict_ssu(chrom, position, strand, reference_genome, Models, '5ss')

print(prediction_3ss, prediction_5ss)


#print(f"Acceptor SSU: {predictions[0]}, Donor SSU: {predictions[1]}")
