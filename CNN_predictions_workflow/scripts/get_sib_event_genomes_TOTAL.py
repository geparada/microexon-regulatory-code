
import csv
from collections import defaultdict 
import glob
import math
from Bio import SeqIO
from Bio.Seq import Seq
from snakemake.utils import min_version

mic_anno = snakemake.input["mic_anno"]
#sib = snakemake.params["sib"]
#output_file = snakemake.output["out_deltasplice"]
#output_file = snakemake.output["out_pangolin"]


Genome = {}

def Genomictabulator(fasta):

    f = open(fasta)

    for chrfa in SeqIO.parse(f, "fasta"):
        Genome[chrfa.id] = chrfa.seq

    f.close()

Genomictabulator(snakemake.input["Genome"])

human_event_coords = dict()

with open(mic_anno) as file:
    
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
        
        human_event_coords[row["event"]] = (C1, C2, ME_start, ME_end-1, chrom, strand, int(lengthDiff))


def apply_variants(chromosome_seq, variants, input_coordinates):
    # variants: set of strings 'chrom:pos:alt:ref'
    # input_coordinates: list of 0-based coordinates, they all must be different numbers
    # chromosome_seq: reference sequence as a string (1-based indexing in VCF)
    #
    # Steps:
    # 1. Parse variants and sort by position.
    # 2. Check no overlaps (now prints a warning and skips them instead of raising an error).
    # 3. Validate ref alleles.
    # 4. Apply variants by slicing strings.
    # 5. Update input_coordinates according to indels.

    input_var_offsets = defaultdict(int)

    # Parse and sort
    parsed = []
    for v in variants:
        chrom, pos, ref, alt = v.split(":")
        pos = int(pos) - 1
        parsed.append((chrom, pos, ref, alt))
    parsed.sort(key=lambda x: x[1])

    # Check overlap:
    filtered_parsed = []
    last_end = -1
    for (chrom, pos, ref, alt) in parsed:
        if "," in alt:
            alt = alt.split(",")[0]
        start = pos
        end = pos + len(ref)
        if start < last_end:
            # Print warning and skip the overlapping variant
            print(f"WARNING: Overlapping variant detected: {chrom}:{pos+1}:{ref}:{alt} . Skipped.")
            continue
        filtered_parsed.append((chrom, pos, ref, alt))
        last_end = end

    # We only work with non-overlapping variants now
    parsed = filtered_parsed

    # Validate and apply
    mutated_seq = chromosome_seq
    offset = 0
    for (chrom, pos, ref, alt) in parsed:
        # Convert pos to 0-based for string indexing
        orig_start = pos
        adj_start = orig_start + offset
        adj_end = adj_start + len(ref)
        # Check ref
        found_ref = mutated_seq[adj_start:adj_end]
        if found_ref != ref:
            raise ValueError(
                f"Ref mismatch at {chrom}:{pos} expecting {ref}, but {found_ref} was found, "
                f"offset: {offset}, context: {mutated_seq[adj_start-3:adj_start+3]}, "
                f"original len: {len(chromosome_seq)}, mutated len: {len(mutated_seq)}"
            )
        # Apply variant
        mutated_seq = mutated_seq[:adj_start] + alt + mutated_seq[adj_end:]
        # Update offsets in input_coordinates
        diff = len(alt) - len(ref)
        if diff != 0:
            for i, in_cord in enumerate(input_coordinates):
                if (pos + 1) <= in_cord:
                    input_var_offsets[in_cord]+=diff 

        offset += diff

    output_coordinates = [  x + input_var_offsets[x] for x in input_coordinates ]

    return mutated_seq, output_coordinates

event_sibs = defaultdict(set)

with open(snakemake.input["total_event_vars"]) as file:

    reader = csv.DictReader(file, delimiter=",")

    for row in reader:

        if row['sample_ID'] == snakemake.params["sib"]:

            for var in row['vars'].split("|"):

                event_sibs[( row['event'], row['sample_ID'])].add(var)


def get_seq(chrom_seq, strand, ss_pos, EL):
    #EL = 30000
    seq_start = ss_pos - EL // 2
    seq_end = seq_start + EL + 1

    seq = chrom_seq[max(seq_start, 0):min(seq_end, len(chrom_seq))]
    if seq_start < 0:
        seq = 'N' * abs(seq_start) + seq
    if seq_end > len(chrom_seq):
        seq = seq + 'N' * (seq_end - len(chrom_seq))

    if strand == '-':
        seq = str(Seq(seq).reverse_complement())

    return seq


with open(snakemake.output["out_deltasplice"], "w") as fa_out_ds, open(snakemake.output["out_pangolin"], "w") as fa_out_p:

    for event, sib in event_sibs:

        vars = event_sibs[(event, sib)]
        C1, C2, ME_start, ME_end, chrom, strand, lengthDiff = human_event_coords[event]
        ME = [ME_start, ME_end]

        mut_chrom, lifted_ME = apply_variants(str(Genome[chrom]).upper(), vars, ME)

        for i, ss in enumerate(lifted_ME):
            if strand == "+":
                if i == 0:
                    out_seq_deltasplice = get_seq(mut_chrom, strand, ss, 30000)
                    out_seq_pangolin = get_seq(mut_chrom, strand, ss, 10000)
                    seq_id = "|".join([sib, event, '3ss'])
                elif i == 1:
                    out_seq_deltasplice = get_seq(mut_chrom, strand, ss, 30000)
                    out_seq_pangolin = get_seq(mut_chrom, strand, ss, 10000)
                    seq_id = "|".join([sib, event, '5ss'])
            else:  # strand == "-"
                if i == 1:
                    out_seq_deltasplice = get_seq(mut_chrom, strand, ss, 30000)
                    out_seq_pangolin = get_seq(mut_chrom, strand, ss, 10000)
                    seq_id = "|".join([sib, event, '3ss'])
                elif i == 0:
                    out_seq_deltasplice = get_seq(mut_chrom, strand, ss, 30000)
                    out_seq_pangolin = get_seq(mut_chrom, strand, ss, 10000)
                    seq_id = "|".join([sib, event, '5ss'])

            fa_out_ds.write(f">{seq_id}\n")
            fa_out_ds.write(f"{out_seq_deltasplice}\n")

            fa_out_p.write(f">{seq_id}\n")
            fa_out_p.write(f"{out_seq_pangolin}\n")

