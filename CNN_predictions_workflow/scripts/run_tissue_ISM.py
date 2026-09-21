from pkg_resources import resource_filename
from pangolin.model import *
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
import argparse
import csv

Genome = {}

def Genomictabulator(fasta):


    f = open(fasta)

    for chrfa in SeqIO.parse(f, "fasta"):
        Genome[chrfa.id] = chrfa.seq


    f.close()




def pangolin_PSI(seqs, strands, seq_names): 

    # Change this to the desired models. The model that each number corresponds to is listed below.
    #model_nums = [4]
    model_nums = [0, 2, 4, 6]
    model_names = ["Heart", "Liver", "Brain", "Testis"]

    # 0 = Heart, P(splice)
    # 1 = Heart, usage
    # 2 = Liver, P(splice)
    # 3 = Liver, usage
    # 4 = Brain, P(splice)
    # 5 = Brain, usage
    # 6 = Testis, P(splice)
    # 7 = Testis, usage

    # Change this to the desired sequences and strand for each sequence. If the sequence is N bases long, Pangolin will
    # return scores for the middle N-10000 bases (so if you are interested in the score for a single site, the input should
    # be: 5000 bases before the site, base at the site, 5000 bases after the site). Sequences < 10001 bases can be padded with 'N'.
    #seqs = [generate_random_dna_sequence(10001, 123)]
    # print(seqs)
    #seqs = [10001*'A']
    #strands = ['+', '+']

    # Load models
    models = []
    for i in model_nums:
        for j in range(1, 6):
            model = Pangolin(L, W, AR)
            if torch.cuda.is_available():
                model.cuda()
                weights = torch.load(resource_filename("pangolin","models/final.%s.%s.3" % (j, i)))
            else:
                weights = torch.load(resource_filename("pangolin","models/final.%s.%s.3" % (j, i)),
                                    map_location=torch.device('cpu'))
            model.load_state_dict(weights)
            model.eval()
            models.append(model)

    # Get scores

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


    results = []

    for i, seq in enumerate(seqs):
        seq = one_hot_encode(seq, strands[i]).T
        seq = torch.from_numpy(np.expand_dims(seq, axis=0)).float()

        if torch.cuda.is_available():
            seq = seq.to(torch.device("cuda"))

        for j, model_num in enumerate(model_nums):
            score = []
            # Average across 5 models
            for model in models[5*j:5*j+5]:
                with torch.no_grad():
                    score.append(model(seq)[0][INDEX_MAP[model_num],:].cpu().numpy())

            results.append([model_names[j], seq_names[i], np.mean(score, axis=0)[0], np.std(score, axis=0)[0]])

    return(results)



def main():

    parser = argparse.ArgumentParser(description='Run tissue ISM analysis.')
    parser.add_argument('--genome', required=True, help='Path to genome file')
    parser.add_argument('--mic-anno', required=True, help='Path to microexon annotation file')
    parser.add_argument('--event', required=True, help='Event identifier')
    parser.add_argument('--output', required=True, help='Output file path')
    args = parser.parse_args()

    # Rest of your script using args.genome, args.mic_anno, args.event, args.output

    Genomictabulator(args.genome)

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

    with open(args.output, 'w') as out:

        writer = csv.writer(out, delimiter="\t")
        writer.writerow(['ss', 'rel_pos', 'ref', 'alt', 'tissue', 'ss_out', 'WT_PSI_mean', 'WT_PSI_sd', 'MUT_PSI_mean', 'MUT_PSI_sd', 'dPSI', 'norm_dPSI'])

        C1, C2, ME_start, ME_end, chrom, strand, lengthDiff = human_event_coords[args.event]

        if strand=="+":
            seqs = [str(Genome[chrom][ME_start-5000:ME_start+5001]).upper(), str(Genome[chrom][ME_end-5000:ME_end+5001]).upper()]
            strands = ['+', '+']

        elif strand=="-":
            seqs = [str(Genome[chrom][ME_end-5000:ME_end+5001].reverse_complement()).upper(), str(Genome[chrom][ME_start-5000:ME_start+5001].reverse_complement()).upper()]
            strands = ['+', '+']

        WT_values = pangolin_PSI(seqs, strands, ["3ss", "5ss"])

        WT_res = dict()

        for row in WT_values:

            tissue, ss, PSI_mean, PSI_sd = row
            WT_res[(tissue, ss)] = (PSI_mean, PSI_sd)


        seq_names = ["3ss", "5ss"]
        #pos_range = [i for i in range(-300, 301) if i != 0]

        w = 300
        #w = 5

        #pos_range = list(range(-w-lengthDiff, w+1+lengthDiff))


        for i, seq in enumerate(seqs):

            ss = seq_names[i]

            if i==0:
                pos_range = list(range(-w, w+1+lengthDiff))
            elif i==1:
                pos_range = list(range(-w-lengthDiff, w+1))

            for rel_pos in pos_range:

                ref = seq[5000+rel_pos]
                alts = set(["A", "G", "C", "T"]) - set(ref)

                for alt in alts: 

                    seq_list = list(seq)
                    seq_list[5000+rel_pos] = alt

                    mut_seq = "".join(seq_list)

                    mut_values = pangolin_PSI([mut_seq], ["+"], [ss])

                    for row in mut_values:

                        tissue, ss_out, MUT_PSI_mean, MUT_PSI_sd = row

                        WT_PSI_mean, WT_PSI_sd = WT_res[(tissue, ss_out)]

                        dPSI = MUT_PSI_mean-WT_PSI_mean
                        norm_dPSI = (MUT_PSI_mean-WT_PSI_mean)/WT_PSI_mean

                        writer.writerow([ss, rel_pos, ref, alt, tissue, ss_out, WT_PSI_mean, WT_PSI_sd, MUT_PSI_mean, MUT_PSI_sd, dPSI, norm_dPSI])


if __name__ == "__main__":
    main()
