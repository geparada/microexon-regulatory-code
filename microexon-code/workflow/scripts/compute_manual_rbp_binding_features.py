"""Computes PWM scores of several, manually defined RBPs"""

#%% Imports
import pandas as pd
from Bio import SeqUtils
from tqdm import tqdm
import twobitreader
from utils import reverse_complement
from variant_genome import VariantGenome, variant_from_str, Interval, get_variant_genome_from_series

#%% Manual RBP binding sites
rbp_motifs = dict(
    nova="YCAY", celf="YGCTYK", star="YTAAY", srrm4="TGC", srsf11=("TCTCT", "CTCTC")
)


def rbp_counts_for_event(event, event_id, genome, rbp_motifs):
    if event.strand == "+":
        up300_seq = genome.get_seq(
            Interval(
                event.chrom,
                event.strand,
                event.upIntEnd - 300,
                event.upIntEnd,
                event.upIntEnd,
                None,
            )
        )
        exon_seq = genome.get_seq(
            Interval(
                event.chrom, event.strand, event.upIntEnd, event.dnIntStart, event.upIntEnd, None
            )
        )
        dn300_seq = genome.get_seq(
            Interval(
                event.chrom,
                event.strand,
                event.dnIntStart,
                event.dnIntStart + 300,
                event.dnIntStart,
                None,
            )
        )
    else:
        up300_seq = genome.get_seq(
            Interval(
                event.chrom,
                event.strand,
                event.upIntStart,
                event.upIntStart + 300,
                event.upIntStart,
                None,
            )
        )
        exon_seq = genome.get_seq(
            Interval(
                event.chrom, event.strand, event.dnIntEnd, event.upIntStart, event.upIntStart, None
            )
        )
        dn300_seq = genome.get_seq(
            Interval(
                event.chrom,
                event.strand,
                event.dnIntEnd - 300,
                event.dnIntEnd,
                event.dnIntEnd,
                None,
            )
        )
    counts = pd.Series(
        name=event_id,
        data={
            "RBP.cov_Nova_exon": len(SeqUtils.nt_search(exon_seq, rbp_motifs["nova"]))
            - 1,
            "RBP.cov_Celf_exon": len(SeqUtils.nt_search(exon_seq, rbp_motifs["celf"]))
            - 1,
            "RBP.cov_Star_exon": len(SeqUtils.nt_search(exon_seq, rbp_motifs["star"]))
            - 1,
            "RBP.cov_Srrm4_exon": len(SeqUtils.nt_search(exon_seq, rbp_motifs["srrm4"]))
            - 1,
            "RBP.cov_Srsf11_exon": sum(
                len(SeqUtils.nt_search(exon_seq, m)) - 1 for m in rbp_motifs["srsf11"]
            ),
            "RBP.cov_Nova_up300": len(SeqUtils.nt_search(up300_seq, rbp_motifs["nova"]))
            - 1,
            "RBP.cov_Celf_up300": len(SeqUtils.nt_search(up300_seq, rbp_motifs["celf"]))
            - 1,
            "RBP.cov_Star_up300": len(SeqUtils.nt_search(up300_seq, rbp_motifs["star"]))
            - 1,
            "RBP.cov_Srrm4_up300": len(
                SeqUtils.nt_search(up300_seq, rbp_motifs["srrm4"])
            )
            - 1,
            "RBP.cov_Srsf11_up300": sum(
                len(SeqUtils.nt_search(up300_seq, m)) - 1 for m in rbp_motifs["srsf11"]
            ),
            "RBP.cov_Nova_dn300": len(SeqUtils.nt_search(dn300_seq, rbp_motifs["nova"]))
            - 1,
            "RBP.cov_Celf_dn300": len(SeqUtils.nt_search(dn300_seq, rbp_motifs["celf"]))
            - 1,
            "RBP.cov_Star_dn300": len(SeqUtils.nt_search(dn300_seq, rbp_motifs["star"]))
            - 1,
            "RBP.cov_Srrm4_dn300": len(
                SeqUtils.nt_search(dn300_seq, rbp_motifs["srrm4"])
            )
            - 1,
            "RBP.cov_Srsf11_dn300": sum(
                len(SeqUtils.nt_search(dn300_seq, m)) - 1 for m in rbp_motifs["srsf11"]
            ),
        },
    )
    return counts


def manual_binding_features(exon_file, genome_file, output_file):
    wt_genome = twobitreader.TwoBitFile(genome_file)
    exons = pd.read_csv(
        exon_file,
        sep="\t",
        index_col="event",
        converters={
            "upIntStart": lambda pos: int(pos) - 1,
            "dnIntStart": lambda pos: int(pos) - 1,
        },
    )
    motif_counts = []

    for event_id, event in tqdm(exons.iterrows(), total=exons.shape[0]):
      mt_genome = get_variant_genome_from_series(event, wt_genome)
      counts = rbp_counts_for_event(event, event_id, mt_genome, rbp_motifs)
      motif_counts.append(counts)
    motif_counts = pd.DataFrame(motif_counts)
    motif_counts.to_csv(output_file, index_label="event")


if __name__ == "__main__":
    manual_binding_features(snakemake.input['exons'], snakemake.input['genome'], snakemake.output[0])