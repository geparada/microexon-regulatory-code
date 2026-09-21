import pandas as pd
import numpy as np


def main(wt_pred_file, mt_pred_file, output_file):
    wt_preds = pd.read_csv(wt_pred_file, index_col="event")
    mt_preds = pd.read_csv(mt_pred_file, index_col="event")

    event, rel_pos, subst = zip(*(i.split("_") for i in mt_preds.index))
    mt_preds.index = event
    ref, alt = zip(*(s.split(":") for s in subst))
    region, pos = zip(*((p[:2], int(p[2:])) for p in rel_pos))

    mt_preds["region"] = region
    mt_preds["rel_pos"] = pos
    mt_preds["ref"] = ref
    mt_preds["alt"] = alt

    mt_preds.rename(dict(neural_high_pred="mt_psi"), axis="columns", inplace=True)
    wt_preds.rename(dict(neural_high_pred="wt_psi"), axis="columns", inplace=True)

    dpsi = mt_preds.join(wt_preds.wt_psi)
    dpsi["dpsi"] = dpsi.mt_psi - dpsi.wt_psi
    dpsi["fold_change"] = np.log2(dpsi.mt_psi) - np.log2(dpsi.wt_psi)

    dpsi.to_csv(output_file, index_label="event")


if __name__ == "__main__":
    main(snakemake.input["wt_pred_file"], snakemake.input["mt_pred_file"], snakemake.output[0])
