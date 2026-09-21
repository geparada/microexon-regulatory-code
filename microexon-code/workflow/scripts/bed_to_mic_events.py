import pandas as pd
import twobitreader


def liftover_coordinates(wt_events, up_bed, up_un_bed, dn_bed, dn_un_bed):
  mics = pd.read_table(
    wt_events,
    index_col="event"
  )

  up_int = pd.read_table(
    up_bed,
    header=None,
    names=["chrom", "upIntStart", "upIntEnd", "event", "score", "strand"],
    index_col="event",
    dtype={'upIntStart': int, 'upIntEnd': int}
  )
  up_un_int = pd.read_table(
    up_un_bed,
    header=None,
    names=["chrom", "upIntStart", "upIntEnd", "event", "score", "strand"],
    index_col="event",
    dtype={'upIntStart': int, 'upIntEnd': int},
    comment="#"
  )
  up_int = pd.concat([up_int, up_un_int], axis=0)

  dn_int = pd.read_table(
    dn_bed,
    header=None,
    names=["chrom", "dnIntStart", "dnIntEnd", "event", "score", "strand"],
    index_col="event",
    dtype={'dnIntStart': int, 'dnIntEnd': int}
  )
  dn_un_int = pd.read_table(
    dn_un_bed,
    header=None,
    names=["chrom", "dnIntStart", "dnIntEnd", "event", "score", "strand"],
    index_col="event",
    dtype={'dnIntStart': int, 'dnIntEnd': int},
    comment="#"
  )
  dn_int = pd.concat([dn_int, dn_un_int], axis=0)

  del mics["upIntStart"], mics["upIntEnd"], mics["dnIntStart"], mics["dnIntEnd"]

  mics = (mics.join(up_int[["upIntStart", "upIntEnd"]])
          .join(dn_int[["dnIntStart", "dnIntEnd"]]))

  return mics


def filter_mics(wt_events, mt_events, wt_genome, mt_genome, intron_window=300, exon_window=100):
  wt_genome = twobitreader.TwoBitFile(wt_genome)
  mt_genome = twobitreader.TwoBitFile(mt_genome)

  wt_events = pd.read_table(
    wt_events,
    index_col="event",
    converters=dict(
      upIntStart=lambda pos: int(pos) - 1,
      dnIntStart=lambda pos: int(pos) - 1
    )
  )

  mt_events_out = []
  skipped_events = []
  for event_id, mt_event in mt_events.iterrows():
    wt_event = wt_events.loc[event_id]


    if wt_event.strand == '+':
      if not mt_event.upIntStart < mt_event.upIntEnd < mt_event.dnIntStart < mt_event.dnIntEnd:
        skipped_events.append((event_id, "invalid event"))
        continue

      wt_c1 = wt_genome[wt_event.chrom][wt_event.upIntStart - exon_window:wt_event.upIntStart + intron_window].upper()
      wt_a = wt_genome[wt_event.chrom][wt_event.upIntEnd - intron_window:wt_event.dnIntStart + intron_window].upper()
      wt_c2 = wt_genome[wt_event.chrom][wt_event.dnIntEnd - intron_window:wt_event.dnIntEnd + exon_window].upper()

      mt_c1 = mt_genome[mt_event.chrom][mt_event.upIntStart - exon_window:mt_event.upIntStart + intron_window].upper()
      mt_a = mt_genome[mt_event.chrom][mt_event.upIntEnd - intron_window:mt_event.dnIntStart + intron_window].upper()
      mt_c2 = mt_genome[mt_event.chrom][mt_event.dnIntEnd - intron_window:mt_event.dnIntEnd + exon_window].upper()

    elif wt_event.strand == '-':
      if not mt_event.dnIntStart < mt_event.dnIntEnd < mt_event.upIntStart < mt_event.upIntEnd:
        skipped_events.append((event_id, "invalid event"))
        continue

      wt_c1 = wt_genome[wt_event.chrom][wt_event.upIntEnd - intron_window:wt_event.upIntEnd + exon_window].upper()
      wt_a = wt_genome[wt_event.chrom][wt_event.dnIntEnd-intron_window:wt_event.upIntStart+intron_window].upper()
      wt_c2 = wt_genome[wt_event.chrom][wt_event.dnIntStart-exon_window:wt_event.dnIntStart+intron_window].upper()

      mt_c1 = mt_genome[mt_event.chrom][mt_event.upIntEnd - intron_window:mt_event.upIntEnd + exon_window].upper()
      mt_a = mt_genome[mt_event.chrom][mt_event.dnIntEnd-intron_window:mt_event.upIntStart+intron_window].upper()
      mt_c2 = mt_genome[mt_event.chrom][mt_event.dnIntStart-exon_window:mt_event.dnIntStart+intron_window].upper()

    else:
      raise ValueError(wt_event.strand)

    if (wt_c1 != mt_c1) or (wt_a != mt_a) or (wt_c2 != mt_c2):
      mt_events_out.append(mt_event)
    else:
      skipped_events.append((event_id, "equal to reference"))

  return pd.DataFrame(mt_events_out), pd.DataFrame(skipped_events, columns=["event", "reason"]).set_index("event")


def main(wt_events, up_bed, up_un_bed, dn_bed, dn_un_bed,
         wt_genome, mt_genome, mt_events, skipped_events, intron_window, exon_window):
  df_mt_events = liftover_coordinates(wt_events, up_bed, up_un_bed, dn_bed, dn_un_bed)
  df_mt_events, df_skipped_events = filter_mics(wt_events, df_mt_events, wt_genome, mt_genome, intron_window, exon_window)

  df_mt_events.upIntStart += 1
  df_mt_events.dnIntStart += 1

  df_mt_events.to_csv(
    mt_events,
    sep="\t",
    index_label="event"
  )

  df_skipped_events.to_csv(
    skipped_events,
    index_label="event"
  )


if __name__ == "__main__":
  main(**snakemake.input, **snakemake.output, **snakemake.params)
