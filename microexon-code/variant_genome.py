import twobitreader
from collections import namedtuple
from typing import List
from utils import reverse_complement

Variant = namedtuple('Variant', ['chromosome', 'pos', 'ref', 'alt'])
Interval = namedtuple("Interval", ["chrom", "strand", "start", "end", "anchor", "anchor_offset"])


def variant_from_str(variant_str):
  chrom, pos, ref, alt = variant_str.split(':')
  pos = int(pos) - 1
  return Variant(chrom, pos, ref, alt)


class _VariantChromosome:
  def __init__(self, genome: twobitreader.TwoBitFile, chromosome: str, variants: List[Variant]):
    self.genome = genome
    self.chromosome = chromosome
    self.variants = variants

  def __getitem__(self, item):
    if isinstance(item, int):
      return self._getseq(item, item+1)
    elif isinstance(item, slice):
      if not (item.step == 1 or item.step is None):
        raise ValueError("Step size of slice must be 1.")
      return self._getseq(item.start, item.stop)

  def _getseq(self, start, stop):
    seq = bytearray(self.genome[self.chromosome][start:stop].upper(), 'utf8')
    len_seq = stop - start

    for variant in self.variants:
      if variant.chromosome != self.chromosome:
        continue

      rel_pos = variant.pos - start

      if rel_pos < 0 or rel_pos >= len_seq:
        continue

      assert seq[rel_pos] == ord(variant.ref)
      seq[rel_pos] = ord(variant.alt)

    return seq.decode()


class VariantGenome:
  def __init__(self, genome: twobitreader.TwoBitFile, variants: List[Variant]):
    self.genome = genome
    self.variants = variants

  def __getitem__(self, chromosome):
    return _VariantChromosome(self.genome, chromosome, self.variants)

  def get_seq(self, interval: Interval):
    return apply_variants(self.genome, self.variants, interval)


def apply_variants(sequence, variants: List[Variant], interval: Interval):
    """Apply variants to a sequence interval.

    This function ignores the interval's strand attribute and always returns
    the sequence on the reference strand.

    Parameters
    ----------
    sequence : sequence-like
        Can be any object that supports the slice interface
        (`sequence[start:end]`), i.e. string, TwoBitFile, bytearray,
        etc.
    variants : list-like
        A list of variants of the form ``chromosome:position:ref:alt``.
    interval : Interval
        An anchored interval/coordinate for which to extract the sequence.
    anchor : int
        Sequence position that should stay fixed when length-changing
        variants are applied.

    Returns
    -------
    str
        Variant sequence
    """

    variants = _preprocess_variants(interval.chrom, variants)

    start = interval.start
    end = interval.end
    anchor = interval.anchor
    anchor_offset = interval.anchor_offset

    if anchor_offset:
      start += anchor_offset
      end += anchor_offset

    if anchor is None:
      var_sequence = _apply_variants_no_anchor(sequence, variants, interval)
    elif start < anchor < end:
      interval_left = Interval(
        interval.chrom, interval.strand,
        start, anchor, None, None
      )

      interval_right = Interval(
        interval.chrom, interval.strand,
        anchor, end, None, None
      )

      var_sequence_left = _apply_variants_right_anchor(
        sequence, variants, interval_left)
      var_sequence_right = _apply_variants_left_anchor(
        sequence, variants, interval_right)

      var_sequence = var_sequence_left + var_sequence_right

    elif anchor == start:
      var_sequence = _apply_variants_left_anchor(
        sequence, variants, interval)

    elif anchor == end:
      var_sequence = _apply_variants_right_anchor(
        sequence, variants, interval)
    elif anchor > end:
      # This method is inefficient since it takes extracts all the sequence
      # up to the anchor.
      start = interval.start
      end = interval.end

      tmp_interval = Interval(
        interval.chrom, interval.strand,
        start, anchor, None, None
      )

      var_sequence = _apply_variants_right_anchor(
        sequence, variants, tmp_interval)

      var_sequence = var_sequence[:end - start]
    elif anchor < start:
      start = interval.start
      end = interval.end

      tmp_interval = Interval(
        interval.chrom, interval.strand,
        anchor, end, None, None
      )

      var_sequence = _apply_variants_left_anchor(
        sequence, variants, tmp_interval
      )

      var_sequence = var_sequence[-(end - start):]
    else:
      raise AssertionError("Should never be reached")  # pragma: no cover
      # raise ValueError("Anchor is outside of interval.")

    if interval.strand == '-':
      var_sequence = reverse_complement(var_sequence)

    return var_sequence


def _preprocess_variants(chromosome, variants):
  variants_split = []
  for variant in variants:
    ref = variant.ref.encode()
    alt = variant.alt.encode()
    if variant.chromosome != chromosome:
      continue

    if ref in b'-.*':
      ref = b""
    if alt in b'-.*':
      alt = b""

    variants_split.append((variant.pos, ref, alt))

  variants_split = sorted(variants_split, reverse=True)
  return variants_split


def _apply_variants_no_anchor(dna, variants, interval):
  start = interval.start
  end = interval.end

  # v_interval = Interval(interval.chrom, '+', start, end, interval.refg, None)
  # variant_sequence = bytearray(dna[v_interval.chrom][v_interval.start:v_interval.end])
  variant_sequence = bytearray(dna[interval.chrom][start:end].upper().encode())

  # Apply the variants
  for v_start, ref, alt in variants:
    v_start_rel = v_start - start
    v_end_rel = v_start_rel + len(ref)

    if v_end_rel < 0 or v_start_rel >= len(variant_sequence):
      continue

    if v_start_rel < 0:
      # Clip left side of variant
      v_offset = start - v_start
      ref = ref[v_offset:]
      alt = alt[v_offset:]
      v_start_rel = 0

    if v_end_rel > len(variant_sequence):
      # Clip right side of variant
      v_offset = end - v_start
      ref = ref[:v_offset]
      alt = alt[:v_offset]
      v_end_rel = len(variant_sequence)

    try:
      assert v_end_rel >= 0
      assert v_start_rel >= 0
    except AssertionError:
      print(f"Assertion error in variant: {interval.chrom}:{v_start+1}:{ref}:{alt}")
      raise
    # assert variant_sequence[v_start_rel:v_end_rel].upper() == ref.upper()
    variant_sequence[v_start_rel:v_end_rel] = alt

  return variant_sequence.decode()


def _apply_variants_right_anchor(dna, variants, interval):
  start = interval.start
  end = interval.end

  variants_processed = []
  for v_start, ref, alt in variants:
    v_end = v_start + len(ref)

    if (v_start >= end and v_end > end) or v_end <= start:
      # Variant falls outside of interval
      continue

    # Clip variant if it overlaps the anchor
    if v_start < end <= v_end:
      v_len = end - v_start
      ref = ref[:v_len]
      alt = alt[:v_len]

    # Clip if variant overlaps the start
    if v_start - len(alt) + len(ref) < start:
      v_offset = start - (v_start - len(alt) + len(ref))
      ref = ref[v_offset:]
      alt = alt[v_offset:]
      v_start = start

    # Move the start to accommodate the variant
    start += min(len(alt) - len(ref), 0)
    # start = min(v_start, start)
    variants_processed.append((v_start, ref, alt))

  # v_interval = _interval.Interval(interval.chrom, '+', start, end, interval.refg)
  # variant_sequence = bytearray(dna(v_interval))
  variant_sequence = bytearray(dna[interval.chrom][start:end].upper().encode())

  # Apply the variants
  for v_start, ref, alt in variants_processed:
    v_start_rel = v_start - start
    v_end_rel = v_start_rel + len(ref)

    assert v_end_rel >= 0
    assert v_start_rel >= 0
    # assert variant_sequence[v_start_rel:v_end_rel].upper() == ref.upper()
    variant_sequence[v_start_rel:v_end_rel] = alt

  # Clip sequence if it's still too long

  len_interval = interval.end - interval.start
  variant_sequence = variant_sequence[-len_interval:]
  assert len(variant_sequence) == len_interval
  return variant_sequence.decode()


def _apply_variants_left_anchor(dna, variants, interval):
  start = interval.start
  end = interval.end

  variants_processed = []
  for v_start, ref, alt in variants:
    v_end = v_start + len(ref)

    if v_end <= start or v_start >= end:
      # Variant falls outside of interval
      continue

    if v_start <= start < v_end:
      # Variant overlaps the anchor
      v_offset = start - v_start
      ref = ref[v_offset:]
      alt = alt[v_offset:]
      v_start = start

    # Clip if variant overlaps the end
    if v_start + len(alt) - len(ref) > end:
      # v_clip = v_start + len(alt) - len(ref) - end
      # ref = ref[:-v_clip]
      # alt = alt[:-v_clip]
      v_clip = end - v_start
      ref = ref[:v_clip]
      alt = alt[:v_clip]

    variants_processed.append((v_start, ref, alt))
    end -= min(len(alt) - len(ref), 0)

  # v_interval = _interval.Interval(interval.chrom, '+', start, end, interval.refg)
  # variant_sequence = bytearray(dna(v_interval))
  variant_sequence = bytearray(dna[interval.chrom][start:end].upper().encode())

  for v_start, ref, alt in variants_processed:
    v_start_rel = v_start - start
    v_end_rel = v_start_rel + len(ref)

    assert v_start_rel <= len(variant_sequence)

    # if v_end_rel > len(variant_sequence):
    #   # Clip variant if it extends beyond the sequence
    #   v_clip = len(variant_sequence) - v_start_rel
    #   ref = ref[:v_clip]
    #   alt = alt[:v_clip]
    #   v_end_rel = len(variant_sequence)

    # assert variant_sequence[v_start_rel:v_end_rel].upper() == ref.upper()
    variant_sequence[v_start_rel:v_end_rel] = alt

  len_interval = interval.end - interval.start
  variant_sequence = variant_sequence[:len_interval]
  assert len(variant_sequence) == len_interval
  return variant_sequence.decode()


def get_variant_genome_from_series(event, wt_genome):
  event.rename({"variants": "variant"}, inplace=True)

  if 'variant' in event:
    variants = event['variant'].lstrip('[').rstrip(']')
    variants = [variant_from_str(v) for v in set(variants.split(',')) if len(v.split(":")) == 4]
  else:
    variants = []
  return VariantGenome(wt_genome, variants)


def get_variant_genome_from_dict(event, wt_genome):
  if 'variants' in event:
    event['variant'] = event['variants']

  if 'variant' in event:
    variants = event['variant'].lstrip('[').rstrip(']')
    variants = [variant_from_str(v) for v in set(variants.split(',')) if len(v.split(":")) == 4]
  else:
    variants = []
  return VariantGenome(wt_genome, variants)
