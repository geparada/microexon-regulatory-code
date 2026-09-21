import numpy as np


def reverse_complement(seq):
  try:
    assert set(seq).issubset(set('ACGTN'))
  except AssertionError:
    raise ValueError(
      "Sequence {} contains invalid values: {}"
        .format(seq, set(seq) - set('ACGTN'))
    )
  return ''.join('TGCAN'['ACGTN'.index(s)] for s in seq.upper()[::-1])


def converter_str(high_val, low_val):
  def f(x):
    if x == low_val:
      return 0
    elif x == 'NA':
      return 1
    elif x == high_val:
      return 2
    else:
      raise ValueError

  return f


def encode_seq(seq):
  """Fast one-hot encoding of DNA sequence using bit-fiddling magic.
  Only works for uppercase letters A, C, G, T, and N.
  N is mapped to 6.
  https://medium.com/@h_76213/efficient-dna-embedding-with-tensorflow-ffce5f499083
  """
  seq_decoded = np.frombuffer(seq.encode(), dtype=np.int8).copy()
  seq_decoded &= ~(1 << 6 | 1 << 4)
  seq_decoded >>= 1
  seq_decoded ^= (seq_decoded & 1 << 1) >> 1
  assert np.all(np.isin(seq_decoded, [0, 1, 2, 3, 6]))
  return seq_decoded


def one_hot(seq, dtype=np.int64, seqweaver=False):
  assert np.all(np.isin(seq, [0, 1, 2, 3, 6]))
  one_hot_seq = np.zeros(seq.shape + (7,), dtype)
  np.put_along_axis(one_hot_seq, np.expand_dims(seq, -1), 1, seq.ndim)

  if seqweaver:
    one_hot_seq[:, [1, 2]] = one_hot_seq[:, [2, 1]]
  return one_hot_seq[:, :4]


def get_sequence_for_event(event, genome, window=150):
  if event.strand == '+':
    up_seq = genome[event.chrom][event.upIntEnd-window:event.upIntEnd].upper()
    ex_seq = genome[event.chrom][event.upIntEnd:event.dnIntStart].upper()
    dn_seq = genome[event.chrom][event.dnIntStart:event.dnIntStart+window].upper()
  else:
    up_seq = reverse_complement(genome[event.chrom][event.upIntStart:event.upIntStart+window].upper())
    ex_seq = reverse_complement(genome[event.chrom][event.dnIntEnd:event.upIntStart].upper())
    dn_seq = reverse_complement(genome[event.chrom][event.dnIntEnd-window:event.dnIntEnd].upper())

  return up_seq, ex_seq, dn_seq


def compute_variant_distance(event, variant):
  pos = int(variant.split(":")[1]) - 1

  if event.strand == "+":
    if pos < event.upIntEnd:
      return (pos - event.upIntEnd, "up")
    elif event.upIntEnd <= pos < event.dnIntStart:
      if (pos - event.upIntEnd) < (event.dnIntStart - pos):
        return (pos - event.upIntEnd, "ex")
      else:
        return (pos - event.dnIntStart, "ex")
    elif event.dnIntStart <= pos:
      return (pos - event.dnIntStart, "dn")
    else:
      raise ValueError
  elif event.strand == "-":
    if pos < event.dnIntEnd:
      return (event.dnIntEnd - pos - 1, "dn")
    elif event.dnIntEnd <= pos < event.upIntStart:
      if (event.upIntStart - pos) < (pos - event.dnIntEnd):
        return (event.upIntStart - pos - 1, "ex")
      else:
        return (event.dnIntEnd - pos - 1, "ex")
    elif event.upIntStart <= pos:
      return (event.upIntStart - pos - 1, "up")
