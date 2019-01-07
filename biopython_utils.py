#!/usr/bin/env python
"""
functions, source code in Biopython
"""

import os
import sys

##----------------------------------------------------------------------------##
## Biopython BEGIN
##
def get_gc(s):
    """Code source: Bio.SeqUtils.GC
    calculate G+C content, return percentage (0-100)
    support ambiguous nucleotide S (G or C)

    >>> from Bio.SeqUtils import GC
    >>> GC('ACTGN')
    40.0
    """
    gc = sum(s.count(x) for x in ['G', 'g', 'C', 'c', 'S', 's'])
    try:
        return gc * 100.0 / len(s)
    except ZeroDivisionError:
        return 0.0


# from Bio.Data.IUPACData
"""Functions to get reverse complement sequence
    >>> from Bio.Seq import Seq
    >>> from Bio.Alphabet import generic_dna
    >>> my_dna = Seq('CCCCCgatA-G', generic_dna)
    >>> my_dna
    Seq('CCCCCgatA-G', DNAAlphabet())
    >>> my_dna.reverse_complement()
    Seq('C-TatcGGGGG', DNAAlphabet())
"""

ambiguous_rna_complement = {
    "A": "U",
    "C": "G",
    "G": "C",
    "U": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N",
    }


ambiguous_dna_complement = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N",
    }


def _maketrans(complement_mapping):
    """Make a python string translation table
    Arguments:
     - complement_mapping - a dictionary such as ambiguous_dna_complement
       and ambiguous_rna_complement from Data.IUPACData.

    Returns a translation table (a string of length 256) for use with the
    python string's translate method to use in a (reverse) complement.

    Compatible with lower case and upper case sequences.

    For internal use only.
    """
    before = ''.join(complement_mapping.keys())
    after = ''.join(complement_mapping.values())
    before += before.lower()
    after += after.lower()
    if sys.version_info[0] == 3:
        return str.maketrans(before, after)
    else:
        return string.maketrans(before, after)


_dna_complement_table = _maketrans(ambiguous_dna_complement)
_rna_complement_table = _maketrans(ambiguous_rna_complement)


def complement(s):
    """
    Return the complement sequence
    use pre-defined table
    """
    if 'U' in s or 'u' in s:
        ttable = _rna_complement_table
    else:
        ttable = _dna_complement_table
    return str(s).translate(ttable)


def reverse_complement(s):
    """
    code source: Bio.Seq.Seq
    Return the reverse complement sequence 
    """
    # Use -1 stride/step to reverse the complement
    return complement(s)[::-1]
