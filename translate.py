#! /usr/bin/env python3

import sys

def translate_sequence(rna_sequence, genetic_code):
    """Translates a sequence of RNA into a sequence of amino acids.

    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the translated amino acids.
    """

    seq = rna_sequence.upper()
    AA = ''
    for i in range(0, len(rna_sequence), 3):
        codon = seq[i:i+3]
        if codon in genetic_code:
            if genetic_code[codon] == '*':
                break
            else:
                AA += genetic_code[codon]
        else:
            break

    return AA

def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    list
        A list of strings; each string is an sequence of amino acids encoded by
        `rna_sequence`.
    """
    translations = []

    for frame in range(3):
        sequence = rna_sequence[frame:]
        if len(sequence) % 3 != 0:
            sequence = sequence[:-(len(sequence) % 3)]
        amino_acids = ''

        for i in range(0, len(sequence), 3):
            codon = sequence[i:i+3].upper()
            if codon == 'AUG':
                amino_acids += genetic_code[codon]
                for b in range(i+3, len(sequence), 3):
                    codon = sequence[b:b+3].upper()
                    if codon in genetic_code and genetic_code[codon] != '*':
                        amino_acids += genetic_code[codon]
                    else:
                        translations.append(amino_acids)
                        amino_acids = ''
                        break
                else:
                    translations.append(amino_acids)
                    amino_acids = ''
    return translations if translations else []

def get_reverse(sequence):
    """Reverse orientation of `sequence`.

    Returns a string with `sequence` in the reverse order.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_reverse('AUGC')
    'CGUA'
    """
    seq = sequence[::-1]
    if not seq:
       return ""
    return seq.upper()

def get_complement(sequence):
    """Get the complement of a `sequence` of nucleotides.

    Returns a string with the complementary sequence of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_complement('AUGC')
    'UACG'
    """
    complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}
    comp = []
    seq = sequence.upper()
    bases = list(seq)
    for i in bases:
        comp.append(complement[i])
        if not seq:
    	    return ""
    return (''.join(comp))


def reverse_and_complement(sequence):
    """Get the reversed and complemented form of a `sequence` of nucleotides.

    Returns a string that is the reversed and complemented sequence
    of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> reverse_and_complement('AUGC')
    'GCAU'
    """
    complement = {'A': 'U', 'C': 'G', 'U': 'A', 'G': 'C'}
    reverse_comp = []
    seq = sequence.upper()
    bases = list(seq)
    bases.reverse()
    for i in bases:
        reverse_comp.append(complement[i])
    if not seq:
        	return ""
    return (''.join(reverse_comp))

def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (the three reading frames of
    `rna_sequence`, and the three reading frames of the reverse and complement
    of `rna_sequence`) and return (as a string) the longest sequence of amino
    acids that it encodes, according to the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the longest sequence of amino acids encoded by
        `rna_sequence`.
    """
    # Make a variable for the forward sequence
    forward_seq = translate_sequence(rna_sequence, genetic_code)

    # Make a variable for the reverse complement sequence
    reverse_complement_seq = reverse_and_complement(rna_sequence)

    # Get all translations on both sequences and append the translations to the same variable
    all_aa_sequences = []
    longest_aa_sequence = ''

    # Explore the three reading frames of the RNA sequence
    for frame in range(3):
        aa_sequence = translate_sequence(forward_seq[frame:], genetic_code)
        all_aa_sequences.append(aa_sequence)
        if len(aa_sequence) > len(longest_aa_sequence):
            longest_aa_sequence = aa_sequence

    # Explore the three reading frames of the reverse complement of the RNA sequence
    for frame in range(3):
        aa_sequence = translate_sequence(reverse_complement_seq[frame:frame+len(rna_sequence)], genetic_code)
        all_aa_sequences.append(aa_sequence)
        if len(aa_sequence) > len(longest_aa_sequence):
            longest_aa_sequence = aa_sequence
    return longest_aa_sequence




if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
 
