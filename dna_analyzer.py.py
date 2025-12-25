#!/usr/bin/env python3
"""
dna_analyzer.py
Simple DNA analyzer for Pydroid 3 (no external libraries).
Usage:
  - Interactive: python3 dna_analyzer.py
  - With file:  python3 dna_analyzer.py -i sample.fasta -k 3
"""
import sys
import argparse
import re
from collections import Counter

RC_MAP = str.maketrans("ACGT", "TGCA")

def read_sequence_from_file(path):
    try:
        with open(path, "r") as f:
            lines = [line.strip() for line in f if line.strip()]
    except Exception as e:
        print(f"Error reading file: {e}")
        return ""
    # If FASTA: skip header lines starting with '>'
    if lines and lines[0].startswith(">"):
        seq = "".join(line for line in lines if not line.startswith(">"))
    else:
        seq = "".join(lines)
    return seq

def clean_seq(seq):
    seq = seq.upper()
    # keep only A,C,G,T
    return re.sub("[^ACGT]", "", seq)

def gc_content(seq):
    if not seq:
        return 0.0
    g = seq.count("G")
    c = seq.count("C")
    return 100.0 * (g + c) / len(seq)

def at_content(seq):
    if not seq:
        return 0.0
    a = seq.count("A")
    t = seq.count("T")
    return 100.0 * (a + t) / len(seq)

def reverse_complement(seq):
    return seq.translate(RC_MAP)[::-1]

def kmer_counts(seq, k=3):
    counts = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if len(kmer) == k:
            counts[kmer] += 1
    return counts

def summary(seq, k=3, topn=5):
    seq_clean = clean_seq(seq)
    length = len(seq_clean)
    print("\n--- DNA Analyzer Summary ---")
    print(f"Sequence length: {length}")
    if length == 0:
        print("No valid A/C/G/T bases found.")
        return
    print(f"GC content: {gc_content(seq_clean):.2f}%")
    print(f"AT content: {at_content(seq_clean):.2f}%")
    print(f"Reverse complement: {reverse_complement(seq_clean)}")
    counts = kmer_counts(seq_clean, k)
    if counts:
        most = counts.most_common(topn)
        print(f"Top {topn} {k}-mers:")
        for kmer, cnt in most:
            print(f"  {kmer}: {cnt}")
    print("-----------------------------\n")

def interactive_mode(default_k=3):
    print("Enter DNA sequence (paste raw sequence or FASTA). Finish with an empty line:")
    lines = []
    try:
        while True:
            line = input()
            if line.strip() == "":
                break
            lines.append(line)
    except EOFError:
        pass
    seq = "\n".join(lines)
    summary(seq, k=default_k)

def main():
    parser = argparse.ArgumentParser(description="Simple DNA Analyzer")
    parser.add_argument("-i", "--input", help="Input file path (plain or FASTA)")
    parser.add_argument("-k", "--kmer", type=int, default=3, help="k-mer size (default 3)")
    args = parser.parse_args()

    if args.input:
        seq = read_sequence_from_file(args.input)
        if not seq:
            print("No sequence found in file or file couldn't be read.")
            return
        summary(seq, k=args.kmer)
    else:
        # interactive/paste mode (good for mobile)
        interactive_mode(default_k=args.kmer)

if __name__ == "__main__":
    main()