#!/usr/bin/env python3

import argparse
import os
import subprocess
import time
from collections import Counter
from math import log2
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
from sklearn.linear_model import LassoCV

def parse_args():
    parser = argparse.ArgumentParser(description="Alignment filtering using LASSO or entropy.")
    parser.add_argument("-s", "--seqfile", required=True, help="FASTA file with aligned sequences")
    parser.add_argument("-l", "--lasso", type=int, choices=[0, 1], default=0, help="Use LASSO analysis (1=yes, 0=no)")
    parser.add_argument("-r", "--replicates", type=int, default=10000, help="Number of replicates for LASSO analysis")
    return parser.parse_args()

def check_alignment(infile):
    seqs = list(SeqIO.parse(infile, "fasta"))
    seqlens = {len(record.seq) for record in seqs}
    if len(seqlens) != 1:
        raise ValueError("Sequences are not aligned (different lengths).")
    return seqs

def run_model_finder(infile):
    print("	(1) Fitting the best substitution model using ModelFinder.")
    subprocess.Popen(["iqtree2", "--quiet", "-s", infile, "-m", "TESTONLY", "-redo"])
    while not os.path.exists(f"{infile}.iqtree"):
        time.sleep(1)
    with open(f"{infile}.iqtree") as f:
        for line in f:
            if "Model of substitution:" in line:
                return line.strip().split(":")[1].strip()
    raise ValueError("Could not find best model in .iqtree file.")

def run_lasso_replicates(infile, model, nseq, replicates):
    outname = f"{infile}.siteloglike.txt"
    print("	(2) Calculating site-wise log-likelihoods of random tree replicates.")
    for _ in range(replicates):
        subprocess.run(["iqtree2", "--quiet", "-s", infile, "-r", str(nseq), "rand.tre", "-redo"])
        subprocess.run(["iqtree2", "--quiet", "-s", infile, "-t", "rand.tre", "--sitelh", "-blfix", "-m", model, "-redo"])
        with open(f"{infile}.sitelh") as f:
            lh = next((line.strip().replace("Site_Lh", "").strip() for line in f if "Site_Lh" in line), "")
        with open(f"{infile}.iqtree") as f:
            lht = next((line.strip().replace("Log-likelihood of the tree:", "").split(" (")[0].strip()
                        for line in f if "Log-likelihood of the tree" in line), "")
        with open(outname, "a") as f:
            f.write(f"{lht} {lh}\n")

def run_python_lasso(infile):
    print("	(3) Running LASSO analysis in Python.")
    dat = pd.read_csv(f"{infile}.siteloglike.txt", sep=" ", header=None)
    y = dat.iloc[:, 0].values
    X = dat.iloc[:, 1:].values
    model = LassoCV(cv=5, random_state=42).fit(X, y)
    coef = model.coef_
    ok = [1 if c != 0 else 0 for c in coef]
    with open("LassoClassification.txt", "w") as out:
        out.write("sitewise_info " + " ".join(map(str, ok)) + "\n")

def write_lasso_filtered_fasta(infile, records):
    with open("LassoClassification.txt") as f:
        line = f.readline().strip().split()
        ok = list(map(int, line[1:]))
    indices = [i for i, v in enumerate(ok) if v == 1]
    with open(f"{infile}.lasso.info.fasta", "w") as out:
        for rec in records:
            new_seq = ''.join([rec.seq[i] for i in indices])
            rec.seq = Seq(new_seq)
            SeqIO.write(rec, out, "fasta")

def shannon_entropy(column):
    counts = Counter([c for c in column if c in "ACGT-"])
    total = sum(counts.values())
    if total == 0:
        return 0
    return -sum((count / total) * log2(count / total) for count in counts.values())

def write_entropy_filtered_fasta(infile, records):
    aln_len = len(records[0].seq)
    indices = []
    for i in range(aln_len):
        column = [rec.seq[i] for rec in records]
        if all(c.upper() not in "N" for c in column):
            ent = shannon_entropy(column)
            if ent >= 0.5:
                indices.append(i)
    with open(f"{infile}.entropy.info.fasta", "w") as out:
        for rec in records:
            new_seq = ''.join([rec.seq[i] for i in indices])
            rec.seq = Seq(new_seq)
            SeqIO.write(rec, out, "fasta")

def main():
    args = parse_args()
    infile = args.seqfile
    use_lasso = args.lasso
    replicates = args.replicates

    records = check_alignment(infile)
    nseq = len(records)

    if use_lasso == 1:
        model = run_model_finder(infile)
        run_lasso_replicates(infile, model, nseq, replicates)
        run_python_lasso(infile)
        write_lasso_filtered_fasta(infile, records)
    else:
        write_entropy_filtered_fasta(infile, records)

if __name__ == "__main__":
    main()
