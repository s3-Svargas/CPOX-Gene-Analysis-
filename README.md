# CPOX Gene Analysis – BIO 312 Final Project

This repository contains the data and results for my BIO 312 term project on the evolution of coproporphyrinogen oxidase (CPOX) in aquatic vs terrestrial vertebrates. It is organized so that someone can see the input files, the main outputs, and the commands used to generate the figures in my paper.

---

## 1. Repository layout

```text
data/
  raw/         # starting files (e.g., NP_000088.3.fa, BLAST output, filtered hit list)
  processed/   # homolog FASTA, alignment, summary tables (net charges, DSSP, etc.)

results/
  trees/       # IQ-TREE output files, midpoint-rooted gene tree
  figures/     # final figures used in the term paper (alignment, net charge, DSSP, tree)
