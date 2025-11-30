# 📌 CPOX Gene Analysis — BIO 312 Final Project

This repository contains the data, results, and a reproducible workflow for my BIO 312 term project. I analyzed the evolution of coproporphyrinogen oxidase (CPOX) across aquatic and terrestrial vertebrates by:

- Identifying homologs using local BLAST
- Creating a multiple sequence alignment (MSA)
- Inferring a maximum-likelihood phylogeny using IQ-TREE
- Analyzing AlphaFold structural predictions for:
  - Net charge
  - Secondary structure
  - Solvent-accessible surface area (ASA)

## 📂 Repository Structure
```text
data/
  raw/
  processed/
results/
  trees/
  figures/
  structures/
code/
  01_make_homologs_fasta.sh
  02_make_alignment.sh
  03_run_iqtree.sh
  04_net_charge_and_plots.sh
  05_dssp_and_plots.sh
  plot_cpox_traits.py
```
🧰 Software
BLAST+, samtools, seqkit, MUSCLE v5, IQ-TREE2, gotree, newick utilities, ImageMagick, Python3 (pandas, matplotlib, mdtraj), DSSP via mdtraj

📊 Figures Used in Term Paper
Figure 1 - results/figures/Figure1_CPOX_MSA.pdf
Figure 2 - results/figures/Figure2_CPOX_NetCharge.png
Figure 3 - results/figures/Figure3_CPOX_DSSP_ASA.png
Figure 4 - results/figures/Figure4_CPOX_Phylogeny.pdf

🧪 Full Reproducible Workflow (Commands)
(All commands needed to regenerate Figures 1–4)

Step 1 — Homolog Search (BLAST)

```bash
samtools faidx ../lab03-s3-Svargas/allprotein.fas 'Hsap|NP_000088.3|CPOX' > data/raw/NP_000088.3.fa
blastp -db ../lab03-s3-Svargas/allprotein.fas -query data/raw/NP_000088.3.fa -evalue 1e-30 -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" -out data/raw/CPOX.blastp.detail.out
awk '{ if ($6 < 1e-30) print $2 }' data/raw/CPOX.blastp.detail.out > data/raw/CPOX.blastp.detail.filtered.out
seqkit grep --pattern-file data/raw/CPOX.blastp.detail.filtered.out ../lab03-s3-Svargas/allprotein.fas > data/processed/CPOX.homologs.fas
```
Step 2 — Multiple Sequence Alignment

```bash
muscle -align data/processed/CPOX.homologs.fas -output data/processed/CPOX.homologs.al.fas
Rscript --vanilla ../lab04-s3-Svargas/plotMSA.R data/processed/CPOX.homologs.al.fas results/figures/Figure1_CPOX_MSA.pdf
```
Step 3 — Phylogeny

```bash
iqtree2 -s data/processed/CPOX.homologs.al.fas -bb 1000 -nt AUTO -pre results/trees/CPOX
gotree reroot midpoint -i results/trees/CPOX.treefile -o results/trees/CPOX.homologs.al.mid.treefile
nw_order -c n results/trees/CPOX.homologs.al.mid.treefile | nw_display -s | convert svg:- results/figures/Figure4_CPOX_Phylogeny.pdf
```
Step 4 — Net Charge

```bash
for pqr in results/structures/*.pqr; do base=$(basename "$pqr" .pqr); awk '/ATOM/ {sum+=$9} END {print base"\t"sum}' base="$base" "$pqr"; done > data/processed/net_charges.tsv
python3 code/plot_cpox_traits.py --net-charges data/processed/net_charges.tsv --species-key data/processed/species_key.csv --outdir results/figures
```
Step 5 — DSSP / ASA

```bash
python3 ../lab06-s3-Svargas/dssp_batch_summary_mdtraj.py --pdb-dir results/structures --species-key data/processed/species_key.csv --refseq-map data/raw/CPOX.blastp.detail.filtered.out --out-csv data/processed/dssp_summary.csv --plots
python3 code/plot_cpox_traits.py --dssp-summary data/processed/dssp_summary.csv --species-key data/processed/species_key.csv --outdir results/figures
```
🔁 Full regeneration
```bash
bash code/01_make_homologs_fasta.sh
bash code/02_make_alignment.sh
bash code/03_run_iqtree.sh
bash code/04_net_charge_and_plots.sh
bash code/05_dssp_and_plots.sh
```
