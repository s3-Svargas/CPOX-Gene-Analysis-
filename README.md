# CPOX Gene Analysis — BIO 312 Final Project

This repository contains the data, results, and a reproducible workflow for my BIO 312 term project. I analyzed the evolution of coproporphyrinogen oxidase (CPOX) across aquatic and terrestrial vertebrates by:

- Identifying homologs using local BLAST
- Creating a multiple sequence alignment (MSA)
- Inferring a maximum-likelihood phylogeny using IQ-TREE
- Analyzing AlphaFold structural predictions for:
  - Net charge
  - Secondary structure
  - Solvent-accessible surface area
 
# CPOX – Lab 3: Identifying Homologs with Local BLAST

## Clone Lab 3 Repository

```bash
cd ~
git clone git@github.com:Bio312/lab03-$MYGIT.git
cd lab03-$MYGIT
```

Create a working directory for CPOX:

```bash
mkdir CPOX
cd CPOX
```

## Make a directory for CPOX work.

```bash
cd ~/lab03-$MYGIT
mkdir CPOX
cd CPOX
```

## Extract the human CPOX sequence as the BLAST query

Use samtools faidx to pull the human CPOX protein from the combined allprotein.fas file.

```bash
samtools faidx ~/lab03-$MYGIT/allprotein.fas 'Hsap|NP_000088.3|CPOX' \
  > ~/lab03-$MYGIT/CPOX/NP_000088.3.fa
```
This creates NP_000088.3.fa containing the human CPOX protein that will be used as the query for BLAST.

## Run BLASTp for CPOX against the course protein database

Run a local BLASTp search of the human CPOX query against the allprotein.fas database you created earlier in Lab 3.

```bash
blastp \
  -db ~/lab03-$MYGIT/allprotein.fas \
  -query ~/lab03-$MYGIT/CPOX/NP_000088.3.fa \
  -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
  -out ~/lab03-$MYGIT/CPOX/CPOX.blastp.detail.out
```
This writes a tab-delimited file CPOX.blastp.detail.out with one line per hit, including sequence IDs, percent identity, alignment length, E-value, bit score, and annotation.

## Filter CPOX hits by E-value to keep high-confidence homologs 

Now keep only strong hits with E-value below a strict cutoff (e.g., 1e-30) and save the matching subject IDs to a new file.

```bash
awk '{ if ($6 < 1e-30) print $2 }' \
  ~/lab03-$MYGIT/CPOX/CPOX.blastp.detail.out \
  > ~/lab03-$MYGIT/CPOX/CPOX.blastp.detail.filtered.out
```

This produces CPOX.blastp.detail.filtered.out, containing the IDs of strong CPOX homologs that you carry forward into Lab 4.

# Lab 4 — CPOX Multiple Sequence Alignment

## Clone Lab 4 Repository

```bash
cd ~
git clone git@github.com:Bio312/lab04-$MYGIT.git
cd lab04-$MYGIT
mkdir CPOX
cd CPOX
```

## Make a FASTA of All Identified CPOX Homologs

Use filtered IDs from Lab 3:

```bash
seqkit grep \
  --pattern-file ~/lab03-$MYGIT/CPOX/CPOX.blastp.detail.filtered.out \
  ~/lab03-$MYGIT/allprotein.fas \
  | seqkit grep -v -p "carpio" \
  > ~/lab04-$MYGIT/CPOX/CPOX.homologs.fas
```

## Align CPOX Homologs with MUSCLE

```bash
muscle \
  -align ~/lab04-$MYGIT/CPOX/CPOX.homologs.fas \
  -output ~/lab04-$MYGIT/CPOX/CPOX.homologs.al.fas
```

## Alignment Visualization

```bash
alv -kli ~/lab04-$MYGIT/CPOX/CPOX.homologs.al.fas \
  | less -RS
```

## Generate MSA PDF

```bash
Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R \
  ~/lab04-$MYGIT/CPOX/CPOX.homologs.al.fas
```

This will create CPOX.homologs.al.fas.pdf. 

# Lab 5 — CPOX Gene Family Phylogeny (IQ-TREE)

## Clone Lab 5 Repository 

```bash
cd ~
git clone git@github.com:Bio312/lab05-$MYGIT.git
cd lab05-$MYGIT
mkdir CPOX
cd CPOX
```

## Copy Aligned CPOX Homologs from Lab 4

```bash
cp ~/lab04-$MYGIT/CPOX/CPOX.homologs.al.fas \
   ~/lab05-$MYGIT/CPOX/
```
## Run IQ-TREE (Maximum Likelihood)

This performs a model selection, an ML tree search, and an Ultrafast bootstrap (1,000 replicates)

```bash
iqtree \
  -s ~/lab05-$MYGIT/CPOX/CPOX.homologs.al.fas \
  -bb 1000 \
  -nt 2 \
  -pre ~/lab05-$MYGIT/CPOX/CPOX
```

You will get three files:
CPOX.treefile, which is the Unrooted ML tree (Newick)
CPOX.iqtree, which has the model info & bootstrap support
CPOX.log, which has the runtime log 

You should now check the model:

```bash
grep -A3 "Best-fit model" ~/lab05-$MYGIT/CPOX/CPOX.iqtree
```

## Display Tree (ASCII, Unrooted) 

```bash
nw_display ~/lab05-$MYGIT/CPOX/CPOX.treefile
```

## Save Unrooted Graphic PDF

```bash
Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R \
  ~/lab05-$MYGIT/CPOX/CPOX.treefile \
  ~/lab05-$MYGIT/CPOX/CPOX.treefile.pdf \
  0.4 15
```

## Midpoint Root the Tree

```bash
gotree reroot midpoint \
  -i ~/lab05-$MYGIT/CPOX/CPOX.treefile \
  -o ~/lab05-$MYGIT/CPOX/CPOX.mid.treefile
```

Display rooted ASCII:

```bash
nw_order -c n ~/lab05-$MYGIT/CPOX/CPOX.mid.treefile | nw_display -
```

## Save Midpoint-Rooted Phylogram PDF
```bash
nw_order -c n ~/lab05-$MYGIT/CPOX/CPOX.mid.treefile \
  | nw_display -s -w 1000 \
      -l 'font-size:10px;font-family:sans' \
      -i 'font-size:8px;font-family:sans' \
      -b 'visibility:hidden' \
      -I l -n -4 -W 4.0 -v 24 - \
  | convert -density 300 svg:- \
      ~/lab05-$MYGIT/CPOX/CPOX.mid.treefile.pdf
```

## Save Midpoint-Rooted Cladogram PDF (Topology-only)

```bash
nw_order -c n ~/lab05-$MYGIT/CPOX/CPOX.mid.treefile \
  | nw_topology - \
  | nw_display -s -w 1000 \
      -l 'font-size:10px;font-family:sans' \
      -i 'font-size:8px;font-family:sans' \
      -b 'visibility:hidden' \
      -S -I l -n -4 -W 4.0 -v 24 - \
  | convert -density 300 svg:- \
      ~/lab05-$MYGIT/CPOX/CPOX.mid.cladogram.pdf
```

## Save Command History + Push to GitHub

```bash
cd ~/lab05-$MYGIT

history > lab5.commandhistory.txt

git add .
git commit -m "CPOX ML Tree, Midpoint Rooted Tree, PDFs, and History"
git pull --no-edit
git push
```

