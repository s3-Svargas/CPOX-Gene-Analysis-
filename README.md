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
# Lab 6 — CPOX Structural Analysis and Net Charge

## 6.0 Clone the Lab 6 repository and create a CPOX directory

```bash
cd ~
git clone git@github.com:Bio312/lab06-"$MYGIT".git
cd lab06-"$MYGIT"
mkdir CPOX
cd CPOX
```
## Map CPOX RefSeq IDs to UniProt accessions

```bash
awk 'NR==FNR {m[$1]=$2; next} {if ($0 in m) print $0"\t"m[$0]; else print $0"\tMISSING"}' \
  ~/lab06-"$MYGIT"/refseq_uniprot.tsv \
  ~/lab06-"$MYGIT"/CPOX/CPOX_refseq_ids.txt \
  > ~/lab06-"$MYGIT"/CPOX/CPOX_refseq_uniprot_subset.tsv
```

This produces a mapping file linking CPOX RefSeq IDs to UniProt accessions.

## Calculate net charges from PQR structures

```bash
out=~/lab06-"$MYGIT"/CPOX/net_charges.tsv
: > "$out"  # Create/truncate output file

for f in ~/lab06-"$MYGIT"/CPOX/*.pqr; do
  base=$(basename "$f")
  Z=$(awk '$1=="ATOM"||$1=="HETATM"{s+=$9}END{printf"%.3f",s}' "$f")
  printf "%s\tNetCharge=%s\n" "$base" "$Z" >> "$out"
done
```

This generates net_charges.tsv, storing net charge values for each protein structure.

## Build helper mapping tables linking structures to habitat
Abbreviation -> status (aquatic/terrestrial)
```bash
awk -F, 'NR>1 {print $2"\t"$5}' \
  ~/lab06-"$MYGIT"/species_key.csv \
  | sort -u \
  > abbr_status.tsv
```
RefSeq -> abbreviation mapping

```bash
awk -F'|' '{print $2"\t"$1}' \
  ~/lab03-"$MYGIT"/CPOX/CPOX.blastp.detail.filtered.out \
  | sort -u \
  > refseq_abbr.tsv
```

RefSeq -> PNG filename mapping

```bash
ls *.png 2>/dev/null \
  | sed 's/.png$//' \
  | awk -F'__' '{print $1"\t"$0".png"}' \
  | sort -u \
  > refseq_png.tsv
```

Join all tables to link PNG to species + habitat

```bash
awk -F'\t' '{print $2"\t"$1}' refseq_abbr.tsv \
  | sort -u \
  > abbr_refseq.tsv

join -t $'\t' -1 1 -2 1 abbr_refseq.tsv abbr_status.tsv \
  > abbr_refseq_status.tsv

awk -F'\t' '{print $2"\t"$1"\t"$3}' abbr_refseq_status.tsv \
  | sort -u \
  > refseq_abbr_status.tsv

join -t $'\t' -1 1 -2 1 refseq_png.tsv refseq_abbr_status.tsv \
  > png_refseq_abbr_status.tsv
```

This produces a table with columns:
```text
RefSeqID PNG Abbr Status
```

## Classify structural images by habitat

```bash
awk -F'\t' '$4=="aquatic" {print $2}' png_refseq_abbr_status.tsv \
  > aquatic_pngs.txt

awk -F'\t' '$4=="terrestrial" {print $2}' png_refseq_abbr_status.tsv \
  > terrestrial_pngs.txt

column -t -s $'\t' png_refseq_abbr_status.tsv \
  | tee PNG_GUIDE.txt
```

This creates habitat-based file lists and a guide table for figures.

## DSSP structural summary and plotting

```bash
python3 ~/lab06-"$MYGIT"/dssp_batch_summary_mdtraj.py \
  --pdb-dir ~/lab06-"$MYGIT"/CPOX \
  --species-key ~/lab06-"$MYGIT"/species_key.csv \
  --refseq-map ~/lab03-"$MYGIT"/CPOX/CPOX.blastp.detail.filtered.out \
  --out-csv ~/lab06-"$MYGIT"/CPOX/dssp_summary.csv \
  --plots
```

This generates dssp_summary.csv and structural comparison plots used in the results.

# Lab 8 - Evaluating sequencing reads and genome assemblies 

## Clone Lab 8 Repository and Set Up CPOX Directory
```bash
cd ~
git clone git@github.com:Bio312/lab08-"$MYGIT".git
cd lab08-"$MYGIT"
mkdir CPOX
cd CPOX
```

## Copy CPOX Multiple Sequence Alignment from Lab 4

```bash
cp ~/lab04-"$MYGIT"/CPOX/CPOX.homologs.al.fas \
   ~/lab08-"$MYGIT"/CPOX/
```

## Infer CPOX Gene Tree with IQ-TREE (ML + Model Selection)

```bash
iqtree \
  -s ~/lab08-"$MYGIT"/CPOX/CPOX.homologs.al.fas \
  -m MFP \
  -bb 1000 \
  -nt 2 \
  -pre ~/lab08-"$MYGIT"/CPOX/CPOX
```

This generates CPOX.treefile, CPOX.log, and model stats using ModelFinder.

## Midpoint-Root the Maximum-Likelihood Tree

```bash
gotree reroot midpoint \
  -i ~/lab08-"$MYGIT"/CPOX/CPOX.treefile \
  -o ~/lab08-"$MYGIT"/CPOX/CPOX.mid.treefile
```

This produces a midpoint-rooted phylogeny for interpretation.

## Save Midpoint-Rooted Phylogram as a PDF

```bash
nw_order -c n ~/lab08-"$MYGIT"/CPOX/CPOX.mid.treefile \
  | nw_display -s -w 1000 \
      -l 'font-size:10px;font-family:sans' \
      -i 'font-size:8px;font-family:sans' \
      -b 'visibility:hidden' \
      -I l -n -4 -W 4.0 -v 24 - \
  | convert -density 300 svg:- \
      ~/lab08-"$MYGIT"/CPOX/CPOX.mid.treefile.pdf
```

This creates the final publication-quality phylogram used in the Results section.

