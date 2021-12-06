# Notebook

```bash
# Creation of a new branch
git checkout -b MSA_error
```

## 2021-12-06 error on msa before consensus

There is a problem in the C-NR-genome, the 'N' may not correspond to mismatches but insteed to the absence of MSA upstream of the EMBOSS `cons`. I will therefore proceed to a multiple alignment and a cunt from the genomes assembled

```bash
# conda create -c bioconda -y -n mafft mafft
# conda create -c bioconda -y -n emboss emboss
# conda create -c bioconda -y -n mummer4 mummer4

# MSA
conda activate mafft
cat *.fasta >> oshv.fna
mafft \
    --thread 4 \
    --auto \
    oshv.fna \
    > oshv.faa
conda deactivate

# Generate consensus
conda activate emboss
cons -sequence oshv.faa \
    -outseq C-NR-genome.fasta \
    -plurality 0.2 -name C-NR-genome.fasta
conda deactivate

# WG comparison
conda activate mummer4
NR_genomes_seq=oshv.fna
Consensus_global=C-NR-genome.fasta
# mumref
nucmer -c 100 -l 15 -f -p nucmer_numref $Consensus_global $NR_genomes_seq
show-coords -r -c -l nucmer_numref.delta > nucmer_numref.coords
show-snps -T nucmer_numref.delta > nucmer_numref.snps

conda activate mafft
mafft --addfragments C-NR-genome.fasta  oshv.faa >> oshv_and_C-NR.faa
```