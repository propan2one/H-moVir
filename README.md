# Genomic epidemiology, genetic diversity and connectivity of the Ostreid Herpesvirus 1 population in France

Scripts used for recreating figures and analyses for [DOI](). Scripts were not coded with the intention of reuse by others. Provided as a courtesy for better transparency and reproducibility of research.

## Abstract

The genetic diversity of viral populations has become a major issue in the understanding of their phylogeographic and dissemination history, but it represents a real challenge when studying whole-genome viral diversity in the natural environment. These molecular ecological approaches to study viral diversity are commonly used for RNA viruses harboring small genomes but have not been applied to DNA viruses with large genomes. In this study we used the “Pacific Oyster Mortality Syndrome” (POMS, a disease that affects oyster farms around the world) as a model to study the genetic diversity of its causative agent, the Ostreid herpes virus 1 (OsHV-1) in the three most important French oyster-farming areas. Using ultra-deep sequencing on individual moribund oysters and new bioinformatics methodology we de novo assembled twenty one OsHV-1 genomes. Then, using a combination of whole-genomes comparison, phylogenetic analysis and quantification of major and minor variants we unveiled the connectivity of OsHV-1 viral populations between the three oyster-farming areas. We propose a scenario in which the Marennes-Oléron basin would be the main source of OsHV-1 diversity which would then be dispersed to other farming areas, a hypothesis consistent with the current knowledge on oyster transfers in France. In conclusion, this study showed that molecular ecological approaches may be applied to large genome virus to unravel the extent of their genetic diversity and better understand the spread of viral populations in natural environments.

![graphical_summary](https://github.com/propan2one/OshV-1-molepidemio/blob/main/image/graphical_summary.png?raw=true)

## Computational platform

- The upstream data treatments were carried out under job scheduler for high-performance computing (PBSpro).

- The downstream data treatments were carried out locally on a Dell computer with 62,8 Go RAM, intel® Core™ i7-7700 CPU @ 3.60GHz × 8 wiith Ubuntu OS 16.04 LTS.

## Retrieval of sequences from databases

### Previously published OsHV genome

The virus genomes used in the study were downloaded directly from NCBI as follows:

| Genome                                                                                                             |  Id accession |
|--------------------------------------------------------------------------------------------------------------------|---------------|
| [Ostreid_herpesvirus_1_strain_microVar_variant_A_complete_genome](https//www.ncbi.nlm.nih.gov/nuccore/KY242785.1/) |    KY242785.1 |
| [Ostreid_herpesvirus_1_complete_genome](https//www.ncbi.nlm.nih.gov/nuccore/AY509253.2/)                           |  AY509253.2   |
| [Ostreid_herpesvirus_1_strain_microVar_variant_B_complete_genome](https//www.ncbi.nlm.nih.gov/nuccore/KY271630.1/) |    KY271630.1 |
| [Ostreid_herpesvirus_1_isolate_ZK0118_complete_genome](https//www.ncbi.nlm.nih.gov/nuccore/MF509813.1/)            |   MF509813.1  |
| [Chlamys_acute_necrobiotic_virus_complete_genome](https//www.ncbi.nlm.nih.gov/nuccore/GQ153938.1/)                 |    GQ153938.1 |
| [Abalone_herpesvirus_Taiwan_2005_complete_genome](https//www.ncbi.nlm.nih.gov/nuccore/KU096999.1/)                 |    KU096999.1 |
| [Ostreid_herpesvirus_1_strain_CDSB2012_complete_genome](https//www.ncbi.nlm.nih.gov/nuccore/KP412538.1/)           |  KP412538.1   |
| [Ostreid_herpesvirus_1_2016_PT_complete_genome](https//www.ncbi.nlm.nih.gov/nuccore/MG561751.2/)                   |  MG561751.2   |

```bash
# Download fasta format
#OsHV_genomes_folder=~/Documents/OshV-1-molepidemio/raw/a-OsHV-1-NCBI-genome # exemple
OsHV_genomes_folder=/PATH/OshV-1-molepidemio/raw/a-OsHV-1-NCBI-genome
cd $OsHV_genomes_folder
while read name id ; do
    ncbi-acc-download -F fasta ${id} -p ${name}
done < $OsHV_genomes_folder/genome_name_and_id.csv

# Download gff3 format
while read name id ; do
    ncbi-acc-download -F gff3 ${id} -p ${name}
done < $OsHV_genomes_folder/genome_name_and_id.csv
```

### Host genome

The Host genome used in the study were downloaded directly from NCBI as follows:

```bash
# Download fasta format
#Host_genome_folder=~/Documents/OshV-1-molepidemio/raw/c-Host-NCBI-genome # exemple
Host_genome_folder=/PATH/OshV-1-molepidemio/raw/a-OsHV-1-NCBI-genome
cd $Host_genome_folder
id=NW_011935992.1
name=oyster.v9
ncbi-acc-download -F fasta ${id} -p ${name}
```

## Upstream analysis

### 00) Metadata analysis

The metadata provided by the sequencing platform was cleaned up using an [00-preliminary_data.Rmd](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/00-preliminary_data.Rmd). This made the file `ID_experiment.csv` synchronizing the filenames and identifiers readable. Subsequently, the file to use can be found in `OshV-1-molepidemio/raw/b-raw_metadatas/ID_experiment.csv` and the insertion size has been added by hand with the help of `libreoffice --calc`.

Note: Login and passwords are missing in the code. The paths used in the code have also been replaced for security reasons.

___

### 01) Data transfer

All data has been transferred and MD5 have been check from the servers of the transfer platform to the datarmor calculation server with : [01-data_transfer.pbs](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/01-data_transfer.pbs)

___

### 02) Genetic diversity rarefaction analysis

The analysis of genetic diversity rarefaction was performed for all samples individually and the command was executed as follows using this script: [02-rarefaction_virus.pbs](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/02-rarefaction_virus.pbs)

```bash
while read h f; do r1=`ls /PATH/OshV-1-molepidemio/raw/${h}*_R1.fastq.gz`; r2=`ls /PATH/OshV-1-molepidemio/raw/${h}*_R2.fastq.gz`; qsub -v "rate=5000, id=${f}, reads1=${r1}, reads2=${r2}, outdir=/home1/scratch/jdelmott/2020-02-12-Rarefaction_Haplofit/, genomefile=/PATH/OshV-1-molepidemio/raw/a-OsHV-1-NCBI-genome/Ostreid_herpesvirus_1_strain_microVar_variant_A_complete_genome_0.fa, gffFile=/PATH/OshV-1-molepidemio/raw/a-OsHV-1-NCBI-genome/Ostreid_herpesvirus_1_strain_microVar_variant_A_complete_genome_0.gff" OshV-1-molepidemio/src/02-rarefaction_virus.pbs; done < OshV-1-molepidemio/raw/b-raw_metadatas/ID_experiment.csv
```

___

### 03) Creation of Non Rredundant-genomes from NCBI

![NR_genomes.png](https://github.com/propan2one/OshV-1-molepidemio/blob/main/image/NR_genomes.png?raw=true)

- In a first step, the sequences were analyzed on Ugene with a dotplot of the sequence against itself to identify the different structures of the genomes. Each repetition was manually annotated as follows: 1-Ul, 2-IRl, 3-X, 4-IRs, 5-Us. 

- In a second step, the annotations were then saved in a csv and concatenate in csv format [NR_genomic_part_coordonate.csv](https://github.com/propan2one/OshV-1-molepidemio/blob/main/raw/a-OsHV-1-NCBI-genome/NR_genomic_part_coordonate.csv) to be cut out using the [seqkit](https://bioinf.shenwei.me/seqkit/) tool and the `seqkit subseq -r` option.

- In a third step, the sequences will be concatenated to create the non-redundant (NR) genomes. 

See [03-NR_genome_genomes_construction.md](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/03-NR_genome_genomes_construction.md) for more details.


Checks were performed using multiple alignment using [mafft](https://mafft.cbrc.jp/alignment/software/) using the script [mafft_MSA.pbs](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/mafft_MSA.pbs). Visualization of the alignment was done with [Aliview](https://github.com/AliView/AliView).

___

### 04) Assembly of the genome of the virus of interest

The construction of a genome for each was carried out for all samples individually and the command was executed as follows using this script: [04-metaviromics.pbs](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/04-metaviromics.pbs)

```bash
while read h f i ; do r1=`ls /PATH/OshV-1-molepidemio/raw/${h}*_R1.fastq.gz`; r2=`ls /PATH/OshV-1-molepidemio/raw/${h}*_R2.fastq.gz`; qsub -v "id=${f}, reads1=${r1},reads2=${r2}, genomefile=/home1/datawork/jdelmott/data_jean/oyster.v9.fa,database=/home1/datawork/jdelmott/data_jean/viral.2.1.genomic.fna,GenomeOsHV1=/home1/datawork/jdelmott/data_jean/OsHV-1_strain_microVar_variant_A.fasta,mincontig=200,minlength=50,outdir=/home1/scratch/jdelmott/2020-03-20-Haplofit_metaviromic,insersize=${i}" OshV-1-molepidemio/src/04-metaviromics.pbs; done < OshV-1-molepidemio/raw/b-raw_metadatas/ID_experiment.csv
```

**Note** that the multiple alignment used to check the sub-structure of the non-redundant genome is not provided in the command line parameters.

___

### 05) Genome cleaning

When necessary, the genomes were cleaned manually. Up to three events were studied: an artifact of a 5' UL sequence, a few small tandem repeat sequences of about 160 bp (polyN) and a location containing a repeat of several Cytosine (polyC). All the analyses were performed by computer and the code can be found in [05-Genome_cleaning.md](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/05-Genome_cleaning.md) file.

The set of assembled non-redundant genomes can be found at this location: [NR-Asm-genome](https://github.com/propan2one/OshV-1-molepidemio/tree/main/results/NR-Asm-genome)

___

### 06) Phylogenetic analysis

Once all the NR genomes were assembled and corrected, an alignment using the NR versions of the redundant genomes was performed. To do this, the genomes were aligned with each other using the [mafft_MSA.pbs](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/mafft_MSA.pbs) script and then the tree construction was performed using the [06-inference_iqtree.pbs](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/06-inference_iqtree.pbs) script (1000 bootstrap). The outgroup has been designated as `NR-genome_Abalone_herpesvirus_Taiwan_2005.fasta`.

```bash
# global
qsub -v "outdir=,MSA=Fig2A_global_phylogeny.msa,outgroup=NR-genome_Abalone_herpesvirus_Taiwan_2005.fasta" 06-inference_iqtree.pbs
```

A few minor modifications were made later as described below:

```python
# changing format
from Bio import AlignIO
AlignIO.convert("NR_genome_corrected.maf", "fasta", "NR_genome_corrected.phy", "phylip-relaxed")
```

```bash
#The best-scoring maximum-likelihood tree
~/Software/iqtree-1.6.7-Linux/bin/iqtree -s NR_genome_corrected.phy -m MFP

#### Reading and visualizing tree files
java -jar ~/Software/FigTree\ v1.4.4/lib/figtree.jar
```

Subsequently the creation and data analysis of the figure was done under R in the [b-figure-02.Rmd](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/b-figure-02.Rmd) file. The `TREEFILE` can be found here: [TreeFile](https://github.com/propan2one/OshV-1-molepidemio/blob/main/results/Phylogenetics_analysis/Fig2A_global_phylogeny.msa.treefile).

___

### 07) Non redundant consensus generation and alignment

To make genome-wide comparisons, we compared all of these genomes individually to the consensus genome of these genomes. To do this we first generated a consensus sequence using the EMBOSS suite with `-plurality 0.2` as follow:

```bash
# Generation of consenus between the 21 NR genomes

# Generation of consenus between the 21 NR genomes
cat PATH/OshV-1-molepidemio/results/NR-Asm-genome/*.fasta >> PATH/OshV-1-molepidemio/results/NR-Asm-genome/NR_genomes_seq.fasta
cons -sequence PATH/OshV-1-molepidemio/results/NR-Asm-genome/NR_genomes_seq.fasta \
    -outseq PATH/OshV-1-molepidemio/results/NR-Asm_consensus/C-NR-genome.fasta \
    -plurality 0.2 -name C-NR-genome.fasta
```

Then the genome-wide comparison is carried out against the consensus genome as described below.

```bash
#conda create -y -n mummer4
source activate mummer4
#conda install -c bioconda mummer4
ln -s PATH/OshV-1-molepidemio/results/NR-Asm_consensus/C-NR-genome.fasta .
cat PATH/OshV-1-molepidemio/results/NR-Asm-genome/*.fasta >> PATH/OshV-1-molepidemio/results/NR-Asm-genome/NR_genomes_seq.fasta # creattion du multifasta de tout les génomes
NR_genomes_seq=PATH/OshV-1-molepidemio/results/NR-Asm-genome/NR_genomes_seq.fasta
Consensus_global=PATH/OshV-1-molepidemio/results/NR-Asm_consensus/C-NR-genome.fasta

# mumref
nucmer -c 100 -l 15 -f -p nucmer_numref $Consensus_global $NR_genomes_seq
show-coords -r -c -l nucmer_numref.delta > nucmer_numref.coords 
show-snps -T nucmer_numref.delta > nucmer_numref.snps
```

Thereafter the analysis is done under R in the [INCOMPLET](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/.Rmd) file.

___

### 08) Variant calling analysis

The variant calling analysis was made against the C-NR-genome (step 07 above) to synchronize the positions with genomic comparison. The script to use is [08-DiVir.pbs](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/08-DiVir.pbs).

```bash

while read h f i ; do r1=`ls /PATH/OshV-1-molepidemio/raw/${h}*_R1.fastq.gz`; r2=`ls /PATH/OshV-1-molepidemio/raw/${h}*_R2.fastq.gz`; qsub -v "name=${f}, reads1=${r1},reads2=${r2},outdir=,viral_genome=OshV-1-molepidemio/results/NR-Asm_consensus/C-NR-genome.fasta" OshV-1-molepidemio/src/08-DiVir.pbs ; done < OshV-1-molepidemio/raw/b-raw_metadatas/ID_experiment.csv
```

Subsequently and to save time the VCF files were transformed into an array using [VCFminer](https://github.com/Steven-N-Hart/vcf-miner) to be analyzed using R in the downstream analysis .

```bash
docker run -d -p 8888:8080 stevenhart/vcf-miner
firefox http://`hostname -I | awk '{print $1}'`:8888/vcf-miner/
```
___

## Downstream analysis

All the cleaning and analysis of the data was done directly on R, or on Graphpadprism. Once the visualizations were done, they were exported in "eps" format and then reworked on Adobe Illustrator to improve the graphic quality of each panel.

- [a-figure-01.Rmd](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/a-figure-01.Rmd) shows that the variability in sequencing depth does not affect the analysis of OsHV-1 diversity. 
- [b-figure-02.Rmd](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/b-figure-02.Rmd) shows that phylogenetic analysis of OsHV-1 non-redundant genomes reveals that the three OsHV-1 populations belong to the µVar genotype.
- 