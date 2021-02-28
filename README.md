# Genomic epidemiology, genetic diversity and connectivity of the Ostreid Herpesvirus 1 population in France

Scripts used for recreating figures and analyses for [DOI](). Scripts were not coded with the intention of reuse by others. Provided as a courtesy for better transparency and reproducibility of research.

## Abstract

The genetic diversity of viral populations has become a major issue in the understanding of their phylogeographic and dissemination history, but it represents a real challenge when studying whole-genome viral diversity in the natural environment. These molecular ecological approaches to study viral diversity are commonly used for RNA viruses harboring small genomes but have not been applied to DNA viruses with large genomes. In this study we used the “Pacific Oyster Mortality Syndrome” (POMS, a disease that affects oyster farms around the world) as a model to study the genetic diversity of its causative agent, the Ostreid herpes virus 1 (OsHV-1) in the three most important French oyster-farming areas. Using ultra-deep sequencing on individual moribund oysters and new bioinformatics methodology we de novo assembled twenty one OsHV-1 genomes. Then, using a combination of whole-genomes comparison, phylogenetic analysis and quantification of major and minor variants we unveiled the connectivity of OsHV-1 viral populations between the three oyster-farming areas. We propose a scenario in which the Marennes-Oléron basin would be the main source of OsHV-1 diversity which would then be dispersed to other farming areas, a hypothesis consistent with the current knowledge on oyster transfers in France. In conclusion, this study showed that molecular ecological approaches may be applied to large genome virus to unravel the extent of their genetic diversity and better understand the spread of viral populations in natural environments.

![graphical_summary](https://github.com/propan2one/OshV-1-molepidemio/blob/main/image/graphical_summary.png?raw=true)

## Ressources

Login and passwords are missing in the code. The paths used in the code have also been replaced for security reasons. They can be found as follows where the variabe `PATH` corresponds to the location where the repository was cloned:

### Previously published OsHV genome

The virus genomes used in the study were downloaded directly from NCBI as follows:

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

### 01) Data transfer

All data has been transferred and MD5 have been check from the servers of the transfer platform to the datarmor calculation server with : [01-data_transfer.pbs](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/01-data_transfer.pbs)

### 02) Genetic diversity rarefaction analysis

The analysis of genetic diversity rarefaction was performed for all samples individually and the command was executed as follows using this script: [02-rarefaction_virus.pbs](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/02-rarefaction_virus.pbs)

```bash
while read h f; do r1=`ls /PATH/OshV-1-molepidemio/raw/${h}*_R1.fastq.gz`; r2=`ls /PATH/OshV-1-molepidemio/raw/${h}*_R2.fastq.gz`; qsub -v "rate=5000, id=${f}, reads1=${r1}, reads2=${r2}, outdir=/home1/scratch/jdelmott/2020-02-12-Rarefaction_Haplofit/, genomefile=/PATH/OshV-1-molepidemio/raw/a-OsHV-1-NCBI-genome/Ostreid_herpesvirus_1_strain_microVar_variant_A_complete_genome_0.fa, gffFile=/PATH/OshV-1-molepidemio/raw/a-OsHV-1-NCBI-genome/Ostreid_herpesvirus_1_strain_microVar_variant_A_complete_genome_0.gff" OshV-1-molepidemio/src/02-rarefaction_virus.pbs; done < OshV-1-molepidemio/raw/b-raw_metadatas/ID_experiment.csv
```

### 03) Creation of NR-genomes from NCBI

- In a first step, the sequences were analyzed on Ugene with a dotplot of the sequence against itself to identify the different structures of the genomes. Each repetition was manually annotated as follows: 1-Ul, 2-IRl, 3-X, 4-IRs, 5-Us. 

- In a second step, the annotations were then saved in a csv and concatenate in csv format [NR_genomic_part_coordonate.csv](https://github.com/propan2one/OshV-1-molepidemio/blob/main/raw/a-OsHV-1-NCBI-genome/NR_genomic_part_coordonate.csv) to be cut out using the [seqkit](https://bioinf.shenwei.me/seqkit/) tool and the `seqkit subseq -r` option.

- In a third step, the sequences will be concatenated to create the non-redundant genomes. 

See [03-NR_genome_genomes_construction.md](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/03-NR_genome_genomes_construction.md) for more details.


Checks were performed using multiple alignment using [mafft](https://mafft.cbrc.jp/alignment/software/) en utilisant le script [mafft_MSA.pbs](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/mafft_MSA.pbs). Visualization of the alignment was done with [Aliview](https://github.com/AliView/AliView).


### 04) Assembly of the genome of the virus of interest

The construction of a genome for each was carried out for all samples individually and the command was executed as follows using this script: [04-metaviromics.pbs](https://github.com/propan2one/OshV-1-molepidemio/blob/main/src/04-metaviromics.pbs)

```bash
while read h f i ; do r1=`ls /PATH/OshV-1-molepidemio/raw/${h}*_R1.fastq.gz`; r2=`ls /PATH/OshV-1-molepidemio/raw/${h}*_R2.fastq.gz`; qsub -v "id=${f}, reads1=${r1},reads2=${r2}, genomefile=/home1/datawork/jdelmott/data_jean/oyster.v9.fa,database=/home1/datawork/jdelmott/data_jean/viral.2.1.genomic.fna,GenomeOsHV1=/home1/datawork/jdelmott/data_jean/OsHV-1_strain_microVar_variant_A.fasta,mincontig=200,minlength=50,outdir=/home1/scratch/jdelmott/2020-03-20-Haplofit_metaviromic,insersize=${i}" OshV-1-molepidemio/src/04-metaviromics.pbs; done < OshV-1-molepidemio/raw/b-raw_metadatas/ID_experiment.csv
```

## Downstream analysis

