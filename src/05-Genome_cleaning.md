# Genome cleaning

During genome assembly, we encountered some poorly assembled regions due to the genomic structure (repeats) of the OsHV-1 virus. We observed additional sequences upstream (5') of UL. Parts of the genome showed N (named polyN). As well as the junction of UL (named polyC). To understand the origin of this scaffolding and to solve them we proceeded as follows.

## Correction of sequences upstream (5') from UL

There are 56 bp upstream corresponding to the extension of the read overlap with UL. These outgoing sequences correspond to the fact that the UL is surrounded by 2 sequences in reverse complement of each other. The assembler will thus extend these sequences during the compaction of the de burjin graph.

![Upstream_corrected](https://github.com/propan2one/OshV-1-molepidemio/blob/main/image/05-Genome_cleaning-upstream_correction.png?raw=true)

Once the UL contiguous orientation is correct the sequences have been corrected using the `sed` function as follows.

```bash
sed -i "s/CTGAGTATCAATTCGAAGTAATCTCCTATACCCAAATCATTACACATCTCGTGCA//" # removed what is out of place in 5'
```

## Resolution of the polyN part

In the structure of the OsHV-1 genome there are small tandem repeated sequences of about 160 bp in the TRs (A, B) and IRs (C) regions. These regions lead to the introduction of 'N' in the sequence assembly in some samples (D).

![polyN_presentation](https://github.com/propan2one/OshV-1-molepidemio/blob/main/image/05-a-resolution_of_the_polyN_part.png?raw=true)
 
### Verification of the assembly by aligning the reads on the de novo assembly genome

First we looked at the alignment of the reads against the assembled scaffold and we were interested in the polyN part.

```bash
while read h f; do r1=`ls /PATH/OshV-1-molepidemio/raw/${h}*_R1.fastq.gz`; r2=`ls /PATH/OshV-1-molepidemio/raw/${h}*_R2.fastq.gz`; qsub -v "id=${f}, reads1=${r1}, reads2=${r2}, outdir=, Genome=, outputName=${f}" mappingfastq_bowtie.pbs ; done < OshV-1-molepidemio/raw/b-raw_metadatas/ID_experiment.csv
```

Where `mappingfastq_bowtie.pbs` correspond to classique analyse of aligning reads against a genome (OsHV-1 µVar A) as follow:

```bash
if [ ! -f ${Genome}.1.bt2 ]
then
    echo -e '\n Building genome index... \n'
    bowtie2-build $Genome $Genome
fi

echo -e "Aligning $(basename $reads1) & $(basename $reads1) against $(basename $Genome)" >> $logfile
bowtie2 -p $NCPUS -x $Genome \
-1 ${reads1} \
-2 ${reads2} | \
samtools view -F 4 -@ $NCPUS -b > ${outputName}.bam \
    2>> ${logfile}
echo -e "Aligning $(basename $reads1) & $(basename $reads1) against $(basename $Genome)" >> $logfile
samtools sort -T ${outputName}.bam -o ${outputName}_sort.bam ${outputName}.bam
samtools flagstat ${outputName}_sort.bam >> ${logfile}
```

We observed (see below) that there were no reads that could cross this ployN.

![polyN_presentation](https://github.com/propan2one/OshV-1-molepidemio/blob/main/image/05-b-resolution_of_the_polyN_part.png?raw=true)

### Verification of the assembly by aligning the reads on the previously assembled genome

Then we looked at how to align the reads in this region (560bp) of the OsHV-1 µVar A genome.

```bash
while read h f; do r1=`ls /PATH/OshV-1-molepidemio/raw/${h}*_R1.fastq.gz`; r2=`ls /PATH/OshV-1-molepidemio/raw/${h}*_R2.fastq.gz`; qsub -v "id=${f}, reads1=${r1}, reads2=${r2}, outdir=, Genome=/PATH/OshV-1-molepidemio/raw/a-OsHV-1-NCBI-genome/Ostreid_herpesvirus_1_strain_microVar_variant_A_complete_genome_0.fa, outputName=${f}" mappingfastq_bowtie.pbs ; done < OshV-1-molepidemio/raw/b-raw_metadatas/ID_experiment.csv
```

The reads are correctly aligned in this region (see below) demonstrating both that the assembled consensus genome is similar to the previously assembled genome in this section and that the assembler failed to resolve this region.

![polyN_presentation_oshv_uVar](https://github.com/propan2one/OshV-1-molepidemio/blob/main/image/05-c-resolution_of_the_polyN_part.png?raw=true)

For the samples `Thau_2018_NSI_milling_ind1_noPCR` and `Thau_2018_NSI_milling_ind3_noPCR`, which have deletions in the consensus on the polyN region present the different configurations, it will have been arbitrarily chosen to take the sequence corresponding to the OsHV-1µVar A.

## Resolution of the polyC part

In some genomes we have been able to observe the presence of incoherence which seems to be due to the presence of a variable length C trac (called polyC).

![polyC_presentation](https://github.com/propan2one/OshV-1-molepidemio/blob/main/image/05-a-resolution_of_the_polyC_part.png?raw=true)

We suspect that the C number in this region is highly variable in the OsHV-1 genome. We therefore made the choice to base the sequence selection on the observation of whole reads overlapping the 'C'.

Example of `Thau_2018_NSI_broyage_ind1_noPCR`

![polyC_presentation](https://github.com/propan2one/OshV-1-molepidemio/blob/main/image/05-b-resolution_of_the_polyC_part.png?raw=true)


Example of `Thau_2018_NSI_broyage_ind3_noPCR`

![polyC_presentation](https://github.com/propan2one/OshV-1-molepidemio/blob/main/image/05-c-resolution_of_the_polyC_part.png?raw=true)

## Modification of the sequences of the different genomes

### For the individual 1 Thau 2018 the sequences that will be modified are

#### Region PolyN

>Thau_2018_NSI_broyage_ind1_quasi-assembly_pilon_pilon
cggaaatgattgctgattggttga-------------------------------------------------------------------------------------------------------------------atgtgatatcatctacctccgcaat


>NR-genome_Ostreid_herpesvirus_1_strain_microVar_variant_A
cggaaatgattgctgattggttgaatgtgatatcatctacctccgcaatgtgttgattgggttgttattattagaccattacttcactctaaaaatataaaccaatgactgtagtcggaaatgattgctgattggttgaatgtgatatcatctacctccgcaat


#### Region Ul-polyC

This region appears to be an artifact

ACTACAGATATATCTCTAGTGTCAGCCTCTTCAATTACATTGGCAAATGTGCAGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGATCTACCATCAATTCCTTGGCCAAAAAATGGTTTAAACCAAACACACA

#### Sequence change code

```bash
export PATH=/appli/anaconda/3.7/bin/:$PATH
source activate seqkit
cd PATH/2020-04-29_final_genome
# Thau ind 1
genome=PATH/2020-04-28-asm_polish_2/Thau_2018_NSI_broyage_ind1_noPCR_asm_polished_contigs.fasta
echo -e ">Thau_2018_NSI_broyage_ind1_noPCR_NR-genome" > Thau_2018_NSI_broyage_ind1_noPCR_NR-genome.fasta
# permet de linéariser la séquence
seqkit seq --seq $genome | sed ':a;N;$!ba;s/\n//g' >> Thau_2018_NSI_broyage_ind1_noPCR_NR-genome.fasta
sed -i 's/ACTACAGATATATCTCTAGTGTCAGCCTCTTCAATTACATTGGCAAATGTGCAGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGATCTACCATCAATTCCTTGGCCAAAAAATGGTTTAAACCAAACACACA//' Thau_2018_NSI_broyage_ind1_noPCR_NR-genome.fasta
sed -i 's/CGGAAATGATTGCTGATTGGTTGAATGTGATATCATCTACCTCCGCAAT/CGGAAATGATTGCTGATTGGTTGAATGTGATATCATCTACCTCCGCAATGTGTTGATTGGGTTGTTATTATTAGACCATTACTTCACTCTAAAAATATAAACCAATGACTGTAGTCGGAAATGATTGCTGATTGGTTGAATGTGATATCATCTACCTCCGCAAT/' Thau_2018_NSI_broyage_ind1_noPCR_NR-genome.fasta
```

### For the individual 3 Thau 2018 the sequences that will be modified

#### Region PolyN

>NR-genome_Ostreid_herpesvirus_1_strain_microVar_variant_B
tcggaaatgattgctgattggttgaatgtgatatcatctacctccgcaatgtgttgattgggttgttattattagaccattacttcactctaaaaatataaacca

>Thau_2018_NSI_broyage_ind3_quasi-assembly_pilon_pilon
tcggaaatgattgctgattggttgaatg----------------------nnnnnnnttgggttgttattattagaccattacttcactctaaaaatataaacca

#### Region Ul-polyC

This region seems to be an artifact, but in the alignment there is a huge "C" stage fright (more important than in the consensus). To reflect this information, we will keep the polyC upstream and the rest will be removed.

ACACACCCCCCCCCCTCTTCCCCAATAAAAGCTTTATGGTGCCGCCTAAAACCTGAGTTTGCTAGTTTTGGACTTTTTTGATATACTATAGTAATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNN

#### Sequence change code

```bash
export PATH=/appli/anaconda/3.7/bin/:$PATH
source activate seqkit
cd PATH/2020-04-29_final_genome
# Thau ind 3
genome=PATH/2020-04-28-asm_polish_2/Thau_2018_NSI_broyage_ind3_noPCR_asm_polished_contigs.fasta
echo -e ">Thau_2018_NSI_broyage_ind3_noPCR_NR-genome" > Thau_2018_NSI_broyage_ind3_noPCR_NR-genome.fasta
# permet de linéariser la séquence
seqkit seq --seq $genome | sed ':a;N;$!ba;s/\n//g' >> Thau_2018_NSI_broyage_ind3_noPCR_NR-genome.fasta
sed -i 's/ACACACCCCCCCCCCTCTTCCCCAATAAAAGCTTTATGGTGCCGCCTAAAACCTGAGTTTGCTAGTTTTGGACTTTTTTGATATACTATAGTAATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNN//' Thau_2018_NSI_broyage_ind3_noPCR_NR-genome.fasta
sed -i 's/TCGGAAATGATTGCTGATTGGTTGAATGNNNNNNNTTGGGTTGTTATTATTAGACCATTACTTCACTCTAAAAATATAAACCA/TCGGAAATGATTGCTGATTGGTTGAATGTGATATCATCTACCTCCGCAATGTGTTGATTGGGTTGTTATTATTAGACCATTACTTCACTCTAAAAATATAAACCA/' Thau_2018_NSI_broyage_ind3_noPCR_NR-genome.fasta
```

### For the individual 5 Thau 2018 the sequences that will be modified

#### Region PolyN

>NR-genome_Ostreid_herpesvirus_1_strain_microVar_variant_B
tgctgattggtctgttattattagaccattacatttattctaaaaatataaaccaatgactgtagtcggaaatgattgctgattggttgaatgtgatatcatctacctccgcaatgtgttgattgggttgttattattagaccattacttcactctaaaaatataaaccaatgactgtagtcggaaatg

>Thau_2018_NSI_broyage_ind5_quasi-assembly_pilon_pilon
tgctgattggtctgttattattagaccattacatttattctaaaaatataa-------------------------------------------------------------------------------------------------------------------accaatgactgtagtcggaaatg

#### Sequence change code

```bash
export PATH=/appli/anaconda/3.7/bin/:$PATH
source activate seqkit
cd PATH/2020-04-29_final_genome
# Thau ind 5
genome=PATH/2020-04-28-asm_polish_2/Thau_2018_NSI_broyage_ind5_noPCR_asm_polished_contigs.fasta
echo -e ">Thau_2018_NSI_broyage_ind5_noPCR_NR-genome" > Thau_2018_NSI_broyage_ind5_noPCR_NR-genome.fasta
# permet de linéariser la séquence
seqkit seq --seq $genome | sed ':a;N;$!ba;s/\n//g' >> Thau_2018_NSI_broyage_ind5_noPCR_NR-genome.fasta
sed -i 's/TGCTGATTGGTCTGTTATTATTAGACCATTACATTTATTCTAAAAATATAAACCAATGACTGTAGTCGGAAATG/TGCTGATTGGTCTGTTATTATTAGACCATTACATTTATTCTAAAAATATAAACCAATGACTGTAGTCGGAAATGATTGCTGATTGGTTGAATGTGATATCATCTACCTCCGCAATGTGTTGATTGGGTTGTTATTATTAGACCATTACTTCACTCTAAAAATATAAACCAATGACTGTAGTCGGAAATG/' Thau_2018_NSI_broyage_ind5_noPCR_NR-genome.fasta
```

### For the individual 6 Thau 2018 the sequences that will be modified

#### Region PolyN

>NR-genome_Ostreid_herpesvirus_1_strain_microVar_variant_B
gaaatgattgctgattggttgaatgtgatatcatctacctccgcaatgtgttgattgggttgttattattagaccattacttcactctaaaaatataaaccaatgactgtagtcggaaatgattgctgattggttgaatgtgatatcatctacctccgcaatgtgttgattgggttgttattattagaccattacttcactcta

>Thau_2018_NSI_broyage_ind6_quasi-assembly_pilon_pilon
gaaatgattgctgattggttga-------------------------------------------------------------------------------------------------------------------atgtgatatcatctacctccgcaatgtgttgattgggttgttattattagaccattacttcactcta

#### Sequence change code

```bash
export PATH=/appli/anaconda/3.7/bin/:$PATH
source activate seqkit
cd PATH/2020-04-29_final_genome
# Thau ind 6
genome=PATH/2020-04-28-asm_polish_2/Thau_2018_NSI_broyage_ind6_noPCR_asm_polished_contigs.fasta
echo -e ">Thau_2018_NSI_broyage_ind6_noPCR_NR-genome" > Thau_2018_NSI_broyage_ind6_noPCR_NR-genome.fasta
# permet de linéariser la séquence
seqkit seq --seq $genome | sed ':a;N;$!ba;s/\n//g' >> Thau_2018_NSI_broyage_ind6_noPCR_NR-genome.fasta
sed -i 's/GAAATGATTGCTGATTGGTTGAATGTGATATCATCTACCTCCGCAATGTGTTGATTGGGTTGTTATTATTAGACCATTACTTCACTCTA/GAAATGATTGCTGATTGGTTGAATGTGATATCATCTACCTCCGCAATGTGTTGATTGGGTTGTTATTATTAGACCATTACTTCACTCTAAAAATATAAACCAATGACTGTAGTCGGAAATGATTGCTGATTGGTTGAATGTGATATCATCTACCTCCGCAATGTGTTGATTGGGTTGTTATTATTAGACCATTACTTCACTCTA/' Thau_2018_NSI_broyage_ind6_noPCR_NR-genome.fasta
```

### Correction of polyN of other genome assemblages

| SAMPLES                            | CORRECTION                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
|------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Brest_2018_NSI_broyage_ind10_noPCR | TTATTAGACCATTACTTCACTCTAAAAATATAAACCAATGACTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTGGGTTGTTA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| Brest_2018_NSI_broyage_ind2_noPCR  | TTATTAGACCATTACTTCACTCTAAAAATATAAACCAATGACTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTGGGTTGTTA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| Brest_2018_NSI_broyage_ind4_noPCR  | TTATTAGACCATTACTTCACTCTAAAAATATAAACCAATGACTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTGGGTTGTTA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| Brest_2018_NSI_broyage_ind6_noPCR  | TCTAAAAATATAAACCAATGACTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| Brest_2018_NSI_broyage_ind9_noPCR  | TTATTAGACCATTACTTCACTCTAAAAATATAAACCAATGACTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTGGGTTGTTA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| LT_2018_NSI_broyage_ind10_noPCR    | GACCATTACTTCACTCTAAAAATATAAACCAATGACTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTGGGTTGTTATTATTA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| LT_2018_NSI_broyage_ind1_noPCR     | AACCAATGACTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCTAAAAATATA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| LT_2018_NSI_broyage_ind3_noPCR     | GNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTGGGTTGTTATTATTAGACCATTACTTCACTCTAAAAATATAAACCAATGACT                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| LT_2018_NSI_broyage_ind4_noPCR     | ATGTGATATCATCTACCTCCGNNNNNNNNNNACCTGCCAAATACAGGATGTGCCCACCTCTTCTGATCAAGGTGCCAATTCAATACAGATTAGTGCCGCCCATAACAATCAGGAAGTGTGATCTGCCAATGGTTGAACTTGATAATGAACTATATATGTTTTTGGACAGTTATACTCCTCCTACCATGGAGGTCAATGGTGAATCTTATATGATAATTGACTCTTGTGATGCACCAACCATGAATATTAATGGAAAAACATATTATGTTATTGATGAATGTTAATTTTTTTTATTGATGAATGTTAATTTGTTGTTGTTTTATACTTGATTAAAAAAAAATGAAAATGAAGTAATATCTTTATTTATATTTTACATTGAAAACACACCCACCAAGGCGCAATCCAATAAAAATTGAGTATGGCCAGGTTAGAGTTCATGAACAATACAAATTGAGTATATTTGTACTTAGAAATTTTTCGCGCAAACCCCACAAACATATACTTTAACTTGAGATTTTTTTCTAAGTTGAAATTTTACATCATTGCGCATGTGTAAGGGAGAAACAAAATGGCGTAACAAAATGGCGGATTGTTCCCGCTATTTTATGAAACCACATTATAACGTATATTTTAAACCAATCATGTATCTACCATTGCTGCTGTTTTGCGTCATTTCCTGTTATGGTGAGCAGATAAACAACCTCGATGACCTCCAAGCAAAATTGGACAGCATGCCGCCCAGTGACTTTATTGACCACAATGGACATAATATATGCCAAGATTGCGACAGGTTATGTCCTTTAATTTCCGATAATCCCACATGTGAAGAAGATTGTTATGGGAGGTGTAACAAGGGGATAACACGGCAATGATTTAATTATTAATTCTATTTTTAAATATAACATCCAATGAAAACAGCCGGAAATTGTACCTGATTGATTTTTACCAATAAAAAAAATCAAAAGATAATGTATATTTTTTAGTATTTTTAATTTAAATTGTCAAATCAATATAGTCATCAATTATAATACAATCATCATCAATAGCATTGTTTACAGCAATATTGTCAATCATATCATTCATAAGTCTATCATTTGACATCTGCGACCTTTCGCTAAGCACCTCCTCCACCAAACTGTCCAATGTATAACGAGGAACATAATTGTCACTTTCATCATCACTAAGTAGGTGGTCTTCTTCCTCAGGTGTAAGTTCATCAATGTCATTATTATTCTCCAAAATCATGATATCATCATTGTCTATCAGCCAATACTTGTAGAGTTTAGCCTCATTTCCATAAATGATCCGCCTCTCCGGCCAGTCCTTTCTGCAGATTGGACACTTGTCTTCCATCTCCTGCCTTATGCAATTGTCACAAATGTAATGAGTTGGGCAGGGCATTGTCGTAACTCTCAATGCCCCTCCCTCAATGTCACATTCAGCCAATGTCTTTTTAGTATAAGTAGGCTCCATGCAGATAACACACACCTTCTCGACAGCATCAACTTTCGCCAACAATGTTTTTAATTTTCTGCCAGCCGCCATTGCTGTCTCTTTCGTCAGCTTCAACTTCTTAATCCCGCCCATGATCCTGGGTACGTCATGGACATGCGCAGTAGTAGTATTTCTGGACTTTACCATTTTCACAGCGTATTTTATATTCCTATATCTGCTAACAGTAGATGTGTACAGTAAGACTCCCAAGCAAGTTATGATTTTAACATGGGCAATAGCGCTGCTTATATAGAGCAGGCCAGGGTGCTGATTGGTCTGTTATTATTAGACCATTACATTTATTCTAAAAATATAAACCAATGACTGTAGTCGGAAATGATTGCTGATTGGTTGA |
| LT_2018_NSI_broyage_ind7_noPCR     | TGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCTAAAAATATAAACCAATGAC                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| LT_2018_NSI_broyage_ind8_noPCR     | CAATGACTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTGGGTTGTTATTATTAGACCATTACTTCACTCTAAAAATATAAAC                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| LT_2018_NSI_broyage_ind9_noPCR     | AACCAATGACTGNNNNNTTGGGTTGTTATTATTAGACCATTACTTCACTCTAAAAATATA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| Thau_2018_NSI_broyage_ind10_noPCR  | GACCATTACTTCACTCTAAAAATATAAACCAATGACTGNNTTGGGTTGTTATTATTA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| Thau_2018_NSI_broyage_ind4_noPCR   | GACCATTACTTCACTCTAAAAATATAAACCAATGACTGNTTGGGTTGTTATTATTA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| Thau_2018_NSI_broyage_ind7_noPCR   | GACCATTACTTCACTCTAAAAATATAAACCAATGACTGNNNTTGGGTTGTTATTATTA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |

#### Sequence change code

```bash
#qsub -I -q omp -l walltime=05:00:00 -l mem=90g -l ncpus=8
export PATH=/appli/anaconda/3.7/bin/:$PATH
source activate seqkit
while read h f; do genome=PATH/2020-04-29_final_genome/${h}_NR-genome.fasta ; echo ">${h}_NR-genome" > $genome ; seqkit seq --seq PATH/2020-04-28-asm_polish_2/${h}_asm_polished_contigs.fasta | sed ':a;N;$!ba;s/\n//g' >> $genome ; sed -i "s/$f//" $genome ; done < OshV-1-molepidemio/raw/d-Manually_corrected_assembly/correction_polyN.csv

## Added des Thau_2018_NSI_broyage_ind8_noPCR Thau_2018_NSI_broyage_ind9_noPCR which was direct well assembled

genome=PATH/2020-04-28-asm_polish_2/Thau_2018_NSI_broyage_ind8_noPCR* 
echo -e ">Thau_2018_NSI_broyage_ind8_noPCR_NR-genome" > Thau_2018_NSI_broyage_ind8_noPCR_NR-genome.fasta
seqkit seq --seq $genome | sed ':a;N;$!ba;s/\n//g' >> Thau_2018_NSI_broyage_ind8_noPCR_NR-genome.fasta

genome=PATH/2020-04-28-asm_polish_2/Thau_2018_NSI_broyage_ind9_noPCR*
echo -e ">Thau_2018_NSI_broyage_ind9_noPCR_NR-genome" > Thau_2018_NSI_broyage_ind9_noPCR_NR-genome.fasta
seqkit seq --seq $genome | sed ':a;N;$!ba;s/\n//g' >> Thau_2018_NSI_broyage_ind9_noPCR_NR-genome.fasta
```