# Genome cleaning

During genome assembly, we encountered some poorly assembled regions due to the genomic structure (repeats) of the OsHV-1 virus. We observed sequences upstream (5') of UL. Parts of the genome showed N (named polyN). As well as the junction of UL named polyC). To understand the origin of this scaffolding and to solve them we proceeded as follows.

## Correction of sequences upstream (5') from UL

```bash
#### En vrai on sait pas dans quelle sens est le Us... à ce niveau la donc j'ai lancer le pipeline sur tout pour les contigs
ls */*/14*_NR_genome_contigs_raw.faa
seqkit seq -r -p Us.fasta > rev_comp.fasta

sed -i "s/CTGAGTATCAATTCGAAGTAATCTCCTATACCCAAATCATTACACATCTCGTGCA//" # enlevé ce qui dépace en 5'

for folder in *_noPCR* ; do mkdir $folder/R01-${folder%_noPCR}_selected_contigs ; done
#for folder in *_noPCR* ; do rm $folder/${folder%_noPCR} ; done
for folder in *_noPCR* ; do cp $folder/output_fasta/* $folder/R01-${folder%_noPCR}_selected_contigs; done
```

## Verification of the assembly by aligning the reads on the de novo assembly genome

```bash
# Mapping on IRl de OsHV-1 µVar A to observe the "NN" of the samples from the thau lagoon
while read h f i ; do r1=`ls /home1/scratch/jdelmott/2020-04-20-Metaviromic_Haplofit/$f/02-trimmomatic/${h}*R1_trim.fastq.gz`; r2=`ls /home1/scratch/jdelmott/2020-04-20-Metaviromic_Haplofit/$f/02-trimmomatic/${h}*R2_trim.fastq.gz`; qsub -v "outdir=/home1/datawork/jdelmott/2020-04-07-Metaviromic_Haplofit_backup/2020-04-29-Thau_pb_v2, reads1=${r1},reads2=${r2}, 
Genome=/home1/datawork/jdelmott/data_jean/OsHV-multiple-region/IRl_Ostreid_herpesvirus_1_strain_microVar_variant_A.fasta, outputName=${f}" ~/mappingfastq_bowtie.pbs ; done < /home1/datawork/jdelmott/2020-04-07-Metaviromic_Haplofit_backup/2020-04-28-Thau_problems/ID_Thau_pb.csv
```