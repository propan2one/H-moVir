# Reconstruction of genomes in full version

The objective will be to redefine the different genome structures by aligning the sequences of OsHV-1 µVar A.







The NR genome structure being 1-Ul, 2-IRl, 3-X, 4-IRs, 5-Us, I will have to identify these 5 regions in each of the genomes, extract them and duplicate them in rev-complement.

https://github.com/propan2one/OshV-1-molepidemio#03-creation-of-non-rredundant-genomes-from-ncbi

```bash
# Working directory
WD=~/Documents/OshV-1-molepidemio/results/Genomes_asm_dir_test
cd $WD
# Download fasta format (exactly like #3)
id=KY242785.1
name=Ostreid_herpesvirus_1_strain_microVar_variant_A_complete_genome
ncbi-acc-download -F fasta ${id} -p ${name}

# Genome slicing
## Recovery of genomic coordinates
grep -e Ul \
    -e IRl \
    -e X \
    -e IRs \
    -e Us -- ~/Documents/OshV-1-molepidemio/raw/a-OsHV-1-NCBI-genome/NR_genomic_part_coordonate.csv | \
    grep -e "Ostreid_herpesvirus_1_strain_microVar_variant_A_complete_genome" > genomic_coordonate_uVarA.csv

## Fragments cutting
while read F N S E L C G ; do cat Ostreid_herpesvirus_1_strain_microVar_variant_A_complete_genome_0.fa | \
    seqkit subseq -r $S:$E >> ${N}_raw.fa ; done < genomic_coordonate_uVarA.csv

for file in *_raw.fa ; do sed -i "s/KY242785.1 Ostreid herpesvirus 1 strain microVar variant A, complete genome/${file%_raw.fa}/" $file ; done

cat *_raw.fa > structure_oshv_uVar.fna

# Visual check
seqkit fx2tab structure_oshv_uVar.fna -l -g -n -i -H
##name	seq	qual	length	GC
#IRl		7338	41.54
#IRs		9776	40.20
#Ul			164268	38.23
#Us			3370	37.12
#X			1510	38.28
```

## Determination of the location of the different genomic structures

In general we can observe using `seqkit subseq -r 1:50 *raw.fa` and `seqkit subseq -r -50:-1 *raw.fa` that:

>Ul_start
GTCTTACCACCGTTATCATCCCTGTCTCTATTTTTTTATACAATTCTTTT
>Ul_end
ACCCATGTATATATGTATCATTGGCGATATACTTAGATTTTCCTTGGTAT

>IRl_start   
TGCACGAGATGTGTAATGATTTGGGTATAGGAGATTACTTCGAATTGATA
>IRl_end
GGAGGGGGAGGAGGTGGGGAGGATGGGGGAGGTGTTGGGGAGGTGGGGGG

>X_start
GCTGTTATATGAATTGAGTAAAGGTAAGAAAATTGCCCTCGCCGTGCCAA
>X_end
TGTTGTTGTGGGATGGGGGAGATTTGAGGCGATTTTGCTGTTGTTTGCTC

>IRs_start
TCCCCAACCCCCCTCTAAACCCCCCTCTCCCCACCACCTCACCACCCCCT
>IRs_end
GTTGTAATCATGGTTTATTTTTAGGCAAGGATGTGTGCGCATTTCACAAT

>Us_start
ACAAAACATTTCAAGAATAACTTTTAAAACTGTTGATATGTCACGTTATT
>Us_end
GAATATGCCCCATTGACGTCAATACTCCTTATATGGCTGGTAAACATCCG


### NR genome Brest 2018 ind2

- `GAAAACGACATA` 12bp between Ul and IRl have been add to UL du to size
- 

```bash
mkdir blast ; cd $WD/blast
ln -s $WD/structure_oshv_uVar.fna .
# mkdb of structures
makeblastdb -in structure_oshv_uVar.fna -parse_seqids -dbtype nucl

db=$WD/blast/structure_oshv_uVar.fna
genome=~/Documents/OshV-1-molepidemio/results/NR-Asm-genome/NR_genome_Brest_2018_NSI_broyage_ind2_noPCR.fasta
blastn \
    -query $genome \
    -db $db \
    -word_size 450 \
    -num_threads 4 \
    -outfmt "6 sseqid qseqid length pident sstart send qstart qend evalue qlen" \
    -evalue 0.0001 \
    -out $(basename ${genome%.fasta}).blastout
# sseqid qseqid                                     length  pident sstart   send    qstart  qend    evalue qlen
# Ul	NR_genome_Brest_2018_NSI_broyage_ind2_noPCR	164285	99.92	1	    164268	1	    164267	0.0	186286
# IRl	NR_genome_Brest_2018_NSI_broyage_ind2_noPCR	7354	99.61	1	    7338	164280	171629	0.0	186286
# X	    NR_genome_Brest_2018_NSI_broyage_ind2_noPCR	1508	99.87	3	    1510	171635	173142	0.0	186286
# IRs	NR_genome_Brest_2018_NSI_broyage_ind2_noPCR	9779	99.83	1	    9776	173143	182915	0.0	186286
# Us	NR_genome_Brest_2018_NSI_broyage_ind2_noPCR	3371	99.85	1	    3369	182916	186286	0.0	186286
```

### NR genome Brest 2018 ind10

- `GAAAACGACATA` 12bp between Ul and IRl have been add to UL du to size
- 

```bash
mkdir blast ; cd $WD/blast
ln -s $WD/structure_oshv_uVar.fna .
# mkdb of structures
makeblastdb -in structure_oshv_uVar.fna -parse_seqids -dbtype nucl

db=$WD/blast/structure_oshv_uVar.fna
genome=~/Documents/OshV-1-molepidemio/results/NR-Asm-genome/NR_genome_Brest_2018_NSI_broyage_ind10_noPCR.fasta
blastn \
    -query $genome \
    -db $db \
    -word_size 450 \
    -num_threads 4 \
    -outfmt "6 sseqid qseqid length pident sstart send qstart qend evalue qlen" \
    -evalue 0.0001 \
    -out $(basename ${genome%.fasta}).blastout
# sseqid qseqid                                           length  pident sstart   send    qstart  qend    evalue qlen
# Ul      NR_genome_Brest_2018_NSI_broyage_ind10_noPCR    164284  99.92   1       164268  1       164266  0.0     186279
# IRs     NR_genome_Brest_2018_NSI_broyage_ind10_noPCR    9779    99.83   1       9776    173136  182908  0.0     186279
# IRl     NR_genome_Brest_2018_NSI_broyage_ind10_noPCR    7349    99.67   1       7338    164279  171623  0.0     186279
# Us      NR_genome_Brest_2018_NSI_broyage_ind10_noPCR    3371    99.85   1       3369    182909  186279  0.0     186279
# X       NR_genome_Brest_2018_NSI_broyage_ind10_noPCR    1508    99.80   3       1510    171628  173135  0.0     186279
```





NR_genome_Brest_2018_NSI_broyage_ind10_noPCR.fasta
NR_genome_Brest_2018_NSI_broyage_ind4_noPCR.fasta
NR_genome_LT_2018_NSI_broyage_ind10_noPCR.fasta
NR_genome_Brest_2018_NSI_broyage_ind9_noPCR.fasta
NR_genome_Brest_2018_NSI_broyage_ind6_noPCR.fasta
NR_genome_LT_2018_NSI_broyage_ind1_noPCR.fasta
NR_genome_LT_2018_NSI_broyage_ind3_noPCR.fasta
NR_genome_LT_2018_NSI_broyage_ind7_noPCR.fasta
NR_genome_LT_2018_NSI_broyage_ind4_noPCR.fasta
NR_genome_LT_2018_NSI_broyage_ind9_noPCR.fasta
NR_genome_LT_2018_NSI_broyage_ind8_noPCR.fasta
NR_genome_Thau_2018_NSI_broyage_ind10_noPCR.fasta
NR_genome_Thau_2018_NSI_broyage_ind4_noPCR.fasta
NR_genome_Thau_2018_NSI_broyage_ind3_noPCR.fasta
NR_genome_Thau_2018_NSI_broyage_ind1_noPCR.fasta
NR_genome_Thau_2018_NSI_broyage_ind7_noPCR.fasta
NR_genome_Thau_2018_NSI_broyage_ind6_noPCR.fasta
NR_genome_Thau_2018_NSI_broyage_ind5_noPCR.fasta
NR_genome_Thau_2018_NSI_broyage_ind9_noPCR.fasta
NR_genome_Thau_2018_NSI_broyage_ind8_noPCR.fasta
