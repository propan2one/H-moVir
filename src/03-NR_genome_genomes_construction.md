# 03) Creation of NR-genomes

- In a first step, the sequences were analyzed on Ugene with a dotplot of the sequence against itself to identify the different structures of the genomes (min length 40, repeats identity 95 %). Each repetition was manually annotated as follows: 1-Ul, 2-IRl, 3-X, 4-IRs, 5-Us. 

- In a second step, the annotations were then saved in a csv and concatenate in csv format [NR_genomic_part_coordonate.csv](https://github.com/propan2one/OshV-1-molepidemio/blob/main/raw/a-OsHV-1-NCBI-genome/NR_genomic_part_coordonate.csv) to be cut out using the [seqkit](https://bioinf.shenwei.me/seqkit/) tool and the `seqkit subseq -r` option as follow:

```bash
# extract seq
while read F N S E L C G ; do cat $G.fasta | seqkit subseq -r $S:$E | seqkit seq --seq | sed ':a;N;$!ba;s/\n//g' > ${N}_$G.fasta ; done < NR_genomic_part_coordonate.csv
# fill a folder of them
while read F N S E L C G ; do if [ ! -d $G ]; then mkdir $G ; fi ; mv ${N}_$G.fasta $G ; done < NR_genomic_part_coordonate.csv
```

- In a third step, the sequences will be concatenated to create the non-redundant genomes.

```bash
cd OshV-1-molepidemio/raw/a-OsHV-1-NCBI-genome

# µVar A
cd Ostreid_herpesvirus_1_strain_microVar_variant_A_complete_genome
output=NR-genome_Ostreid_herpesvirus_1_strain_microVar_variant_A.fa
for seq in {Ul,IRl,X,IRs,Us} ; do cat $seq* >> ${output}_transi; done
sed -i ':a;N;$!ba;s/\n//g' ${output}_transi
echo -e ">${output%.fa}" > $output
cat ${output}_transi >> $output
rm ${output}_transi
cd ..


# µVar B
cd Ostreid_herpesvirus_1_strain_microVar_variant_B_complete_genome
output=NR-genome_Ostreid_herpesvirus_1_strain_microVar_variant_B.fa
for seq in {Ul,IRl,X,IRs,Us} ; do cat $seq* >> ${output}_transi; done
sed -i ':a;N;$!ba;s/\n//g' ${output}_transi
echo -e ">${output%.fa}" > $output
cat ${output}_transi >> $output
rm ${output}_transi
cd ..

# Davison
cd Ostreid_herpesvirus_1_complete_genome
output=NR-genome_Ostreid_herpesvirus_1.fa
for seq in {Ul,IRl,X,IRs,Us} ; do cat $seq* >> ${output}_transi; done
sed -i ':a;N;$!ba;s/\n//g' ${output}_transi
echo -e ">${output%.fa}" > $output
cat ${output}_transi >> $output
rm ${output}_transi
cd ..

# PT
cd Ostreid_herpesvirus_1_2016_PT
output=NR-genome_Ostreid_herpesvirus_1_2016_PT.fa
for seq in {Ul,IRl,IRs,Us} ; do cat $seq* >> ${output}_transi; done # modif pas de X labelisé ici
sed -i ':a;N;$!ba;s/\n//g' ${output}_transi
echo -e ">${output%.fa}" > $output
cat ${output}_transi >> $output
rm ${output}_transi
cd ..

# ZK0118
cd Ostreid_herpesvirus_1_isolate_ZK0118_complete_genome
output=NR-genome_Ostreid_herpesvirus_1_isolate_ZK0118.fa
for seq in {Ul,IRl,X,IRs,Us} ; do cat $seq* >> ${output}_transi; done
sed -i ':a;N;$!ba;s/\n//g' ${output}_transi
echo -e ">${output%.fa}" > $output
cat ${output}_transi >> $output
rm ${output}_transi
cd ..

## J I C : the verification of pairwaise is ok and the beginning of the sequence on ugene corresponds well to the beginning of the UL
#stretcher \
#    -asequence Ostreid_herpesvirus_1_complete_genome.fasta \
#    -bsequence NR-genome_Ostreid_herpesvirus_1.fa \
#    -outfile VERIFICATION_NR

# Creating Multiple Alignment
cd OshV-1-molepidemio/raw/a-OsHV-1-NCBI-genome
cat OshV-1-molepidemio/raw/a-OsHV-1-NCBI-genome/*/NR-genome_*.fa > multiple_fasta_NR-genome_OsHV.faa
```