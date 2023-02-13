# Analyses

# Genome
## PacBio
### Preliminary
BAM --> FASTQ
cat *.fastq > all_pacbio_libraries.fq

### LRBinner
```bash
LRBinner -r all_pacbio_libraries.fq -o lrb_output_all_PB_libs --ae-epochs 200 --resume -mbs 1000 -bit 0 -bs 10 -bc 10 --threads 56

separate_reads.py --reads all_pacbio_libraries.fq --bins lrb_output_all_PB_libs/binning_result.pkl --outpath seperated_reads
```
#### Assembly of Bins
flye assembly on separate bins

#### BUSCO Lineage Assignment
busco auto-lineage on all bins

##### Alveolata
* bin-0
* bin-1
* bin-2

##### Bacterial
* bin-3
* bin-5
* bin-7

##### Unclassified
* bin-4
* bin-6
* bin-8 to bin-22

### Assembly of Combined 'Alveolata' Identifed Bins
```bash
flye --threads 56 --meta --pacbio-raw bin-012reads.fastq.gz -o bin-012reads_assembly_raw
```

#### Rough Clean & Mask
```bash
funannotate clean -i assembly.fasta -o assembly_cleaned.fasta --exhaustive --cpus 56
funannotate sort -i assembly_cleaned.fasta -o assembly_cleaned_sorted.fasta

RepeatMasker -e hmmer -pa 28 -qq -xsmall -gff --species Paramecium -dir .  assembly_cleaned_sorted.fasta
mkdir repeatmasker
mv assembly_cleaned_sorted.fasta.* repeatmasker/

ln -s repeatmasker/assembly_cleaned_sorted.fasta.masked assembly_cleaned_sorted_masked.fasta

ln -s assembly_cleaned_sorted_masked.fasta assembly.fasta
```

### Pilon
Two rounds of Pilon with Illumina Nova-Seq Libraries

### Funannotate
Had to modify the funnanotate scripts (train/library) to force PASA + Trinity/Trinotate scripts to enable -G Ciliate mode. Also edited the funannoate library to enable --gcode 6 for GeneMark-ES v4.71. Very hacky, would be nice to fix, add CLI option and add a pull request. Update - it's too messy I think. I had to edit a bunch of other scripts EVM and P2G and the main "translate" code (not using biopython) too! Too many options with different names: gcode 6, --stops ATG, -G Ciliate. Eugh. 

Also had to add this to the headers so that tbl2asn knows to translate stop codons properly at the end of 'predict'. Annoying this can't be passed as an option - wtf NCBI. This also messes up GeneMark. Sigh. 

```bash
sed -i 's/>.*/& [gcode=6]/' genome_no_bac.fasta
```

```bash
funannotate train -i genome_no_bac.fasta -o training --left 1A_S28_R1_001.fastq.gz  1B_S29_R1_001.fastq.gz  1C_S30_R1_001.fastq.gz  1E_S31_R1_001.fastq.gz  3A_S32_R1_001.fastq.gz  3C_S33_R1_001.fastq.gz  3D_S34_R1_001.fastq.gz  3E_S35_R1_001.fastq.gz  5B_S36_R1_001.fastq.gz  5C_S37_R1_001.fastq.gz  5D_S38_R1_001.fastq.gz    --right 1A_S28_R2_001.fastq.gz  1B_S29_R2_001.fastq.gz  1C_S30_R2_001.fastq.gz  1E_S31_R2_001.fastq.gz  3A_S32_R2_001.fastq.gz  3C_S33_R2_001.fastq.gz  3D_S34_R2_001.fastq.gz  3E_S35_R2_001.fastq.gz  5B_S36_R2_001.fastq.gz  5C_S37_R2_001.fastq.gz  5D_S38_R2_001.fastq.gz --stranded RF --pacbio_isoseq clustered.hq.fasta.gz --jaccard_clip --memory 1000G -c 25 --species "Paramecium bursaria" --strain "186b" --cpus 56

funannotate predict -i genome_no_bac.fasta -o training -s "Paramecium bursaria" --strain "186b" --cpus 56 --busco_db alveolata_stramenophiles --transcript_evidence ~/guy/pb_isoseq/final/3_final_set/pb_isoseq_collapsed_isoforms_all.rep.fa --genemark_gtf ../gmes_petap_P.bur_186b/genemark.gtf --repeats2evm -d /databases/funannotate/1.8.13/ --SeqCenter "Exeter Sequencing Service" --name "XXXXXX" --organism other --keep_no_stops --busco_seed_species tetrahymena  2>&1 | tee 4_predict.out
```

### Telomeres
```bash
telomeric_repeats.fasta
>5_prime
CCCCAACCCCAA
>3_prime
TTGGGGTTGGGG
```

```bash
blastn -task blastn-short -db assembly.fasta -query telomeric_repeats.fasta -outfmt '6 std slen' -out telomere_scaffolds.blastn

# this finds matches that occur within the first 12bp of a contig
awk '$9>=1 && $9<=50{print $0}' telomere_scaffolds.blastn | cut -f 2 | sort | uniq >identified_telomere_starts.list
# sam for the ends
awk '$9<=$13 && $9>=($13-50){print $0}' telomere_scaffolds.blastn | cut -f 2 | sort | uniq >identified_telomere_ends.list

comm -12 identified_telomere_starts.list identified_telomere_ends.list >identified_telomeres_complete.list
comm -13 identified_telomere_starts.list identified_telomere_ends.list >identified_telomere_start_only.list
comm -23 identified_telomere_starts.list identified_telomere_ends.list >identified_telomere_end_only.list

faSomeRecords v4_new_pilon2_assembly_cleaned_sorted_masked.fasta identified_telomeres_complete.list identified_telomeres_complete.fasta
faSomeRecords v4_new_pilon2_assembly_cleaned_sorted_masked.fasta identified_telomere_start_only.list identified_telomere_start_only.fasta
faSomeRecords v4_new_pilon2_assembly_cleaned_sorted_masked.fasta identified_telomere_end_only.list identified_telomere_end_only.fasta

telomeres=(scaffold_149 scaffold_170 scaffold_201 scaffold_230 scaffold_372 scaffold_421 scaffold_477 scaffold_493 scaffold_541 scaffold_551 scaffold_561 scaffold_578 scaffold_607 scaffold_621 scaffold_631 scaf>
for i in "${telomeres[@]}"; do   echo "${i}";   sed -i "s/^>${i}$/>${i}_putative_chromosome/" assembly.fasta; done

telomeres_start=(scaffold_10 scaffold_1001 scaffold_1030 scaffold_1038 scaffold_1096 scaffold_11 scaffold_116 scaffold_137 scaffold_159 scaffold_168 scaffold_174 scaffold_177 scaffold_186 scaffold_19 scaffold_2>
for i in "${telomeres_start[@]}"; do   echo "${i}";   sed -i "s/^>${i}$/>${i}_partial_chr_start/" assembly.fasta; done

telomeres_end=(scaffold_1046 scaffold_112 scaffold_113 scaffold_118 scaffold_126 scaffold_138 scaffold_144 scaffold_145 scaffold_146 scaffold_15 scaffold_157 scaffold_16 scaffold_163 scaffold_18 scaffold_181 sc>
for i in "${telomeres_end[@]}"; do   echo "${i}";   sed -i "s/^>${i}$/>${i}_partial_chr_end/" assembly.fasta; done
```

### Mitochondria
```bash
mitofinder -a assembly.fasta -p 56 -r p_caudatum_mitochondria.gb --new-genes --allow-intron --numt --intron-size 35 --max-contig-size 60000 -o 6 -j pb_mito


mitochondrial=(scaffold_350 scaffold_862 scaffold_698)
for i in "${mitochondrial[@]}"; do   echo "${i}";   sed -i "s/^>${i}$/>${i}_putative_mito/" assembly.fasta; done
```

## Illumina Novaseq

### Adapter Trimming and QC
```bash
fastp -i 10024_Paramecium_r1.fq.gz -I 10024_Paramecium_r2.fq.gz -o 10024_Paramecium_r1_trimmed.fq.gz -O 10024_Paramecium_r2_trimmed.fq.gz --unpaired1 10024_Paramecium_unpaired_trimmed.fq.gz --unpaired2 10024_Paramecium_unpaired_trimmed.fq.gz --detect_adapter_for_pe --trim_poly_g -c -x -w 16
```

# Iso-Seq

## isoseq3
```bash
# identify and remove the 5'/3' cDNA primers
lima --isoseq --peek-guess --dump-clips -j 112 ../m64176e_220706_095424.hifi_reads.bam ../primers.fasta hifi_reads.bam

# remove polyA tails and artificial concatemers
isoseq3 refine --require-polya -j 112 hifi_reads.consensusreadset.xml ../primers.fasta flnc.bam

# cluster Full Length Non Chimeric/Concatemer reads
isoseq3 cluster -j 112 flnc.bam clustered.bam --verbose --use-qvs
# outputs: clustered.hq.fasta.gz = ALL transcripts, but we don't want to use this file yet...

# Map polished transcripts to v5 genome assembly
pbmm2 align -j 112 --preset ISOSEQ --sort clustered.bam pb186b_assembly_v5.fasta clustered_mapped_186b_v5.bam

# collapse isoforms
isoseq3 collapse -j 112 --do-not-collapse-extra-5exons clustered_mapped_186b_v5.bam collapsed.gff
```

## Transdecoder
```bash
# extract the long open reading frames
TransDecoder.LongOrfs -t ../3_final_set/pb_isoseq_collapsed_isoforms_all.rep.fa -G Ciliate

# run blast
blastp -query pb_isoseq_collapsed_isoforms_all.rep.fa.transdecoder_dir/longest_orfs.pep -db /databases/uniprot/uniprot_sprot.fasta -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 >longest_peps_vs_uniprot_sprot.tab

# run pfam
hmmscan --cpu 56 --domtblout pfam.domtblout /databases/pfam/35/Pfam-A.hmm ../pb_isoseq_collapsed_isoforms_all.rep.fa.transdecoder_pb_isoseq_collapsed_isoforms_all.rep.fa.transdecoder_dir/longest_orfs.pep

# predict coding regions
TransDecoder.Predict -t ../3_final_set/pb_isoseq_collapsed_isoforms_all.rep.fa -G Ciliate --single_best_only --retain_blastp_hits blast/longest_peps_vs_uniprot_sprot.tab --retain_pfam_hits pfam/pfam.domtblout

# headers are a mess
seqkit replace -p '(.+)' -r '{nr}' pb_isoseq_collapsed_isoforms_all.rep.fa.transdecoder.pep >pb_isoseq_collapsed_isoforms_all.rep.fa.transdecoder.renamed.pep
```

```bash
pblat P.bur_186b.fas pb_isoseq_collapsed_isoforms_all.rep.fa.transdecoder.cds -out=blast8 cds_vs_iso.out -threads=112
awk '$3>90{print $0}' cds_vs_iso.out | cut -f1 | sort | uniq | wc -l
```

## BUSCO

### peptides
```bash
busco -m prot -i ./pb_isoseq_collapsed_isoforms_all.rep.fa.transdecoder.renamed.pep --metaeuk_parameters="--translation-table=6" --metaeuk_rerun_parameters="--translation-table=6" -c 56 -o busco_prot_auto_isoseq --auto-lineage
# eukaryota_odb10: C:55.3%[S:32.9%,D:22.4%],F:3.1%,M:41.6%,n:255
# alveolata_odb10: C:86.0%[S:42.1%,D:43.9%],F:0.6%,M:13.4%,n:171
# archaea_odb10:   C:25.8%[S:10.5%,D:15.3%],F:8.9%,M:65.3%,n:124
# bacteria_odb10:  C:46.9%[S:24.7%,D:22.2%],F:5.2%,M:47.9%,n:1
```
### transcriptome
```bash
# headers are a mess
seqkit replace -p '(.+)' -r '{nr}' pb_isoseq_collapsed_isoforms_all.rep.fa >pb_isoseq_collapsed_isoforms_all.rep.renamed.fa

# busco
busco -m tran -i ./pb_isoseq_collapsed_isoforms_all.rep.renamed.fa --metaeuk_parameters="--translation-table=6" --metaeuk_rerun_parameters="--translation-table=6" -c 56 -o busco_tran_auto_isoseq --auto-lineage
# eukaryota_odb10: C:55.7%[S:28.6%,D:27.1%],F:3.5%,M:40.8%,n:255
# alveolata_odb10: C:87.1%[S:36.8%,D:50.3%],F:0.6%,M:12.3%,n:171
# archaea_odb10:   C:09.8%[S:04.1%,D:05.7%],F:0.0%,M:90.2%,n:194
# bacteria_odb10:  C:04.0%[S:00.8%,D:03.2%],F:1.6%,M:94.4%,n:12
```
