# Analyses

# Genome
## Preliminary
BAM --> FASTQ
cat *.fastq > all_pacbio_libraries.fq

## LRBinner
```bash
LRBinner -r all_pacbio_libraries.fq -o lrb_output_all_PB_libs --ae-epochs 200 --resume -mbs 1000 -bit 0 -bs 10 -bc 10 --threads 56

separate_reads.py --reads all_pacbio_libraries.fq --bins lrb_output_all_PB_libs/binning_result.pkl --outpath seperated_reads
```
### Assembly of Bins
flye assembly on separate bins

### BUSCO Lineage Assignment
busco auto-lineage on all bins

#### Alveolata
* bin-0
* bin-1
* bin-2

#### Bacterial
* bin-3
* bin-5
* bin-7

#### Unclassified
* bin-4
* bin-6
* bin-8 to bin-22

## Assembly of Combined 'Alveolata' Identifed Bins
```bash
flye --threads 56 --meta --pacbio-raw bin-012reads.fastq.gz -o bin-012reads_assembly_raw
```

### Rough Clean & Mask
```bash
funannotate clean -i assembly.fasta -o assembly_cleaned.fasta --exhaustive --cpus 56
funannotate sort -i assembly_cleaned.fasta -o assembly_cleaned_sorted.fasta

RepeatMasker -e hmmer -pa 28 -qq -xsmall -gff --species Paramecium -dir .  assembly_cleaned_sorted.fasta
mkdir repeatmasker
mv assembly_cleaned_sorted.fasta.* repeatmasker/

ln -s repeatmasker/assembly_cleaned_sorted.fasta.masked assembly_cleaned_sorted_masked.fasta

ln -s assembly_cleaned_sorted_masked.fasta assembly.fasta
```

## Pilon
Two rounds of Pilon with Illumina Nova-Seq Libraries

## Telomeres
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

## Mitochondria
```bash
mitofinder -a assembly.fasta -p 56 -r p_caudatum_mitochondria.gb --new-genes --allow-intron --numt --intron-size 35 --max-contig-size 60000 -o 6 -j pb_mito


mitochondrial=(scaffold_350 scaffold_862 scaffold_698)
for i in "${mitochondrial[@]}"; do   echo "${i}";   sed -i "s/^>${i}$/>${i}_putative_mito/" assembly.fasta; done
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
