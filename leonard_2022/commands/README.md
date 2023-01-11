# Analyses

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
