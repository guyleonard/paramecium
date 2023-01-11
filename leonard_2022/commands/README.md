# Analyses

## Preliminary
BAM --> FASTQ
cat *.fastq > all_pacbio_libraries.fq

## LRBinner
```bash
LRBinner -r all_pacbio_libraries.fq -o lrb_output_all_PB_libs --ae-epochs 200 --resume -mbs 1000 -bit 0 -bs 10 -bc 10 --threads 56
```
