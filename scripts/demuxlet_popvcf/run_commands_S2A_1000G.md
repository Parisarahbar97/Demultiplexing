Commands used for the S2A 1000G population-VCF test run
========================================================

All commands run from `/home/pr422/RDS/live/Users/Parisa/Demultiplexing`.

1. dsc-pileup using the reheadered population VCF (OMP threads = 40)

```
OMP_NUM_THREADS=40 bash /home/pr422/RDS/live/Users/Parisa/Demultiplexing/scripts/demuxlet_popvcf/01_run_pileup_1000G.sh \
  --sam /home/pr422/RDS/live/Users/Parisa/alex_output/epilep_cellranger_outputs/S2A_mapped/outs/possorted_genome_bam.bam \
  --barcodes /home/pr422/RDS/live/Users/Parisa/demux_manual/S2A/run1/S2A.barcodes.txt \
  --vcf /home/pr422/RDS/live/Users/Parisa/demux_manual/S2A/1000G_test/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.bamorder.vcf.gz \
  --out /home/pr422/RDS/live/Users/Parisa/demux_manual/S2A/1000G_test/pileup_1000G
```

2. demuxlet genotype-error sweep on the new pileup (OMP threads = 40)

```
OMP_NUM_THREADS=40 bash /home/pr422/RDS/live/Users/Parisa/Demultiplexing/scripts/demuxlet_popvcf/02_run_demuxlet_sweep.sh \
  --plp /home/pr422/RDS/live/Users/Parisa/demux_manual/S2A/1000G_test/pileup_1000G \
  --vcf /home/pr422/RDS/live/Users/Parisa/demux_manual/S2A/run1/S2A.bamorder.vcf.gz \
  --barcodes /home/pr422/RDS/live/Users/Parisa/demux_manual/S2A/run1/S2A.barcodes.txt \
  --outdir /home/pr422/RDS/live/Users/Parisa/demux_manual/S2A/1000G_test/demuxlet_GT_1000G
```
