D2A freemuxlet with 1000G site panel (reproducible test)

Purpose
- Run freemuxlet on D2A using a general SNP site panel (1000G) instead of the MIS pool VCF and compare results to the MIS-based run.

Outputs will be written to: `/home/pr422/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test`.

Prerequisites
- Docker images: `parisa/genotype:2.6` and `parisa/demux:2.1`.
- D2A BAM: `/home/pr422/RDS/live/Users/Parisa/EPILEP/healthy/mapping/output/D2A_mapped/outs/possorted_genome_bam.bam`
- D2A CellBender CSV: `/home/pr422/RDS/live/Users/Parisa/EPILEP/healthy/qc/output_latest/D2A/D2A_cellbender_output/cellbender_out_cell_barcodes.csv`
- 1000G VCF (sites): `/home/pr422/RDS/live/Users/Parisa/alex_output/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf`

Run order (from repo root)
1. `bash tests/d2a_1000g_freemuxlet/00_prepare_panel.sh`
2. `bash tests/d2a_1000g_freemuxlet/01_make_barcodes.sh`
3. `THREADS=40 bash tests/d2a_1000g_freemuxlet/02_dsc_pileup.sh`
4. `bash tests/d2a_1000g_freemuxlet/03_run_freemuxlet.sh`
5. `bash tests/d2a_1000g_freemuxlet/04_compare_to_mis.sh`

Notes
- The scripts do not modify original inputs. The 1000G VCF is copied, bgzipped, and reheadered to match the BAM’s contig order, then indexed.
- If the 1000G contig naming does not match the BAM (e.g., `1` vs `chr1`), adjust the panel before running pileup.

Manual run commands (copy/paste)

```bash
# Prepare 1000G panel (bgzip + index)
docker run --rm -u $(id -u):$(id -g) \
  -v /home/pr422:/home/pr422 \
  --entrypoint bash parisa/genotype:2.6 -lc "
    set -euo pipefail
    bgzip -c /home/pr422/RDS/live/Users/Parisa/alex_output/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf > \
      /home/pr422/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test/1000G_sites.vcf.gz
    tabix -f -p vcf /home/pr422/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test/1000G_sites.vcf.gz
    bcftools index -f /home/pr422/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test/1000G_sites.vcf.gz
  "

# Reheader 1000G panel to BAM contig order
docker run --rm -u $(id -u):$(id -g) \
  -v /home/pr422:/home/pr422 \
  --entrypoint bash parisa/genotype:2.6 -lc "
    set -euo pipefail
    bash /home/pr422/RDS/live/Users/Parisa/Demultiplexing/tests/d2a_1000g_freemuxlet/reorder_vcf_to_bam.sh \
      /home/pr422/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test/1000G_sites.vcf.gz \
      /home/pr422/RDS/live/Users/Parisa/EPILEP/healthy/mapping/output/D2A_mapped/outs/possorted_genome_bam.bam \
      /home/pr422/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test/1000G_sites.bamorder.vcf.gz
  "

# Build barcode whitelist
awk -F'[,	]' 'NR==1 && tolower($1) ~ /barcode/ {next} {print $1}' \
  /home/pr422/RDS/live/Users/Parisa/EPILEP/healthy/qc/output_latest/D2A/D2A_cellbender_output/cellbender_out_cell_barcodes.csv \
  | LC_ALL=C sort -u > /home/pr422/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test/D2A.barcodes.txt

# dsc-pileup with 40 threads (OMP_NUM_THREADS)
docker run --rm -u $(id -u):$(id -g) \
  -v /home/pr422:/home/pr422 \
  -e OMP_NUM_THREADS=40 \
  --entrypoint /opt/conda/bin/popscle \
  parisa/demux:2.1 \
  dsc-pileup \
    --sam /home/pr422/RDS/live/Users/Parisa/EPILEP/healthy/mapping/output/D2A_mapped/outs/possorted_genome_bam.bam \
    --vcf /home/pr422/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test/1000G_sites.bamorder.vcf.gz \
    --group-list /home/pr422/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test/D2A.barcodes.txt \
    --tag-group CB \
    --tag-UMI UB \
    --min-MQ 30 \
    --min-BQ 20 \
    --cap-BQ 40 \
    --out /home/pr422/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test/D2A_1000G_pileup

# freemuxlet for K = 3, 4, 5
for K in 3 4 5; do
  docker run --rm -u $(id -u):$(id -g) \
    -v /home/pr422:/home/pr422 \
    --entrypoint /opt/conda/bin/popscle \
    parisa/demux:2.1 \
    freemuxlet \
      --plp /home/pr422/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test/D2A_1000G_pileup \
      --nsample $K \
      --group-list /home/pr422/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test/D2A.barcodes.txt \
      --out /home/pr422/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test/freemux1000G_K${K}
  echo
  gzip -cd /home/pr422/RDS/live/Users/Parisa/demux_manual/D2A_1000G_freemuxlet_test/freemux1000G_K${K}.clust1.samples.gz \
    | awk -F'\t' 'NR>1{c[$5]++} END{for(k in c) printf "%s\t%d\n",k,c[k]}' \
    | sort
  echo "====="
done

# Compare against MIS-based K=4 baseline
echo "=== MIS-based freemux (K=4) ==="
gzip -cd /home/pr422/RDS/live/Users/Parisa/demux_manual/D2A/misVCF_run/freemux_K4.clust1.samples.gz \
  | awk -F'\t' 'NR>1{c[$5]++} END{for(k in c) printf "%s\t%d\n",k,c[k]}' \
  | sort
```
