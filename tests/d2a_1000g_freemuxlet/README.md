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
3. `bash tests/d2a_1000g_freemuxlet/02_dsc_pileup.sh`
4. `bash tests/d2a_1000g_freemuxlet/03_run_freemuxlet.sh`
5. `bash tests/d2a_1000g_freemuxlet/04_compare_to_mis.sh`

Notes
- The scripts do not modify original inputs. The 1000G VCF is copied, bgzipped, and indexed into the test output directory.
- If the 1000G contig naming does not match the BAM (e.g., `1` vs `chr1`), stop and adjust the panel before running pileup.
