nextflow.enable.dsl = 2

include { HARMONIZE } from './modules/harmonize_vcf_bam.nf'
include { PILEUP    } from './modules/dsc_pileup.nf'
include { DEMUXLET  } from './modules/run_demuxlet.nf'
include { QC        } from './modules/qc_demuxlet.nf'

/** Read manifest */
Channel
  .fromPath(params.manifest)
  .splitCsv(header:true)
  .map { row ->
    // tuple(sample_id, bam, barcodes, vcf, sm_list)
    tuple( row.sample_id as String,
           file(row.bam),
           file(row.barcodes),
           file(row.vcf),
           file(row.sm_list) )
  }
  .set { SAMPLES }

/** 1) Harmonize + fill-tags + polymorphic filter → clean VCF */
CLEAN_VCFS = HARMONIZE(
  SAMPLES.map{ sample_id, bam, barcodes, vcf, sm_list -> tuple(sample_id, bam, vcf) }
)

/** Join back bam/barcodes/sm_list with clean vcf */
JOINED = SAMPLES
  .map{ sample_id, bam, barcodes, vcf, sm_list -> tuple(sample_id, bam, barcodes, sm_list) }
  .join( CLEAN_VCFS )  // yields (sample_id, bam, barcodes, sm_list, clean_vcf)

/** 2) Pileup */
PILEUPS = PILEUP(
  JOINED.map{ sample_id, bam, barcodes, sm_list, clean_vcf ->
    tuple(sample_id, bam, barcodes, clean_vcf, sm_list)
  }
)

/** 3) Demuxlet */
DEMUX_BEST = DEMUXLET( PILEUPS )

/** 4) QC summary (optional but handy to report) */
QC( DEMUX_BEST )
