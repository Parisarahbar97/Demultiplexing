nextflow.enable.dsl = 2

include { harmonize_vcf_bam } from './modules/harmonize_vcf_bam.nf'
include { dsc_pileup        } from './modules/dsc_pileup.nf'
include { run_demuxlet      } from './modules/run_demuxlet.nf'
include { qc_demuxlet       } from './modules/qc_demuxlet.nf'

workflow {
  // 1) Load sample sheet
  ch_samples = Channel
    .fromPath(params.samples)
    .splitCsv(header:true)
    .map { row ->
      tuple( row.sample_id as String,
             file(row.bam),
             file(row.barcodes),
             file(row.vcf),
             file(row.sm_list) )
    }

  // 2) Harmonize + fill-tags + clean
  harm = ch_samples
    .map { sid, bam, barcodes, vcf, sm -> tuple(sid, bam, vcf) }
    | harmonize_vcf_bam

  // 3) Join harmonized VCFs back to sample metadata
  joined = ch_samples
    .map { sid, bam, barcodes, vcf, sm -> tuple(sid, bam, barcodes, sm) }
    .join(harm, by: 0)
    .map { sid, bam, barcodes, sm,
           sid2, bamlike, bamlike_tbi, af, af_tbi, clean, clean_tbi ->
      assert sid == sid2
      // choose files for downstream
      def vcf_for_pileup = (params.pileup_vcf == 'af') ? af : clean
      def vcf_for_demux  = (params.demux_vcf  == 'af') ? af : clean
      tuple(sid, bam, barcodes, sm, vcf_for_pileup, vcf_for_demux)
    }

  // 4) dsc-pileup
  pile = joined
    .map { sid, bam, barcodes, sm, vcf_p, vcf_d ->
      tuple(sid, bam, barcodes, vcf_p)
    }
    | dsc_pileup

  // 5) demuxlet
  demux = pile
    .join(joined, by: 0)
    .map { sid,
           plp_plp, plp_var, plp_cel,
           sid2, bam, barcodes, sm, vcf_p, vcf_d ->
      assert sid == sid2
      tuple(sid, plp_plp, plp_var, plp_cel, vcf_d, barcodes, sm)
    }
    | run_demuxlet

  // 6) QC
  demux
    .join(ch_samples, by: 0)
    .map { sid, best, sid2, bam, barcodes, vcf, sm ->
      assert sid == sid2
      tuple(sid, best, barcodes, sm)
    }
    | qc_demuxlet
}
