/*
Nextflow subworkflow to prepare refflat and ribosomal intervals files for Picard
Tools (e.g., CollectRnaSeqMetrics).

Requires:
  -- GTF_RRNA module
  -- GTF_REFFLAT module
  -- gtf_rrna.py and gtf_reflat.py

Input:
  -- gtf: a gzipped GTF file
  -- star_genome: (optional) STAR index to read contig info from
  -- genome_faidx (optional*): genomic fa index to read contig info from
  * at least one of genome or faidx is required

Output:
  -- a ribosomal intervals file, emitted as `rrna_intervals`
  -- a refflat file, emitted as `refflat`
  -- a channel containing software version info in yaml format, emitted as `versions`
*/

nextflow.enable.dsl = 2

include { GTF_RRNA } from '../modules/custom/gtf_rrna/main'
include { GTF_REFFLAT } from '../modules/custom/gtf_refflat/main'

workflow MAKE_PICARD_FILES {
    take:
    gtf
    star_genome

    main:
    GTF_RRNA(gtf, star_genome)
    GTF_REFFLAT(gtf)

    emit:
    rrna_intervals = GTF_RRNA.out.rrna_intervals
    refflat = GTF_REFFLAT.out.refflat
    versions = GTF_RRNA.out.version.mix(GTF_REFFLAT.out.version)
}
