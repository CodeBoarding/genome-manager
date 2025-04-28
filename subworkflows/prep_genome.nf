/*
Subworkflow to prepare all genome-related files needed by the pipeline.
*/

import groovy.json.JsonSlurper

include { STAR_GENOMEGENERATE } from '../modules/star/genomegenerate/main'
include { GFFREAD_TO_TRANSCRIPTOME } from '../modules/gffread/to_transcriptome/main'
include { MAKE_PICARD_FILES } from './make_picard_files.nf'


workflow PREP_GENOME {
    take:
    genome_fasta
    gtf
    genome_metadata
    junctions

    main:

    ch_versions = Channel.empty()

    // extract genome_metadata provided by download_genome and attach it to both fasta and gtf
    // genome_metadata file is in a channel so need to map to JsonSlurper and collect to read it
    def metadata = genome_metadata.map{
      it -> new JsonSlurper().parse(it)
      }.collect{it}

    // metadata is now a DataflowVariable so use .get() to extract the value and index to first position
    ch_genome_fasta = genome_fasta.map{
        it -> [metadata.get()[0], it[1]]
        }

    ch_gtf = gtf.map{
        it -> [metadata.get()[0], it[1]]
        }

    STAR_GENOMEGENERATE(ch_genome_fasta, ch_gtf, junctions)
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.version)

    GFFREAD_TO_TRANSCRIPTOME(ch_gtf, ch_genome_fasta)
    ch_versions = ch_versions.mix(GFFREAD_TO_TRANSCRIPTOME.out.version)

    MAKE_PICARD_FILES(ch_gtf, STAR_GENOMEGENERATE.out.index)
    ch_versions = ch_versions.mix(MAKE_PICARD_FILES.out.versions.first())

    emit:
    star_index = STAR_GENOMEGENERATE.out.index
    transcriptome_fasta_gzip = GFFREAD_TO_TRANSCRIPTOME.out.fasta_gzip
    transcriptome_fasta = GFFREAD_TO_TRANSCRIPTOME.out.fasta
    rrna_intervals = MAKE_PICARD_FILES.out.rrna_intervals
    refflat = MAKE_PICARD_FILES.out.refflat
    versions = ch_versions
}
