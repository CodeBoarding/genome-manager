/*
Subworkflow to load a genome. Works in tandem with genome registries created by genome_manager.py.
*/

import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import org.yaml.snakeyaml.Yaml

nextflow.enable.dsl=2

include { UNPIGZ as UNPIGZ_GENOME_FA; UNPIGZ as UNPIGZ_TRANSCRIPTOME_FA; UNPIGZ as UNPIGZ_GTF } from '../modules/unpigz/unpigz'
include { GTF_TX2GENE } from '../modules/custom/gtf_tx2gene/gtf_tx2gene'
include { GTF_METADATA } from '../modules/custom/gtf_metadata/gtf_metadata'
include { CONCAT as CONCAT_GTF; CONCAT as CONCAT_FASTA } from '../modules/custom/concat/concat'
include { LOAD_USER_GENES } from '../modules/genome_manager/load_user_genes/load_user_genes'
include { PREP_GENOME } from './prep_genome.nf'
include { GTF2GFF3 } from '../modules/custom/gtf2gff3/gtf2gff3'

workflow LOAD_GENOME {
    main:

    ch_versions = Channel.empty()

    def genome_release = params.genome_id.tokenize(':')[-1].padLeft(3, "0")
    def genome_config = new File("${params.genome_registry}/.conf/genome-registry/${genome_release}.json")
    def genome_registry = new JsonSlurper().parse(genome_config)
    def selected_genome = genome_registry.genomes[params.genome_id.toLowerCase()]

    // genome metadata is expected as a file by PREP_GENOME so pass it as an InputStream
    def metadata_json = groovy.json.JsonOutput.toJson(selected_genome.base.metadata)
    def genome_metadata = new ByteArrayInputStream(metadata_json.getBytes())
    ch_genome_metadata = Channel.from(genome_metadata)

    ch_genome_fasta_in = Channel.fromList([
        [selected_genome.base['genome_fasta'].plus(selected_genome.base.metadata),
        file(selected_genome.base['genome_fasta']['path'][params.system_name])]
        ]).collect()
    UNPIGZ_GENOME_FA(ch_genome_fasta_in)
    ch_versions = ch_versions.mix(UNPIGZ_GENOME_FA.out.version.first())

    ch_gtf_in = Channel.fromList([
        [selected_genome.base['gtf'].plus(selected_genome.base.metadata),
        file(selected_genome.base['gtf']['path'][params.system_name])]
        ]).collect()
    UNPIGZ_GTF(ch_gtf_in)
    ch_versions = ch_versions.mix(UNPIGZ_GTF.out.version.first())

    if( params.add_genes ) { 
        LOAD_USER_GENES(params.add_genes)
        ch_versions = ch_versions.mix(LOAD_USER_GENES.out.version.first())

        CONCAT_GTF(UNPIGZ_GTF.out.decompressed, LOAD_USER_GENES.out.gtf)
        // CONCAT only uses Unix cat so no version to report
        ch_gtf = CONCAT_GTF.out.merged

        CONCAT_FASTA(UNPIGZ_GENOME_FA.out.decompressed, LOAD_USER_GENES.out.fasta)
        // CONCAT only uses Unix cat so no version to report
        ch_genome_fasta = CONCAT_FASTA.out.merged

        ch_junctions = file(params.junctions)
        PREP_GENOME( ch_genome_fasta, ch_gtf, ch_genome_metadata, ch_junctions )
        ch_transcriptome_fasta = PREP_GENOME.out.transcriptome_fasta
        ch_star_index = PREP_GENOME.out.star_index
        ch_refflat = PREP_GENOME.out.refflat
        ch_rrna_intervals = PREP_GENOME.out.rrna_intervals
        ch_versions = ch_versions.mix(PREP_GENOME.out.versions.first())
    }

    else {
        ch_genome_fasta = UNPIGZ_GENOME_FA.out.decompressed
        ch_gtf = UNPIGZ_GTF.out.decompressed

        ch_transcriptome_fasta_in = Channel.fromList([
            [selected_genome['transcriptome_fasta'].plus(base_info),
            file(selected_genome['transcriptome_fasta']['path'][params.system_name])]
            ]).collect()
    
        UNPIGZ_TRANSCRIPTOME_FA(ch_transcriptome_fasta_in)
        ch_versions = ch_versions.mix(UNPIGZ_TRANSCRIPTOME_FA.out.version.first())
        ch_transcriptome_fasta = UNPIGZ_TRANSCRIPTOME_FA.out.decompressed

        ch_star_index = Channel.fromList([
            [selected_genome['star_index'].plus(base_info),
            file(selected_genome['star_index']['path'][params.system_name])]
            ]).collect()

        ch_refflat = Channel.fromList([
            [selected_genome['refflat'].plus(base_info),
            file(selected_genome['refflat']['path'][params.system_name])]
            ]).collect()

        ch_rrna_intervals = Channel.fromList([
            [selected_genome['rrna_interval_list'].plus(base_info),
            file(selected_genome['rrna_interval_list']['path'][params.system_name])]
            ]).collect()
    }

    GTF2GFF3( ch_gtf )
    ch_versions = ch_versions.mix(GTF2GFF3.out.version.first())

    GTF_TX2GENE( ch_gtf )
    ch_versions = ch_versions.mix(GTF_TX2GENE.out.version.first())
    
    GTF_METADATA( ch_gtf )
    ch_versions = ch_versions.mix(GTF_METADATA.out.version.first())

    emit:
    genome_fasta_gzip = ch_genome_fasta_in
    gtf_gzip = ch_gtf_in
    genome_fasta = ch_genome_fasta
    gtf = ch_gtf
    gff3 = GTF2GFF3.out.gff3
    transcriptome_fasta = ch_transcriptome_fasta
    star_index = ch_star_index
    refflat = ch_refflat
    rrna_intervals = ch_rrna_intervals
    versions = ch_versions
}
