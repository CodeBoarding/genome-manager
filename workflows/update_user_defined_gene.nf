nextflow.enable.dsl=2

include { UPDATE_GENE } from '../modules/genome_manager/update_gene/update_gene'

workflow UPDATE_USER_DEFINED_GENE {
    ch_gene_model = Channel.fromPath(params.gene_model)
    UPDATE_GENE(ch_gene_model)
}
