nextflow.enable.dsl=2

include { REGISTER_GENE } from '../modules/genome_manager/register_gene/register_gene'

workflow REGISTER_USER_DEFINED_GENE {
    ch_gene_fasta = Channel.fromPath(params.gene_fasta)
    ch_gene_model = Channel.fromPath(params.gene_model)
    REGISTER_GENE(ch_gene_fasta, ch_gene_model)
}
