/*
Nextflow module to specifically make a transcriptome fasta using gffread

Requires: none

Params: none

Input:
  -- tuple of metadata and gtf, path to a GTF file
  -- tuple of metadata and fasta, path to a genome fasta file

Output:
  -- tuple of metadata and path to the transcriptome fasta file, emitted as `fasta`
  -- tuple of metadata and path to the gzipped transcriptome fasta file, emitted as `fasta_gzip`
  -- tuple of metadata and path to the input GTF file, emitted as `gtf`
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in module.config to add optional arguments
    publishDir in module.config

*/

process GFFREAD_TO_TRANSCRIPTOME {
    tag "${meta_gtf.id}"

    label 'process_2cpu_4gb_30min'

    input:
    tuple val(meta_gtf), path(gtf)
    tuple val(meta_fasta), path(fasta)

    output:
    // use meta_gtf for transcriptome fasta output as the gtf actually defines the transcriptome
    // and not the genome fasta
    tuple val(meta_gtf), path("$outname"), emit: fasta
    tuple val(meta_gtf), path('*transcriptome*.fa.gz'), emit: fasta_gzip
    tuple val(meta_gtf), path('*gtf', includeInputs: true), emit: gtf
    path "_version.yaml", emit: version

    script:
    def args = task.ext.args ? task.ext.args.join(' \\\n      ').trim(): ''
    def fields = meta_gtf.species.toLowerCase().tokenize("_")
    def basename =  fields[0].substring(0,1) + fields[-1].substring(0,3)
    // genome customization can be identified by the presence of a meta_gtf.id2 (or meta_gtf.id1) field
    // use this to build a more informative output filename
    outname = meta_gtf.id2 ? "${basename}.${meta_fasta.assembly}.${meta_fasta.id2}.transcriptome.${meta_fasta.assembly_type}.fa" :
      "${basename}.${meta_fasta.assembly}.transcriptome.${meta_fasta.assembly_type}.fa"
    """
    # pigz -dc $fasta > genome.fa
    gffread \\
        $args \\
        -g $fasta \\
        -w $outname \\
        $gtf
    pigz -k $outname

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        gffread: \$(gffread --version)
        pigz: \$(echo \$(pigz --version 2>&1 | sed "s/^pigz //"))
    END_VERSION
    """
}
