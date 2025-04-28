/*
Nextflow DSL2 module to produce a tx2gene file from a GTF file.

Params: none

Requires: gtf_tx2gene.py

Input:
  -- gtf: path to a decompressed GTF file

Output:
  -- the tx2gene file, emitted as `tx2gene`
  -- _version.yaml containing software version info, emit: `version`

Customization:
    publishDir in module.config

*/

process GTF_TX2GENE {
    tag "${meta.id}"

    label 'process_2cpu_4gb_30min'

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path('*.tx2gene.csv'), emit: tx2gene
    path "_version.yaml", emit: version

    script:
    def outname = "${meta.assembly}.${meta.version}.gtf.tx2gene.csv"
    """
    gtf_tx2gene.py --input $gtf --output $outname

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        gtf_tx2gene.py: \$(echo \$(gtf_tx2gene.py --version) | sed 's/gtf_tx2gene.py //g')
        python: \$(python --version 2>&1 | sed -nre 's/^[^0-9]*([0-9]+\\.[0-9]+(\\.[0-9]+)?)/\\1/p')
    END_VERSION
    """
}
