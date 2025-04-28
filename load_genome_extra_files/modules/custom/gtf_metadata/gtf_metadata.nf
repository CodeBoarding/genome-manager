/*
Nextflow DSL2 module to produce a gene or transcript metadata file from a GTF file.

Params: none

Requires: gtf_metadata.py

Input:
  -- gtf: path to a decompressed GTF file

Output:
  -- the metadata file, emitted as `metadata`
  -- _version.yaml containing software version info, emitted as `version`

Customization:
    task.ext.args in module.config
    publishDir in module.config

*/

process GTF_METADATA {
    tag "${meta.id}"

    label 'process_2cpu_4gb_30min'

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path('*.metadata.csv'), emit: metadata
    path "_version.yaml", emit: version

    script:
    def outname = "${meta.assembly}.${meta.version}.gtf.metadata.csv"
    def args = task.ext.args ? task.ext.args.join(' \\\n      ').trim : ''
    """
    gtf_metadata.py \\
      --input $gtf \\
      --output $outname \\
      $args

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        gtf_metadata.py: \$(echo \$(gtf_metadata.py --version) | sed 's/gtf_metadata.py //g')
        python: \$(python --version 2>&1 | sed -nre 's/^[^0-9]*([0-9]+\\.[0-9]+(\\.[0-9]+)?)/\\1/p')
    END_VERSION
    """
}
