/*
Nextflow DSL2 module to convert a GTF file to refflat format, as used by Picard Tools.

Params: none

Requires: gtf_refflat.py

Input:
  -- gtf: path to a GTF file

Output:
  -- the refflat file, emitted as `refflat`
  -- _version.yaml containing the software version info, emitted as `version`

Customization:
    publishDir in module.config

*/

process GTF_REFFLAT {
    tag "${meta.id}"

    label 'process_2cpu_4gb_30min'

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path('*.refflat'), emit: refflat
    path "_version.yaml", emit: version

    script:
    // def outname = meta.derived_basename + '.refflat'
    def outname = "${meta.species_short}.${meta.assembly}.ens${meta.release}.refflat"
    """
    gtf_refflat.py --input $gtf --output $outname

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        gtf_refflat.py: \$(echo \$(gtf_refflat.py --version) | sed 's/gtf_refflat.py //g')
        python: \$(python --version 2>&1 | sed -nre 's/^[^0-9]*([0-9]+\\.[0-9]+(\\.[0-9]+)?)/\\1/p')
    END_VERSION
    """
}
