/*
Nextflow DSL2 module to create an rRNA intervals file for a genome, as used by Picard Tools.

Params: none

Requires: gtf_rrna.py

Input:
  -- gtf: path to a GTF file
  -- star_index: STAR index to read contig info from

Output:
  -- the rrna intervals file, emitted as `rrna_intervals`
  -- _version.yaml containing software version info, emitted as `version`

Customization:
    publishDir in module.config

*/

process GTF_RRNA {
    tag "${meta_gtf.id}"

    label 'process_2cpu_4gb_30min'

    input:
    tuple val(meta_gtf), path(gtf)
    tuple val(meta_index), path(star_index)

    output:
    tuple val(meta_gtf), path('*.rrna'), emit: rrna_intervals
    path "_version.yaml", emit: version

    script:
    // def outname = meta_gtf.derived_basename + '.rrna'
    def outname = "${meta_gtf.species_short}.${meta_gtf.assembly}.ens${meta_gtf.release}.rrna"
    """
    gtf_rrna.py \\
    --input $gtf \\
    --genome $star_index \\
    --output $outname

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        gtf_rrna.py: \$(echo \$(gtf_rrna.py --version) | sed 's/gtf_rrna.py //g')
        python: \$(python --version 2>&1 | sed -nre 's/^[^0-9]*([0-9]+\\.[0-9]+(\\.[0-9]+)?)/\\1/p')
    END_VERSION
    """
}
