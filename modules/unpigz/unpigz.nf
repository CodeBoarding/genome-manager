/*
Nextflow DSL2 module to decompress files using pigz.

Input:
  -- tuple of metadata and path to gzipped file

Output:
  -- tuple of metadata and unzipped file, emitted as `decompressed`
  -- tuple of metadata and the input gzip file, emitted as `gzip`
  -- _version.yaml file containing tool version info, emitted as `version`
*/

process UNPIGZ {
    tag "$gzip"

    label 'process_8cpu_32mb_10min'

    input:
    tuple val(meta), path(gzip)

    output:
    tuple val(meta), path("$outfile"), emit: decompressed
    tuple val(meta), path("*.gz", includeInputs: true), emit: gzip
    path "_version.yaml", emit: version

    script:
    outfile = gzip.toString() - '.gz'
    def args = task.ext.args ? task.ext.args.join(' ').trim() : ''
    """
    pigz -dkc $args $gzip > $outfile

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        pigz: \$(echo \$(pigz --version 2>&1 | sed "s/^pigz //"))
    END_VERSION
    """
}
