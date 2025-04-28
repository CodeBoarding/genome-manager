/*
Nextflow module to run gtf2gff3.pl script from MAJIQ
*/

process GTF2GFF3 {
    tag "${meta.id}"

    label 'process_8cpu_16gb_30min'

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path('*.gff3'), emit: gff3
    path "_version.yaml", emit: version

    script:
    def config = task.ext.config == null ? '--cfg /usr/bin/gtf2gff3.cfg' : task.ext.config
    def outname = gtf.name.take(gtf.name.lastIndexOf('.')) + '.gff3'
    """
    gtf2gff3.pl $config $gtf > $outname
 
    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        gtf2gff3.pl: \$(gtf2gff3.pl --help 2>&1 | grep version | sed -nre 's/^[^0-9]*([0-9]+\\.[0-9]+)/\\1/p')
        perl: \$(perl --version 2>&1 | grep -Eo 'v[0-9]+(\\.[0-9]+)*' | sed -nre 's/v//p')
    END_VERSION
   """
}
