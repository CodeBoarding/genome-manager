process CUSTOM_DUMPSOFTWAREVERSIONS {
   // multiqc sif image contains python 3.9 with pyyaml 6.0
    label 'process_1cpu_512mb_30min'

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    template 'dumpsoftwareversions.py'
}
