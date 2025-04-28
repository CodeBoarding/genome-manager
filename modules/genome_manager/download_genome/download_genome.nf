/*
Nextflow module to download the fasta and gtf source files from Ensembl with genome_manager.py

Requires: genome_manager.py

Params:
  -- genome_registry: path to the genome registry to register this genome with
  -- species: string, the scientific name of the species with an underscore, e.g., mus_musculus
  -- release: integer, the Ensembl release version number
  -- assembly: (optional) string, official name of the genome build as used by Ensembl, e.g., GRCh38

Inputs: None

Outputs:
  -- tuple val(meta_gtf), Path(gtf): gzipped GTF file and associated metadata, emitted as `gtf`
  -- tuple val(meta_fasta), Path(fasta): gzipped fasta file and associated metadata, emitted as `genome_fasta`
  -- path to metadata.json, emitted as `metadata`
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization: none
*/

params.assembly = null

def abbreviate_species(species) {
    fields = species.toLowerCase().tokenize("_")
    return fields[0].substring(0,1) + fields[-1].substring(0,3)
}

process DOWNLOAD_GENOME {
    tag "${params.species}:${params.release}"

    label 'process_1cpu_64mb_30min_local'

    output:
    tuple val(meta), path('*.gtf.gz'), emit: gtf
    tuple val(meta), path('*.{fa,fasta}.gz'), emit: genome_fasta
    path('metadata.json'), emit: metadata
    path "_version.yaml", emit: version

    script:
    meta = [
        id: "${params.species}:${params.release}",
        species: "${params.species}",
        species_short: abbreviate_species(params.species),
        release: "${params.release}",
        assembly: "${params.assembly}"
    ]
    def assembly_arg = params.assembly ? "--assembly-name ${params.assembly}" : ""

    """
    genome_manager.py \\
      download-genome \\
      --registry-path ${params.genome_registry} \\
      --species ${params.species} \\
      --release ${params.release} \\
      --use_cwd \\
      $assembly_arg

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        genome_manager.py: \$(echo \$(genome_manager.py --version) | sed 's/genome_manager.py //g')
        python: \$(python --version 2>&1 | sed -nre 's/^[^0-9]*([0-9]+\\.[0-9]+(\\.[0-9]+)?)/\\1/p')
    END_VERSION
    """
}