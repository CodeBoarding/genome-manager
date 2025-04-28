/*
Nextflow module to register a genome with genome_manager.py

Requires: genome_manager.py

Params:
  -- genome_registry: path to the genome registry to register this genome with
  -- system_name: the system where genome_registry will be accessible (e.g., hpc, AWS, etc.)

Inputs:
  -- tuple val(meta_gtf), path(gtf): gzipped GTF file and associated metadata
  -- tuple val(meta_genome), path(genome_fasta): gzipped genome fasta file and associated metadata
  -- tuple val(meta_transcriptome), path(transcriptome_fasta): gzipped transcriptome fasta and associated metadata
  -- tuple val(meta_star), path(star_index): STAR index directory and associated metadata
  -- tuple val(meta_refflat), path(refflat): refflat file and associated metadata
  -- tuple val(meta_rrna), path(rrna_interval_list): rRNA intervals_list format file and associated metadata
  -- path(genome_metadata): JSON file containing genome metadata following the BaseGenomeMetadata schema

Outputs:
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization: none
*/

import groovy.json.JsonBuilder

process REGISTER_GENOME {
    tag "${meta_genome.id}"

    label 'process_1cpu_64mb_30min_local'

    input:
    tuple val(meta_gtf), path(gtf)
    tuple val(meta_genome), path(genome_fasta)
    tuple val(meta_transcriptome), path(transcriptome_fasta)
    tuple val(meta_star), path(star_index)
    tuple val(meta_refflat), path(refflat)
    tuple val(meta_rrna), path(rrna_interval_list)
    path(genome_metadata)
    
    output:
    path "_version.yaml", emit: version

    script:
    // this json must match the expected input format declared in genome_manager.py
    genome_params = [
        base: [
            id: "${meta_genome.id}",
            species: "${meta_genome.species}",
            release: "${meta_genome.release}",
            assembly: "${meta_genome.assembly}",
            assembly_type: "${meta_genome.assembly_type}",
            sequence_type: "${meta_genome.type}"
            ],
        genome_fasta: genome_fasta.toString(),
        gtf: gtf.toString(),
        transcriptome_fasta: transcriptome_fasta.toString(),
        refflat: refflat.toString(),
        rrna_interval_list: rrna_interval_list.toString(),
        star_index: star_index.toString()
    ]

    json_string = new JsonBuilder(genome_params).toString()

    """
    # printf '$json_string' > args.json
    genome_manager.py \\
      register-genome \\
      --registry-path ${params.genome_registry} \\
      --system-name ${params.system_name} \\
      --genome-metadata ${genome_metadata} \\
      --input-dir \$PWD

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        genome_manager.py: \$(echo \$(genome_manager.py --version) | sed 's/genome_manager.py //g')
        python: \$(python --version 2>&1 | sed -nre 's/^[^0-9]*([0-9]+\\.[0-9]+(\\.[0-9]+)?)/\\1/p')
    END_VERSION
    """
}
