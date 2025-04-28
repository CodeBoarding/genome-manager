/*
Nextflow DSL2 module to load user-defined gene(s) from a genome registry created with genome_manager.py

Requires: genome_manager.py

Params:
  -- genome_registry: path to the genome registry top-level directory
  -- system_name: the system to lookup filepaths on

Input:
  -- gene_ids: a list of gene-ids corresponding to the entries available in the genome registry

Output:
  -- tuple of metadata and path to the fasta file
  -- tuple of metadata and path to the gtf file
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`
*/

process LOAD_USER_GENES {
    tag "$gene_ids"

    label 'process_1cpu_64mb_30min_local'

    input:
    val gene_ids
    // path genome_registry    

    output:
    tuple val(meta), path('*.fa'), emit: fasta
    tuple val(meta), path('*.gtf'), emit: gtf
    path "_version.yaml", emit: version

    script:
    def id_list = gene_ids.replaceAll('\\s', '').tokenize(',')
    meta = [id: id_list.join('.')]
    def ids = id_list.join(' ')
    """
    genome_manager.py get-genes \\
      --registry-path ${params.genome_registry} \\
      --system-name ${params.system_name} \\
      --gene-ids $ids

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        genome_manager.py: \$(echo \$(genome_manager.py --version) | sed 's/genome_manager.py //g')
        generate_gtf_entry.py: \$(generate_gtf_entry.py --version 2>&1 | sed 's/generate_gtf_entry.py //g')
        python: \$(python --version 2>&1 | sed -nre 's/^[^0-9]*([0-9]+\\.[0-9]+(\\.[0-9]+)?)/\\1/p')
    END_VERSION
    """
}
