/*
Nextflow module to register a user-defined gene with genome_manager.py

Requires: genome_manager.py

Params:
  -- genome_registry: path to the genome registry to register this genome with
  -- system_name: the system where genome_registry will be accessible (e.g., hpc, AWS, etc.)
  -- gene_fasta: fasta file (uncompressed) containing the full sequence of the gene/chromosome/plasmid/etc.
  -- gene_model: yaml file describing the gene/transcripts 

Inputs: None

Outputs:
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization: none
*/

process REGISTER_GENE {
    tag "${new File(params.gene_fasta).getName()}"

    label 'process_1cpu_64mb_30min_local'

    input:
    path gene_fasta
    path gene_model

    output:
    path "_version.yaml", emit: version

    script:
    """
    genome_manager.py \\
      register-gene \\
      --registry-path ${params.genome_registry} \\
      --system-name ${params.system_name} \\
      --fasta $gene_fasta \\
      --yaml-file $gene_model

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        genome_manager.py: \$(echo \$(genome_manager.py --version) | sed 's/genome_manager.py //g')
        python: \$(python --version 2>&1 | sed -nre 's/^[^0-9]*([0-9]+\\.[0-9]+(\\.[0-9]+)?)/\\1/p')
    END_VERSION
    """
}
