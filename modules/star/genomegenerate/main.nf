/*
Nextflow module to run STAR genomeGenerate on an input fasta file.

Requires: none

Params:
  -- readlength (default = 101), the length of the reads to be aligned
       [see Additional Details below]

Input:
  -- (meta_fasta, fasta): tuple of metadata and uncompressed fasta file of the genome to index
     Note that while STAR does allow multiple fasta files as input, this module
     requires a single fasta file due to the file size calculation it performs.
  -- (meta_gtf, gtf): tuple of gtf metadata and file for the genome
  -- junctions (optional) path to a junctions files to include in genome index, if none use
       special value "NO_FILE" and pass as a file

Output:
  -- tuple of meta_fasta and path to the indexed genome files, emitted as `star_index`
  -- path to _version.yaml that contains the version of the tool used in the run, emitted as `version`

Customization:
    ext.args in modules.config to add optional arguments
    publishDir in modules.config

Additional Details:
--sjdbOverhang is calculated from params.readlength as readlength - 1 as
  recommended in the STAR manual (default: 101 leading to --sjdbOverhang 100 which
  is the STAR default recommended to work in most situations)
--genomeSAindexNbases is estimated from the size of the input fasta file using the formula
  recommended in the STAR manual (this should be a good enough approximation of the genome size)

*/

process STAR_GENOMEGENERATE {
    container 'quay.io/biocontainers/star:2.7.1a--0'

    tag "${meta_fasta.id}"
    
    label 'process_20cpu_48gb_1h'

    input:
    tuple val(meta_fasta), path(fasta)
    tuple val(meta_gtf), path(gtf)
    path junctions

    output:
    tuple val(meta_fasta), path("star-index*"), emit: index
    path "_version.yaml", emit: version

    script:
    if (meta_fasta.species != meta_gtf.species) {
      throw new Exception('fasta and GTF species do not match')
    }
    if (meta_fasta.version != meta_gtf.version) {
      log.warn('fasta and GTF versions do not match')
    }
    def overhang = params.readlength ? params.readlength - 1 : 100 // if not provided, use the default recommended by STAR documentation
    def junctions_arg = junctions.name != 'NO_FILE' ? "--sjdbFileChrStartEnd $junctions" : ''
    def genome_length_est = fasta.size()
    def compare = Math.floor( Math.log(genome_length_est) / Math.log(2) * 0.5 - 1 )
    def saindex = (int)Math.min(14, Math.floor(compare))
    def args = task.ext.args ? task.ext.args.join(' \\\n      ').trim() : ''
    """
    VERSIONGENOME=\$(STAR --help | grep versionGenome | awk '{print \$2}' | sed -e 's/\\.//g')
    mkdir "star-index-v\${VERSIONGENOME}"
    STAR \\
      --runThreadN ${task.cpus} \\
      --runMode genomeGenerate \\
      --genomeDir "star-index-v\${VERSIONGENOME}/" \\
      --genomeFastaFiles $fasta \\
      --sjdbGTFfile $gtf \\
      $junctions_arg  \\
      --sjdbOverhang ${overhang} \\
      --genomeSAindexNbases ${saindex} \\
      $args

    cat <<-END_VERSION > _version.yaml
    "${task.process}":
        star: \$(STAR --version)
    END_VERSION
    """
}
