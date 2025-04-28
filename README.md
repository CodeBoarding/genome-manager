# Genome Manager
### A Python tool for managing data files representing genomes
Base genomes and their annotations, generally represented as fasta files and associated GTF/GFF annotation files, are fundamental inputs for NGS analysis. Various analysis tools require additional files with specific formats that are derived from the base genome, for example STAR genome indexing or PicardTools modules. Therefore, it is convenient to pre-process these files, store them, and serve them to applications that need them. This tool is designed to register these files along with custom user-defined genes and to make them available for pipelines. It can be particularly helpful in serving Nextflow pipelines.

## Tracked File Types
  - genome fasta file
  - GTF file
  - transcriptome fasta file
  - STAR index (directory containing output of `STAR genomegenerate`)
  - refflat format file (required by [Picard CollectRnaSeqMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360036884731-CollectRnaSeqMetrics-Picard-))
  - rRNA interval list file (required by [Picard CollectRnaSeqMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360036884731-CollectRnaSeqMetrics-Picard-))

New file types can be requested, although genome registries will likely need to be completely rebuilt upon adding them.

## Basic Usage
```
usage: genome_manager.py [-h] [--version]
                         {init,register-genome,register-gene,update-gene,get-genes,add-mountpoint,remove-mountpoint,list-mountpoints,list-genomes,list-genes,download-genome,clean}
                         ...

positional arguments:
  {init,register-genome,register-gene,update-gene,get-genes,add-mountpoint,remove-mountpoint,list-mountpoints,list-genomes,list-genes,download-genome,clean}
    init                initialize a new genome registry
    register-genome     register a genome
    register-gene       register a user-defined gene
    update-gene         update an existing user-defined gene with a new YAML gene model
    get-genes           retrieve fasta and YAML genes models for selected user-defined genes
    add-mountpoint      add absolute paths for a new mount point when a single physical genome registry is shared  across a network
    remove-mountpoint   remove a previously registered mountpoint by system-name
    list-mountpoints    list all registered mount points in a two column format of <system-name> <mount point>
    list-genomes        list all registered genomes arranged by species
    list-genes          list all registered user-defined genes by ID
    download-genome     download the source genome files from Ensembl FTP server
    clean               delete any temporary files/directories to recovery space

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
```

## Initializing a new genome registry
```
usage: genome_manager.py init [-h] --registry-path REGISTRY_PATH --system-name SYSTEM_NAME [--group-name GROUP_NAME]

options:
  -h, --help            show this help message and exit
  --registry-path REGISTRY_PATH
                        path to the genome registry
  --system-name SYSTEM_NAME
                        string identifying the system where paths will exist (e.g., HPC name, workstation name, etc.)
  --group-name GROUP_NAME
                        optional permission group name for write access to non-restricted subdirectories (e.g., user_defined_genes)
```

**IMPORTANT**: Setting the optional `group-name` argument is highly recommended as it will restrict write access to only the users beloging to that group. If this argument is not provided, then the registry paths will belong to the creating user's default group. Any user in that group will be able to perform potentially very destructive actions on the genome registry. Therefore, it is best practice to set `group-name` to a permission group comprised only of those users who should have admin access to the registry.

## Adding a new genome
Adding a genome requires download of large source files and generation of additional files derived from them. Therefore, the actual process of doing this has been implemented using Nextflow ([see below](#nextflow)). For completeness of the documentation, here is the actual genome_manager.py command for registering a new genome, although it is strongly recommended to not do this directly:

```
usage: genome_manager.py register-genome [-h] --registry-path REGISTRY_PATH --system-name SYSTEM_NAME --genome-metadata-file
                                         GENOME_METADATA_FILE --input-dir INPUT_DIR

options:
  -h, --help            show this help message and exit
  --registry-path REGISTRY_PATH
                        path to the genome registry
  --system-name SYSTEM_NAME
                        string identifying the system where paths will exist (e.g., HPC name, workstation name, etc.)
  --genome-metadata-file GENOME_METADATA_FILE
                        path to JSON file containing genome metadata
  --input-dir INPUT_DIR
                        path to input directory containing all required genome files
```

Currently, the project is designed around the genomes and annotations hosted by Ensembl. The process of adding a new genome to a registry is automated via a Nextflow script that requires a `species` (e.g., `homo_sapiens`) and an Ensembl `release` version number (e.g., `110`), along with the path to the genome registry and the name of the system you are adding to. If the genome to be added is a release that uses the most recent assembly for **human**, **mouse**, **rat**, **cynomolgus**, **rhesus macaque**, **pig**, or **Chinese hamster ovary (CHO) cells**, then the `assembly` does not need to be provided. However, when adding an older release that used a previous assembly (e.g., GRCm38, used by Ensembl release versions <= 102) or an organism other than those listed, the `assembly` needs to also be provided in the call to Nextflow. The map of releases to corresponding assemblies for all organsisms can be found [here](https://www.ensembl.org/info/website/archives/assembly.html). Note that if there are multiple mount points registered for the genome registry (see ["Adding a new mount point"](#adding-a-new-mount-point) below), then paths will be built for all of those systems as well.

Example adding a recent version of a common genome:
```
nextflow run pfizer-opensource/genome-manager --command register-genome --genome-registry /path/to/my/registry --system-name my-system-name --species homo_sapiens --release 110
```
Note: You should also use the `-profile` flag to select an appropriate profile, either from `nextflow.config` or from a custom config file you provide via the `-c` flag. This profile should enable use of Docker or Singularity at the very least.

Example adding a more esoteric genome:
```
nextflow run pfizer-opensource/genome-manager --command register-genome --genome-registry /path/to/my/registry --system-name my-system-name --species cairina_moschata_domestica --release 110 --assembly 'CaiMos1.0'
```

## Adding a custom user-defined gene
The genome_manager.py script also has a command for adding a custom gene defined by the user. This is useful for when one needs to supplement a base organism genome with custom sequence(s) to align to, such as a gene therapy vector. This requires a fasta file containing the sequence of the gene, which can also include non-expressed surrounding sequences like flanking regions, and a yaml file that describes the transcript(s) produced from the sequence in the fasta file. See below for an example of the expected yaml file format.

```
usage: genome_manager.py register-gene [-h] --registry-path REGISTRY_PATH --yaml-file
                                       YAML_FILE --fasta FASTA --system-name SYSTEM_NAME

options:
  -h, --help            show this help message and exit
  --registry-path REGISTRY_PATH
                        path to the genome registry
  --yaml-file YAML_FILE
                        path to YAML file containing gene model for custom fasta
  --fasta FASTA         path to fasta file containing sequence of a genome modification
  --system-name SYSTEM_NAME
                        string identifying the system where paths will exist (e.g., HPC
                        name, workstation name, etc.)
```

### YAML format for a user-defined gene
```
# My-Custom-Gene
- gene_id: "my-gene-id"
  gene_name: "my-gene-name"
  strand: '+'
  gene_biotype: "exogenous_protein_coding"
  transcripts:
    - name: "my-gene-01"
      exons:
        - start: 600
          end: 2663
      cds_start: 848
      cds_end: 2433
    - name: "my-gene-02"
      exons:
        - start:812
          end: 2354
      cds_start: 908
      cds_end: 2130
```

If it becomes necessary to update the annotation for any reason, this can be done using the `update-gene` command along with an updated YAML file containing the new information.
```
usage: genome_manager.py update-gene [-h] --registry-path REGISTRY_PATH --yaml-file
                                     YAML_FILE --system-name SYSTEM_NAME

options:
  -h, --help            show this help message and exit
  --registry-path REGISTRY_PATH
                        path to the genome registry
  --yaml-file YAML_FILE
                        path to YAML file containing gene model for custom fasta
  --system-name SYSTEM_NAME
                        string identifying the system where paths will exist (e.g., HPC
                        name, workstation name, etc.)
```
The gene version number will be auto-incremented meaning that the previous version will still be stored and can be called by its version number for reproducibility. By default, when accessing user-defined genes, genome_manager.py will return the most recent version.

## Working with user-defined genes
To get a user-defined gene back out of the registry, use the `get-genes` command:
```
usage: genome_manager.py get-genes [-h] --registry-path REGISTRY_PATH --gene-ids
                                   GENE_IDS [GENE_IDS ...] --system-name SYSTEM_NAME
                                   [--outdir OUTDIR] [--version-delim VERSION_DELIM]

options:
  -h, --help            show this help message and exit
  --registry-path REGISTRY_PATH
                        path to the genome registry
  --gene-ids GENE_IDS [GENE_IDS ...]
                        comma-separated list of gene IDs to retrieve from registry
  --system-name SYSTEM_NAME
                        string identifying the system where paths will exist (e.g., HPC
                        name, workstation name, etc.)
  --outdir OUTDIR       target directory to write output to
  --version-delim VERSION_DELIM
                        delimiter when specifying a specific gene model version number
                        after a gene-id (default = '.')
```

This will return the fasta file along with a GTF file containing all the necessary entries to describe the gene and any transcripts defined in the YAML file. In the context of a Nextflow pipeline, these can then be added to channels, combined with a base organism genome, and routed to the appropriate modules to rebuild the necessary genome files for your pipeline (e.g., STAR genomegenerate or making a refflat format file for PicardTools, etc.).

## Adding a new mount point
The genome registry is designed to support multiple mount points for the same physical volume (e.g., NFS mounts). If the fileshare is mounted to a different path on a system, use the `add-mountpoint` command to register that mount point along with a new `system-name` that can be used to reference it.

```
usage: genome_manager.py add-mountpoint [-h] --registry-path REGISTRY_PATH --mountpoint
                                        MOUNTPOINT --system-name SYSTEM_NAME

options:
  -h, --help            show this help message and exit
  --registry-path REGISTRY_PATH
                        path to the genome registry
  --mountpoint MOUNTPOINT
                        path to the new mount point for the registry (i.e., the top-
                        level genome registry directory)
  --system-name SYSTEM_NAME
                        string identifying the system where paths will exist (e.g., HPC
                        name, workstation name, etc.)
```

## Removing a mount point
Mountpoints can also be removed using the `remove-mountpoint` command. Note, however, that the default system, i.e., the one that was created with the `init` command, cannot be removed.

```
usage: genome_manager.py remove-mountpoint [-h] --registry-path REGISTRY_PATH --remove-system-name REMOVE_SYSTEM_NAME

options:
  -h, --help            show this help message and exit
  --registry-path REGISTRY_PATH
                        path to the genome registry
  --remove-system-name REMOVE_SYSTEM_NAME
                        the system-name string of the mount point to be removed
```

## List functions

### List all currently registered genomes
```
usage: genome_manager.py list-genomes [-h] --registry-path REGISTRY_PATH [--species SPECIES]

options:
  -h, --help            show this help message and exit
  --registry-path REGISTRY_PATH
                        path to the genome registry
  --species SPECIES     optional species to search for
```
If `species` is not provided, this will output a list of all registered genomes to stdout. Otherwise, it will output only those genomes matching `species`. It is arranged by species and the first column provides an ID that can be used in Nextflow pipelines via the load_genome workflow.

### List all currently registerd user-defined genes
```
usage: genome_manager.py list-genes [-h] --registry-path REGISTRY_PATH

options:
  -h, --help            show this help message and exit
  --registry-path REGISTRY_PATH
                        path to the genome registry
```

This command will print an alphanumerically sorted listing of IDs for all user-defined genes that have been added to the genome registry.

### List all current mount points
```
usage: genome_manager.py list-mountpoints [-h] --registry-path REGISTRY_PATH

options:
  -h, --help            show this help message and exit
  --registry-path REGISTRY_PATH
                        path to the genome registry
```
This will output a list of all registered mount points to stdout. The first column is the `system_name` to use for selecting the path in column two.

## Clean
Some commands will store temporary files or backup copies of files into a .tmp directory inside the registry. Generally, these are automatically cleaned up. However, in the event they ever are not, the `clean` command can be used to remove them.
```
usage: genome_manager.py clean [-h] --registry-path REGISTRY_PATH

options:
  -h, --help            show this help message and exit
  --registry-path REGISTRY_PATH
                        path to the genome registry
```

# Docker
The scripts for genome-manager are located in the /genome_manager directory, along with a Dockerfile and `env.yaml`. A Docker image can be built directly from this directory using the Dockerfile. Note, however, that the base image is a custom micromamba Docker image that must also be built. The Dockerfile to build this can be found at [TODO: add link to pfizer-opensource repo]

# Nextflow
The main commands for adding new genomes and user-defined genes, as well as updating user-defined genes, are also provided as a Nextflow pipeline for increased reproducibility.

Usage:
```
nextflow run pfizer-opensource/genome-manager --command [register-genome, register-gene, update-gene]
```

Each command has additional parameters:
 
`register-genome`: genome_registry, system_name, species, release, assembly (optional)  
`register-gene`: genome_registry, system_name, gene_fasta, gene_model  
`udpate-gene`: genome_registry, system_name, gene_model  

Additionally, a profile should be chosen from `profiles.config` file to provide the necessary settings for Singularity/Docker, job manager, and other parameters. This file has settings for the major Pfizer HPC environments. A custom profiles.config file will need to be substituted for use on other systems.

## Incorporating genome-manager in a Nextflow pipeline
This repo contains everything necessary to initialize, populate, and maintain a genome registry. It can then be accessed from any Nextflow pipeline by adding the appropriate modules and subworkflows to it. The `load_genome_extra_files` directory contains a collection of additional modules that need to be added to your pipeline. The key file is `load_genome_extra_files/subworkflows/load_genome.nf` as this is what will ultimately be called directly from within the Nextflow pipeline to import a specified genome +/- additional genes.

The required files and corresponding directory structure (relative to the top-level directory of your pipeline) are listed below. These relative paths are also where the files are found in this repository, with the exception of those found in `/load_genome_extra_files` which is indicated parenthetically where appropriate:

/modules/unpigz/unpigz.nf  
/modules/custom/gtf_tx2gene/gtf_tx2gene.nf (/load_genome_extra_files)  
/modules/custom/gtf_metadata/gtf_metadata.nf (/load_genome_extra_files)  
/modules/custom/concat/concat.nf (/load_genome_extra_files)  
/modules/genome_manager/load_user_genes/load_user_genes.nf  
/modules/custom/gtf2gff3/gtf2gff3.nf (/load_genome_extra_files)  
/modules/star/genomegenerate/main.nf  
/modules/gffread/to_transcriptome/main.nf  
/subworkflows/make_picard_files.nf  
/subworkflows/prep_genome.nf  
/subworkflows/load_genome.nf (/load_genome_extra_files)

## Author
Rob Moccia
