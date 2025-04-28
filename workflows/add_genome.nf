nextflow.enable.dsl=2

include { UNPIGZ as UNZIP_FASTA; UNPIGZ as UNZIP_GTF } from '../modules/unpigz/unpigz'
include { PREP_GENOME } from '../subworkflows/prep_genome.nf'
include { DOWNLOAD_GENOME } from '../modules/genome_manager/download_genome/download_genome'
include { REGISTER_GENOME } from '../modules/genome_manager/register_genome/register_genome'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

workflow ADD_GENOME {
    // hard-code ch_junctions as NO_FILE for now; only enable optional junctions file if genome_manager.py is
    // updated to track the junctions file and ensure that genome is readily discernable from the base unmodified
    // version
    ch_junctions = file('NO_FILE')
    ch_versions = Channel.empty()

    DOWNLOAD_GENOME()
    ch_versions = ch_versions.mix(DOWNLOAD_GENOME.out.version)

    UNZIP_FASTA(DOWNLOAD_GENOME.out.genome_fasta)
    ch_versions = ch_versions.mix(UNZIP_FASTA.out.version)

    UNZIP_GTF(DOWNLOAD_GENOME.out.gtf)
    ch_versions = ch_versions.mix(UNZIP_GTF.out.version)

    PREP_GENOME(
        UNZIP_FASTA.out.decompressed,
        UNZIP_GTF.out.decompressed,
        DOWNLOAD_GENOME.out.metadata,
        ch_junctions)
    ch_versions = ch_versions.mix(PREP_GENOME.out.versions)

    REGISTER_GENOME(
        DOWNLOAD_GENOME.out.gtf,
        DOWNLOAD_GENOME.out.genome_fasta,
        PREP_GENOME.out.transcriptome_fasta_gzip,
        PREP_GENOME.out.star_index,
        PREP_GENOME.out.refflat,
        PREP_GENOME.out.rrna_intervals,
        DOWNLOAD_GENOME.out.metadata
    )
    ch_versions = ch_versions.mix(REGISTER_GENOME.out.version)

    CUSTOM_DUMPSOFTWAREVERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
}

// function to convert full species string in <genus>_<species> format to an abbreviated string
// def abbreviate_species(species) {
//     fields = species.toLowerCase().tokenize("_")
//     return fields[0].substring(0,1) + fields[-1].substring(0,3)
// }

// // function to reformat assembly names for cleaner filenames
// def format_assembly(assembly) {
//     if (assembly.toLowerCase().startsWith("macaca")) {
//         fields = assembly.toLowerCase().replaceAll("\\.", "").tokenize("_")
//         formatted = fields[0].substring(0,1) + fields[1].substring(0,3) + fields[2]
//     }
//     else {
//         formatted = assembly.replaceAll("\\.", "")
//     }
//     return formatted
// }

// assembly_type_map = [
//     'toplevel': 'tl',
//     'primary_assembly': 'pa'
// ]

// def (_,release,filename) = (params.genome_fasta =~ /^.+(?:release|ensembl)-(\d+).*\/(.+)$/)[0]
// def fasta_fields = (filename =~ /^(?<species>[^.]+)\.(?<assembly>[^.]+(?:\.\d+)?)\.(?<type>[^.]+)\.(?<assemblytype>[^.]+)\.(?:fa|fasta)\.gz$/)
// assert fasta_fields.matches()
// def fasta_assembly = format_assembly(fasta_fields.group('assembly'))

// def fasta_meta = [
//     id: fasta_assembly + ':' + release.padLeft(3, "0"),
//     species: fasta_fields.group('species').toLowerCase(),
//     species_short: abbreviate_species(fasta_fields.group('species')),
//     assembly: fasta_assembly,
//     release: release.padLeft(3, "0"),
//     type: fasta_fields.group('type'),
//     assembly_type: assembly_type_map.containsKey(fasta_fields.group('assemblytype')) ?
//         assembly_type_map[fasta_fields.group('assemblytype')] : fasta_fields.group('assemblytype'),
//     ]

// File gtf = new File(params.gtf)
// def gtf_fields = (gtf.getName() =~ /^(?<species>[^.]+)\.(?<assembly>[^.]+(?:\.\d+)?)\.(?<release>\d+)\.gtf\.gz$/)
// assert gtf_fields.matches()
// def gtf_species_short = abbreviate_species(gtf_fields.group('species'))
// def gtf_assembly = format_assembly(gtf_fields.group('assembly'))
// def gtf_version = gtf_fields.group('release').padLeft(3, "0")
// def gtf_meta = [
//     id: gtf_assembly + ':' + gtf_fields.group('release'),
//     species: gtf_fields.group('species').toLowerCase(),
//     species_short: gtf_species_short,
//     assembly: gtf_assembly,
//     release: gtf_version
// ]

// workflow {
//     ch_genome_fasta = Channel.fromList( [
//             [ fasta_meta, file(params.genome_fasta) ]
//             ] ).collect()

//     ch_gtf = Channel.fromList( [
//             [ gtf_meta, file(params.gtf) ]
//             ] ).collect()

//     ch_junctions = file('NO_FILE')

//     PREP_GENOME_FA(ch_genome_fasta)
//     PREP_GTF(ch_gtf)
//     PREP_GENOME(PREP_GENOME_FA.out.decompressed, PREP_GTF.out.decompressed, ch_junctions)
//     REGISTER_GENOME(
//         ch_gtf,
//         ch_genome_fasta,
//         PREP_GENOME.out.transcriptome_fasta_gzip,
//         PREP_GENOME.out.star_index,
//         PREP_GENOME.out.refflat,
//         PREP_GENOME.out.rrna_intervals
//     )
// }
