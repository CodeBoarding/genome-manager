//
// This file holds additional check(s) specific to the genome-manager workflow
//

class WorkflowAddGenome {

    //
    // Check and validate parameters
    //
    public static void initialize(params, log) {
        // Verify that provided genome_registry path is a valid registry
        genomeRegistryExistsError(params, log)
    }

    private static void genomeRegistryExistsError(params, log) {
        // def registry_paths = [params.genome_registry, 'genomes']
        // def genome_registry = new File(registry_paths.join(File.separator))
        def genome_registry = new File(params.genome_registry)
        if (!genome_registry.exists()) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome registry not found at '${params.genome_registry}'.\n" +
                "  Check path to genome registry.\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }
    }
}