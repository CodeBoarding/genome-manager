#!/usr/bin/env nextflow

nextflow.enable.dsl=2

if (params.command == 'register-genome') {
    include { ADD_GENOME } from './workflows/add_genome'
}
else if (params.command == 'register-gene') {
    include { REGISTER_USER_DEFINED_GENE } from './workflows/register_user_defined_gene'
}
else if (params.command == 'update-gene') {
    include { UPDATE_USER_DEFINED_GENE } from './workflows/update_user_defined_gene'
}

workflow {
    WorkflowMain.initialize(workflow, params, log)
    WorkflowAddGenome.initialize(params, log)
    
    if (params.command == 'register-genome') {
        ADD_GENOME()
    }
    else if (params.command == 'register-gene') {
        REGISTER_USER_DEFINED_GENE()
    }
    else if (params.command == 'update-gene') {
        UPDATE_USER_DEFINED_GENE()
    }
}

// workflow {
//     GENOME_MANAGER()
// }