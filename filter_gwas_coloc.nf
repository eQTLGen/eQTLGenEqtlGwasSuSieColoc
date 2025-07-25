#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2


def helpmessage() {

log.info"""

eQtlGwasSuSieColoc v${workflow.manifest.version}"
===========================================================
Pipeline for running fast eQTL-GWAS colocalisation analyses on the fine-mapped eQTLs.

Usage:

nextflow run main.nf 
--coloc
--ld_reference
--variant_reference
--OutputDir

Mandatory arguments:
--coloc                     Output files from coloc pipeline
--ld_reference              Parquet dataset with permutation-based LD information
--variant_reference         Parquet file with variant information

Optional arguments:
--OutputDir                 Output directory, defaults to 'results'
""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

// Default parameters
params.OutputDir = 'results'

//Show parameter values
log.info """===============================================================
eQTLGen GWAS-eQTL fine-mapped coloc pipeline v${workflow.manifest.version}"
==========================================================================="""
def summary = [:]
summary['Pipeline Version']                         = workflow.manifest.version
summary['Current user']                             = "$USER"
summary['Current home']                             = "$HOME"
summary['Current path']                             = "$PWD"
summary['Working dir']                              = workflow.workDir
summary['Script dir']                               = workflow.projectDir
summary['Config Profile']                           = workflow.profile
summary['Container Engine']                         = workflow.containerEngine
if(workflow.containerEngine) summary['Container']   = workflow.container
summary['Output directory']                         = params.OutputDir
summary['Coloc files']                              = params.coloc
summary['LD reference']                             = params.ld_reference
summary['Variant reference']                        = params.variant_reference

// import modules
include { COLOCFILTER; FilterColoc} from './modules/eQtlGwasColoc.nf'

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

// Define argument channels
coloc_file_ch = Channel.fromPath(params.coloc)

variant_reference = Channel.fromPath(params.variant_reference)
ld_ch = Channel.fromPath(params.ld_reference, type: 'dir')

workflow {
       COLOCFILTER(coloc_file_ch, variant_reference, ld_ch)
       COLOCFILTER.out.flatten().collectFile(name: 'EqtlGwasColocSusieResults_withLdRColumn.txt', keepHeader: true, sort: true, storeDir: "${params.OutputDir}")
    }

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
