#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

process ListGwasSignals {
    input:
        path(finemapping_loci)
    output:
        path "finemapping_summary.txt"
    script:
        """
        ListGwasSignals.R \
        --loci ${finemapping_loci} \
        --min_lbf 2 \
        --output "finemapping_summary.txt"
        """
}

def helpmessage() {

log.info"""

List GWAS signals v${workflow.manifest.version}"
===========================================================

""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

// Default parameters
params.OutputDir = 'results'

gwas_loci_ch = Channel.fromPath(params.loci).collate(100, remainder = true )

workflow {
    gwas_signals = ListGwasSignals(gwas_loci_ch)
    gwas_signals.collectFile(name: 'gwas_signals.txt', keepHeader: true, skip: 1, sort: true, storeDir: "${params.OutputDir}")
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
