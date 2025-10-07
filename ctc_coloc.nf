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
--finemapped_signals \
--ld_reference \
--allele_info \
--gtf \
--OutputDir

Mandatory arguments:
--allele_info               Parquet file with alleles and SNP positions for eQTL dataset.
--finemapped_gwas           Folder with fine-mapped GWAS loci.
--gwas_info_file            File with GWAS file name and GWAS type (quant/cc).
--gtf                       GTF file for gene annotation.

Optional arguments:
--coloc_threshold           COLOC PP4 threshold to declare colocalisation. Defaults to 0.8.
--gwas_window               GWAS locus window size for fine-mapping. Defaults to 1000000 (+/-1Mb from lead variant).
--full_output               Whether to output all coloc results, including those which show no colocalisation (true/false). Defaults to false.
--mrlink2                   Whether to run also MRLink2 analysis (true/false). Defaults to false.
--genefilter                File with optional gene filter to confine. Need to include ENSG IDs and no need for header.
--OutputDir                 Output directory. Defaults to "results".
""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

// Default parameters
params.OutputDir = 'results'
params.coloc_threshold = 0.8
params.full_output = false
params.mrlink2 = false
params.gwas_window = 1000000

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
summary['Fine-mapped eQTL files']                   = params.finemapped_eqtl
summary['Fine-mapped GWAS files']                   = params.finemapped_gwas
summary['GWAS info file']                           = params.gwas_info_file
summary['LD reference']                             = params.ld_reference
summary['Allele info file']                         = params.allele_info
summary['GTF file']                                 = params.gtf
summary['COLOC PP4 threshold']                      = params.coloc_threshold
summary['GWAS window size']                         = params.gwas_window
summary['COLOC full output']                        = params.full_output
summary['MRLink2']                                  = params.mrlink2
summary['Gene Filter']                              = params.genefilter
// import modules
include { UNTAR; PARSELOCI; MAKEANNOTATIONTABLE; COLOCLBF; COLOCLBFLEAN; PARSEGENENAMES; PARSEGWAS; FINEMAPGWAS; Untar; ParseLoci; MakeAnnotationTable; ColocLbf; ColocLbfLean; ParseGeneNames; ParseGwas; FinemapGwas; MakeAnnotationTableGwasPairs; ColocLbfGwas} from './modules/eQtlGwasColoc.nf'

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

// Define argument channels
allele_ch = Channel.fromPath(params.allele_info)
gtf_ch = Channel.fromPath(params.gtf)

gwas_info_file_ch = Channel.fromPath(params.gwas_info_file)
    .splitCsv(header: true, sep: '\t')
    .branch { row ->
         ctc: row.bctc_trait == "TRUE"
         nonctc: row.bctc_trait == "FALSE"
    }

ctc_trait_ch = gwas_info_file_ch.ctc
    .map { row -> row.pheno}.collect()
nonctc_trait_ch = gwas_info_file_ch.nonctc
    .map { row -> row.pheno}.collect()

ctc_gwas_traits = gwas_info_file_ch.ctc.map { row -> row[0] }
nonctc_gwas_traits = gwas_info_file_ch.nonctc.map { row -> row[0] }

finemapped_gwas_ch = Channel.fromPath(params.finemapped_gwas, type: 'dir')
finemapped_gwas_files_ch = Channel.fromPath("${params.finemapped_gwas}/*.gz")

ld_ch = Channel.fromPath(params.ld_reference, type: 'dir')

coloc_pp4_ch = Channel.value(params.coloc_threshold)
gwas_window_ch = Channel.value(params.gwas_window)
full_output_ch = Channel.value(params.full_output)
mrlink2_ch = Channel.value(params.mrlink2)

workflow {
    overlap_ch = MakeAnnotationTableGwasPairs(
        finemapped_gwas_ch.collect(), ctc_trait_ch, nonctc_trait_ch
        )

    overlap_ch_buffered = overlap_ch.splitText(by: 500, keepHeader: true, file: true)

    ColocLbfGwas(
        overlap_ch_buffered,
        finemapped_gwas_ch.collect(),
        coloc_pp4_ch,
        full_output_ch,
        allele_ch.collect())

    annotate_input_ch = ColocLbfGwas.out.colocalizations
        .flatten().filter { it.name.endsWith('_coloc_results.txt') }
        .collectFile(name: 'CtcGwasColocSusieResults.txt', keepHeader: true, sort: true, storeDir: "${params.OutputDir}")

    coloc_summary_output_ch = ColocLbfGwas.out.summary
        .collectFile(name: 'CtcGwasColocSummary.txt', keepHeader: true, skip: 1, sort: true, storeDir: "${params.OutputDir}")
}


workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
