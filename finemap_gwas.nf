#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2


def helpmessage() {

log.info"""

FinemapGwas v${workflow.manifest.version}"
===========================================================
Pipeline for SuSiE fine-mapping on GWAS summary statistics.

Usage:

nextflow run main.nf 
--finemapped_signals \
--gwas_database \
--ld_reference \
--allele_info \
--maf_file \
--gtf \
--OutputDir

Mandatory arguments:
--finemapped_eqtl           Folder with fine-mapped eQTL loci.
--allele_info               Parquet file with alleles and SNP positions for eQTL dataset.
--gwas_database             Folder with fine-mapped GWAS loci.
--gwas_info_file            File with GWAS file name, GWAS type (quant/cc) and effective N.
--maf_file                  File with variant index and MAF.
--ld_reference              Folder with permuted LD reference files.
--gtf                       GTF file for gene annotation.

Optional arguments:
--coloc_threshold           COLOC PP4 threshold to declare colocalisation. Defaults to 0.8.
--gwas_window               GWAS locus window size for fine-mapping. Defaults to 1000000 (+/-1Mb from lead variant).
--full_output               Whether to output all coloc results, including those which show no colocalisation (true/false). Defaults to false.
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
params.gwas_window = 1000000

//Show parameter values
log.info """===============================================================
eQTLGen cis-trans fine-mapped coloc pipeline v${workflow.manifest.version}"
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
summary['eQTL files']                               = params.finemapped_signals
summary['GWAS files']                               = params.gwas_database
summary['GWAS info file']                           = params.gwas_info_file
summary['MAF file']                                 = params.maf_file
summary['LD reference']                             = params.ld_reference
summary['Allele info file']                         = params.allele_info
summary['GTF file']                                 = params.gtf
summary['COLOC PP4 threshold']                      = params.coloc_threshold
summary['GWAS window size']                         = params.gwas_window
summary['COLOC full output']                        = params.full_output
// import modules
include { UNTAR; PARSELOCI; MAKEANNOTATIONTABLE; COLOCLBF; PARSEGENENAMES; PARSEGWAS; FINEMAPGWAS; Untar; ParseLoci; MakeAnnotationTable; ColocLbf; ParseGeneNames; ParseGwas; FinemapGwas} from './modules/eQtlGwasColoc.nf'

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

// Define argument channels
finemapped_ch = Channel.fromPath(params.finemapped_signals, type: 'dir')
allele_ch = Channel.fromPath(params.allele_info)
gtf_ch = Channel.fromPath(params.gtf)
gwas_ch = Channel.fromPath(params.gwas_database, type: 'dir')
gwas_info_file_ch = Channel.fromPath(params.gwas_info_file)
maf_file_ch = Channel.fromPath(params.maf_file)
ld_ch = Channel.fromPath(params.ld_reference, type: 'dir')
coloc_pp4_ch = Channel.value(params.coloc_threshold)
gwas_window_ch = Channel.value(params.gwas_window)
full_output_ch = Channel.value(params.full_output)

workflow {
        
        gwas_info_file_ch
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def gwas_name = row.pheno
            def gwas_type = row.type
            return [gwas_name, gwas_type]
        }
        .set { gwas_info_tuple_ch }

       PARSEGWAS(gwas_ch.combine(gwas_info_tuple_ch).combine(gwas_window_ch))

       PARSEGWAS.out.
       collect().
       flatten().
       buffer( size: 10, remainder: true ).
       map { items -> [items] }.  
       set{ finemap_input_ch }
       
       ld_ch.
       combine(gwas_info_file_ch).
       combine(maf_file_ch).
       combine(finemap_input_ch).
       set{ finemap_input_ch }
       
       FINEMAPGWAS(finemap_input_ch)

       gwas_loci_ch = FINEMAPGWAS.out[0].flatten()
       gwas_log_ch  = FINEMAPGWAS.out[1].flatten()

       gwas_log_ch.collectFile(name: 'gwas_combined_analysis.log', keepHeader: true, skip: 1).set { gwas_final_log_ch }
    
    }


workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
