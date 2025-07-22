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
--finemapped_eqtl           Folder with fine-mapped eQTL loci.
--allele_info               Parquet file with alleles and SNP positions for eQTL dataset.
--finemapped_gwas           Folder with fine-mapped GWAS loci.
--gwas_info_file            File with GWAS file name and GWAS type (quant/cc).
--ld_reference              Folder with permuted LD reference files.
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
params.genefilter = 'data/help_input.txt'

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
include { UNTAR; PARSELOCI; MAKEANNOTATIONTABLE; COLOCLBF; COLOCLBFLEAN; PARSEGENENAMES; PARSEGWAS; FINEMAPGWAS; Untar; ParseLoci; MakeAnnotationTable; ColocLbf; ColocLbfLean; ParseGeneNames; ParseGwas; FinemapGwas} from './modules/eQtlGwasColoc.nf'

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

// Define argument channels
finemapped_eqtl_ch = Channel.fromPath(params.finemapped_eqtl, type: 'dir')
finemapped_eqtl_files_ch = Channel.fromPath("${params.finemapped_eqtl}/*.gz")

allele_ch = Channel.fromPath(params.allele_info)
gtf_ch = Channel.fromPath(params.gtf)

finemapped_gwas_ch = Channel.fromPath(params.finemapped_gwas, type: 'dir')
finemapped_gwas_files_ch = Channel.fromPath("${params.finemapped_gwas}/*.gz")

gwas_info_file_ch = Channel.fromPath(params.gwas_info_file)
ld_ch = Channel.fromPath(params.ld_reference, type: 'dir')

genefilter_ch = Channel.fromPath(params.genefilter)

coloc_pp4_ch = Channel.value(params.coloc_threshold)
gwas_window_ch = Channel.value(params.gwas_window)
full_output_ch = Channel.value(params.full_output)
mrlink2_ch = Channel.value(params.mrlink2)

workflow {
       MAKEANNOTATIONTABLE(finemapped_eqtl_ch.combine(finemapped_gwas_ch).combine(genefilter_ch))

       MAKEANNOTATIONTABLE.out
         .splitCsv(header: true, sep: '\t')
         .map { row ->
             def eqtl_gene = row.eqtl_gene
             def eqtl_file = row.eqtl_file
             def gwas_file = row.gwas_file
             return [eqtl_gene, eqtl_file, gwas_file]
         }
         .set { gene_file_tuples_ch }


        eqtl_loci_ch = finemapped_eqtl_files_ch.map { file ->
            def name = file.getName()
            tuple(name, file)
        }

        gwas_loci_ch = finemapped_gwas_files_ch.map { file ->
            def name = file.getName()
            tuple(name, file)
        }

        gene_file_tuples_ch
        .map { row ->
            def (gene, eqtl_file, gwas_file) = row
            tuple(eqtl_file, gene, gwas_file)
        }
        .combine(eqtl_loci_ch, by: 0)
        .map { row ->
             def (eqtl_file, gene, gwas_file, eqtl_path) = row
             tuple(gwas_file, eqtl_file, gene, eqtl_path)
        }
        .combine(gwas_loci_ch, by: 0)
        .map { row ->
            def (gwas_file, eqtl_file, gene, eqtl_path, gwas_path) = row
            tuple(gene, eqtl_path, gwas_path)
        }
        .unique()
        .groupTuple()
        .map { row ->
        def (gene, eqtl_path, gwas_path) = row
            def unique_eqtl = eqtl_path*.toString().unique().collect { file(it) }
            def unique_gwas = gwas_path*.toString().unique().collect { file(it) }
            tuple(gene, unique_eqtl, unique_gwas)
        }
        .collate(1)
        .map { batch ->
        def genes      = batch.collect { it[0] }
        def eqtl_paths  = batch.collect { it[1] }
        def gwas_paths = batch.collect { it[2] }
        tuple(eqtl_paths, gwas_paths)
        }
        .map { eqtl_list, gwas_list ->
        def flat_eqtl = eqtl_list.flatten().unique()
        def flat_gwas = gwas_list.flatten().unique()
        return [flat_eqtl, flat_gwas]
        }
        .set { gene_file_tuples_ch2 }
        gene_file_tuples_ch2

        PARSEGENENAMES(gtf_ch)

        coloc_input_ch = MAKEANNOTATIONTABLE.out
        .combine(finemapped_gwas_ch)
        .combine(coloc_pp4_ch)
        .combine(full_output_ch)
        .combine(PARSEGENENAMES.out)
        .combine(allele_ch)
        .combine(gene_file_tuples_ch2).view()

        COLOCLBFLEAN(coloc_input_ch)
        annotate_input_ch = COLOCLBFLEAN.out.flatten().filter { it.name.endsWith('_coloc_results.txt') }.collectFile(name: 'EqtlGwasColocSusieResults.txt', keepHeader: true, sort: true, storeDir: "${params.OutputDir}")

    
    }


workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
