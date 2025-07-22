#!/bin/bash nextflow

process Untar {
    scratch '$TMPDIR'

    input:
        path(tar_archive)

    output:
        path("finemapped.result.*")

    script:
    """
    tar -xf ${tar_archive}
    """
}

process ParseLoci {
    scratch '$TMPDIR'

    input:
        tuple path(raw_file), path(ref), path(gtf)

    output:
       path("ENSG*txt.gz")
       path("*_AnnotationFile.txt")

    script:
        """
        ParseLoci.R \
        --input_file ${raw_file} \
        --reference ${ref} \
        --gtf ${gtf}
        """
}

process ParseGwas {
    scratch '$TMPDIR'

    input:
        tuple path(gwas_folder), val(gwas_id), val(gwas_type), val(gwas_window)

    output:
       path("*txt.gz")

    script:
        """
        ParseGwas.R \
        --gwas_folder ${gwas_folder} \
        --gwas_id ${gwas_id} \
        --win_size ${gwas_window}
        """
}

process FinemapGwas {
    scratch '$TMPDIR'

    publishDir "${params.OutputDir}/GWAS_finemap", mode: 'copy', overwrite: true, pattern: "*.txt.gz"
    
    input:
        tuple path(ld), path(gwas_info_file), path(maf_file), path(loci) 

    output:
       path("*gwas.txt.gz"), optional: true

    script:
        """
        FinemapGwas.R \
        --gwas_info_file ${gwas_info_file} \
        --maf_file ${maf_file} \
        --ld_folder ${ld}
        """
}

process MakeAnnotationTable {
    scratch '$TMPDIR'
    
    publishDir "${params.OutputDir}", mode: 'copy', overwrite: true, pattern: "GwasEqtlOverlaps.txt"

    input:  
        tuple path(finemapped_eqtl, stageAs: 'finemapped_eqtl'), path(finemapped_gwas, stageAs: 'finemapped_gwas'), path(genefilter)

    output:
        path("EqtlGwasOverlaps.txt")

    script:
    """
    ConstructAndAnnotateTable.R \
    --finemapped_eqtl finemapped_eqtl \
    --finemapped_gwas finemapped_gwas \
    --genefilter ${genefilter}
    """

}

process ParseGeneNames {
    scratch '$TMPDIR'

    tag "${gene}"
    input:
        path(gtf)
    output:
        path("GeneNames.txt")

    shell:
        """
        ParseGeneNames.R \
        --gtf ${gtf} 
        """
}

process ColocLbf {
    scratch '$TMPDIR'

    tag "${gene}"
    input:
        tuple path(overlapfile), path(gwas_ch), path(ld), val(coloc_th), val(full_output), val(mrlink2), path(gene_names), path(ref), path(eqtl_paths), path(gwas_paths)
    output:
        path("*_coloc_results.txt")

    shell:
        """
        Coloc_LBF.R \
        --overlap_file EqtlGwasOverlaps.txt \
        --coloc_threshold ${coloc_th} \
        --gene_names ${gene_names} \
        --reference ${ref} \
        --ld_folder ${ld} \
        --full_results ${full_output} \
        --mrlink2 ${mrlink2}
        """
}


process ColocLbfLean {
    scratch '$TMPDIR'

    tag "${gene}"
    input:
        tuple path(overlapfile), path(gwas_ch), val(coloc_th), val(full_output), path(gene_names), path(ref), path(eqtl_paths), path(gwas_paths)
    output:
        path("*_coloc_results.txt")

    shell:
        """
        Coloc_LBF_lean.R \
        --overlap_file EqtlGwasOverlaps.txt \
        --coloc_threshold ${coloc_th} \
        --gene_names ${gene_names} \
        --reference ${ref} \
        --full_results ${full_output}
        """
}



workflow UNTAR {
    take:
        data

    main:
        untar_ch = Untar(data)
        
    emit:
        untar_ch

}

workflow PARSELOCI {
    take:
        data

    main:
        ParseLoci_output_ch = ParseLoci(data)

    emit:
        loci_ch = ParseLoci.out[0]
        log_ch  = ParseLoci.out[1]
}

workflow PARSEGWAS {
    take:
        data

    main:
        ParseGwas_output_ch = ParseGwas(data)

    emit:
        ParseGwas_output_ch
}

workflow FINEMAPGWAS {
    take:
        data

    main:
        FinemapGwas_output_ch = FinemapGwas(data)

    emit:
        loci_ch = FinemapGwas.out[0]
        log_ch  = FinemapGwas.out[1]
}

workflow PARSEGENENAMES {
    take:
        data

    main:
        ParseLoci_output_ch = ParseGeneNames(data)

    emit:
        ParseLoci_output_ch
}


workflow MAKEANNOTATIONTABLE {
    take:
        data

    main:
        MakeAnnotationTable_output_ch = MakeAnnotationTable(data)

    emit:
        MakeAnnotationTable_output_ch
}

workflow COLOCLBF {
    take:
        data

    main:
        ColocLbf_output_ch = ColocLbf(data)

    emit:
        ColocLbf_output_ch
}


workflow COLOCLBFLEAN {
    take:
        data

    main:
        ColocLbf_output_ch = ColocLbfLean(data)

    emit:
        ColocLbf_output_ch
}
