process prepare_orfquant {

    label "orfquant_prep"
    publishDir "${outdir}/orfquant", mode: 'copy'

    input:
    val psites_file
    val orfquant_prefix
    val outdir

    output:
    path "${orfquant_prefix}_for_ORFquant", emit: psites_merged

    script:

    """
    merge_psites.R \
    ${psites_file} \
    ${orfquant_prefix}
    """
}

process orfquant {

    label "orfquant"
    publishDir "${outdir}/orfquant", mode: 'copy'

    input:
    val psites_merged
    val orfquant_prefix
    val rannot
    val pandoc_dir
    val orfquant_annot_package
    val package_install_loc
    val outdir

    output:
    path "${orfquant_prefix}_final_ORFquant_results", emit: orfquant_results_file
    path "${orfquant_prefix}_Detected_ORFs.gtf", emit: orfquant_gtf
    path "${orfquant_prefix}_Protein_sequences.fasta", emit: orfquant_fasta
    path "${orfquant_prefix}_*"

    script:
    """
    run_ORFquant.R \
    ${psites_merged} \
    ${orfquant_prefix} \
    ${rannot} \
    $task.cpus \
    ${pandoc_dir} \
    ${orfquant_annot_package} \
    ${package_install_loc} \
    ${outdir}
    """
}
