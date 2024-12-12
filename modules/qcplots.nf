process riboseqc_plots {

    label "riboseqc_plots"
    publishDir "${outdir}/qc", mode: 'copy'

    input:
    val data_files
    val outdir
    val pandoc_dir
    val render_file
    val orfquant_prefix

    output:
    "*.html"

    script:
    """
    Rscript riboseqc_html.R \
    ${input_files}
    ${pandoc_dir} \
    ${render_file} \
    ${orfquant_prefix}
    """
}

process create_qc_plots {

    label "qc_plots"
    publishDir "${outdir}/qc", mode: 'copy'

    input:
    val multiqc_input
    val contaminants_input
    val riboseqc_input
    val outdir

    output:
    "*.pdf"

    script:

    """
    """
}
