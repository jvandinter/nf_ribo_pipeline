process multiqc {

    label "multiqc"
    publishDir "${outdir}/multiqc", mode: 'copy'

    input:
    val samtools_input
    val contaminants_input
    val fastp_input
    val outdir

    output:


    script:

    """
    multiqc
    """

}
