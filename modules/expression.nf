process intersect_psites {

    // Intersect a reference BED with p-site positions with p-sites from a sample

    tag "${meta.sample_id}"
    label "intersect_psites"
    publishDir "${outdir}/bedfiles", mode: 'copy'

    input:
    tuple val(meta), val(sample_psite_bed)
    val ref_psite_bed
    val outdir

    output:
    tuple val(meta), path("${meta.sample_id}_intersect.bed"), emit: sample_intersect_bed

    script:
    """
    bedtools intersect \
      -a ${sample_psite_bed} \
      -b ${ref_psite_bed} \
      -wa \
      -wb \
      -header \
      -f 1.00 \
      -s \
      -sorted > "${meta.sample_id}_intersect.bed"
      """
}

process ppm_matrix {

    // Create a matrix object for raw P-sites and PPM for each ORF and
    // each sample included in the cohort

    label "calculate_matrix"
    publishDir "${outdir}/orf_expression", mode: 'copy'

    input:
    val ref_psite_bed
    val sample_intersect_bed
    val orfcaller_name
    val outdir

    output:
    path("${orfcaller_name}_psites_permillion.txt"), emit: ppm_matrix
    path("${orfcaller_name}_psites.txt"), emit: psite_matrix

    script:

    """
    Rscript psite_matrix.R \
    "${ref_psite_bed}" \
    "${sample_intersect_bed}" \
    "${orfcaller_name}"
    """
}
