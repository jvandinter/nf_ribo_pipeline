process intersect_psites {

    // Intersect a reference BED with p-site positions with p-sites from a sample

    tag "${meta.sample_id}"
    label "sample_psites"
    publishDir "${outdir}/bedfiles", mode: 'copy'

    input:
    sample_psite_bed
    ref_psite_bed
    outdir

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

    input:
    ref_psite_bed
    sample_intersect_bed
    ref_name

    output:

    script:

    """
    Rscript calculate_psites_percodon_p_all.R \
    "${ref_psite_bed}" \
    "${sample_intersect_bed}" \
    "${ref_name}"
    """
}
