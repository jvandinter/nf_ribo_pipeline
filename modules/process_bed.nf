process ref_psites {

    // Create reference in-frame P-sites file from GTF

    label "ref_psites"
    publishDir "${outdir}/annotation", mode: 'copy'

    input:
    val gtf

    output:
    val(), emit:
    path(${ref_base}.gtf_psites_p0.sorted.bed), emit: ref_psite_bed

    script:
    def ref_base = gtf.baseName.replaceAll(/\.gtf$/, '')
    """
    python3 calculate_psite_coords_final.py \
    -i ${gtf} \
    -a "${is_annotation}" \
    -o "${bedfiles_dir}" \
    -t "${id_type}"

    sort -k1,1 -k2,2n "$${ref_base}.gtf_psites_plus_partial.bed" > "${ref_base}.gtf_psites_p0.sorted.bed"
    """

}

process sample_psites {

    // Create per sample a BED file with sorted in-frame P-sites found in the sample

    tag "${meta.sample_id}"
    label "sample_psites"
    publishDir "${outdir}/bedfiles", mode: 'copy'

  input:
  tuple val(meta),val()
  outdir

  output:
  tuple val(meta), path("${meta.sample_id}_psites.sorted.bed"), emit: sample_psite_bed
  
  script:
    """
    Rscript psites_bed_from_riboseqc.R \
    ${}
    ${meta.sample_id} \
    ${ref_psite_bed} \
    sort -k1,1 -k2,2n "${meta.sample_id}_psites.bed" > "${meta.sample_id}_psites.sorted.bed"
    """
}