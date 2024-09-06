// used FASTP NFcore module as inspiration
process trimming {

    // Remove low-quality and short reads from fastq file

    tag "${sample_id}"
    label "trimming"
    publishDir "${outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample_id),path(reads)    // Input FASTQ reads
    val outdir                          // Output directory

    output:
    tuple val(sample_id), path('*.fastp.fastq.gz') , emit: reads
    path('*.json')
    path('*.html')
    path('*.log')

    script:
    """
    fastp \
    --thread $tasks.cpu \
    --in1 "${reads}" \
    --out1 "${sample_id}.fastp.fastq.gz" \
    --json ${sample_id}.fastp.json \
    --html ${sample_id}.fastp.html \
    --length_required 25 \
    2> >(tee ${sample_id}.fastp.log >&2)
    """
}
