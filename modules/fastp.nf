// used FASTP NFcore module heavily as inspiration
process trimming {

    // Remove low-quality and short reads from fastq file

    tag "${sample_id}"
    label "trimming"
    publishDir "${outdir}/fastp", mode: 'copy'

    cpus 4
    time '12h'
    memory '4 GB'

    input:
    tuple(val(sample_id),path(reads))

    output:
    tuple val(meta), path('*.fastp.fastq.gz') , emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.log')            , emit: log

    script:
    """
    fastp \
    --thread ${cpu} \
    --in1 "${reads}" \
    --out1 "${sample_id}.fastp.fastq.gz" \
    --json ${sample_id}.fastp.json \
    --html ${sample_id}.fastp.html \
    --length_required 25 \
    2> >(tee ${sample_id}.fastp.log >&2)
    """
}
