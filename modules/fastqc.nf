// used fastqc NFcore module heavily as inspiration
process fastqc {
    // generated QC files

    tag "${sample_id}"
    label "fastqc"
    publishDir "${outdir}/fastqc", mode: 'copy'

    cpus 1
    time '12h'
    memory '2 GB'

    input:
    tuple val(sample_id, path(reads))

    output:
    tuple val(sample_id), path("*.html"), emit: fastqc_html
    tuple val(sample_id), path("*.zip") , emit: fastqc_zip

    script:
    // The total amount of allocated RAM by FastQC is equal to the number of threads defined (--threads) time the amount of RAM defined (--memory)
    // https://github.com/s-andrews/FastQC/blob/1faeea0412093224d7f6a07f777fad60a5650795/fastqc#L211-L222
    // Dividing the task.memory by task.cpu allows to stick to requested amount of RAM in the label
    def memory_in_mb = MemoryUnit.of("${task.memory}").toUnit('MB') / task.cpus
    // FastQC memory value allowed range (100 - 10000)
    def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb)
    
    """
    fastqc \
        --threads ${cpus} \
        --memory ${fastqc_memory} \
        --outdir "${sample_id}"
    """

}
