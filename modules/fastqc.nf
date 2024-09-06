// Used fastqc NFcore module as inspiration
process fastqc {

    // Generate QC files

    tag "${sample_id}"
    label "fastqc"
    publishDir "${outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id, path(reads)) // Trimmed FASTQ read from fastp
    val outdir                        // Output directory

    output:
    path("*.html") // Output QC summary 
    path("*.zip")  // QC files

    script:
    //TODO: put in config file
    // The total amount of allocated RAM by FastQC is equal to the number of threads defined (--threads) time the amount of RAM defined (--memory)
    // https://github.com/s-andrews/FastQC/blob/1faeea0412093224d7f6a07f777fad60a5650795/fastqc#L211-L222
    // Dividing the task.memory by task.cpu allows to stick to requested amount of RAM in the label
    def memory_in_mb = MemoryUnit.of("${task.memory}").toUnit('MB') / task.cpus
    // FastQC memory value allowed range (100 - 10000)
    def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb)

    """
    fastqc \
        --threads $tasks.cpu \
        --memory $tasks.memory \
        --outdir "${sample_id}"
    """

}
