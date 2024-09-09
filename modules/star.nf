process star_index {

    // Create index for STAR

    label "STAR index"
    publishDir "${outdir}/star_index/", mode: 'copy'

    input: 
    val fasta       // Reference genome fasta file
    val outdir      // Output directory
    val gtf         // Transcriptome GTF file

    output:
    path 'star_index', emit: star_index_path

    script:
    """
    STAR \
    --runMode genomeGenerate \
    --runThreadN $task.cpu \
    --sjdbGTFfile ${gtf} \
    --sjdbOverhang 29 \
    --genomeDir "star_index" \
    --genomeFastaFiles ${fasta}
    """
}

process star {

    // Aligns RPF reads to the 

    tag ${sample_id}
    label "alignment"
    publishDir "${outdir}/star/", mode: 'copy'

    input: 
    tuple val(sample_id), path(reads)   // Trimmed RPF reads
    val outdir                          // Output directory
    val gtf                             // Transcriptome GTF file
    val star_index                      // STAR index
    val run_price                       // Do we need to generate the PRICE bam

    output:
    path("${sample_id}/${sample_id}.*")
    tuple val("${sample_id}"), path("${sample_id}/${sample_id}.*.bam"), emit: bam

    script:
    // GROOVY START
    def star_params = "--"

    def star_params_price = ""
    // GROOVY END

    // Check which BAM files need to be generated
    if (run_price == true) {

        """
        # ORFquant BAM
        STAR \
        --genomeDir ${star_index} \
        --sjdbGTFfile ${gtf} \
        --readFilesIn ${reads} \
        --outSAMattrRGline ID:${sample_id} LB:${sample_id} PL:IllUMINA SM:${sample_id}
        --outFileNamePrefix "${sample_id}/${sample_id}." \
        --runThreadN $task.cpus \
        ${star_params}

        # PRICE BAM
        STAR \
        --genomeDir ${star_index} \
        --sjdbGTFfile ${gtf} \
        --readFilesIn ${reads} \
        --outSAMattrRGline ID:${sample_id} LB:${sample_id} PL:IllUMINA SM:${sample_id}
        --outFileNamePrefix "${sample_id}/${sample_id}.end2end." \
        --runThreadN $task.cpus \
        ${star_params_price}
        """

    } if (run_price == true && run_orfquant == false) {

        """
        # PRICE BAM
        STAR \
        --genomeDir ${star_index} \
        --sjdbGTFfile ${gtf} \
        --readFilesIn ${reads} \
        --outSAMattrRGline ID:${sample_id} LB:${sample_id} PL:IllUMINA SM:${sample_id}
        --outFileNamePrefix "${sample_id}/${sample_id}.end2end." \
        --runThreadN $task.cpus \
        ${star_params_price}
        """

    } else {

        """
        # ORFquant BAM
        STAR \
        --genomeDir ${star_index} \
        --sjdbGTFfile ${gtf} \
        --readFilesIn ${reads} \
        --outSAMattrRGline ID:${sample_id} LB:${sample_id} PL:IllUMINA SM:${sample_id}
        --outFileNamePrefix "${sample_id}/${sample_id}." \
        --runThreadN $task.cpus \
        ${star_params}
        """

    }

}

process samtools_index {

    tag ${sample_id}
    label "samtools"
    publishDir "${outdir}/star/", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    path("${sample_id}/${sample_id}.*")

    script:
    """
    samtools index \
        -@ $task.cpus \
        ${bam}
    """
}
