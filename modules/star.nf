process star_index {

    // Create index for STAR

    label "STAR_index"
    publishDir "${outdir}/star_index/", mode: 'copy'

    input: 
    val fasta       // Reference genome fasta file
    val outdir      // Output directory
    val gtf         // Transcriptome GTF file

    output:
    val "star_index", emit: star_index_path
    path "star_index/*"

    script:
    """
    STAR \
    --runMode genomeGenerate \
    --runThreadN $task.cpus \
    --sjdbGTFfile ${gtf} \
    --sjdbOverhang 29 \
    --genomeDir "star_index" \
    --genomeFastaFiles ${fasta}
    """
}

process align {

    // Aligns RPF reads to the 

    tag "${meta.sample_id}"
    label "alignment"
    publishDir "${outdir}/star/", mode: 'copy'

    input: 
    tuple val(meta), path(reads)   // Trimmed RPF reads
    val outdir                     // Output directory
    val gtf                        // Transcriptome GTF file
    val star_index_path            // STAR index
    val run_price                  // Do we need to generate the PRICE bam
    val run_orfquant               // Do we need to generate the normal bam

    output:
    path("${meta.sample_id}/${meta.sample_id}.*")
    tuple val(meta), path("${meta.sample_id}/${meta.sample_id}.local.*.bam"), optional: true, emit: bams
    tuple val(meta), path("${meta.sample_id}/${meta.sample_id}.end2end.*.bam"), optional: true, emit: bams_end2end

    script:
    def sample_id = meta.sample_id
    def star_params = "--readFilesCommand zcat " +
                      "--runDirPerm All_RWX " +
                      "--twopassMode Basic " +
                      "--outFilterMismatchNmax 2" +
                      "--outFilterMultimapNmax 20" +
                      "--outSAMattributes All" +
                      "--limitOutSJcollapsed 10000000" +
                      "--limitIObufferSize=300000000" +
                      "--outFilterType BySJout" +
                      "--alignSJoverhangMin 1000" +
                      "--outTmpKeep None"

    def star_params_end2end = "--readFilesCommand zcat " +
                              "--runDirPerm All_RWX" +
                              "--twopassMode Basic" +
                              "--outFilterMismatchNmax 2" +
                              "--outFilterMultimapNmax 20" +
                              "--outSAMattributes All" +
                              "--limitOutSJcollapsed 10000000" +
                              "--limitIObufferSize=300000000" +
                              "--outFilterType BySJout" +
                              "--alignSJoverhangMin 1000" +
                              "--alignEndsType EndToEnd" +
                              "--outTmpKeep None"

    // Check which BAM files need to be generated
    if (run_price == true) {

        """
        # ORFquant BAM
        STAR \
        --genomeDir ${star_index_path} \
        --sjdbGTFfile ${gtf} \
        --readFilesIn ${reads} \
        --outSAMattrRGline ID:${sample_id} LB:${sample_id} PL:IllUMINA SM:${sample_id}
        --outFileNamePrefix "${sample_id}/${sample_id}.local." \
        --runThreadN $task.cpus \
        ${star_params}

        # PRICE BAM
        STAR \
        --genomeDir ${star_index_path} \
        --sjdbGTFfile ${gtf} \
        --readFilesIn ${reads} \
        --outSAMattrRGline ID:${sample_id} LB:${sample_id} PL:IllUMINA SM:${sample_id}
        --outFileNamePrefix "${sample_id}/${sample_id}.end2end." \
        --runThreadN $task.cpus \
        ${star_params_end2end}
        """

    } else if (run_price == true && run_orfquant == false) {

        """
        # PRICE BAM
        STAR \
        --genomeDir ${star_index_path} \
        --sjdbGTFfile ${gtf} \
        --readFilesIn ${reads} \
        --outSAMattrRGline ID:${sample_id} LB:${sample_id} PL:IllUMINA SM:${sample_id}
        --outFileNamePrefix "${sample_id}/${sample_id}.end2end." \
        --runThreadN $task.cpus \
        ${star_params_end2end}
        """

    } else {

        """
        # ORFquant BAM
        STAR \
        --genomeDir ${star_index_path} \
        --sjdbGTFfile ${gtf} \
        --readFilesIn ${reads} \
        --outSAMattrRGline ID:${sample_id} LB:${sample_id} PL:IllUMINA SM:${sample_id}
        --outFileNamePrefix "${sample_id}/${sample_id}." \
        --runThreadN $task.cpus \
        ${star_params}
        """

    }

}
