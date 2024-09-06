process bowtie_index {

    // Create index for bowtie2 alignment

    label "bowtie2 index"
    publishDir "${outdir}", mode: 'copy'

    input:
    val contaminants_fasta // FASTA with unwanted sequences
    val outdir             // Output directory

    output:
    path 'bowtie2_index', emit: bowtie2_index

    script:
    """
    bowtie2-build \
    -f \
    --seed 24 \
    --threads $tasks.cpu \
    ${contaminants_fasta} \
    bowtie2_index

    """
}

process bowtie {

    // Use BOWTIE2 to align against unwanted sequences such
    // as rRNA, tRNA, snRNA, snoRNA, mtDNA and keep true
    // ribosome-protected fragments for mapping and ORF calling

    tag "${sample_id}"
    label "bowtie2"
    publishDir "${outdir}/bowtie2", mode: 'copy'

    input:
    val bowtie2_index                   // Bowtie2 reference index
    tuple val(sample_id), path(reads)   // Trimmed fastp reads
    val keep_sam                        // Boolean, keep big SAM file for debugging
    val outdir                          // Output directory

    output:
    tuple(sample_id, path("*_filtered.fastq.gz")), emit: trimmed_reads
    tuple(sample_id, path("*_contaminants.sam")), emit: samfile // TODO: check whether this can be empty

    script:
    """
    bowtie2 \
    --seedlen=25 \
    --threads $tasks.cpu \
    --time \
    --un-gz "${sample_id}_filtered.fastq.gz" \
    -x ${bowtie2_index} \
    -U "${reads}" \
    -S "${sample_id}_contaminants.sam"
    if(${keep_sam} == FALSE) {
        rm "${sample_id}_contaminants.sam"
        # Create harmless empty file to not crash the workflow
        touch "${sample_id}_empty_contaminants.sam"
    } fi
    """
}
