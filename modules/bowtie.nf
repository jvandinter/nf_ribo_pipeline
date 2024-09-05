process bowtieIndex {

    // Create index for bowtie2 alignment

    tag "${sample_id}"
    label "bowtie2 index"

    cpus 4
    time '12h'
    memory '4 GB'

    input:

    output:

    script:
    """

    """
}

process bowtie {

    // Use BOWTIE2 to align against 

    tag "${sample_id}"
    label "bowtie2"

    cpus 4
    time '12h'
    memory '8 GB'

    input:
    bowtie2_index
    reads
    TMPDIR

    output:
    trimmed_reads

    script:
    """
    bowtie2 \
    --seedlen=25 \
    --threads ${threads} \
    --time \
    --un-gz "${sample_id}_filtered.fastq.gz" \
    -x ${bowtie2_index} \
    -U "${reads}" \
    -S "${TMPDIR}/${sample_id}_contaminants.sam"
    """
}
