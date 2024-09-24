include { trimming } from "../modules/fastp.nf"
include { fastqc } from "../modules/fastqc.nf"

workflow SELECTION {

    // Initial selection of RPF reads for mapping
    // Also outputs QC based on fastq information

    take:
    reads               // Riboseq reads with associated sample ID in tuple format
    bowtie2_index       // Precomputed contaminants index for bowtie2
    contaminants_fasta  // Fasta file with rRNA, tRNA, and other contaminants
    outdir              // Output directory

    main:
    // Run FASTP
    trimming(reads, outdir)
    // Gather output tuples into single list

    // Run FASTQC

    // Check if bowtie2 index needs to run
    if (!path(bowtie2_index).exists()) {
        // Run bowtie2 - index
        bowtie_index(contaminants_fasta, outdir)
        bowtie2_index_ch = bowtie_index.out.bowtie2_index
    } else {
    bowtie2_index_ch = value(path(params.bowtie2_index))
    }

    // Run bowtie2
    // Gather output tuples into single list
    rpf_reads = bowtie2.out.reads.collect()
    samfiles = bowtie2.out.samfile.collect()

    emit:
    rpf_reads   // Selected riboseq reads
    samfiles    // Reads considered to be non-RPF
}