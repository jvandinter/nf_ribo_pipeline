// trimming_test.nf

nextflow.enable.dsl=2

// Load the trimming process from the main script
include { trimming } from '../modules/fastp.nf'
include { fastqc } from '../modules/fastqc.nf'

// Test parameters

workflow {

    // Mock input
    ch_input = Channel.fromPath(params.input)
    .splitCsv(header:true)
    .map { row -> 
        def meta = [:]
        meta.sample_id = row.sample_id
        meta.subject_id = row.subject_id
        meta.sample_type = row.sample_type
        meta.sequence_type = row.sequence_type
        meta.file_type = row.file_type
        // Use 'realpath' to resolve the symlink
        def resolvedPath = "realpath ${row.filename_1}".execute().text.trim()
        
        [ meta, file(resolvedPath) ]
    }

    // Filter for ribo-seq fastq files
    ch_reads = ch_input.filter { meta, fastq ->
        meta.file_type == "fastq" && meta.sequence_type == "ribo"
    }
    ch_reads.view()
    def outdir = file(params.outdir)

    // Run the trimming process
    trimming(ch_reads, outdir)
    trimming_output = trimming.out.reads.collect()

    // Run FASTQC
    fastqc(trimming_output, outdir)
}
