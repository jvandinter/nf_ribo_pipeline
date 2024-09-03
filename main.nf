include { printStartLogs; check_params} from "./modules/helper_functions.nf"
include { reads_qc } from './subworkflows/reads_qc'
include { rpf_filtering } from './subworkflows/rpf_filtering'
include { riboseq_alignment } from './subworkflows/riboseq_alignment'
include { riboseq_qc } from './subworkflows/riboseq_qc'
include { orf_expression } from './subworkflows/orf_expression'
include { orf_annotation } from './subworkflows/orf_annotation'

workflow {
    // Initialise workflow
    printStartLogs()
    check_params()

    // Create reads input channel
    reads =  channel.fromPath(params.reads_path, , checkIfExists: true)
    // assign R1 and resolve symlinks
    r1 = reads.map { keys, files -> new File(files[0].toString()).canonicalPath }

   //write R1 file
    r1.collectFile(
        name: 'r1_files.txt',
        storeDir: "${params.project_folder}/documentation/",
        newLine: true, sort: true)

    //write samples file
    reads
    .map { keys, files -> keys }
    .set { sample_ids }

     sample_ids
    .collectFile(
        name: 'sample_ids.txt',
        storeDir: "${params.project_folder}/documentation/",
        newLine: true, sort: true)

    // Step 01: QC
    if (params.qc) {
        rnaseq_qc(reads, params.outdir)
        star_input = rnaseq_qc.out.trimmed_reads
    } else {
        // Look for trimmed reads at the usual location
        default_trimmed_reads = "${params.outdir}/trimgalore/**/*R1_trimmed*.{fastq.gz,fq.gz}"
        if (file(default_trimmed_reads).isEmpty()) {
            star_input = reads
            println  "No trimmed reads found in path: ${default_trimmed_reads}, using ${params.reads_path}"
        } else {
            star_input = Channel
            .fromPath("${default_trimmed_reads}")
        }

        // Step 04: Align
    if (params.align){
        rnaseq_alignment(star_input, params.paired_end, params.outdir)
        bam = rnaseq_alignment.out.bam
    } else {
         // Look for alignment files
        default_bams = "${params.outdir}/star/**/*.Aligned.sortedByCoord.out.bam"
        end2end_bams = "${params.outdir}/star/**/*.Aligned.sortedByCoord.out.bam"
        if (params.bam_files) {
            Channel.fromFilePairs(params.bam_files, size: 1, checkIfExists: true)
            .set { bam }
        } else if (!file(default_bams).isEmpty()) {
            Channel.fromFilePairs(default_bams, size: 1, checkIfExists: true)
            .set { bam }
        } else {
            bam = null
        }
    }
}
