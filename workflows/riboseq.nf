include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { printStartLogs} from "./modules/helper_functions.nf"
include { SELECTION } from './subworkflows/selection.nf'
include { ALIGNMENT } from './subworkflows/alignment.nf'
include { RIBOQC } from './subworkflows/riboqc.nf'
include { ORFQUANT } from './subworkflows/orfquant.nf'
include { PRICE } from './subworkflows/price.nf'
include { EXPRESSION } from './subworkflows/expression.nf'
include { ANNOTATION } from './subworkflows/annotation.nf'

workflow RIBOSEQ {

    // Validate input parameters
    validateParameters()

    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    // Validate files
    checkInputFiles()

    //Create input channel from samplesheet
    ch_input = fromList(samplesheetToList(params.input, "assets/schema_input.json"))

    //Get channel for reads
    ch_fastq = ch_input
        .filter { meta, filename_1
            meta.filetype == "fastq" && meta.sequence == "ribo"
        }
        .map { meta, fastq_1, ->
                meta, [ fastq_1 ]
        }
    ch_reads = ch_fastq
        .multiMap { meta, reads ->
            reads: tuple(meta.id, reads)
            paired_end:   meta.paired_end
        }
    reads = ch_reads.reads
    // TODO: handle annotations from samplesheet

}
