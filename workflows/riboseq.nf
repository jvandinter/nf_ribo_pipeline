/*
Riboseq Nextflow Pipeline

This pipeline processes ribosome profiling (Riboseq) data, including quality control,
alignment, ORF quantification, and expression analysis.

Main steps:
1. SELECTION: Initial read processing and contaminant removal
2. ALIGNMENT: Align reads to the reference genome
3. RIBOQC: Quality control of Riboseq data
4. ORFQUANT: ORF quantification (optional)
5. PRICE: Alternative ORF calling method (optional)
6. ANNOTATION: Annotate identified ORFs
7. EXPRESSION: Analyze ORF expression
*/

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

    Channel.fromPath(,checkIfExists = true)
    // Validate input parameters
    validateParameters()

    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    // Validate files

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

    reads = ch_reads.reads
    // TODO: handle annotations from samplesheet

    SELECTION(reads,
              params.contaminants_fasta,
              params.outdir)
    
    ALIGNMENT(SELECTION.out.rpf_reads,
              params.genome,
              params.gtf,
              params.star_index_path,
              params.run_price,
              params.run_orfquant,
              params.outdir)
    
    RIBOQC(ALIGNMENT.out.bam_files, gtf_ch)
    
    if (params.run_orfquant) {
        ORFQUANT(ALIGNMENT.out.bam_files, gtf_ch)
    }
    
    if (params.run_price) {
        PRICE(ALIGNMENT.out.end2end_bam_files, params.price_index)
    }
    
    ANNOTATION(
        ORFQUANT.out.orfquant_results,
        PRICE.out.price_results,
        gtf_ch,
        params.run_orfquant,
        params.run_price
    )
    
    EXPRESSION(
        ANNOTATION.out.annotated_orfs,
        ALIGNMENT.out.bam_files,
        params.run_orfquant,
        params.run_price
    )
}
}
