include { star_index; align; samtools } from '../modules/star.nf'

workflow ALIGNMENT {

/*
Aligns the data in two STAR modes for ORFquant and
PRICE respectively and creates an index for the ORFquant bam
*/

    take:
    rpf_reads       // Output from SELECTION subworkflow
    genome              // Reference genome used for STAR index
    star_index_path     // Location of precomputed STAR index
    gtf                 // Transcriptome used for STAR
    run_price           // Boolean: whether to run PRICE
    run_orfquant    // Boolean: whether to run ORFquant
    outdir              // Output directory

    main:

    // List of files to check
    def filesToCheck = [
        'chrLength.txt',
        'chrStart.txt',
        'geneInfo.tab',
        'SA',
        'sjdbList.fromGTF.out.tab',
        'chrNameLength.txt',
        'exonGeTrInfo.tab',
        'Genome',
        'SAindex',
        'sjdbList.out.tab',
        'chrName.txt',
        'exonInfo.tab',
        'genomeParameters.txt',
        'sjdbInfo.txt',
        'transcriptInfo.tab'
    ]

    // Check if STAR index exists
    def fileFound = false
    for (file in filesToCheck) {
        def filePath = file(star_index_path, file)
        if (filePath.exists()) {
            fileFound = true
            star_index_ch = star_index_path
            break
        } else {
            star_index(genome, outdir)
            star_index_ch = star_index.out.star_index_path
        }
    }

    align(rpf_reads, outdir, gtf, star_index_ch, run_price, run_orfquant)

    // Wait untill all STAR runs are completed
    bam_list_end2end = align.out.bams_end2end.collect()

    if (params.run_orfquant = true) {
        samtools(align.out.bams, outdir)
        bam_list = samtools.out.bams.collect()
    }
    
    emit:
    bam_list         // bam files for ORFquant
    bam_list_end2end // bam files for PRICE
}
