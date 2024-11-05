include { riboseqc; create_annotation } from "../modules/riboseqc.nf"
include { multiqc } from "../modules/multiqc.nf"
include { create_qc_plots } from '../modules/qcplots.nf'

workflow RIBOQC {

/*
Take output from the previous processes and generate QC figures
using multiQC and riboseQC
*/

    take:
    orfquant_annotation // ORFquant annotation file
    orfquant_annot_package // ORFquant annotation R package
    package_install_loc // Location where R package is installed
    pandoc_dir // Location of pandoc for R HTML creation
    orfquant_bams // Output from ALIGNMENT subworkflow
    orfquant_annotation_exists // Boolean: whether to generate annotation
    gtf // Transcriptome used for STAR
    twobit // UCSC file format for the fasta
    outdir // Output directory
    orfquant_prefix // Naming of ORFquant files
    //contaminants // RPF filtering statistics
    //star_output // STAR statistics per file
    //samtools_output // SAMTOOLS statistics per file

    main:

    if(!orfquant_annotation_exists) {
        // Create riboseqc annotation
        create_annotation(gtf,
                          twobit,
                          package_install_loc,
                          orfquant_prefix)
        rannot_ch = create_annotation.out.orfquant_annotation
        package_ch = create.annotation.out.annotation_package
    } else {
        rannot_ch = orfquant_annotation
        package_ch = orfquant_annot_package
    }

    // Create riboseqc files
    riboseqc(orfquant_bams,
             outdir,
             rannot_ch,
             pandoc_dir,
             package_ch,
             package_install_loc)

    // Check annotation status
    if(!orfquant_annotation_exists) {
        // Create riboseqc annotation
        create_annotation(gtf,
                          twobit,
                          package_install_loc,
                          orfquant_prefix)
        rannot_ch = create_annotation.out.orfquant_annotation
        package_ch = create.annotation.out.annotation_package
    } else {
        rannot_ch = orfquant_annotation
        package_ch = orfquant_annot_package
    }

    // Create riboseqc files
    riboseqc(orfquant_bams,
             outdir,
             rannot_ch,
             pandoc_dir,
             package_ch,
             package_install_loc)

    /*
    // run multiqc
    multiqc(samtools,
            star,
            outdir)

    // Create figures
    create_qc_plots(multiqc.out.data_files,
            riboseqc.out.data_files,
            contaminants,
            outdir)
    */
    // Gather output tuples into single list
    for_orfquant_files = riboseqc.out.orfquant_psites

    emit:
    rannot_ch
    package_ch
    for_orfquant_files // files for ORFquant

}
