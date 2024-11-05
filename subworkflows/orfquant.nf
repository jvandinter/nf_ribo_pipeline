include { prepare_orfquant; orfquant } from '../modules/orfquant.nf'

process writePsitesPaths {
    input:
    val collected_paths

    output:
    path "psites_paths.txt", emit: psites_file_channel

    script:
    """
    printf "%s\n" "${collected_paths.join('\n')}" > psites_paths.txt
    """
}

workflow ORFQUANT {

    take:
    orfquant_psites
    orfquant_annotation
    orfquant_annot_package
    package_install_loc
    pandoc_dir
    orfquant_prefix
    outdir

    main:

    orfquant_psites
    .map { meta, path -> path } // Extract only the paths
    .collect()                   // Collect all paths into a single list
    .set { collected_paths }     // Set this list into a new channel
    writePsitesPaths(collected_paths)

    prepare_orfquant(writePsitesPaths.out.psites_file_channel,
                     orfquant_prefix,
                     outdir)

    orfquant(prepare_orfquant.out.psites_merged,
             orfquant_prefix,
             orfquant_annotation,
             pandoc_dir,
             orfquant_annot_package,
             package_install_loc,
             outdir)

    emit:
    orfquant.out.gtf

}
