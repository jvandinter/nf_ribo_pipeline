def printHeader () {
    def logMessage =  """
        .-------------------------------------------------------.
        |                _____                 _      _     _   |
        | _ _ ___ ___   |  |  |___ ___ ___ ___| |_   | |___| |_ |
        || | | .'|   |  |     | -_| -_|_ -|  _|   |  | | .'| . ||
        | \\_/|__,|_|_|  |__|__|___|___|___|___|_|_|  |_|__,|___||
        '-------------------------------------------------------'

        ${workflow.manifest.name} ${workflow.manifest.version}
        ==========================
        """
        log.info logMessage.stripIndent()
}

//TODO: make this work for the riboseq pipeline
def check_files(name, path, type) {
    if(!path) {
        error "When running merge without assembly you must provide `${name}`."
    } else {
        if (type == "dir") {
            file_to_check = file(path, type: "dir")
        } else {
            file_to_check = file(path, type: "file")
        }
        if (file_to_check != null) {
            // Check if it's a list of files
            if (file_to_check instanceof List) {
                // Check each file in the list
                file_to_check.flatten().each { file ->
                    if (!file.exists()) {
                        error "--${name}: ${type} doesn't exist, check path ${file}"
                    }
                }
            } else {
                // Check the single file
                if (!file_to_check.exists()) {
                    error "--${name}: ${type} doesn't exist, check path ${path}"
                }
            }
        } else {
            error "--${name}: No files found at ${path}"
        }
    }
}

def checkInputFiles() {
    //Check inputs

    // Locate bams
    default_bams = "${params.outdir}/star/**/*.Aligned.sortedByCoord.out.bam"
    bam_avail = true
    if (!params.align) {
        if (!params.bam_files && !file(default_bams).isEmpty()) {
            log.info "bam files     : ${default_bams}".stripIndent()
        } else if (params.bam_files && !file(params.bam_files).isEmpty()) {
            log.info "bam files     : ${params.bam_files}".stripIndent()
        } else {
            bam_avail = null
        }
    }

    // Check reference_gtf
    check_files("reference_gtf", params.reference_gtf, "file")

    // Check contaminants fasta

    // Check bowtie2 index
    if (params.qc){
        check_files("kallisto_index", params.kallisto_index, "file")
    }

    // Check references for build_annotaiton
    if (params.build_annotation) {
        check_files("twobit", "${params.twobit}*", "file")
    }

    log.info "\n==========================\n"
}

process write_psites_paths {
    input:
    val collected_paths

    output:
    path "file_paths.txt", emit: psites_file_channel

    script:
    """
    printf "%s\n" "${collected_paths.join('\n')}" > file_paths.txt
    """
}
