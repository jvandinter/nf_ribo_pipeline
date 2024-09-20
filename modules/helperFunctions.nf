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
    default_strand =  "${params.outdir}/check_strandedness/strandedness_all.txt"
    if (!params.qc && ( params.align || params.assembly )) {
        if (!params.strand_info && file(default_strand, type : "file").exists()) {
            log.info "strand info   : ${default_strand}".stripIndent()
        } else if (params.strand_info && file(params.strand_info, type : "file").exists()) {
            log.info "strand info   : ${params.strand_info}".stripIndent()
        } else {
            error  """
                No strand information found in ${default_strand}
                If `qc = false`, define your strandedness info file with `strand_info` in `params.config`
                """.stripIndent()
        }
    }

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

    // Check star_index_basedir
    if (params.align){
        check_files("star_index_basedir", params.star_index_basedir, "dir")
    }

    // Check references for build_annotaiton
    if (params.build_annotation) {
        check_files("twobit", "${params.twobit}*", "file")
    }

    // Check expression parameters
    if (params.expression) {
        assert params.expression_mode in ["sq", "sa", "sqfc", "safc"], "`expression_mode` must be one of the following: sq, sa, sqfc, safc"

        if (!bam_avail & (params.expression_mode != "sq")) {
            log.info "expression mode  : ${params.expression_mode} --> sq [forced to quasi-mapping, no bam avail]"
        }
    }

    log.info "\n==========================\n"
}
