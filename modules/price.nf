process price_index {

    label "price_index"

    input:
    path fasta // Genome fasta used for alingment
    path gtf // Transcriptome GTF used for alignment
    val price_prefix // name of price files
    val gedi_exec_loc // Location of gedi installation, until containerisation works

    output:
    path "${price_prefix}.oml", emit: price_index
    path "${gtf.baseName}.*"
    path "${fasta.baseName}.*"

    script:
    """
    ${gedi_exec_loc}/gedi -e IndexGenome \
        -s "${fasta}" \
        -a "${gtf}" \
        -f "." \
        -o "${price_prefix}.oml" \
        -nobowtie \
        -nostar \
        -nokallisto
    """

}

process price {

    label "price"
    publishDir "${outdir}/price", mode: 'copy'

    input:
    path bamlist // File that lists all BAM files
    path price_index // Index for PRICE
    val price_prefix // Name of price files
    val gedi_exec_loc // Location of gedi installation, until containerisation works
    val outdir // Output directory

    output:
    path "${price_prefix}.orfs.cit.bed", emit: price_orfs
    path "${price_prefix}.*"

    script:

    """
    ${gedi_exec_loc}/gedi -e Price \
        -reads ${bamlist} \
        -genomic ${price_index} \
        -prefix "${price_prefix}"
    
    ${gedi_exec_loc}/gedi Nashorn -e \
        'load("'${price_prefix}.orfs.cit'").ei().map(function(o) new BedEntry(o.data.getStartStop(o,true).toMutable().setData(new NameAnnotation(o.data.getGeneId()+"__"+o.data.getTranscript()+"__"+o.data.getType()+"__"+o.data.getOrfid()+"__"+o.data.getStartCodon())))).print()' \
        > "${price_prefix}.orfs.cit.bed"
    """

}
