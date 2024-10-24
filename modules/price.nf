process price_index {

    label "price_index"

    input:
    path fasta
    path gtf
    val price_prefix
    val gedi_exec_loc

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
    path bamlist
    path price_index
    val price_prefix
    val gedi_exec_loc
    val outdir

    output:
    path "${price_prefix}.cit.bed", emit: price_orfs
    path "${price_prefix}.*"

    script:

    """
    ${gedi_exec_loc}/gedi -e Price \
    -reads ${bamlist} \
    -genomic ${price_index} \
    -prefix "${price_prefix}"

    ${gedi_exec_loc}/gedi Nashorn \
    -e 'load("'${price_prefix}.cit'").ei().map(function(o) new BedEntry(o.data.getStartStop(o,true).toMutable().setData(new NameAnnotation(o.data.getGeneId()+"__"+o.data.getTranscript()+"__"+o.data.getType()+"__"+o.data.getOrfid()+"__"+o.data.getStartCodon())))).print()' > "${price_prefix}.cit.bed"
    """

}
