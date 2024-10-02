process GENERATEBARCODES {
    time '48h'
    cpus 8
    memory '24 GB'
    label 'process_high'


    input:
    path(bamfile)
    path(baifile)
    val(filename)
    output:
    path(filename)
    script:
    """
        snv_genotyping_utils generate-cell-barcodes --bamfile ${bamfile} --output ${filename}
    """
}