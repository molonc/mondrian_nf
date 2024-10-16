process SPLITBAM {
    time '48h'
    cpus 16
    memory '100 GB'
    label 'process_high'

  input:
    path(bamfile)
    path(baifile)
    val(chromosomes)
    val(num_threads)
  output:
    path("outdir/*bam"), emit: bams
  script:
    def chromosomes = "--chromosomes " + chromosomes.join(" --chromosomes ")
    """
        io_utils split-bam-by-barcode \
          --infile ${bamfile} \
          --outdir outdir --tempdir tempdir \
          ${chromosomes} \
          --ncores ${num_threads}
    """
}
