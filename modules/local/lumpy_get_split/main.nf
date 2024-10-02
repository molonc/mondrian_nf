process GETSPLIT {
    time '48h'
    cpus 16
    memory '32 GB'
    label 'process_high'

  input:
    path(bam)
  output:
    path("split.sorted.bam"), emit: bam
  script:
    """
        samtools view -h ${bam} | lumpy_extractSplitReads_BwaMem -i stdin | samtools view -Sb - > split.bam
        samtools sort split.bam -o split.sorted.bam
    """

}
