process GETDISCORDANT {
    time '48h'
    cpus 1
    memory '6 GB'
    label 'process_high'

  input:
    path(bam)
  output:
    path("discordant.sorted.bam"), emit: bam
  script:
    """
        samtools view -bh -F 1294 ${bam} > discordant.bam
        samtools sort discordant.bam -o discordant.sorted.bam
    """

}
