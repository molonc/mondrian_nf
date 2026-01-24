process MERGEBAMS {
    time '48h'
    cpus 32
    memory '240 GB'
    label 'process_high'

  input:
    path(bams, stageAs: "?/*")
    val(filename)
  output:
    path("${filename}.bam"), emit: bam
    path("${filename}.bam.bai"), emit: bai
  script:
    def input_bam = ''
    if (bams instanceof nextflow.util.BlankSeparatedList){
        input_bam = '--inputs ' + bams.join(' --inputs ')
    } else {
        input_bam = '--inputs ' + bams
    }
    """
        variant_utils merge-bams \
        $input_bam \
        --output ${filename}.bam \
        --tempdir temp \
        --threads ${task.cpus}
        samtools index ${filename}.bam
    """
}
