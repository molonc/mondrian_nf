process RECOPY {
    time '48h'
    cpus 8
    memory '24 GB'
    label 'process_high'

  input:
    path(input_path, stageAs: "?/*")
    val(filename)
  output:
    path("${filename}"), emit: output_path
  script:
    """
        cp ${input_path} ${filename}
    """
}
