process INFERHAPSMETADATA {
    time '48h'
    cpus 1
    memory '6 GB'
    label 'process_high'

  input:
    path(csv)
    path(yaml)
    path(metadata_input)
  output:
    path("metadata.yaml"), emit: metadata
  script:
    """
        haplotype_utils generate-infer-haps-metadata \
        --csv ${csv} --yaml ${yaml} --metadata_yaml ${metadata_input} \
        --metadata_output metadata.yaml
    """
}
