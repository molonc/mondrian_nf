process ALIGNMETRICS {
    time '24h'
    cpus 2
    memory '8 GB'
    label 'process_high'

  input:
    tuple(
        val(cell_id), val(lanes), val(flowcells),
        path(primary_reference), val(primary_reference_version), val(primary_reference_name),
        path(primary_reference_fai), path(primary_reference_amb),path(primary_reference_ann),
        path(primary_reference_bwt),path(primary_reference_pac),path(primary_reference_sa),
        path(metadata),
        path("fastqscreen/${cell_id}"),
        path("cleanalignsort/${cell_id}"),
        path("aligned.bam"),
        path("aligned.bam.bai")
    )
  output:
    tuple(
        val(cell_id),
        path("aligned.bam"),
        path("aligned.bam.bai"),
        path("metrics.csv.gz"),
        path("metrics.csv.gz.yaml"),
        path("${cell_id}_gc_metrics.csv.gz"),
        path("${cell_id}_gc_metrics.csv.gz.yaml"),
        path("${cell_id}.tar.gz")
    )
  script:
    
    // alignment_metrics(metadata_yaml, reference, tempdir, cell_id, wgs_metrics_mqual, wgs_metrics_bqual, wgs_metrics_count_unpaired,bam_output, metrics_output,metrics_gc_output,tar_output, num_threads
    

    """
        alignment_utils alignmentmetrics \
        --metadata_yaml ${metadata} \
        --reference ${primary_reference_name},${primary_reference_version},${primary_reference} \
        --tempdir tempdir \
        --cell_id $cell_id \
        --wgs_metrics_mqual 20 \
        --wgs_metrics_bqual 20 \
        --wgs_metrics_count_unpaired false \
        --bam_output aligned.bam \
        --metrics_output metrics.csv.gz \
        --metrics_gc_output ${cell_id}_gc_metrics.csv.gz \
        --tar_output ${cell_id}.tar.gz \
        --num_threads ${task.cpus}

        rm -rf tempdir
    """
}
