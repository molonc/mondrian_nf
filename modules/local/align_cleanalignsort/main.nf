process ALIGNCLEANALIGNSORT {
    time '24h'
    cpus 8
    memory '8 GB'
    label 'process_high'

  input:
    tuple(
      val(cell_id), val(lanes), val(flowcells),
      path(primary_reference), val(primary_reference_version), val(primary_reference_name),
      path(primary_reference_fai), path(primary_reference_amb),path(primary_reference_ann),
      path(primary_reference_bwt),path(primary_reference_pac),path(primary_reference_sa),
      path(metadata),
      path("fastqscreen/${cell_id}")
    )
  output:
    tuple(
        val(cell_id), val(lanes), val(flowcells),
        path(metadata),
        path("fastqscreen/${cell_id}"),
        path("cleanalignsort/${cell_id}"),
        path("aligned.bam"),
        path("aligned.bam.bai")
    )




  script:
    def lanes = lanes.join(' ')
    def flowcells = flowcells.join(' ')
    

    // alignment_cleanalignsort(fastq_pairs, metadata_yaml, reference, tempdir, adapter1, adapter2, cell_id, bam_output, num_threads, run_fastqc=False):

    """

        fastqs_cmd=`python -c 'x=["${lanes}","${flowcells}"];x=[v.split() for v in x];x=[",".join(v) for v in zip(*x)];x=" --fastq_pairs ".join(x);print(x)'`

        alignment_utils alignmentcleanalignsort \
        --fastq_pairs \${fastqs_cmd} \
        --metadata_yaml ${metadata} \
        --reference ${primary_reference_name},${primary_reference_version},${primary_reference} \
        --tempdir tempdir \
        --adapter1 CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
        --adapter2 CTGTCTCTTATACACATCTGACGCTGCCGACGA \
        --cell_id $cell_id \
        --bam_output aligned.bam \
        --num_threads ${task.cpus}


        mkdir -p cleanalignsort/${cell_id}
        cp -r tempdir/${cell_id}/* cleanalignsort/${cell_id}/

        rm -rf tempdir
    """
}
