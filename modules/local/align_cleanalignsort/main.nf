process ALIGNCLEANALIGNSORT {
    time '24h'
    cpus 8
    memory '8 GB'
    label 'process_high'

  input:
    tuple(
      val(cell_id), val(lanes), val(flowcells), path(fastqs1), path(fastqs2),
      path(primary_reference), val(primary_reference_version), val(primary_reference_name),
      path(primary_reference_fai), path(primary_reference_amb),path(primary_reference_ann),
      path(primary_reference_bwt),path(primary_reference_pac),path(primary_reference_sa),
      path(secondary_reference_1), val(secondary_reference_1_version), val(secondary_reference_1_name),
      path(secondary_reference_1_fai), path(secondary_reference_1_amb),path(secondary_reference_1_ann),
      path(secondary_reference_1_bwt),path(secondary_reference_1_pac),path(secondary_reference_1_sa),
      path(secondary_reference_2), val(secondary_reference_2_version), val(secondary_reference_2_name),
      path(secondary_reference_2_fai), path(secondary_reference_2_amb),path(secondary_reference_2_ann),
      path(secondary_reference_2_bwt),path(secondary_reference_2_pac),path(secondary_reference_2_sa),
      path(secondary_reference_3), val(secondary_reference_3_version), val(secondary_reference_3_name),
      path(secondary_reference_3_fai), path(secondary_reference_3_amb),path(secondary_reference_3_ann),
      path(secondary_reference_3_bwt),path(secondary_reference_3_pac),path(secondary_reference_3_sa),
      path(secondary_reference_4), val(secondary_reference_4_version), val(secondary_reference_4_name),
      path(secondary_reference_4_fai), path(secondary_reference_4_amb),path(secondary_reference_4_ann),
      path(secondary_reference_4_bwt),path(secondary_reference_4_pac),path(secondary_reference_4_sa),
      path(secondary_reference_5), val(secondary_reference_5_version), val(secondary_reference_5_name),
      path(secondary_reference_5_fai), path(secondary_reference_5_amb),path(secondary_reference_5_ann),
      path(secondary_reference_5_bwt),path(secondary_reference_5_pac),path(secondary_reference_5_sa),
      path(secondary_reference_6), val(secondary_reference_6_version), val(secondary_reference_6_name),
      path(secondary_reference_6_fai), path(secondary_reference_6_amb),path(secondary_reference_6_ann),
      path(secondary_reference_6_bwt),path(secondary_reference_6_pac),path(secondary_reference_6_sa),
      path(metadata),
      path("fastqscreen/${cell_id}")
    )
  output:
    tuple(
        val(cell_id), val(lanes), val(flowcells), path(fastqs1), path(fastqs2),
        path(primary_reference), val(primary_reference_version), val(primary_reference_name),
        path(primary_reference_fai), path(primary_reference_amb),path(primary_reference_ann),
        path(primary_reference_bwt),path(primary_reference_pac),path(primary_reference_sa),
        path(secondary_reference_1), val(secondary_reference_1_version), val(secondary_reference_1_name),
        path(secondary_reference_1_fai), path(secondary_reference_1_amb),path(secondary_reference_1_ann),
        path(secondary_reference_1_bwt),path(secondary_reference_1_pac),path(secondary_reference_1_sa),
        path(secondary_reference_2), val(secondary_reference_2_version), val(secondary_reference_2_name),
        path(secondary_reference_2_fai), path(secondary_reference_2_amb),path(secondary_reference_2_ann),
        path(secondary_reference_2_bwt),path(secondary_reference_2_pac),path(secondary_reference_2_sa),
        path(secondary_reference_3), val(secondary_reference_3_version), val(secondary_reference_3_name),
        path(secondary_reference_3_fai), path(secondary_reference_3_amb),path(secondary_reference_3_ann),
        path(secondary_reference_3_bwt),path(secondary_reference_3_pac),path(secondary_reference_3_sa),
        path(secondary_reference_4), val(secondary_reference_4_version), val(secondary_reference_4_name),
        path(secondary_reference_4_fai), path(secondary_reference_4_amb),path(secondary_reference_4_ann),
        path(secondary_reference_4_bwt),path(secondary_reference_4_pac),path(secondary_reference_4_sa),
        path(secondary_reference_5), val(secondary_reference_5_version), val(secondary_reference_5_name),
        path(secondary_reference_5_fai), path(secondary_reference_5_amb),path(secondary_reference_5_ann),
        path(secondary_reference_5_bwt),path(secondary_reference_5_pac),path(secondary_reference_5_sa),
        path(secondary_reference_6), val(secondary_reference_6_version), val(secondary_reference_6_name),
        path(secondary_reference_6_fai), path(secondary_reference_6_amb),path(secondary_reference_6_ann),
        path(secondary_reference_6_bwt),path(secondary_reference_6_pac),path(secondary_reference_6_sa),
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

        fastqs_cmd=`python -c 'x=["${lanes}","${flowcells}","${fastqs1}","${fastqs2}"];x=[v.split() for v in x];x=[",".join(v) for v in zip(*x)];x=" --fastq_pairs ".join(x);print(x)'`

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
