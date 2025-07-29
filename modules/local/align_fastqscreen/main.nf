process ALIGNFASTQSCREEN {
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
      path(metadata)
    )
  output:
    tuple(
      val(cell_id), val(lanes), val(flowcells), path(fastqs1), path(fastqs2),
      path(primary_reference), val(primary_reference_version), val(primary_reference_name),
      path(primary_reference_fai), path(primary_reference_amb),path(primary_reference_ann),
      path(primary_reference_bwt),path(primary_reference_pac),path(primary_reference_sa),
      path(metadata),
      path("fastqscreen/${cell_id}")
    )

  script:
    def lanes = lanes.join(' ')
    def flowcells = flowcells.join(' ')
    def supplementary_2 = ''
    if(secondary_reference_2_name) {
        supplementary_2 = '--supplementary_references ' + secondary_reference_2_name + ',' + secondary_reference_2_version + ',' + secondary_reference_2
    } else {
        supplementary_2 = ''
    }
    def supplementary_3 = ''
    if(secondary_reference_3_name) {
        supplementary_3 = '--supplementary_references ' + secondary_reference_3_name + ',' + secondary_reference_3_version + ',' + secondary_reference_3
    } else {
        supplementary_3 = ''
    }
    def supplementary_4 = ''
    if(secondary_reference_4_name) {
        supplementary_4 = '--supplementary_references ' + secondary_reference_4_name + ',' + secondary_reference_4_version + ',' + secondary_reference_4
    } else {
        supplementary_4 = ''
    }
    def supplementary_5 = ''
    if(secondary_reference_5_name) {
        supplementary_5 = '--supplementary_references ' + secondary_reference_5_name + ',' + secondary_reference_5_version + ',' + secondary_reference_5
    } else {
        supplementary_5 = ''
    }
    def supplementary_6 = ''
    if(secondary_reference_6_name) {
        supplementary_6 = '--supplementary_references ' + secondary_reference_6_name + ',' + secondary_reference_6_version + ',' + secondary_reference_6
    } else {
        supplementary_6 = ''
    }


    //alignment_fastqscreen(fastq_pairs, reference, supplementary_references, tempdir, cell_id, num_threads):


    """

        fastqs_cmd=`python -c 'x=["${lanes}","${flowcells}","${fastqs1}","${fastqs2}"];x=[v.split() for v in x];x=[",".join(v) for v in zip(*x)];x=" --fastq_pairs ".join(x);print(x)'`

        alignment_utils alignmentfastqscreen \
        --fastq_pairs \${fastqs_cmd} \
        --reference ${primary_reference_name},${primary_reference_version},${primary_reference} \
        --supplementary_references ${secondary_reference_1_name},${secondary_reference_1_version},${secondary_reference_1} \
        ${supplementary_2} \
        ${supplementary_3} \
        ${supplementary_4} \
        ${supplementary_5} \
        ${supplementary_6} \
        --tempdir tempdir \
        --cell_id $cell_id \
        --num_threads ${task.cpus}

        mkdir -p fastqscreen/${cell_id}
        cp -r tempdir/${cell_id}/* fastqscreen/${cell_id}/

        rm -rf tempdir
    """
}
