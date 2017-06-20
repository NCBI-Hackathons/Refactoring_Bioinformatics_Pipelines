
input_paired_1 = params.input_paired_1
input_paired_2 = params.input_paired_2
cuffmerge_wrapper = params.cm_wrapper
annotationFile = file(params.input_gtf)
bam2sam_out1 = "${params.output_dir}/bam2sam_out/1.sam"
bam2sam_out2 = "${params.output_dir}/bam2sam_out/2.sam"

mate_std_dev = params.mate_std_dev
anchor_length = params.anchor_length
segment_length = params.segment_length

index_file = file(params.index)
index_name = index_file.getFileName()
index_dir = index_file.getParent()
genomeFile = file("${index_dir}/hg38.fa")

log.info "\n\n\t\t\t =====================================================\n\t\t\t ==  R N A - S e q   T U X E D O   P I P E L I N E  ==\n\t\t\t =====================================================\n\n"
log.info "\t===================================== I N F O M A T I O N ===================================="
log.info "\tInput pair 1         : ${file(input_paired_1)}"
log.info "\tInput pair 2         : ${file(input_paired_2)}"
log.info "\tAnnotation           : ${params.input_gtf}"
log.info "\tGenome file          : ${genomeFile}"
log.info "\tCm_wrapper           : ${params.cm_wrapper}"
log.info "\tOutput path          : ${params.output_dir}"
log.info "\tIndex_name           : ${index_name}"
log.info "\tIndex_dir            : ${index_dir}"
log.info "\tIndex_prefix         : ${index_file}"
log.info "\tMate inner distance  : ${mate_std_dev}"
log.info "\tAnchor length        : ${anchor_length}"
log.info "\tRead segments length : ${segment_length}"
log.info "\t==============================================================================================\n\n"

Channel
    .fromPath(genomeFile)
    .into {genomes}

Channel
    .fromPath(annotationFile)
    .into {annotation_1; annotation_2; annotation_3; annotation_4}

Channel
    .fromPath(params.cm_wrapper)
    .into {cuffmerge_wrapper}

read_paired_1 = Channel
        .fromPath(params.input_paired_1)
        .map{file(it).toAbsolutePath()}

read_paired_2 = Channel
        .fromPath(params.input_paired_2)
        .map{file(it).toAbsolutePath()}


/*
process Index_setting {
    input:
        file(index_dir)
        val(index_name)

    output:
        set val(index_name), file("genomeIndex") into genome_index

    script:
    """
    mkdir -p genomeIndex
    cp ${index_dir}/${index_name}.* genomeIndex/.
    """
}
*/

process Get_userName {
    executor 'local'

    output:
        stdout user_name
        'printf $(id -u -n)'
}

process Get_userGroup {
    executor 'local'
    output:
        stdout user_group
        'printf $(id -g -n)'
}

user_name.println()
user_group.println()

process Index_setting {
    executor 'local'

    input:
        val(u_name) from user_name
        val(u_group) from user_group
        file(index_dir)
        val(index_name)

    output:
        file("genomeIndex/*") into genome_index

    script:
    """
    mkdir -p genomeIndex
    cp ${index_dir}/${index_name}.* genomeIndex/
    chown -R ${u_name}:${u_group} genomeIndex/*
    """
}

process Input_setting_1 {
    executor 'local'

    input:
        val(u_name) from user_name
        val(u_group) from user_group
        val (input_paired_1)
        val (reads) from read_paired_1.collect()

    output:
        file("input_1/*") into read_paired_1_1

    script:
    input_path = reads[0].getParent()
    """
    mkdir -p input_1
    cp $input_path/* input_1/
    chown -R ${u_name}:${u_group} input_1/*
    """
}

process Input_setting_2 {
    executor 'local'

    input:
        val(u_name) from user_name
        val(u_group) from user_group
        val (input_paired_2)
        val (reads) from read_paired_2.collect()

    output:
        file("input_2/*") into read_paired_2_1

    script:
    input_path = reads[0].getParent()
    """
    mkdir -p input_2
    cp $input_path/* input_2/
    chown -R ${u_name}:${u_group} input_2/*
    """
}

/*
process Show_genomeIndex{
    input:
        val(reads) from read_paired_1_2.collect()
        val(index_list) from genome_index
    script:
        index_path = index_list[0].getParent()
        """
        echo "tophat2 -p 8 -r ${mate_std_dev} -a ${anchor_length} --segment-length=${segment_length} $index_path/${index_name} ${reads[0]} ${reads[1]}"
        """
}
*/

process Tophat_1 {
        input:
            val (index_name)
            set val (index_name), file(index_dir) from genome_index.first()
            val (reads) from read_paired_1_1.collect()

        output:
            file "tophat_out/accepted_hits_1.bam" into tophat_out_1_1, tophat_out_1_2

        script:
        indexPath = index_name.getParent()
        """
        mkdir -p ${baseDir}/output/tophat
        tophat2 -p 8 -r ${mate_std_dev} -a ${anchor_length} --segment-length=${segment_length} $indexPath/hg38 ${reads[0]} ${reads[1]}
        cp tophat_out/align_summary.txt ${baseDir}/output/tophat/align_summary_1.txt
        mv tophat_out/accepted_hits.bam tophat_out/accepted_hits_1.bam
        """
}

process Tophat_2 {
        input:
            val (index_name)
            set val (index_name), file(index_dir) from genome_index.first()
            val (reads) from read_paired_2_1.collect()

        output:
            file "tophat_out/accepted_hits_2.bam" into tophat_out_2_1, tophat_out_2_2

        script:
        indexPath = index_name.getParent()
        """
        mkdir -p ${baseDir}/output/tophat
        tophat2 -p 8 -r ${mate_std_dev} -a ${anchor_length} --segment-length=${segment_length} $indexPath/hg38 ${reads[0]} ${reads[1]}
        cp tophat_out/align_summary.txt ${baseDir}/output/tophat/align_summary_2.txt
        mv tophat_out/accepted_hits.bam tophat_out/accepted_hits_2.bam
        """
}

/*
process Tophat_1 {
        input:
            val (index_name)
            val (index_dir)
            set val(name), file(reads) from read_paired_1

        output:
            file "tophat_out/accepted_hits_1.bam" into tophat_out_1_1, tophat_out_1_2

        script:
        """
        mkdir -p ${baseDir}/output/tophat
        tophat2 -p 8 -r ${mate_std_dev} -a ${anchor_length} --segment-length=${segment_length} ${index_dir}/${index_base} ${reads}
        cp tophat_out/align_summary.txt ${baseDir}/output/tophat/align_summary_1.txt
        mv tophat_out/accepted_hits.bam tophat_out/accepted_hits_1.bam
        """
}

*/
/*
process Tophat_1 {
        input:
            val (index_name)
            file(index_file) from genome_index
            set val(name), file(reads) from read_paired_1

        output:
            file "tophat_out/accepted_hits_1.bam" into tophat_out_1_1, tophat_out_1_2

        script:
        index_base = "${index_file[0]}.getParent()/${index_name}"
        """
        echo "$index_base"
        mkdir -p ${baseDir}/output/tophat
        tophat2 -p 8 -r ${mate_std_dev} -a ${anchor_length} --segment-length=${segment_length} $index_base ${reads}
        cp tophat_out/align_summary.txt ${baseDir}/output/tophat/align_summary_1.txt
        mv tophat_out/accepted_hits.bam tophat_out/accepted_hits_1.bam
        """
}


process Tophat_2 {
        input:
            file bw2_indices from genome_index.collect()
            set val(name), file(reads) from read_paired_2

        output:
            file "tophat_out/accepted_hits_2.bam" into tophat_out_2_1, tophat_out_2_2

        script:
        index_base = bt2_indices[0].toString() - ~/.\d.bt2/
        """
        mkdir -p ${baseDir}/output/tophat
        tophat2 -p 8 -r 300 -r ${mate_std_dev} -a ${anchor_length} --segment-length=${segment_length} $index_base ${reads}
        cp tophat_out/align_summary.txt ${baseDir}/output/tophat/align_summary_2.txt
        mv tophat_out/accepted_hits.bam tophat_out/accepted_hits_2.bam
        """
}

*/


process BamToSam_1 {
    publishDir "$baseDir/output/bam2sam", mode: 'copy', overwrite: false

    input:
        file tophat_output1 from tophat_out_1_1

    output:
        file "BamToSam_1.sam" into bam2sam_out_1

    """
    mkdir -p ${baseDir}/output/bam2sam
    samtools view -o BamToSam_1.sam -h ${tophat_output1}
    """
}

process BamToSam_2 {
    publishDir "$baseDir/output/bam2sam", mode: 'copy', overwrite: false

    input:
        file tophat_output2 from tophat_out_2_1

    output:
        file "BamToSam_2.sam" into bam2sam_out_2

    """
    mkdir -p ${baseDir}/output/bam2sam
    samtools view -o BamToSam_2.sam -h ${tophat_output2}
    """
}

process Cufflinks_1 {
    publishDir "$baseDir/output/cufflinks", mode: 'copy', overwrite: false

    input:
        file tophat_output1 from tophat_out_1_2
        file anno_file from annotation_1

    output:
        file "transcripts_1.gtf" into cufflinks_gtf_out_1
    """
    mkdir -p ${baseDir}/output/cufflinks
    cufflinks -q --no-update-check -p 8 -I 300000 -F 0.1 -j 0.15 -G ${anno_file} ${tophat_output1}
    mv transcripts.gtf transcripts_1.gtf
    """
}

process Cufflinks_2 {
    publishDir "$baseDir/output/cufflinks", mode: 'copy', overwrite: false

    input:
        file tophat_output2 from tophat_out_2_2
        file anno_file from annotation_2

    output:
        file "transcripts_2.gtf" into cufflinks_gtf_out_2
    """
    mkdir -p ${baseDir}/output/cufflinks
    cufflinks -q --no-update-check -p 8 -I 300000 -F 0.1 -j 0.15 -G ${anno_file} ${tophat_output2}
    mv transcripts.gtf transcripts_2.gtf
    """
}

process Cuffmerge {
    publishDir "$baseDir/output/cuffmerge", mode: 'copy', overwrite: false

    input:
        file cufflinks_gtf_1 from cufflinks_gtf_out_1
        file cufflinks_gtf_2 from cufflinks_gtf_out_2
        file anno_file from annotation_3
        file cuffmergeWrapper from cuffmerge_wrapper

    output:
        file "merged_transcripts.gtf" into cuffmerge_gtf_out

    """
    mkdir -p ${baseDir}/output/cuffmerge
    python ${cuffmergeWrapper} -p 8 -g ${anno_file} --min-isoform-fraction="0.05" --merged-transcripts="merged_transcripts.gtf" ${cufflinks_gtf_1} ${cufflinks_gtf_2}
    """
}

process Cuffcompare {
    publishDir "$baseDir/output/cuffcompare", mode: 'copy', overwrite: false

    input:
        file genome from genomes
        file anno_file from annotation_4
        file cuffmerge_gtf from cuffmerge_gtf_out

    output:
        file "cuffcmp.combined.gtf" into cuffcompare_out

    """
    mkdir -p ${baseDir}/output/cuffcompare
    cuffcompare -r ${anno_file} -s ${genome} -e 100 -d 100 ${cuffmerge_gtf}
    """
}

process Cuffdiff {
    publishDir "$baseDir/output/cuffdiff", mode: 'copy', overwrite: false

    input:
        file cuffcompare_gtf from cuffcompare_out
        file bam2sam_sam_1 from bam2sam_out_1
        file bam2sam_sam_2 from bam2sam_out_2

    output:
        file "*.diff" into cuffdiff_out
    """
    mkdir -p ${baseDir}/output/cuffdiff
    cuffdiff --no-update-check --FDR=0.05 -p 8 --min-alignment-count=10 --library-norm-method=geometric --dispersion-method=pooled --labels "," ${cuffcompare_gtf} ${bam2sam_sam_1} ${bam2sam_sam_2}
    """
}

workflow.onComplete {
    println(workflow.success ? "\n\t\t  Done!\n\n": "\n\t\t  Oops .. something went wrong\n\n")
}
