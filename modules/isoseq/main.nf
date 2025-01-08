process lima {
    /* primer removal and demultiplex primers */
    tag "$sample_id"
    
    input:

    tuple val(sample_id), path(bam), path(pbi_file)
    path(primer)
    
    output:

    path "${sample_id}.output.*.bam", emit: lima_bam
    path "${sample_id}.output.*.bam.pbi", emit: lima_bam_pbi
    path "${sample_id}.output.json", emit: lima_json
    path "${sample_id}.output.lima.clips", emit: lima_clips
    path "${sample_id}.output.lima.counts", emit: lima_counts
    path "${sample_id}.output.lima.guess", emit: lima_guess
    path "${sample_id}.output.lima.report", emit: lima_report
    path "${sample_id}.output.lima.summary", emit: lima_summary
    
    script:
    """
    lima --version
    lima $bam $primer ${sample_id}.output.bam --isoseq --peek-guess
    """
}

process refine {

    /* generate full length non-concatemer reads (FLNC) by trimming polyA tails remove concatmers */

    publishDir "${params.output_dir}", mode: 'copy'
    tag "$sample_id"

    input:
    val sample_id
    path lima_bam 
    path lima_pbi
    path primer

    output:
    tuple path("${sample_id}.flnc.bam"), path("${sample_id}.flnc.bam.pbi"), emit: flnc_bam
    path "${sample_id}.flnc.*.xml", emit: flnc_xml
    path "${sample_id}.*.report.json", emit: flnc_json
    path "${sample_id}.*.report.csv", emit:  flnc_report

    script:

    """
    isoseq refine --version
    isoseq refine $lima_bam $primer ${sample_id}.flnc.bam  --require-polya
    """
}

process cluster2 {

    /* clusters reads together - output is unaligned BAM of clustered reads */

    publishDir "${params.output_dir}", mode: 'copy'
    tag "$sample_id"

    input:
    val sample_id
    tuple  path(flnc_bam ), path(flnc_pbi)
    

    output:

    tuple path("${sample_id}.transcript.bam"), path("${sample_id}.transcript.bam.pbi"), emit: transcript_bam
    path  "${sample_id}.transcript.cluster_report.csv", emit: transcript_report

    script:
    """
    
    isoseq cluster2 --version

    isoseq cluster2 $flnc_bam ${sample_id}.transcript.bam
    
    """
  

}


process pbmm2_align {

    /* read alignment to reference of un-aligned, clustered reads */
    label 'high_memory'
    publishDir "${params.output_dir}", mode: 'copy'
    tag "$sample_id"

    input:
    path reference
    val(sample_id)
    tuple path(bam), path(pbi)
    path (mmi)
    val threads
    val sort_threads

    output:
    tuple path("${sample_id}.pbmm2.mapped.bam"), path("${sample_id}.pbmm2.mapped.bam.bai"), emit: mapped_bam
    //path "${$sample_id}.read_length_and_quality.tsv", emit: bam_rl_qual
    
    script:
    """
    pbmm2 --version
    pbmm2 align \\
        --sort \\
        -j $threads \\
        -J $sort_threads \\
        --preset ISOSEQ \\
        --sample ${sample_id} \\
        --log-level INFO \\
        --unmapped \\
        --bam-index BAI \\
        $reference \\
        $bam \\
        ${sample_id}.pbmm2.mapped.bam

    
    """

}

process collapse {

    /* after reads are aligned, collapse them into unique isoforms */
    publishDir "${params.output_dir}", mode: 'copy'
    tag "$sample_id"

    input:
    val (sample_id)
    tuple path(mapped_bam), path(mapped_bai)
    tuple  path(flnc_bam ), path(flnc_pbi)

    output:
    path("${sample_id}.collapsed.gff"), emit: collapse_gff
    path("${sample_id}.collapsed.flnc_count.txt"), emit: flnc_count
    path ("${sample_id}.collapsed.group.txt"), emit: collapsed_group
    path ("${sample_id}.collapsed.fasta"), emit: collapsed_fasta


    script:
    
    """
    isoseq collapse --do-not-collapse-extra-5exons $mapped_bam $flnc_bam ${sample_id}.collapsed.gff   
    """


}

process pigeon_prepare {

    /* sort and index GFF file of collapsed reads */

    publishDir "${params.output_dir}", mode: 'copy'
    tag "$sample_id"

    input:
    val (sample_id)
    path (gff)

    output:
     path("${sample_id}.collapsed.sorted.gff"), emit: sorted_gff
     path("${sample_id}.collapsed.sorted.gff.pgi"), emit: sorted_gff_idx

    script:

    """
    pigeon --version
    pigeon prepare  $gff
    
    """

}

process pigeon_classify {

    /* pigeon classify annotates isoforms */

    publishDir "${params.output_dir}", mode: 'copy'
    tag "$sample_id"

    input:
    val (sample_id)
    path (sorted_gff)
    path (sorted_gff_idx)
    path (gtf_annotations)
    path (gtf_annotations_idx)
    path (reference)
    path (reference_idx)
    path (tss_bed)
    path (tss_bed_idx)
    path (polyA)
    path (flnc_counts)


    output:

    path ("${sample_id}.report.json"), emit: classify_report_json
    path ("${sample_id}.summary.txt"), emit: classify_summary
    path ("${sample_id}_classification.txt"), emit: classification_result
    path ("${sample_id}_junctions.txt"), emit: classify_junctions

    script:
    """
    pigeon --version
    pigeon classify $sorted_gff $gtf_annotations $reference --cage-peak $tss_bed --poly-a $polyA --fl $flnc_counts
    """

}




process pigeon_filter {

    /* pigeon filters out isoforms */

    publishDir "${params.output_dir}", mode: 'copy'
    tag "$sample_id"


    input:
    val(sample_id)
    path(pigeon_classification)
    path (pigeon_junctions)
    path (sorted_gff)

    output:

    path("${sample_id}_classification.filtered_lite_classification.txt"), emit: pigeon_filter_lite_classification
    path("${sample_id}_classification.filtered_lite_junctions.txt"), emit: pigeon_filter_lite_junctions
    path("${sample_id}_classification.filtered_lite_reasons.txt"), emit: pigeon_filer_lite_reasons
    path("${sample_id}_classification.filtered.summary.txt"), emit: pigeon_filter_summary
    path("${sample_id}_classification.filtered.report.json"), emit: pigeon_filter_report_json
    
    script:

    """
    pigeon --version

    pigeon filter $pigeon_classification --isoforms $sorted_gff
    
    """


}

