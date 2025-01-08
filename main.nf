#!/usr/local/bin/nextflow

nextflow.enable.dsl=2

include {lima; refine; cluster2; pbmm2_align; collapse; pigeon_prepare; pigeon_classify; pigeon_filter} from './modules/isoseq'

def required_parms = ['samplesheet_dir', 'samplesheet', 'output_dir']
required_parms.each { param ->
    if (!params.containsKey(param)) {
        error "Parameter '$param' is required!"
    }
}

Channel.fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> 
        def sample_id = row.sample_id
        def bam_file = file(row.bam_file)
        def pbi_file = file(row.pbi_file)
        if (!bam_file.exists()) {
            error "BAM file not found: ${bam_file}"
        }
        return tuple(sample_id, bam_file, pbi_file)
    }
    .set { input_bams_ch }

    input_bams_ch
    .map { sample_id, bam_file, pbi_file ->
        sample_id
    }
    .set { sample_id_ch }

    workflow { 

    /* lima */

    lima( input_bams_ch, params.primer_file)

    lima.out.lima_bam
        .flatten()
        .set{individual_bams}
    lima.out.lima_bam_pbi
        .flatten()
        .set{individual_pbis}
    
    refine( sample_id_ch, individual_bams, individual_pbis, params.primer_file)

    cluster2(sample_id_ch, refine.out.flnc_bam)

    pbmm2_align( params.reference, sample_id_ch, cluster2.out.transcript_bam, params.mmi_index, params.cpu, params.sort_threads )

    collapse( sample_id_ch, pbmm2_align.out.mapped_bam, refine.out.flnc_bam)

    pigeon_prepare( sample_id_ch, collapse.out.collapse_gff )

    pigeon_classify (sample_id_ch, 
                    pigeon_prepare.out.sorted_gff,
                    pigeon_prepare.out.sorted_gff_idx,
                    params.pigeon_annotations_gtf,
                    params.pigeon_annotations_gtf_idx,
                    params.reference, 
                    params.reference_index,
                    params.tss_file,
                    params.tss_file_idx,
                    params.polyA_file,
                    collapse.out.flnc_count
                )

        pigeon_filter ( sample_id_ch, pigeon_classify.out.classification_result, pigeon_classify.out.classify_junctions, pigeon_prepare.out.sorted_gff )

    
    }

