//
// Align RNA reads from existing alignments, realigning via STAR
//
// NOTE(KO): This subworkflow mirrors READ_ALIGNMENT_RNA but sources reads from an existing BAM/CRAM
// rather than FASTQ, emitting the same channel shapes so downstream stages consume it unchanged.
// There is exactly one alignment per sample, so no Sambamba merge step is required.
//

include { GATK4_MARKDUPLICATES } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { SAMTOOLS_SORT        } from '../../../modules/nf-core/samtools/sort/main'
include { STAR_ALIGN_FROM_BAM  } from '../../../modules/local/star/align_from_bam/main'

workflow READ_ALIGNMENT_RNA_FROM_BAM {
    take:
    // Sample data
    ch_inputs         // channel: [mandatory] [ meta ]

    // Reference data
    genome_star_index // channel: [mandatory] /path/to/genome_star_index/

    main:
    //
    // STEP: Handle inputs
    //
    // Sort inputs
    // runnable: channel: [ meta ]
    // skip: channel: [ meta ]
    ch_inputs_sorted = ch_inputs
        .branch { meta ->
            runnable: Utils.hasTumorRnaBam(meta)
            skip: true
        }

    //
    // MODULE: STAR alignment
    //
    // Create process input channel
    // channel: [ meta_star, [ rg_line ], aln, idx ]
    ch_star_inputs = ch_inputs_sorted.runnable
        .map { meta ->

            def meta_sample = Utils.getTumorRnaSample(meta)
            def sample_id = Utils.getTumorRnaSampleName(meta)
            def library_id = meta_sample.library_id

            // NOTE(KO): as for the DNA path, a single synthetic read group is constructed per sample since
            // lane and flowcell are not available for alignment inputs. STAR expects each field quoted and
            // space delimited rather than the tab delimited form used by bwa-mem2.
            def rg_id = "${sample_id}.${library_id}.realign"
            def rg_entries = [ID: rg_id, SM: sample_id, LB: library_id]
            def rg_line = rg_entries.collect { k, v -> "'${k}:${v}'" }.join(' ')

            def meta_star = [
                key: meta.group_id,
                id: "${meta.group_id}_${sample_id}",
                sample_id: sample_id,
                library_id: library_id,
            ]

            return [meta_star, [rg_line], Utils.getTumorRnaBam(meta), Utils.getTumorRnaBai(meta)]
        }

    // Run process
    STAR_ALIGN_FROM_BAM(
        ch_star_inputs,
        genome_star_index,
    )

    //
    // MODULE: SAMtools sort
    //
    // Create process input channel
    // channel: [ meta_sort, aln ]
    ch_sort_inputs = channel.topic('star_align_bam')
        .map { meta_star, aln ->
            def meta_sort = meta_star + [prefix: meta_star.sample_id]
            return [meta_sort, aln]
        }

    // Run process
    SAMTOOLS_SORT(
        ch_sort_inputs,
    )

    //
    // MODULE: GATK4 markduplicates
    //
    // Create process input channel
    // channel: [ meta_markdups, aln ]
    ch_markdups_inputs = WorkflowOncoanalyser.restoreMeta(channel.topic('samtools_sort_bam'), ch_inputs)
        .map { meta, aln ->
            def meta_markdups = [
                key: meta.group_id,
                id: meta.group_id,
                sample_id: Utils.getTumorRnaSampleName(meta),
            ]
            return [meta_markdups, aln]
        }

    // Run process
    GATK4_MARKDUPLICATES(
        ch_markdups_inputs,
        [],
        [],
    )

    //
    // STEP: Handle outputs
    //
    // Combine BAMs and BAIs
    // channel: [ meta, aln, idx ]
    ch_alns_ready = WorkflowOncoanalyser.groupByMeta(
        WorkflowOncoanalyser.restoreMeta(channel.topic('gatk4_markduplicates_bam'), ch_inputs),
        WorkflowOncoanalyser.restoreMeta(channel.topic('gatk4_markduplicates_bai'), ch_inputs),
    )

    // Combine STAR log with QC and MarkDuplicates metrics
    // channel: [ meta, star_log, md_metrics ]
    ch_qc_files_ready = WorkflowOncoanalyser.groupByMeta(
        WorkflowOncoanalyser.restoreMeta(channel.topic('star_align_qc_log'), ch_inputs),
        WorkflowOncoanalyser.restoreMeta(channel.topic('gatk4_markduplicates_metrics'), ch_inputs),
    )

    // Set outputs
    // channel: [ meta, aln, idx ]
    ch_outputs_aln = channel.empty()
        .mix(
            ch_alns_ready,
            ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
        )

    // channel: [ meta, star_log, md_metrics ]
    ch_outputs_qc_files = channel.empty()
        .mix(
            ch_qc_files_ready,
            ch_inputs_sorted.skip.map { meta -> [meta, [], []] },
        )

    emit:
    tumor    = ch_outputs_aln      // channel: [ meta, aln, idx ]
    qc_files = ch_outputs_qc_files // channel: [ meta, star_log, md_metrics ]
}
