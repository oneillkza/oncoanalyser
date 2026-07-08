//
// Align RNA reads from an existing BAM/CRAM by streaming through samtools fastq into STAR
//
// NOTE(KO): This subworkflow mirrors READ_ALIGNMENT_RNA but sources reads from an existing BAM/CRAM
// (collated to disk, then fed to STAR via samtools fastq process substitutions) instead of FASTQ.
// It intentionally emits rna_tumor as [ meta, bam, bai ] — identical to READ_ALIGNMENT_RNA — so that
// downstream stages (isofox, lilac, cider, neo) consume its output unchanged.
//
// Unlike the FASTQ path, there is only ever a single input BAM per sample, so no Sambamba merge step
// is needed.
//

include { GATK4_MARKDUPLICATES } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { SAMTOOLS_SORT         } from '../../../modules/nf-core/samtools/sort/main'
include { STAR_ALIGN_FROM_BAM   } from '../../../modules/local/star/align_from_bam/main'

workflow READ_ALIGNMENT_RNA_FROM_BAM {
    take:
    // Sample data
    ch_inputs         // channel: [mandatory] [ meta ]

    // Reference data
    genome_star_index // channel: [mandatory] /path/to/genome_star_index/

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Sort inputs
    // channel: runnable: [ meta ]
    // channel: skip: [ meta ]
    ch_inputs_sorted = ch_inputs
        .branch { meta ->
            runnable: sample.Inputs.hasTumorRnaBam(meta)
            skip: true
        }

    // Create BAM/CRAM input channel for STAR_ALIGN_FROM_BAM
    // channel: [ meta_star, bam ]
    ch_star_inputs = ch_inputs_sorted.runnable
        .map { meta ->
            def meta_sample = sample.Inputs.getTumorRnaSample(meta)
            def meta_star = [
                key:        meta.group_id,
                id:         "${meta.group_id}_${meta_sample.sample_id}",
                sample_id:  meta_sample.sample_id,
                read_group: "${meta_sample.sample_id}.realign",
            ]
            return [meta_star, sample.Inputs.get(meta, sample.FileKey.BAM_RNA_TUMOR)]
        }

    // Run STAR alignment
    STAR_ALIGN_FROM_BAM(
        ch_star_inputs,
        genome_star_index,
    )

    ch_versions = ch_versions.mix(STAR_ALIGN_FROM_BAM.out.versions)

    // Sort the unsorted BAM output from STAR
    // channel: [ meta_sort, bam ]
    ch_sort_inputs = STAR_ALIGN_FROM_BAM.out.bam
        .map { meta_star, bam ->
            [[ *:meta_star, prefix: meta_star.read_group ], bam]
        }

    SAMTOOLS_SORT(ch_sort_inputs)

    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    // Markduplicates
    // Single BAM per sample — no merge step required (unlike the multi-lane FASTQ path)
    // channel: [ meta_markdups, bam ]
    ch_markdups_inputs = channels.WorkflowChannels.restoreMeta(SAMTOOLS_SORT.out.bam, ch_inputs)
        .map { meta, bam ->
            def meta_markdups = [
                key:       meta.group_id,
                id:        meta.group_id,
                sample_id: sample.Inputs.getTumorRnaSampleName(meta),
            ]
            return [meta_markdups, bam]
        }

    GATK4_MARKDUPLICATES(ch_markdups_inputs, [], [])

    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)

    // Combine BAMs and BAIs — identical to READ_ALIGNMENT_RNA
    // channel: [ meta, bam, bai ]
    ch_bams_ready = channels.WorkflowChannels.groupByMeta(
        channels.WorkflowChannels.restoreMeta(GATK4_MARKDUPLICATES.out.bam, ch_inputs),
        channels.WorkflowChannels.restoreMeta(GATK4_MARKDUPLICATES.out.bai, ch_inputs),
    )

    // Set outputs
    // channel: [ meta, bam, bai ]
    ch_bam_out = Channel.empty()
        .mix(
            ch_bams_ready,
            channels.PlaceholderChannels.bamBai(ch_inputs_sorted.skip),
        )

    emit:
    rna_tumor = ch_bam_out   // channel: [ meta, bam, bai ]

    versions  = ch_versions  // channel: [ versions.yml ]
}
