//
// Align DNA reads from an existing BAM/CRAM by streaming through samtools fastq
//
// NOTE(KO): This subworkflow mirrors READ_ALIGNMENT_DNA but sources reads from an existing BAM/CRAM
// (streamed via samtools fastq) instead of FASTQ. It intentionally emits the same channel shapes as
// READ_ALIGNMENT_DNA (dna_tumor/dna_normal/dna_donor as [ meta, [bam, ...], [bai, ...] ]) so that the
// downstream REDUX and analysis stages consume its outputs unchanged.
//
// Unlike the FASTQ path, there is no option to chunk up the BAM (i.e. no max_fastq_records splitting).
//

include { BWAMEM2_ALIGN_FROM_BAM } from '../../../modules/local/bwa-mem2/mem_from_bam/main'

workflow READ_ALIGNMENT_DNA_FROM_BAM {
    take:
    // Sample data
    ch_inputs            // channel: [mandatory] [ meta ]

    // Reference data
    genome_fasta         // channel: [mandatory] /path/to/genome_fasta
    genome_bwamem2_index // channel: [mandatory] /path/to/genome_bwa-mem2_index_dir/

    main:
    // Channel for version.yml files
    // channel: [ versions.yml ]
    ch_versions = Channel.empty()

    // Sort inputs, separate by tumor, normal and donor
    // channel: runnable: [ meta ]
    // channel: skip: [ meta ]
    ch_inputs_tumor_sorted = ch_inputs
        .branch { meta ->
            runnable: sample.Inputs.hasTumorDnaBam(meta)
            skip: true
        }

    ch_inputs_normal_sorted = ch_inputs
        .branch { meta ->
            runnable: sample.Inputs.hasNormalDnaBam(meta)
            skip: true
        }

    ch_inputs_donor_sorted = ch_inputs
        .branch { meta ->
            runnable: sample.Inputs.hasDonorDnaBam(meta)
            skip: true
        }

    // Build a flat channel of [ meta_bwamem2, bam ] from input BAM/CRAM
    // channel: [ meta_bwamem2, bam ]
    ch_bwamem2_inputs = Channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { meta ->
                def meta_sample = sample.Inputs.getTumorDnaSample(meta)
                def sample_id   = meta_sample.getOrDefault('longitudinal_sample_id', meta_sample['sample_id'])
                def meta_bwamem2 = [
                    key:         meta.group_id,
                    id:          "${meta.group_id}_${sample_id}",
                    sample_id:   sample_id,
                    sample_type: 'tumor',
                    read_group:  "${sample_id}.realign",
                ]
                def bam = sample.Inputs.get(meta, sample.FileKey.BAM_DNA_TUMOR)
                assert bam != null && bam != [] : "Missing tumor DNA BAM/CRAM for ${meta.group_id}"
                return [meta_bwamem2, bam]
            },
            ch_inputs_normal_sorted.runnable.map { meta ->
                def meta_sample = sample.Inputs.getNormalDnaSample(meta)
                def sample_id   = meta_sample['sample_id']
                def meta_bwamem2 = [
                    key:         meta.group_id,
                    id:          "${meta.group_id}_${sample_id}",
                    sample_id:   sample_id,
                    sample_type: 'normal',
                    read_group:  "${sample_id}.realign",
                ]
                def bam = sample.Inputs.get(meta, sample.FileKey.BAM_DNA_NORMAL)
                assert bam != null && bam != [] : "Missing normal DNA BAM/CRAM for ${meta.group_id}"
                return [meta_bwamem2, bam]
            },
            ch_inputs_donor_sorted.runnable.map { meta ->
                def meta_sample = sample.Inputs.getDonorDnaSample(meta)
                def sample_id   = meta_sample['sample_id']
                def meta_bwamem2 = [
                    key:         meta.group_id,
                    id:          "${meta.group_id}_${sample_id}",
                    sample_id:   sample_id,
                    sample_type: 'donor',
                    read_group:  "${sample_id}.realign",
                ]
                def bam = sample.Inputs.get(meta, sample.FileKey.BAM_DNA_DONOR)
                assert bam != null && bam != [] : "Missing donor DNA BAM/CRAM for ${meta.group_id}"
                return [meta_bwamem2, bam]
            },
        )
        .filter { meta_bwamem2, bam -> bam != null && bam != [] }

    // Run process
    BWAMEM2_ALIGN_FROM_BAM(
        ch_bwamem2_inputs,
        genome_fasta,
        genome_bwamem2_index,
    )

    ch_versions = ch_versions.mix(BWAMEM2_ALIGN_FROM_BAM.out.versions)

    // Reunite BAMs
    // First, count expected BAMs per sample for non-blocking groupTuple op (one per sample for the simple
    // single-BAM-in case, but written generically to match the FASTQ path's reunite pattern)
    // channel: [ meta_count, group_size ]
    ch_sample_counts = ch_bwamem2_inputs
        .map { meta_bwamem2, bam -> [[key: meta_bwamem2.key, sample_type: meta_bwamem2.sample_type], meta_bwamem2] }
        .groupTuple()
        .map { meta_count, metas -> [meta_count, metas.size()] }

    // Now, group with expected size then sort into tumor, normal and donor channels
    // channel: [ meta_group, [bam, ...], [bai, ...] ]
    ch_bams_united = ch_sample_counts
        .cross(
            // First element to match meta_count above for `cross`
            BWAMEM2_ALIGN_FROM_BAM.out.bam
                .map { meta_bwamem2, bam, bai -> [[key: meta_bwamem2.key, sample_type: meta_bwamem2.sample_type], bam, bai] }
        )
        .map { count_tuple, bam_tuple ->
            def group_size = count_tuple[1]
            def (meta_bam, bam, bai) = bam_tuple
            return tuple(groupKey([*:meta_bam], group_size), bam, bai)
        }
        .groupTuple()
        .branch { meta_group, bams, bais ->
            assert ['tumor', 'normal', 'donor'].contains(meta_group.sample_type)
            tumor:  meta_group.sample_type == 'tumor'
            normal: meta_group.sample_type == 'normal'
            donor:  meta_group.sample_type == 'donor'
            placeholder: true
        }

    // Set outputs, restoring original meta
    // channel: [ meta, [bam, ...], [bai, ...] ]
    ch_bam_tumor_out = Channel.empty()
        .mix(
            channels.WorkflowChannels.restoreMeta(ch_bams_united.tumor,  ch_inputs),
            channels.PlaceholderChannels.bamBai(ch_inputs_tumor_sorted.skip),
        )

    ch_bam_normal_out = Channel.empty()
        .mix(
            channels.WorkflowChannels.restoreMeta(ch_bams_united.normal, ch_inputs),
            channels.PlaceholderChannels.bamBai(ch_inputs_normal_sorted.skip),
        )

    ch_bam_donor_out = Channel.empty()
        .mix(
            channels.WorkflowChannels.restoreMeta(ch_bams_united.donor,  ch_inputs),
            channels.PlaceholderChannels.bamBai(ch_inputs_donor_sorted.skip),
        )

    emit:
    dna_tumor  = ch_bam_tumor_out   // channel: [ meta, [bam, ...], [bai, ...] ]
    dna_normal = ch_bam_normal_out  // channel: [ meta, [bam, ...], [bai, ...] ]
    dna_donor  = ch_bam_donor_out   // channel: [ meta, [bam, ...], [bai, ...] ]

    versions   = ch_versions        // channel: [ versions.yml ]
}
