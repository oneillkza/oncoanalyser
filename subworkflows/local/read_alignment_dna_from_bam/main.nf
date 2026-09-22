//
// Align DNA reads from existing alignments, realigning via bwa-mem2
//
// NOTE(KO): This subworkflow mirrors READ_ALIGNMENT_DNA but sources reads from an existing BAM/CRAM
// rather than FASTQ. It emits the same channel shapes as READ_ALIGNMENT_DNA (tumor/normal/donor as
// [ meta, [aln, ...], [idx, ...] ]) so that the downstream REDUX and analysis stages consume its
// outputs unchanged.
//
// Unlike the FASTQ path there is no lane or flowcell information available, and read group overrides
// are rejected for non-FASTQ inputs, so a single synthetic read group is constructed per sample from
// the sample and library identifiers. Consequently there is also no FASTQ splitting, and exactly one
// alignment is produced per sample.
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
    //
    // STEP: Handle inputs
    //
    // Sort inputs, separating by sample type
    // runnable: channel: [ meta ]
    // skip: channel: [ meta ]
    ch_inputs_tumor_sorted = ch_inputs
        .branch { meta ->
            runnable: Utils.hasTumorDnaBam(meta)
            skip: true
        }

    ch_inputs_normal_sorted = ch_inputs
        .branch { meta ->
            runnable: Utils.hasNormalDnaBam(meta)
            skip: true
        }

    ch_inputs_donor_sorted = ch_inputs
        .branch { meta ->
            runnable: Utils.hasDonorDnaBam(meta)
            skip: true
        }

    //
    // MODULE: BWA-MEM2
    //
    // Create process input channel
    // channel: [ meta_bwamem2, aln, idx ]
    ch_bwamem2_inputs = channel.empty()
        .mix(
            ch_inputs_tumor_sorted.runnable.map { meta ->
                [meta, Utils.getTumorDnaSample(meta), 'tumor', Utils.getTumorDnaSampleName(meta), Utils.getTumorDnaBam(meta), Utils.getTumorDnaBai(meta)]
            },
            ch_inputs_normal_sorted.runnable.map { meta ->
                [meta, Utils.getNormalDnaSample(meta), 'normal', Utils.getNormalDnaSampleName(meta), Utils.getNormalDnaBam(meta), Utils.getNormalDnaBai(meta)]
            },
            ch_inputs_donor_sorted.runnable.map { meta ->
                [meta, Utils.getDonorDnaSample(meta), 'donor', Utils.getDonorDnaSampleName(meta), Utils.getDonorDnaBam(meta), Utils.getDonorDnaBai(meta)]
            },
        )
        .map { meta, meta_sample, sample_type, sample_id, aln, idx ->

            // NOTE(KO): the FASTQ path builds the read group identifier from sample, library, lane, and
            // flowcell. Only sample and library are meaningful when collapsing an entire alignment into a
            // single read group, so the identifier is formed from those with a suffix marking the realignment.
            def library_id = meta_sample.library_id
            def rg_id = "${sample_id}.${library_id}.realign"
            def rg_entries = [ID: rg_id, SM: sample_id, LB: library_id]
            def rg_line = '@RG\\t' + rg_entries.collect { k, v -> "${k}:${v}" }.join('\\t')

            def meta_bwamem2 = [
                key: meta.group_id,
                id: "${meta.group_id}_${sample_id}",
                rg_line: rg_line,
                sample_id: sample_id,
                library_id: library_id,
                output_file_id: rg_id,
                sample_type: sample_type,
            ]

            return [meta_bwamem2, aln, idx]
        }

    // Run process
    BWAMEM2_ALIGN_FROM_BAM(
        ch_bwamem2_inputs,
        genome_fasta,
        genome_bwamem2_index,
    )

    // Collect alignments, grouping by sample
    // NOTE(KO): exactly one alignment is produced per sample here, but the counting and grouping pattern of
    // READ_ALIGNMENT_DNA is retained so that the emitted channel shape is identical and downstream stages
    // require no special handling.
    // channel: [ meta_group, group_size ]
    ch_sample_counts = ch_bwamem2_inputs
        .map { meta_bwamem2, _aln, _idx ->

            def meta_group = [
                key: meta_bwamem2.key,
                sample_type: meta_bwamem2.sample_type,
            ]

            return [meta_group, meta_bwamem2]
        }
        .groupTuple()
        .map { meta_group, metas_bwamem2 -> return [meta_group, metas_bwamem2.size()] }

    // Now, group with expected size then sort into tumor, normal and donor channels
    // channel: [ meta_group, [aln, ...], [idx, ...] ]
    ch_alns_united = ch_sample_counts
        // channel: [ [ meta_group, count ], [ meta_group, aln, idx ] ]
        .cross(
            // First element to match meta_group above for `cross`
            channel.topic('bwamem2_align_bam').map { meta_bwamem2, aln, idx -> [[key: meta_bwamem2.key, sample_type: meta_bwamem2.sample_type], aln, idx] }
        )
        .map { count_tuple, inputs_tuple ->
            def group_size = count_tuple[1]
            def (meta_group, aln, idx) = inputs_tuple

            return tuple(groupKey(meta_group, group_size), aln, idx)
        }
        .groupTuple()
        .branch { meta_group, alns, idxs ->
            assert ['tumor', 'normal', 'donor'].contains(meta_group.sample_type)
            tumor: meta_group.sample_type == 'tumor'
            normal: meta_group.sample_type == 'normal'
            donor: meta_group.sample_type == 'donor'
            placeholder: true
        }

    //
    // STEP: Handle outputs
    //
    // Set outputs, restoring original meta
    // channel: [ meta, [aln, ...], [idx, ...] ]
    ch_outputs_tumor = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_alns_united.tumor, ch_inputs),
            ch_inputs_tumor_sorted.skip.map { meta -> [meta, [], []] },
        )

    // channel: [ meta, [aln, ...], [idx, ...] ]
    ch_outputs_normal = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_alns_united.normal, ch_inputs),
            ch_inputs_normal_sorted.skip.map { meta -> [meta, [], []] },
        )

    // channel: [ meta, [aln, ...], [idx, ...] ]
    ch_outputs_donor = channel.empty()
        .mix(
            WorkflowOncoanalyser.restoreMeta(ch_alns_united.donor, ch_inputs),
            ch_inputs_donor_sorted.skip.map { meta -> [meta, [], []] },
        )

    emit:
    tumor  = ch_outputs_tumor  // channel: [ meta, [aln, ...], [idx, ...] ]
    normal = ch_outputs_normal // channel: [ meta, [aln, ...], [idx, ...] ]
    donor  = ch_outputs_donor  // channel: [ meta, [aln, ...], [idx, ...] ]
}
