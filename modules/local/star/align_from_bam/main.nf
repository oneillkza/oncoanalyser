process STAR_ALIGN_FROM_BAM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-69a6f67cb46e41b4c393f71634a9956d5e31f3e9:46745d95bbdc75d9503849416a66ac6555567ff0-0' :
        'biocontainers/mulled-v2-69a6f67cb46e41b4c393f71634a9956d5e31f3e9:46745d95bbdc75d9503849416a66ac6555567ff0-0' }"

    input:
    tuple val(meta), val(rg_lines), path(aln_input), path(idx_input)   // BAM or CRAM to be realigned, plus index
    path genome_star_index

    output:
    tuple val(meta), path('*bam')                         , topic: star_align_bam
    tuple val(meta), path('*Log.final.out')               , topic: star_align_qc_log
    tuple val(meta), val('star_align'), path('.command.*'), topic: command_files
    path 'versions.yml'                                   , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def rg_lines_str = rg_lines.join(' , ')

    // NOTE(KO): reads are recovered from the input alignment rather than read from FASTQ. samtools collate
    // restores read-name grouping, then a single samtools fastq pass writes matched R1/R2 files. Writing the
    // pair to disk (rather than streaming both mates via process substitution) is required for correctness:
    // STAR needs R1 and R2 in strictly corresponding order, and a single samtools fastq invocation guarantees
    // this while discarding singleton and orphaned mates. Secondary and supplementary records are dropped
    // with -F 0x900 so that only primary reads are realigned.

    """
    samtools collate \\
        -O \\
        -u \\
        -@ ${task.cpus} \\
        ${aln_input} | \\
        \\
        samtools fastq \\
            -@ ${task.cpus} \\
            -n \\
            -F 0x900 \\
            -1 reads_R1.fastq.gz \\
            -2 reads_R2.fastq.gz \\
            -0 /dev/null \\
            -s /dev/null \\
            -

    STAR \\
        ${args} \\
        --readFilesIn reads_R1.fastq.gz reads_R2.fastq.gz \\
        --outSAMattrRGline ${rg_lines_str} \\
        --genomeDir ${genome_star_index} \\
        --runThreadN ${task.cpus} \\
        --readFilesCommand zcat \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --alignSplicedMateMapLmin 35 \\
        --alignSplicedMateMapLminOverLmate 0.33 \\
        --chimJunctionOverhangMin 10 \\
        --chimOutType WithinBAM SoftClip \\
        --chimScoreDropMax 70 \\
        --chimScoreJunctionNonGTAG 0 \\
        --chimScoreMin 1 \\
        --chimScoreSeparation 1 \\
        --chimSegmentMin 10 \\
        --chimSegmentReadGapMax 3 \\
        --limitOutSJcollapsed 3000000 \\
        --outBAMcompression 0 \\
        --outFilterMatchNmin 35 \\
        --outFilterMatchNminOverLread 0.33 \\
        --outFilterMismatchNmax 3 \\
        --outFilterMultimapNmax 10 \\
        --outFilterScoreMinOverLread 0.33 \\
        --outSAMattributes All \\
        --outSAMtype BAM Unsorted \\
        --outSAMunmapped Within \\
        --peOverlapNbasesMin 10 \\
        --runRNGseed 0

    rm reads_R1.fastq.gz reads_R2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
        samtools: \$(samtools --version | sed -n '/^samtools / { s/^.* //p }')
    END_VERSIONS
    """

    stub:
    """
    touch Aligned.out.bam
    touch Log.final.out

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}
