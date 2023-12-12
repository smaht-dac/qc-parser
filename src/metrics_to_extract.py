# Tools
SAMTOOLS_STATS = 'samtools_stats'
BAMSTATS = 'bamstats'
RNASEQQC = 'rnaseqqc'
PICARD_COLLECT_ALIGNMENT_SUMMARY_METRICS = 'picard_CollectAlignmentSummaryMetrics'
PICARD_COLLECT_INSERT_SIZE_METRICS = 'picard_CollectInsertSizeMetrics'
PICARD_COLLECT_WGS_METRICS = 'picard_CollectWgsMetrics'
FASTQC = 'fastqc'


########################################################################
# Metrics to extract for each tool
########################################################################
samtools_stats_metrics = {
    'raw total sequences': {
        'key': 'Total Sequences [Samtools]',
        'tooltip': 'Total number of reads in a file, excluding supplementary and secondary reads',
        'derived_from': 'samtools_stats:raw_total_sequences',
        'type': int
    },
    'average length': {
        'key': 'Average Length [Samtools]',
        'tooltip': 'Ratio between total length and sequences',
        'derived_from': 'samtools_stats:average_length',
        'type': float
    },
    'maximum length': {
        'key': 'Maximum Length [Samtools]',
        'tooltip': 'Length of the longest read (includes hard-clipped bases).',
        'derived_from': 'samtools_stats:maximum_length',
        'type': int
    },
    'bases mapped': {
        'key': 'Bases Mapped [Samtools]',
        'tooltip': 'Number of processed bases that belong to reads mapped',
        'derived_from': 'samtools_stats:bases_mapped',
        'type': int
    },
    'sequences': {
        'key': 'Processed Sequences [Samtools]',
        'tooltip': 'Number of processed reads',
        'derived_from': 'samtools_stats:sequences',
        'type': int
    },
    '1st fragments': {
        'key': '1st Fragments [Samtools]',
        'tooltip': 'Number of first fragment reads',
        'derived_from': 'samtools_stats:1st_fragments',
        'type': int
    },
    'last fragments': {
        'key': 'Last Fragments [Samtools]',
        'tooltip': 'Number of last fragment reads',
        'derived_from': 'samtools_stats:last_fragments',
        'type': int
    },
    'reads paired': {
        'key': 'Reads Paired [Samtools]',
        'tooltip': 'Number of paired reads, mapped or unmapped, that are neither secondary nor supplementary (paired-end technology bit set)',
        'derived_from': 'samtools_stats:reads_paired',
        'type': int
    },
    'reads mapped': {
        'key': 'Reads Mapped [Samtools]',
        'tooltip': 'Number of reads, paired or single, that are mapped',
        'derived_from': 'samtools_stats:reads_mapped',
        'type': int
    },
    'reads mapped and paired': {
        'key': 'Reads Mapped and Paired [Samtools]',
        'tooltip': 'Number of mapped paired reads (paired-end technology bit set with both mates mapped)',
        'derived_from': 'samtools_stats:reads_mapped_and_paired',
        'type': int
    },
    'reads properly paired': {
        'key': 'Reads Properly Paired [Samtools]',
        'tooltip': 'Number of mapped paired reads (proper-pair bit set, both mates mapped within the expected distance)',
        'derived_from': 'samtools_stats:reads_properly_paired',
        'type': int
    },
    'reads unmapped': {
        'key': 'Reads Unmapped [Samtools]',
        'tooltip': 'Number of unmapped reads',
        'derived_from': 'samtools_stats:reads_unmapped',
        'type': int
    },
    'reads duplicated': {
        'key': 'Reads Duplicated [Samtools]',
        'tooltip': 'Number of duplicate reads',
        'derived_from': 'samtools_stats:reads_duplicated',
        'type': int
    },
    'reads MQ0': {
        'key': 'Reads MQ0 [Samtools]',
        'tooltip': 'Number of mapped reads with mapping quality 0',
        'derived_from': 'samtools_stats:reads_MQ0',
        'type': int
    },
    'reads QC failed': {
        'key': 'Reads QC-Failed [Samtools]',
        'tooltip': 'Number of reads that failed the quality checks',
        'derived_from': 'samtools_stats:reads_QC_failed',
        'type': int
    },
    'non-primary alignments': {
        'key': 'Non-Primary Alignments [Samtools]',
        'tooltip': 'Number of secondary reads',
        'derived_from': 'samtools_stats:non_primary_alignments',
        'type': int
    },
    'supplementary alignments': {
        'key': 'Supplementary Alignments [Samtools]',
        'tooltip': 'Number of supplementary reads',
        'derived_from': 'samtools_stats:supplementary_alignments',
        'type': int
    },
    'pairs on different chromosomes': {
        'key': 'Pairs on Different Chromosomes [Samtools]',
        'tooltip': 'Number of pairs where one read is on one chromosome and the mate read is on a different chromosome',
        'derived_from': 'samtools_stats:pairs_on_different_chromosomes',
        'type': int
    },
    'percentage of properly paired reads (%)': {
        'key': 'Percentage of Properly Paired Reads [Samtools]',
        'tooltip': None,
        'derived_from': 'samtools_stats:percentage_of_properly_paired_reads',
        'type': float
    },
}

picard_CollectAlignmentSummaryMetrics_metrics = {
    'PF_ALIGNED_BASES': {
        'key': 'Aligned Bases [Picard]',
        'tooltip': 'The total number of aligned bases',
        'derived_from': 'picard_collect_alignment_summary_metrics:pf_aligned_bases',
        'type': int
    },
    'PF_HQ_ALIGNED_BASES': {
        'key': 'Aligned Bases (High Quality) [Picard]',
        'tooltip': 'The number of aligned bases in reads that were mapped at high quality',
        'derived_from': 'picard_collect_alignment_summary_metrics:pf_hq_aligned_bases',
        'type': int
    },
    'PF_MISMATCH_RATE': {
        'key': 'Aligned Bases Mismatch Rate [Picard]',
        'tooltip': 'The fraction of bases mismatching the reference for all aligned bases',
        'derived_from': 'picard_collect_alignment_summary_metrics:pf_mismatch_rate',
        'type': float
    },
    'PF_HQ_ERROR_RATE': {
        'key': 'Aligned Bases Mismatch Rate (High Quality) [Picard]',
        'tooltip': 'The fraction of bases mismatching the reference in reads that were mapped at high quality',
        'derived_from': 'picard_collect_alignment_summary_metrics:pf_hq_error_rate',
        'type': float
    },
    'PF_INDEL_RATE': {
        'key': 'Indel Rate [Picard]',
        'tooltip':  'The number of insertion and deletion events per 100 aligned bases',
        'derived_from': 'picard_collect_alignment_summary_metrics:pf_indel_rate',
        'type': float
    },
    'MEAN_READ_LENGTH': {
        'key': 'Mean Read Length [Picard]',
        'tooltip': 'The mean length of the set of reads examined',
        'derived_from': 'picard_collect_alignment_summary_metrics:mean_read_length',
        'type': int
    },
    'SD_READ_LENGTH': {
        'key': 'Read Length Standard Deviation [Picard]',
        'tooltip': 'The standard deviation for the length of the set of reads examined',
        'derived_from': 'picard_collect_alignment_summary_metrics:sd_read_length',
        'type': int
    },
}

picard_CollectInsertSizeMetrics_metrics = {
    # [5]
    'MEAN_INSERT_SIZE': {
        'key': 'Mean Insert Size',
        'tooltip': 'The mean insert size for the pair orientation',
        'derived_from': 'picard_collect_insert_size_metrics:mean_insert_size',
        'type': float
    },
    # [6]
    'STANDARD_DEVIATION': {
        'key': 'Insert Size Standard Deviation',
        'tooltip': 'Standard deviation of insert sizes for the pair orientation',
        'derived_from': 'picard_collect_insert_size_metrics:standard_deviation',
        'type': float
    },
    # [7]
    'READ_PAIRS': {
        'key': 'Total Number of Read Pairs',
        'tooltip': 'The total number of read pairs that were examined for the pair orientation',
        'derived_from': 'picard_collect_insert_size_metrics:read_pairs',
        'type': int
    },
    # 'PAIR_ORIENTATION' [8]
}

picard_CollectWgsMetrics_metrics = {
    'GENOME_TERRITORY': {
        'key': 'Effective Genome Size [Picard]',
        'tooltip': 'The number of non-N bases in the genome',
        'derived_from': 'picard_collect_wgs_metrics:genome_territory',
        'type': int
    },
    'MEAN_COVERAGE': {
        'key': 'Mean Coverage [Picard]',
        'tooltip': 'The mean coverage of the genome',
        'derived_from': 'picard_collect_wgs_metrics:mean_coverage',
        'type': float
    },
    'SD_COVERAGE': {
        'key': 'Coverage Standard Deviation [Picard]',
        'tooltip': 'The standard deviation for the coverage',
        'derived_from': 'picard_collect_wgs_metrics:sd_coverage',
        'type': float
    },
}

bamstats_metrics = {
    'Estimate_Average_Coverage': {
        'key': 'Average coverage (estimated)',
        'tooltip': 'Estimated average coverage',
        'derived_from': 'bamstats:estimate_average_coverage',
        'type': float
    },
    'Number_of_Uniquely_Aligned_Reads_with_Q_>=_10': {
        'key': 'Uniquely Aligned Reads with Q>=10 [bamstats]',
        'tooltip': None,
        'derived_from': 'bamstats:number_of_uniquely_aligned_reads_with_q_>=_10',
        'type': float
    },
    'Number_of_Uniquely_Aligned_Reads_without_Dups_and_Q_>=_10': {
        'key': 'Uniquely Aligned Reads without Duplicates with Q>=10 [bamstats]',
        'tooltip': None,
        'derived_from': 'bamstats:number_of_uniquely_aligned_reads_without_dups_and_q_>=_10',
        'type': float
    }
}

rnaseqqc_metrics = {
    "3' bias MAD_Std": {
        'key': "3' bias MAD_Std",
        'tooltip': "3' Bias statistics (Mean, Median, Std Deviation, Median Absolute Deviation, 25th percentile, 75th percentile): These aggregate statistics are based on the total coverage in 100 bp windows on both the 3' and 5' ends of a gene. The windows are both offset 150 bases into the gene. This computation is only performed on genes at least 600bp long and with at least 5 unambiguous reads. A gene with even coverage in both it's 3' and 5' windows would have a bias of 0.5; bias near 1 or 0 may indicate degredation",
        'derived_from': 'rnaseqqc:3p_bias_mad_std',
        'type': float
    },
    "3' bias Std": {
        'key': "3' bias Std",
        'tooltip': "3' Bias statistics (Mean, Median, Std Deviation, Median Absolute Deviation, 25th percentile, 75th percentile): These aggregate statistics are based on the total coverage in 100 bp windows on both the 3' and 5' ends of a gene. The windows are both offset 150 bases into the gene. This computation is only performed on genes at least 600bp long and with at least 5 unambiguous reads. A gene with even coverage in both it's 3' and 5' windows would have a bias of 0.5; bias near 1 or 0 may indicate degredation",
        'derived_from': 'rnaseqqc:3p_bias_std',
        'type': float
    },
    "3' Bias, 25th Percentile": {
        'key': "3' Bias, 25th Percentile",
        'tooltip': "3' Bias statistics (Mean, Median, Std Deviation, Median Absolute Deviation, 25th percentile, 75th percentile): These aggregate statistics are based on the total coverage in 100 bp windows on both the 3' and 5' ends of a gene. The windows are both offset 150 bases into the gene. This computation is only performed on genes at least 600bp long and with at least 5 unambiguous reads. A gene with even coverage in both it's 3' and 5' windows would have a bias of 0.5; bias near 1 or 0 may indicate degredation",
        'derived_from': 'rnaseqqc:3p_bias_25th_percentile',
        'type': float
    },
    "3' Bias, 75th Percentile": {
        'key': "3' Bias, 75th Percentile",
        'tooltip': "3' Bias statistics (Mean, Median, Std Deviation, Median Absolute Deviation, 25th percentile, 75th percentile): These aggregate statistics are based on the total coverage in 100 bp windows on both the 3' and 5' ends of a gene. The windows are both offset 150 bases into the gene. This computation is only performed on genes at least 600bp long and with at least 5 unambiguous reads. A gene with even coverage in both it's 3' and 5' windows would have a bias of 0.5; bias near 1 or 0 may indicate degredation",
        'derived_from': 'rnaseqqc:3p_bias_75th_percentile',
        'type': float
    },
    "Genes used in 3' bias": {
        'key': "Genes used in 3' bias",
        'tooltip': "",
        'derived_from': 'rnaseqqc:genes_used_in_3p_bias',
        'type': int
    },
    "Mean 3' bias": {
        'key': "Mean 3' bias",
        'tooltip': "",
        'derived_from': 'rnaseqqc:mean_3p_bias',
        'type': float
    },
    "Median 3' bias": {
        'key': "Median 3' bias",
        'tooltip': "",
        'derived_from': 'rnaseqqc:median_3p_bias',
        'type': float
    },
    'Alternative Alignments': {
        'key': 'Alternative Alignments',
        'tooltip': 'Duplicate read entries providing alternative coordinates',
        'derived_from': 'rnaseqqc:alternative_alignments',
        'type': int
    },
    'Base Mismatch': {
        'key': 'Base Mismatch',
        'tooltip': 'The total number of mismatched bases (as determined by the "NM" tag) of all "Mapped Reads" (as defined above) divided by the total aligned length of all "Mapped Reads".',
        'derived_from': 'rnaseqqc:base_mismatch',
        'type': float
    },
    'Chimeric Reads': {
        'key': 'Chimeric Reads',
        'tooltip': 'Chimeric Reads',
        'derived_from': 'rnaseqqc:chimeric_reads',
        'type': int
    },
    'Duplicate Rate of Mapped': {
        'key': 'Duplicate Rate of Mapped',
        'tooltip': 'This is the proportion of all reads which were marked as PCR/Optical Duplicates out of all "Mapped Reads" (as defined above; excludes Secondary and Vendor QC Failed reads). This is complementary to the "Unique Rate of Mapped".',
        'derived_from': 'rnaseqqc:duplicate_rate_of_mapped',
        'type': int
    },
    'End 1 Antisense': {
        'key': 'End 1 Antisense',
        'tooltip': 'Number of reads that were sequenced in the antisense direction',
        'derived_from': 'rnaseqqc:end_1_antisense',
        'type': int
    },
    'End 1 Mapping Rate': {
        'key': 'End 1 Mapping Rate',
        'tooltip': 'The proportion of Paired reads which were marked as First or Second in the pair, respectively, out of all "Mapped Reads" (above).',
        'derived_from': 'rnaseqqc:end_1_mapping_rate',
        'type': float
    },
    'End 1 Mismatch Rate': {
        'key': 'End 1 Mismatch Rate',
        'tooltip': 'The proportion of mismatched bases (as determined by the "NM" tag) belonging to First or Second mates, divided by the total aligned length of all "Mapped" (above) First or Second mates, respectively.',
        'derived_from': 'rnaseqqc:end_1_mismatch_rate',
        'type': float
    },
    'End 1 Sense': {
        'key': 'End 1 Sense',
        'tooltip': 'Number of End 1 reads that were sequenced in the sense direction',
        'derived_from': 'rnaseqqc:end_1_sense',
        'type': int
    },
    'End 1 Sense Rate': {
        'key': 'End 1 Sense Rate',
        'tooltip': 'The proportion of First or Second Mate reads which intersected with a Sense Strand feature out of all First or Second Mate reads which intersected with any features, respectively.',
        'derived_from': 'rnaseqqc:end_1_sense_rate',
        'type': float
    },
    'End 2 Antisense': {
        'key': 'End 2 Antisense',
        'tooltip': 'Number of reads that were sequenced in the antisense direction',
        'derived_from': 'rnaseqqc:end_2_antisense',
        'type': int
    },
    'End 2 Mapping Rate': {
        'key': 'End 2 Mapping Rate',
        'tooltip': 'The proportion of Paired reads which were marked as First or Second in the pair, respectively, out of all "Mapped Reads" (above).',
        'derived_from': 'rnaseqqc:end_2_mapping_rate',
        'type': float
    },
    'End 2 Mismatch Rate': {
        'key': 'End 2 Mismatch Rate',
        'tooltip': 'The proportion of mismatched bases (as determined by the "NM" tag) belonging to First or Second mates, divided by the total aligned length of all "Mapped" (above) First or Second mates, respectively.',
        'derived_from': 'rnaseqqc:end_2_mismatch_rate',
        'type': float
    },
    'End 2 Sense': {
        'key': 'End 2 Sense',
        'tooltip': 'Number of End 2 reads that were sequenced in the sense direction',
        'derived_from': 'rnaseqqc:end_2_sense',
        'type': int
    },
    'End 2 Sense Rate': {
        'key': 'End 2 Sense Rate',
        'tooltip': 'The proportion of First or Second Mate reads which intersected with a Sense Strand feature out of all First or Second Mate reads which intersected with any features, respectively.',
        'derived_from': 'rnaseqqc:end_2_sense_rate',
        'type': float
    },
    'Estimated Library Complexity': {
        'key': 'Estimated Library Complexity',
        'tooltip': 'An estimation of the number of unique cDNA fragments present in the library. This computation follows the same formula as Picard EstimateLibraryComplexity',
        'derived_from': 'rnaseqqc:estimated_library_complexity',
        'type': int
    },
    'Exonic Rate': {
        'key': 'Exonic Rate',
        'tooltip': 'The proportion of "Mapped Reads" (above) for which all aligned segments unambiguously aligned to exons of the same gene.',
        'derived_from': 'rnaseqqc:exonic_rate',
        'type': float
    },
    'Exonic/Intron ratio': {
        'key': 'Exonic/Intron ratio',
        'tooltip': 'Exonic Reads/Intronic Reads',
        'derived_from': 'rnaseqqc:exonic_intron_ratio',
        'type': float
    },
    'Expression Profiling Efficiency': {
        'key': 'Expression Profiling Efficiency',
        'tooltip': 'The proportion of "Exonic Reads" (see "Exonic Rate", below) out of all reads which were not Secondary Alignments or Platform/Vendor QC Failing reads.',
        'derived_from': 'rnaseqqc:expression_profiling_efficiency',
        'type': float
    },
    'Exons with >0 reads': {
        'key': 'Exons with >0 reads',
        'tooltip': 'Exons with >0 reads',
        'derived_from': 'rnaseqqc:Exons_with_>0_reads',
        'type': int
    },
    'Exons with >=2 reads': {
        'key': 'Exons with >=2 reads',
        'tooltip': 'Exons with >=2 reads',
        'derived_from': 'rnaseqqc:Exons_with_>2_reads',
        'type': int
    },
    'Exons with >=10 reads': {
        'key': 'Exons with >=10 reads',
        'tooltip': 'Exons with >=10 reads',
        'derived_from': 'rnaseqqc:Exons_with_>10_reads',
        'type': int
    },
    'Failed Vendor QC': {
        'key': 'Failed Vendor QC',
        'tooltip': 'Failed Vendor QC',
        'derived_from': 'rnaseqqc:failed_vendor_qc',
        'type': int
    },
    'Genes Detected': {
        'key': 'Genes Detected',
        'tooltip': 'The number of genes which had at least 5 unambiguous reads. The detection threshold can be changed with --detection-threshold',
        'derived_from': 'rnaseqqc:genes_detected',
        'type': int
    },
    'Genes with >0 reads': {
        'key': 'Genes with >0 reads',
        'tooltip': 'Genes with >0 reads',
        'derived_from': 'rnaseqqc:genes_with_>0_reads',
        'type': int
    },
    'Genes with >=2 reads': {
        'key': 'Genes with >=2 reads',
        'tooltip': 'Genes with >=2 reads',
        'derived_from': 'rnaseqqc:genes_with_>2_reads',
        'type': int
    },
    'Genes with >=10 reads': {
        'key': 'Genes with >=10 reads',
        'tooltip': 'Genes with >=10 reads',
        'derived_from': 'rnaseqqc:genes_with_>10_reads',
        'type': int
    },
    'Intergenic Rate': {
        'key': 'Intergenic Rate',
        'tooltip': 'The proportion of "Mapped Reads" (above) for which none of the aligned segments intersected any genes.',
        'derived_from': 'rnaseqqc:intergenic_rate',
        'type': float
    },
    'Intragenic Rate': {
        'key': 'Intragenic Rate',
        'tooltip': 'The sum of "Exonic" and "Intronic" rates (see "Exonic Rate" and "Intronic Rate" above)',
        'derived_from': 'rnaseqqc:intragenic_rate',
        'type': float
    },
    'Intronic Rate': {
        'key': 'Intronic Rate',
        'tooltip': 'The proportion of "Mapped Reads" (above) for which all aligned segments unambiguously aligned to the same gene, but none of which intersected any exons of the gene.',
        'derived_from': 'rnaseqqc:intronic_rate',
        'type': float
    },
    'Mapped Reads': {
        'key': 'Mapped Reads',
        'tooltip': 'Mapped Reads',
        'derived_from': 'rnaseqqc:mapped_reads',
        'type': int
    },
    'Mapped Unique Reads': {
        'key': 'Mapped Unique Reads',
        'tooltip': 'Mapped Unique Reads',
        'derived_from': 'rnaseqqc:mapped_unique_reads',
        'type': int
    },
    'Mapping Rate': {
        'key': 'Mapping Rate',
        'tooltip': 'The proportion of all reads in the Bam which were Mapped, and not Secondary Alignments or Platform/Vendor QC Failing reads ("Mapped Reads").',
        'derived_from': 'rnaseqqc:mapping_rate',
        'type': float
    },
    'Read Length': {
        'key': 'Read Length',
        'tooltip': 'The longest aligned length observed in any read',
        'derived_from': 'rnaseqqc:read_length',
        'type': int
    },
    'rRNA Rate': {
        'key': 'rRNA Rate',
        'tooltip': 'The proportion of "Mapped Reads" (above) which at least partially intersected with an annotated rRNA gene. This is not complementary to any other rates.',
        'derived_from': 'rnaseqqc:rrna_rate',
        'type': float
    },
    'rRNA Reads': {
        'key': 'rRNA Reads',
        'tooltip': 'rRNA Reads',
        'derived_from': 'rnaseqqc:rrna_reads',
        'type': int
    },
    'Total Mapped Pairs': {
        'key': 'Total Mapped Pairs',
        'tooltip': 'Total Mapped Pairs',
        'derived_from': 'rnaseqqc:total_mapped_pairs',
        'type': int
    },
    'Total Reads': {
        'key': 'Total Reads',
        'tooltip': 'Total Reads',
        'derived_from': 'rnaseqqc:total_reads',
        'type': int
    },
    'Unique Rate of Mapped': {
        'key': 'Unique Rate of Mapped',
        'tooltip': 'This is the proportion of reads which were not marked as PCR/Optical Duplicates out of all "Mapped Reads"',
        'derived_from': 'rnaseqqc:unique_rate_of_mapped',
        'type': float
    },
    'Unpaired Reads': {
        'key': 'Unpaired Reads',
        'tooltip': 'Unpaired Reads',
        'derived_from': 'rnaseqqc:unpaired_reads',
        'type': int
    },

}

fastqc_metrics = {
    'Basic Statistics': {
        'key': 'Basic Statistics',
        'tooltip': 'Basic Statistics',
        'derived_from': 'fastqc:basic_statistics',
        'type': str
    },
    'Per base sequence quality': {
        'key': 'Per base sequence quality',
        'tooltip': 'Per base sequence quality',
        'derived_from': 'fastqc:per_base_sequence_quality',
        'type': str
    },
    'Per tile sequence quality': {
        'key': 'Per tile sequence quality',
        'tooltip': 'Per tile sequence quality',
        'derived_from': 'fastqc:per_tile_sequence_quality',
        'type': str
    },
    'Per sequence quality scores': {
        'key': 'Per sequence quality scores',
        'tooltip': 'Per sequence quality scores',
        'derived_from': 'fastqc:per_sequence_quality_scores',
        'type': str
    },
    'Per base sequence content': {
        'key': 'Per base sequence content',
        'tooltip': 'Per base sequence content',
        'derived_from': 'fastqc:per_base_sequence_content',
        'type': str
    },
    'Per sequence GC content': {
        'key': 'Per sequence GC content',
        'tooltip': 'Per sequence GC content',
        'derived_from': 'fastqc:per_sequence_gc_content',
        'type': str
    },
    'Per base N content': {
        'key': 'Per base N content',
        'tooltip': 'Per base N content',
        'derived_from': 'fastqc:per_base_n_content',
        'type': str
    },
    'Sequence Length Distribution': {
        'key': 'Sequence Length Distribution',
        'tooltip': 'Sequence Length Distribution',
        'derived_from': 'fastqc:sequence_length_distribution',
        'type': str
    },
    'Sequence Duplication Levels': {
        'key': 'Sequence Duplication Levels',
        'tooltip': 'Sequence Duplication Levels',
        'derived_from': 'fastqc:sequence_duplication_levels',
        'type': str
    },
    'Overrepresented sequences': {
        'key': 'Overrepresented sequences',
        'tooltip': 'Overrepresented sequences',
        'derived_from': 'fastqc:overrepresented_sequences',
        'type': str
    },
    'Adapter Content': {
        'key': 'Adapter Content',
        'tooltip': 'Adapter Content',
        'derived_from': 'fastqc:adapter_content',
        'type': str
    },
    'Total Sequences': {
        'key': 'Total Sequences',
        'tooltip': 'Total Sequences',
        'derived_from': 'fastqc:total_sequences',
        'type': int
    },
    'Sequence Length': {
        'key': 'Sequence Length',
        'tooltip': 'Sequence Length',
        'derived_from': 'fastqc:sequence_length',
        'type': int
    }
}


metrics = {
    SAMTOOLS_STATS: samtools_stats_metrics,
    BAMSTATS: bamstats_metrics,
    RNASEQQC: rnaseqqc_metrics,
    PICARD_COLLECT_ALIGNMENT_SUMMARY_METRICS: picard_CollectAlignmentSummaryMetrics_metrics,
    PICARD_COLLECT_INSERT_SIZE_METRICS: picard_CollectInsertSizeMetrics_metrics,
    PICARD_COLLECT_WGS_METRICS: picard_CollectWgsMetrics_metrics,
    FASTQC: fastqc_metrics
}
