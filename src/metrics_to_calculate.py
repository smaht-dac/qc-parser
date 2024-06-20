########################################################################
# Metrics to calculate in a postprocessing step from the extracted metrics
########################################################################

samtools_stats_calculated_metrics = [
    {
        "key": "Percentage of Reads Paired [Samtools]",
        "tooltip": "Percentage of paired reads out of processed reads (sequences)",
        "derived_from": "samtools_stats_postprocessed:percentage_paired_reads",
        "type": float,
        "formula": "{samtools_stats:reads_paired}/{samtools_stats:sequences}*100",
        "visible": True,
    },
    {
        "key": "Percentage of Reads Mapped [Samtools]",
        "tooltip": "Percentage of mapped reads out of processed reads (sequences)",
        "derived_from": "samtools_stats_postprocessed:percentage_mapped_reads",
        "type": float,
        "formula": "{samtools_stats:reads_mapped}/{samtools_stats:sequences}*100",
        "visible": True,
    },
    {
        "key": "Percentage of Reads Mapped and Paired [Samtools]",
        "tooltip": "Percentage of mapped and paired reads out of processed reads (sequences)",
        "derived_from": "samtools_stats_postprocessed:percentage_reads_mapped_and_paired",
        "type": float,
        "formula": "{samtools_stats:reads_mapped_and_paired}/{samtools_stats:sequences}*100",
        "visible": True,
    },
    {
        "key": "Percentage of Reads Unmapped [Samtools]",
        "tooltip": "Percentage of unmapped reads out of processed reads (sequences)",
        "derived_from": "samtools_stats_postprocessed:percentage_unmapped_reads",
        "type": float,
        "formula": "{samtools_stats:reads_unmapped}/{samtools_stats:sequences}*100",
    },
    {
        "key": "Percentage of Reads Duplicated [Samtools]",
        "tooltip": "Percentage of duplicated reads out of processed reads (sequences)",
        "derived_from": "samtools_stats_postprocessed:percentage_duplicated_reads",
        "type": float,
        "formula": "{samtools_stats:reads_duplicated}/{samtools_stats:sequences}*100",
        "visible": True,
    },
    {
        "key": "Percentage of Reads MQ0 [Samtools]",
        "tooltip": "Percentage of MQ0 reads out of processed reads (sequences)",
        "derived_from": "samtools_stats_postprocessed:percentage_reads_MQ0",
        "type": float,
        "formula": "{samtools_stats:reads_MQ0}/{samtools_stats:sequences}*100",
        "visible": True,
    },
]

bamstats_calculated_metrics = [
    {
        "key": "Percentage of Uniquely Aligned Reads with Q>=10 [bamstats]",
        "tooltip": "Percentage of uniquely aligned reads with Q>=10",
        "derived_from": "bamstats_postprocessed:percentage_number_of_uniquely_aligned_reads_with_q_>=_10",
        "type": float,
        "formula": "{bamstats:number_of_uniquely_aligned_reads_with_q_>=_10}/{bamstats:Total_Number_Of_Reads}*100",
        "visible": True,
    },
    {
        "key": "Percentage of Uniquely Aligned Reads without Duplicates with Q>=10 [bamstats]",
        "tooltip": "Percentage of uniquely aligned reads without duplicates with Q>=10",
        "derived_from": "bamstats_postprocessed:percentage_number_of_uniquely_aligned_reads_without_dups_and_q_>=_10",
        "type": float,
        "formula": "{bamstats:number_of_uniquely_aligned_reads_without_dups_and_q_>=_10}/{bamstats:Total_Number_Of_Reads}*100",
        "visible": True,
    },
]
