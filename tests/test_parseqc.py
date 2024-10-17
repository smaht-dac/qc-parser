import pytest
import os, platform
import json
import pprint


def test_fastqc():
    pv = platform.python_version()
    fastqc_stats = "./tests/data/fastqc_summary.txt"
    qc_values = f"tmp.fastqc.qc_values.{pv}.json"
    metrics_zip = f"tmp.fastqc.metrics.{pv}.zip"
    cmd = (
        f"parse-qc -n 'FastQC' "
        f"--metrics fastqc {fastqc_stats} "
        f"--output-zip {metrics_zip} "
        f"--output-json {qc_values}"
    )
    os.system(cmd)

    assert os.path.exists(metrics_zip) == True
    assert os.path.exists(qc_values) == True

    qc_file = open(qc_values)
    data = json.load(qc_file)
    qc_file.close()

    assert len(data["qc_values"]) == 13

    seq_len = next(
        item
        for item in data["qc_values"]
        if item["derived_from"] == "fastqc:sequence_length"
    )
    assert seq_len["value"] == 151

    os.system(f"rm -f {qc_values} {metrics_zip}")


def test_fastqc_2():
    pv = platform.python_version()
    fastqc_stats = "./tests/data/fastqc_summary_2.txt"
    qc_values = f"tmp.fastqc2.qc_values.{pv}.json"
    metrics_zip = f"tmp.fastqc2.metrics.{pv}.zip"
    cmd = (
        f"parse-qc -n 'FastQC2' "
        f"--metrics fastqc {fastqc_stats} "
        f"--output-zip {metrics_zip} "
        f"--output-json {qc_values}"
    )
    os.system(cmd)

    assert os.path.exists(qc_values) == True

    qc_file = open(qc_values)
    data = json.load(qc_file)
    qc_file.close()

    assert len(data["qc_values"]) == 13

    seq_len = next(
        item
        for item in data["qc_values"]
        if item["derived_from"] == "fastqc:sequence_length"
    )
    assert seq_len["value"] == '10-150'

    os.system(f"rm -f {qc_values} {metrics_zip}")


def test_samtools():
    pv = platform.python_version()
    samtools_stats = "./tests/data/samtool_stats_res.txt"
    qc_values = f"tmp.samtool.qc_values.{pv}.json"
    metrics_zip = f"tmp.samtool.metrics.{pv}.zip"
    cmd = (
        f"parse-qc -n 'BAM Quality Metrics' "
        f"--metrics samtools_stats {samtools_stats} "
        f"--output-zip {metrics_zip} "
        f"--output-json {qc_values}"
    )
    os.system(cmd)

    assert os.path.exists(metrics_zip) == True
    assert os.path.exists(qc_values) == True

    qc_file = open(qc_values)
    data = json.load(qc_file)
    qc_file.close()

    # Test postprocessed metrics
    percentage_reads_paired = next(
        item
        for item in data["qc_values"]
        if item["derived_from"]
        == "samtools_stats_postprocessed:percentage_reads_paired"
    )
    assert percentage_reads_paired["value"] == 99.94050974320187

    # assert data["name"] == "BAM Quality Metrics"
    pprint.pprint(data["qc_values"])
    assert len(data["qc_values"]) == 24

    os.system(f"rm -f {qc_values} {metrics_zip}")


def test_picard_1():
    pv = platform.python_version()
    metrics = "./tests/data/picard_collectAlignmentSummaryMetrics.txt"
    qc_values = f"tmp.picard1.qc_values.{pv}.json"
    metrics_zip = f"tmp.picard1.metrics.{pv}.zip"
    cmd = (
        f"parse-qc -n 'picard_collectAlignmentSummaryMetrics' "
        f"--metrics picard_CollectAlignmentSummaryMetrics {metrics} "
        f"--output-zip {metrics_zip} "
        f"--output-json {qc_values}"
    )
    os.system(cmd)

    assert os.path.exists(metrics_zip) == True
    assert os.path.exists(qc_values) == True

    qc_file = open(qc_values)
    data = json.load(qc_file)
    qc_file.close()

    # assert data["name"] == "picard_collectAlignmentSummaryMetrics"
    assert len(data["qc_values"]) == 7

    os.system(f"rm -f {qc_values} {metrics_zip}")


def test_samtools_picard():
    pv = platform.python_version()
    samtools_stats = "./tests/data/samtool_stats_res.txt"
    metrics = "./tests/data/picard_collectAlignmentSummaryMetrics.txt"
    qc_values = f"tmp.picard1.qc_values.{pv}.json"
    metrics_zip = f"tmp.picard1.metrics.{pv}.zip"
    cmd = (
        f"parse-qc -n 'samtools picard' "
        f"--metrics picard_CollectAlignmentSummaryMetrics {metrics} "
        f"--metrics samtools_stats {samtools_stats} "
        f"--output-zip {metrics_zip} "
        f"--output-json {qc_values}"
    )
    os.system(cmd)

    assert os.path.exists(metrics_zip) == True
    assert os.path.exists(qc_values) == True

    qc_file = open(qc_values)
    data = json.load(qc_file)
    qc_file.close()

    # assert data["name"] == "samtools picard"
    assert len(data["qc_values"]) == 31

    os.system(f"rm -f {qc_values} {metrics_zip}")


def test_bamstats():
    pv = platform.python_version()
    bamstats_res = "./tests/data/bamstats_res.txt"
    qc_values = f"tmp.bamstats.qc_values.{pv}.json"
    metrics_zip = f"tmp.bamstats.metrics.{pv}.zip"
    cmd = (
        f"parse-qc -n 'bamstats' "
        f"--metrics bamstats {bamstats_res} "
        f"--output-zip {metrics_zip} "
        f"--output-json {qc_values}"
    )
    os.system(cmd)

    assert os.path.exists(metrics_zip) == True
    assert os.path.exists(qc_values) == True

    qc_file = open(qc_values)
    data = json.load(qc_file)
    qc_file.close()

    assert len(data["qc_values"]) == 6

    os.system(f"rm -f {qc_values} {metrics_zip}")


def test_nanoplot():
    pv = platform.python_version()
    metrics = "./tests/data/NanoPlot.txt"
    qc_values = f"tmp.nanoplot.qc_values.{pv}.json"
    metrics_zip = f"tmp.nanoplot.metrics.{pv}.zip"
    cmd = (
        f"parse-qc -n 'nanoplot' "
        f"--metrics nanoplot {metrics} "
        f"--output-zip {metrics_zip} "
        f"--output-json {qc_values}"
    )
    os.system(cmd)

    assert os.path.exists(metrics_zip) == True
    assert os.path.exists(qc_values) == True

    qc_file = open(qc_values)
    data = json.load(qc_file)
    qc_file.close()

    assert data["qc_values"][0]["key"] == "Mean Read Length [NanoPlot]"
    assert data["qc_values"][0]["value"] == 59364.0
    assert data["qc_values"][11]["key"] == "Reads >Q12 [NanoPlot]"
    assert data["qc_values"][11]["value"] == 86.8
    assert len(data["qc_values"]) == 13

    os.system(f"rm -f {qc_values} {metrics_zip}")


def test_verifybamid():
    pv = platform.python_version()
    metrics = "./tests/data/verifybamid2.out"
    qc_values = f"tmp.verifybamid2.qc_values.{pv}.json"
    metrics_zip = f"tmp.verifybamid2.metrics.{pv}.zip"
    cmd = (
        f"parse-qc -n 'verifybamid2' "
        f"--metrics verifybamid2 {metrics} "
        f"--output-zip {metrics_zip} "
        f"--output-json {qc_values}"
    )
    os.system(cmd)

    assert os.path.exists(metrics_zip) == True
    assert os.path.exists(qc_values) == True

    qc_file = open(qc_values)
    data = json.load(qc_file)
    qc_file.close()

    assert data["qc_values"][0]["key"] == "Estimate of Contamination [VerifyBamID2]"
    assert data["qc_values"][0]["value"] == 0.185126
    assert len(data["qc_values"]) == 1

    os.system(f"rm -f {qc_values} {metrics_zip}")


def test_kraken():
    pv = platform.python_version()
    metrics = "./tests/data/kraken2_report.txt"
    qc_values = f"tmp.kraken2.qc_values.{pv}.json"
    metrics_zip = f"tmp.kraken2.metrics.{pv}.zip"
    cmd = (
        f"parse-qc -n 'kraken2' "
        f"--metrics kraken2 {metrics} "
        f"--output-zip {metrics_zip} "
        f"--output-json {qc_values}"
    )
    os.system(cmd)

    assert os.path.exists(metrics_zip) == True
    assert os.path.exists(qc_values) == True

    qc_file = open(qc_values)
    data = json.load(qc_file)
    qc_file.close()

    assert len(data["qc_values"]) == 3
    assert data["qc_values"][0]["derived_from"] == "kraken2:taxonomic_id_9606"
    assert data["qc_values"][0]["value"] == 99.66
    assert data["qc_values"][1]["derived_from"] == "kraken2:taxonomic_id_2"
    assert data["qc_values"][1]["value"] == 0.00
    assert data["qc_values"][2]["derived_from"] == "kraken2:taxonomic_id_10239"
    assert data["qc_values"][2]["value"] == 0.00

    os.system(f"rm -f {qc_values} {metrics_zip}")


def test_rnaseqc():
    pv = platform.python_version()
    metrics = "./tests/data/RNA-SeQC.txt"
    qc_values = f"tmp.rnaseqc.qc_values.{pv}.json"
    metrics_zip = f"tmp.rnaseqc.metrics.{pv}.zip"
    cmd = (
        f"parse-qc -n 'rnaseqc' "
        f"--metrics rnaseqc {metrics} "
        f"--output-zip {metrics_zip} "
        f"--output-json {qc_values}"
    )
    os.system(cmd)

    assert os.path.exists(metrics_zip) == True
    assert os.path.exists(qc_values) == True

    qc_file = open(qc_values)
    data = json.load(qc_file)
    qc_file.close()

    assert data["qc_values"][0]["key"] == "Mapping Rate [RNA-SeQC]"
    assert data["qc_values"][0]["value"] == 0.937137
    assert data["qc_values"][16]["key"] == "Alternative Alignments [RNA-SeQC]"
    assert data["qc_values"][16]["value"] == 52179522
    assert len(data["qc_values"]) == 46

    os.system(f"rm -f {qc_values} {metrics_zip}")


def test_mosdepth():
    pv = platform.python_version()
    metrics = "./tests/data/mosdepth.summary.txt"
    qc_values = f"tmp.mosdepth.qc_values.{pv}.json"
    metrics_zip = f"tmp.mosdepth.metrics.{pv}.zip"
    cmd = (
        f"parse-qc -n 'mosdepth' "
        f"--metrics mosdepth {metrics} "
        f"--output-zip {metrics_zip} "
        f"--output-json {qc_values}"
    )
    os.system(cmd)

    assert os.path.exists(metrics_zip) == True
    assert os.path.exists(qc_values) == True

    qc_file = open(qc_values)
    data = json.load(qc_file)
    qc_file.close()

    assert data["qc_values"][0]["key"] == "Estimated Average Coverage [mosdepth]"
    assert data["qc_values"][0]["value"] == 58.57
    assert len(data["qc_values"]) == 1

    os.system(f"rm -f {qc_values} {metrics_zip}")