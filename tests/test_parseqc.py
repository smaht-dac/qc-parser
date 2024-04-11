import pytest
import os, platform
import json


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

    assert len(data["qc_values"]) == 11

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

    #assert data["name"] == "BAM Quality Metrics"
    assert len(data["qc_values"]) == 18

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

    #assert data["name"] == "picard_collectAlignmentSummaryMetrics"
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

    #assert data["name"] == "samtools picard"
    assert len(data["qc_values"]) == 25

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

    assert data["qc_values"][0]["key"] == 'Mean Read Length'
    assert data["qc_values"][0]["value"] == 59364.0
    assert data["qc_values"][11]["key"] == 'Reads >Q12'
    assert data["qc_values"][11]["value"] == 86.8
    assert len(data["qc_values"]) == 13

    os.system(f"rm -f {qc_values} {metrics_zip}")
