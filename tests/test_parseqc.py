import pytest
import os, platform
import json


def test_samtools():
    pv = platform.python_version()
    samtools_stats = "./tests/data/samtool_stats_res.txt"
    qc_values = f"tmp.samtool.qc_values.{pv}.json"
    metrics_zip = f"tmp.samtool.metrics.{pv}.zip"
    cmd = (
        f"parse-qc -n 'BAM Quality Metrics' "
        f"--metrics samtools {samtools_stats} "
        f"--output-zip {metrics_zip} "
        f"--output-json {qc_values}"
    )
    os.system(cmd)
    
    assert os.path.exists(metrics_zip) == True
    assert os.path.exists(qc_values) == True

    qc_file = open(qc_values)
    data = json.load(qc_file)
    qc_file.close()

    assert data["name"] == "BAM Quality Metrics"
    assert len(data["qc_values"]) == 15

    os.system(f"rm -f {qc_values} {metrics_zip}")
