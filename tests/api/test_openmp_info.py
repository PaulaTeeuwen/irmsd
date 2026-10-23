"""Build-time OpenMP state of the compiled library.

Release wheels are built with IRMSD_REQUIRE_OPENMP=ON and tested with
IRMSD_EXPECT_OPENMP=1 (see [tool.cibuildwheel] in pyproject.toml), so a wheel
whose OpenMP runtime went missing fails here instead of shipping serial.
"""

import json
import os
import subprocess
import sys

import pytest

from irmsd import openmp_info


def test_openmp_info_values():
    info = openmp_info()
    assert set(info) == {"enabled", "max_threads"}
    assert isinstance(info["enabled"], bool)
    assert info["max_threads"] >= 1
    if not info["enabled"]:
        assert info["max_threads"] == 1


def test_openmp_expected_in_wheels():
    if os.environ.get("IRMSD_EXPECT_OPENMP") != "1":
        pytest.skip("OpenMP only required where IRMSD_EXPECT_OPENMP=1 (wheel CI)")
    assert openmp_info()["enabled"], "library was built without OpenMP"


@pytest.mark.skipif(not openmp_info()["enabled"], reason="library built without OpenMP")
def test_openmp_num_threads_is_respected():
    script = "import json; from irmsd import openmp_info; print(json.dumps(openmp_info()))"
    env = dict(os.environ, OMP_NUM_THREADS="2")
    out = subprocess.run(
        [sys.executable, "-c", script], env=env, capture_output=True, text=True, check=True
    )
    assert json.loads(out.stdout.strip().splitlines()[-1])["max_threads"] == 2
