import numpy as np
import xarray as xr
import os


def test_section():
    ds = xr.open_dataset("tests/dummydataset.nc")

    os.system(
        "shipspy sections -i tests/dummydataset.nc -o tests/testsection_output_rn.nc -s tests/dummysections.txt"
    )

    out = xr.open_dataset("tests/testsection_output_rn.nc")

    assert set(ds.data_vars) == set(out.data_vars)
    assert "section" in out.coords
    assert set(ds.coords) == set(out.coords) - {"section"}
    assert ds.time[0] < out.time[0]

    os.system("rm tests/testsection_output_rn.nc")
