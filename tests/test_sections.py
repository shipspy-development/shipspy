import numpy as np
import xarray as xr
import os


def test_section():

    ds = xr.open_dataset("tests/dummydataset.nc")

    os.system("shippy sections -i tests/dummydataset.nc -o tests/testsection_output.nc -s tests/dummysections.txt")
    
    out = xr.open_dataset("tests/testsection_output.nc")
    
    assert set(ds.data_vars) == set(out.data_vars)
    assert "section" in out.coords
    assert set(ds.coords) == set(out.coords) - {'section'}
    
    os.system("rm tests/testsection_output.nc")
