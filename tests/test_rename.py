import numpy as np
import xarray as xr
import pint_xarray
import os

from metpy.units import units
import numpy.testing as npt


def test_rename():
    ds = xr.open_dataset("tests/dummydataset.nc")

    os.system(
        "shipspy rename -i tests/dummydataset.nc -o tests/testrename_output.nc -a tests/dummyvariables.yaml -d test"
    )

    out = xr.open_dataset("tests/testrename_output.nc")

    assert "lat" in list(out.coords)
    assert "lon" in list(out.coords)
    npt.assert_array_almost_equal(
        (ds.t * units(ds.t.attrs["units"])).pint.to("K").values, out.t_air.values
    )
    npt.assert_array_almost_equal(
        (ds.p * units(ds.p.attrs["units"])).pint.to("Pa").values, out.p_air.values
    )

    os.system("rm tests/testrename_output.nc")
