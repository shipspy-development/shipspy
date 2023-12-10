import xarray as xr
import os
import numpy as np
import numpy.testing as npt


def test_dship():
    os.system(
        "shipspy dship -i tests/dummydshipdata.txt -o tests/testdship_output.nc -a tests/dummyvariables.yaml -s Test"
    )

    out = xr.open_dataset("tests/testdship_output.nc")

    assert (out.time.values[1:] - out.time.values[:-1] == np.timedelta64(1, "m")).all()
    assert out.t_air.values[0] == 274.15
    npt.assert_almost_equal(
        (np.linspace(10, 16, 61)[:-1]).mean() + 273.15, out.t_air.values[-2]
    )
    assert "lat" in list(out.coords)
    assert "lon" in list(out.coords)

    os.system("rm tests/testdship_output.nc")
