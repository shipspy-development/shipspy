import xarray as xr
import os
import numpy as np


def test_pamos():
    os.system(
        "shipspy pamos -i tests/dummypamosdir -o tests/testpamos_output.nc -a tests/dummypamosvariables.yaml -c tests/dummypamosheader"
    )

    out = xr.open_dataset("tests/testpamos_output.nc")

    assert (out.time.values[1:] - out.time.values[:-1] == np.timedelta64(1, "m")).all()
    assert out.t_air.values[0] == 274.15

    assert "lat" in list(out.coords)
    assert "lon" in list(out.coords)

    assert out.pump_flag.values[0] == 0
    assert np.isnan(out.number_conc.values)[0] == True
    assert out.pump_flag.values[-1] == 1
    assert np.isnan(out.number_conc.values)[-1] == False

    os.system("rm tests/testpamos_output.nc")
