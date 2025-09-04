import numpy as np
import xarray as xr
import pint_xarray
import pandas as pd
import yaml
import re
from scipy.stats import circmean
from metpy.units import units
import glob


def configure_pamos_parser(parser):
    parser.add_argument(
        "-i",
        "--inputdir",
        metavar="INPUT_DIR",
        help="Directory with input raw files.",
        default=None,
        required=True,
    )

    parser.add_argument(
        "-o",
        "--outputfile",
        metavar="OUTPUT_FILE",
        help="Output file (netCDF)",
        default=None,
        required=True,
    )

    parser.add_argument(
        "-a",
        "--attributes",
        metavar="ATTRIBUTE_DICT",
        help="Dictionary with variable names and attributes (yaml), some are provided in folder rv_information",
        default=None,
        required=True,
    )

    parser.add_argument(
        "-c",
        "--columns",
        metavar="Header for raw data files",
        help="The raw data files need a header file which is given here.",
        required=True,
    )

    parser.add_argument(
        "-f",
        "--flags",
        metavar="quality_flags",
        help="Dictionary for quality flags.",
        required=False,
    )

    parser.add_argument(
        "-e",
        "--extra",
        metavar="additonal_attributes",
        help="Additional attributes for added variables.",
        required=False,
    )

    parser.add_argument(
        "-v",
        "--verbose",
        metavar="DEBUG",
        help='Set the level of verbosity [DEBUG, INFO," " WARNING, ERROR]',
        required=False,
        default="INFO",
    )

    parser.set_defaults(func=run)


def CircMean(data):
    return circmean(data, high=360, low=0)


def circ_mean(da):
    return xr.apply_ufunc(
        CircMean, da, input_core_dims=[["time"]], vectorize=True, dask="parallelized"
    )


def min_means(ds, var):
    """
    Calculate minutely means
    For wind speed and heading, circular means are calculated
    """
    attrs = ds[var].attrs
    if var in list(["wdir_min", "wdir_max", "wdir_mean", "ship_heading"]):
        var_mean = (
            circ_mean(
                ds[var].where(ds[var] >= 0, other=np.nan)
                #                .dropna("time")
                .resample(time="min")
            )
            .rename({"__resample_dim__": "time"})
            .to_dataset()
        )
    elif var in [
        "quality_position",
        "quality_meteorology",
        "quality_aerosols",
        "quality_tracegases",
        "quality_housekeeping",
        "precip_acc",
    ]:
        var_mean = ds[[var]].resample(time="min").max()
    elif var in ["pump_flag"]:
        var_mean = ds[[var]].resample(time="min").max()
    else:
        var_mean = ds[[var]].resample(time="min").mean()

    var_mean[var].attrs = attrs
    return var_mean


def clean_dataarray_to_float(da: xr.DataArray) -> xr.DataArray:
    """
    Converts a DataArray with dtype=object or mixed types to float,
    replacing any non-numeric values with NaN.

    Parameters:
        da (xr.DataArray): The input DataArray (possibly object dtype)

    Returns:
        xr.DataArray: A float DataArray with invalid values replaced by NaN
    """
    return xr.apply_ufunc(
        pd.to_numeric,
        da,
        kwargs={"errors": "coerce"},
        output_dtypes=[float],
        dask="allowed",
    )


def fix_unit(ds, dict, v):
    if not v in ["lat_or", "lon_or"]:
        unit_old = dict[v]["units_old"]
        unit_new = dict[v]["units"]
        ds[v] = (clean_dataarray_to_float(ds[v]) * units(unit_old)).pint.to(unit_new)
    return ds


def add_attrs(ds, dict, v):
    attrs = dict[v]
    del attrs["varname"]
    del attrs["units_old"]
    ds[v].attrs = attrs
    return ds


def add_extra_attrs(ds, dict, v):
    attrs = dict[v]
    ds[v].attrs = attrs
    return ds


def get_input_data(header_file, indir):

    header_ls = pd.read_csv(header_file, delimiter=";", header=None)
    header_vals = header_ls.values[0]
    header = []
    num = 0
    for h in header_vals:
        if type(h) == str:
            header.append(h)
        else:
            header.append(num)
            num = num + 1

    files = glob.glob(f"{indir}/*.emb*")

    f = files[0]
    data = pd.read_csv(f, delimiter=";", header=None, names=header)
    data = data.set_index("timestamp").to_xarray()

    for f in files[1:]:
        next = pd.read_csv(f, delimiter=";", header=None, names=header)
        next = next.set_index("timestamp").to_xarray()
        data = xr.concat([data, next], dim="timestamp")

    data = data.sortby("timestamp")
    datetime = np.datetime64("1970-01-01T00:00") + data[
        "timestamp"
    ].values * np.timedelta64(1, "s")
    data = data.assign_coords(time=("timestamp", datetime))
    data = data.swap_dims({"timestamp": "time"})

    return data


def get_pumpflag(ds):
    flags = np.array(ds["plc_reg10_11"].values, dtype="int")
    pump_flag = []
    for i in flags:
        flag = int(bin(i)[-1])
        pump_flag.append(flag)
    pump_flag = np.array(pump_flag)
    return pump_flag


def get_eraser(ds):
    pump_flag = get_pumpflag(ds)
    eraser = []
    for i in np.arange(len(pump_flag)):
        if pump_flag[i] == 1:
            eraser.append(1.0)
        else:
            eraser.append(np.nan)

    eraser = np.array(eraser)
    return eraser


def erasing(ds, v, eraser):
    if v in [
        "uv_bc1",
        "blue_bc1",
        "green_bc1",
        "red_bc1",
        "ir_bc1",
        "opcr2_pm10",
        "opcr2_pm2p5",
        "opcr2_pm1",
        "number",
        "diameter",
        "ldsa",
        "co2_conc",
        "co_conc",
        "ch4_conc",
        "so2_conc",
    ]:
        ds[v] = (clean_dataarray_to_float(ds[v])) * eraser
    return ds


def repair_position(ds):
    signs = ((ds.lat_or == "N").values.astype(int) - 0.5) * 2
    vals = pd.to_numeric(ds.lat.values, errors="coerce")
    lat = vals * signs
    signs = ((ds.lon_or == "E").values.astype(int) - 0.5) * 2
    vals = pd.to_numeric(ds.lon.values, errors="coerce")
    lon = vals * signs
    ds = ds.assign(lat=("time", lat))
    ds = ds.assign(lon=("time", lon))
    return ds


def quality_flag(ds, vars, flag_dict, deltat):
    flag = np.zeros(len(ds.time))
    for v in vars:
        range_max = flag_dict[v]["max"]
        range_min = flag_dict[v]["min"]
        var = flag_dict[v]["var"]
        too_large = (np.argwhere(ds[v].values > range_max)).flatten()
        too_small = (np.argwhere(ds[v].values < range_min)).flatten()
        too_noisy = (
            np.argwhere(ds[v].rolling(time=deltat, center=True).var().values > var)
        ).flatten()
        nans = (np.argwhere(np.isnan(ds[v].values))).flatten()
        alert = np.concatenate([too_large, too_small, too_noisy, nans])
        alert = np.unique(np.sort(alert))
        flag[alert] = 1
    return flag


def replace_non_floats_with_nan(value):
    return value if isinstance(value, float) else np.nan


def run(args):
    indir = args.inputdir
    outname = args.outputfile
    attrs_dict_name = args.attributes
    with open(attrs_dict_name, "r") as stream:
        vardict = yaml.safe_load(stream)
    varlist = list(vardict.keys())
    header = args.columns

    ds = get_input_data(header, indir)

    pump_flag = get_pumpflag(ds)
    ds = ds.assign(pump_flag=("time", pump_flag))
    eraser = get_eraser(ds)

    ds = repair_position(ds)

    ds = ds[varlist + ["pump_flag"]]
    rename_dict = {}
    for v in varlist:
        if not v in ["lat_or", "lon_or"]:
            rename_dict[v] = vardict[v]["varname"]
            ds = erasing(ds, v, eraser)
            ds = fix_unit(ds, vardict, v)
            ds = add_attrs(ds, vardict, v)
            minutely_data = min_means(ds, v)
            if v == list(ds.keys())[0]:
                pamos = minutely_data
            else:
                pamos = xr.merge([pamos, minutely_data])
    pamos = pamos.rename(rename_dict)
    minutely_data = min_means(ds, "pump_flag")
    pamos = pamos.assign(pump_flag=("time", minutely_data.pump_flag.values))

    precip_rate = pamos["precip_acc"].resample(time="10min").max().diff("time")
    precip_rate = precip_rate.resample(time="1min").interpolate("nearest")
    pamos = pamos.assign(precip_rate=precip_rate / 600)  # pro 10 min in 1 s

    try:
        flags = args.flags
        with open(flags, "r") as stream:
            flag_checks = yaml.safe_load(stream)
        tmean = 10  # in min
        meteo_vars = ["t_air", "p_air", "rh", "wspd_mean", "wdir_mean"]
        meteo_flag = quality_flag(pamos, meteo_vars, flag_checks, tmean)
        aero_vars = [
            "number_conc",
            "pm10",
            "pm2p5",
            "pm1",
            "bc_375nm",
            "bc_470nm",
            "bc_528nm",
            "bc_625nm",
            "bc_880nm",
        ]
        aero_flag = quality_flag(pamos, aero_vars, flag_checks, tmean)
        gas_vars = ["co2_conc", "co_conc", "ch4_conc", "so2_conc"]
        gas_flag = quality_flag(pamos, gas_vars, flag_checks, tmean)
        position_vars = ["lat", "lon"]
        pos_flag = quality_flag(pamos, position_vars, flag_checks, tmean)
        housekeeping_vars = ["t_cabinet"]
        housekeeping_flag = quality_flag(pamos, housekeeping_vars, flag_checks, tmean)

        pamos = pamos.assign(quality_position=("time", pos_flag))
        pamos = pamos.assign(quality_meteorology=("time", meteo_flag))
        pamos = pamos.assign(quality_aerosols=("time", aero_flag))
        pamos = pamos.assign(quality_tracegases=("time", gas_flag))
        pamos = pamos.assign(quality_housekeeping=("time", housekeeping_flag))

    except:
        print("No quality flag dictionary available. Skip quality flags.")

    try:
        extra_attrs = args.extra
        with open(extra_attrs, "r") as stream:
            additions = yaml.safe_load(stream)
        for v in additions.keys():
            pamos = add_extra_attrs(pamos, additions, v)
    except:
        print("No additional attributes available")

    var_list = list(pamos.keys())
    var_list = sorted(var_list)
    pamos = pamos[var_list]
    pamos = pamos.set_coords(["lat", "lon"])
    pamos.time.attrs = {"time zone": "UTC"}
    pamos = pamos.pint.dequantify()

    pamos.to_netcdf(outname)
