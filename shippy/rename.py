import numpy as np
import xarray as xr
import pint_xarray
import pandas as pd
import yaml
import os

from metpy.units import units
from datetime import datetime as dt


def configure_rename_parser(parser):
    parser.add_argument(
        "-i",
        "--inputfile",
        metavar="INPUT_FILE",
        help="Postprocessed input file",
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
        help="Dictionary with variable names and attributes (yaml)",
        default=None,
        required=True,
    )

    parser.add_argument(
        "-d",
        "--instrument",
        metavar="INSTRUMENT_NAME",
        help="Instrument/device name, options are: calitoo, microtops, ctd, radiosondes, uav, hatpro, ceilometer, test",
        required=True,
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


def import_data(filename, instrument):
    """
    Import data and convert to xarray dataset
    """
    if instrument == "calitoo":
        ds = xr.open_dataset(filename, group="Calitoo", decode_times=False)
    elif instrument == "microtops":
        series = pd.read_csv(filename, delimiter=",", header=4)
        ds = series.to_xarray()
    elif instrument == "dusttrak":
        ds = xr.open_dataset(filename, group="DustTrak", decode_times=False)
    elif instrument == "radiosondes":
        if os.path.isfile(filename):
            ds = xr.open_dataset(filename)
        else:
            ds = xr.open_mfdataset(f"{filename}/*.nc")
    else:
        ds = xr.open_dataset(filename)
    return ds


def fix_times(ds, instrument):
    """
    Fix timestamps as datetime64[ns] objects
    """
    if instrument in ["calitoo", "dusttrak"]:
        timestamps = ds.time * np.timedelta64(1, "s").astype(
            "timedelta64[ns]"
        ) + np.datetime64("1970-01-01T00:00:00").astype("datetime64[ns]")
        ds["time"] = timestamps
    elif instrument == "microtops":
        dates = ds["Date(dd:mm:yyyy)"].values
        times = ds["Time(hh:mm:ss)"].values
        timestamps = np.array([], dtype="datetime64[s]")
        for d, t in zip(dates, times):
            tstamp = np.array(
                [np.datetime64(dt.strptime(d + " " + t, "%d:%m:%Y %H:%M:%S"), "s")]
            )
            timestamps = np.append(timestamps, tstamp)
        timestamps = timestamps.astype("datetime64[ns]")
        ds = ds.assign_coords(time=("index", timestamps))
        ds = ds.swap_dims({"index": "time"})
        ds = ds.drop_vars(["Date(dd:mm:yyyy)", "Time(hh:mm:ss)", "index"])
    return ds


def keep_old_attrs(ds):
    old_attrs = {}
    for v in list(ds.data_vars):
        old_attrs[v] = ds[v].attrs
    for v in list(ds.coords):
        old_attrs[v] = ds[v].attrs
    return old_attrs


def get_units_old(ds, instrument):
    """
    Returns dict of the units of the input file
    """
    units_old = {}
    if instrument == "calitoo":
        for v in list(ds.keys()):
            if ds[v].attrs["units"] == "-":
                units_old[v] = "dimensionless"
            elif ds[v].attrs["units"] == "Dobson Units":
                units_old[v] = "0.44615 mmol/m^2"
            elif ds[v].attrs["units"] == "degrees C":
                units_old[v] = "degC"
            else:
                units_old[v] = ds[v].attrs["units"]
    elif instrument == "microtops":
        units_old = {}
        for v in list(ds.keys()):
            units_old[v] = "dimensionless"
        units_old["Latitude"] = "degree_north"
        units_old["Longitude"] = "degree_east"
        units_old["Water Vapor(cm)"] = "cm"
        units_old["STD_Water_Vapor(cm)"] = "cm"
    elif instrument == "ctd":
        for v in list(ds.keys()):
            unit = ds[v].long_name[
                ds[v].long_name.find("[") + 1 : ds[v].long_name.find("]")
            ]
            unit = ds[v].long_name[
                ds[v].long_name.find("[") + 1 : ds[v].long_name.find("]")
            ]
            if not unit.find(",") == -1:
                unit = unit[unit.find(",") + 2 :]
            if unit == "deg C":
                unit = "degC"
            elif unit == "db":
                unit = "dbar"
            elif unit == "PSU":
                unit = "g/kg"
            elif unit == "% saturation":
                unit = "percent"
            elif unit == "NTU":
                unit = "dimensionless"
            units_old[v] = unit
        units_old["flag"] = "dimensionless"
    elif instrument == "radiosondes":
        for v in list(ds.keys()):
            try:
                unit = ds[v].attrs["units"]
                if unit == "1":
                    unit = "dimensionless"
            except KeyError:
                unit = "dimensionless"
                if v == "alt_bnds":
                    unit = "m"
            units_old[v] = unit
    elif instrument == "uav":
        for v in list(ds.keys()):
            try:
                unit = ds[v].attrs["units"]
                if unit == "1":
                    unit = "dimensionless"
            except KeyError:
                unit = "dimensionless"
            units_old[v] = unit
    elif instrument == "ceilometer":
        for v in list(ds.keys()):
            try:
                unit = ds[v].attrs["units"]
                if unit == "1":
                    unit = "dimensionless"
                units_old[v] = unit
            except KeyError:
                continue
        units_old["sci"] = "dimensionless"
    else:
        for v in list(ds.keys()):
            try:
                unit = ds[v].attrs["units"]
                if unit == "1":
                    unit = "dimensionless"
                if unit == "%":
                    unit = "percent"
                units_old[v] = unit
            except KeyError:
                continue
    return units_old


def get_units_new(variables_dict):
    """
    Returns dict of units of the final dataset
    """
    units_new = {}
    for v in list(variables_dict.keys()):
        try:
            units_new[v] = variables_dict[v]["units"]
        except KeyError:
            continue
    return units_new


def fix_units(ds, units_old, units_new):
    """
    Convert old to new units
    """
    for v in list(ds.keys()):
        try:
            u = units(units_new[v])
            if units_old[v] == "1/8":
                ds[v] = (ds[v] * units("dimensionless") / 8).pint.to(str(u))
            else:
                ds[v] = (ds[v] * units(units_old[v])).pint.to(str(u))
        except KeyError:
            ds = ds.drop_vars([v])
        ds = ds.pint.dequantify()
    return ds


def add_extras(ds, instrument):
    """
    Add extra coordinates like identifiers or altitude = -depth for CTD
    """
    if instrument == "ctd":
        ds = ds.assign_coords(
            alt=("DEPTH", -ds.DEPTH.values)
        )  # assign altitude as auxiliary coordinate
    if instrument == "radiosondes":
        ds = ds.sortby("launch_time").swap_dims({"sounding": "launch_time"})

        identifiers = []
        number = 0
        for t in ds.launch_time.values:
            dz = ds.sel(launch_time=t).dz.mean().values
            if dz > 0:
                number = number + 1
                ident = f"RS{int(number):03}_up"
            elif dz < 0:
                if dz * dz_old > 0:
                    number = number + 1
                ident = f"RS{int(number):03}_down"
            identifiers.append(ident)
            dz_old = dz

        ds = ds.assign_coords(identifier=("launch_time", np.array(identifiers)))
    if instrument == "uav":
        identifiers = []
        for number in np.arange(len(ds.start_time)):
            ident = f"UAV{int(number):03}"
            identifiers.append(ident)
        ds = ds.assign_coords(identifier=("start_time", np.array(identifiers)))
    return ds


def put_attrs(ds, variable_dict, old_attrs, instrument):
    """
    Add attributes to variables
    """
    varname_dic = {
        key: variable_dict[key]["varname"] for key in list(variable_dict.keys())
    }
    for v in variable_dict.keys():
        attrs = variable_dict[v]
        try:
            if instrument == "uav":
                del old_attrs[v]["standard_name"]
            attrs = old_attrs[v] | attrs
        except KeyError:
            continue
        del attrs["varname"]
        ds[v].attrs = attrs
    ds = ds.rename(varname_dic)
    return ds


def fix(ds, instrument):
    """
    If needed: fixes and removing unphysical values
    """
    if instrument == "microtops":
        ds["angstrom_exp"] = ds.angstrom_exp.where(
            ds.angstrom_exp > 0, other=np.nan
        )  # remove unphysical (negative) values
    if instrument in ["uav", "hatpro", "radiosondes"]:
        ds["rh"] = ds.rh.where(ds.rh < 1.0, other=np.nan)  # remove unphysical rh > 1
        ds["rh"] = ds.rh.where(
            ds.rh >= 0.0, other=np.nan
        )  # remove unphysical nergative rh
    return ds


def clean_up(ds, instrument):
    """
    Final cleaning up
    Sets latitude and longitude as coordinates
    """
    ds = ds.set_coords(["lat", "lon"])
    if instrument == "radiosondes":
        ds = ds.drop_vars(["sounding"])
    return ds


def save_nc(ds, outfile, instrument):
    """
    Save dataset to netCDF file with time encoding
    """
    if instrument == "ctd":
        for v in list(ds.coords):
            ds[v].encoding = {}
        for v in list(ds.data_vars):
            ds[v].encoding = {}
        ds.to_netcdf(outfile)
    else:
        ds.to_netcdf(outfile)


def run(args):
    attrs_file = args.attributes
    infilename = args.inputfile
    instrument = args.instrument
    outfile = args.outputfile

    with open(attrs_file, "r") as stream:
        variables_dict = yaml.safe_load(stream)

    ds = import_data(infilename, instrument)
    ds = fix_times(ds, instrument)

    old_attrs = keep_old_attrs(ds)
    units_old = get_units_old(ds, instrument)
    units_new = get_units_new(variables_dict)

    ds = fix_units(ds, units_old, units_new)
    ds = add_extras(ds, instrument)
    ds = put_attrs(ds, variables_dict, old_attrs, instrument)
    ds = fix(ds, instrument)
    ds = clean_up(ds, instrument)

    save_nc(ds, outfile, instrument)
