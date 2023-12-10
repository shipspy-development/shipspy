import numpy as np
import xarray as xr
import pint_xarray
import pandas as pd
import yaml
import re
from scipy.stats import circmean
from metpy.units import units


def configure_dship_parser(parser):
    parser.add_argument(
        "-i",
        "--inputfile",
        metavar="INPUT_FILE",
        help="Postprocessed file (netCDF)",
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
        "-s",
        "--ship",
        metavar="Merian, Sonne, or Meteor",
        help="Ship name where the data comes from. Options are Merian, Sonne, Meteor, Test.",
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
    if var in list(["wdir", "ship_heading"]):
        var_mean = (
            circ_mean(
                ds[var]
                .where(ds[var] > 0, other=np.nan)
                .dropna("time")
                .resample(time="min")
            )
            .rename({"__resample_dim__": "time"})
            .to_dataset()
        )
    else:
        var_mean = ds[[var]].resample(time="min").mean()
    return var_mean


def min_std(ds, var):
    """
    Calculate minutely standard deviations
    """
    var_std = ds[[var]].resample(time="min").std()
    return var_std


def remove_fill_vals(ds, ship, var):
    """
    Remove fill values. For Sonne, fill value is 9.
    Also remove unphysical values.
    """
    if ship in ["Sonne", "Meteor", "Test"]:
        ds[var] = ds[var].where(ds[var] != 9.0, other=np.nan)
    if (ship == "Merian") & (var == "sst_7m"):
        ds[var] = ds[var].where(ds[var] != 0, other=np.nan)
    if var in ["wave_height", "wave_length", "wave_period"]:
        ds[var] = ds[var].where(ds[var] > 0, other=np.nan)
    if var == "wave_dir":
        ds[var] = ds[var].where((ds[var] > -0.01) & (ds[var] < 360.1), other=np.nan)
    return ds


def fix(ds, ship, var):
    """
    Fixes on variables
    """
    if var in ["t_air", "sst_2m", "sst_7m", "sst_4m", "sst_3m"]:
        ds[var] = ds[var].pint.to("kelvin")  # convert deg Celcius to K
    if var == "p_air":
        ds[var] = ds[var].pint.to("pascal")  # convert hPa to Pa
    if var == "rh":
        ds[var] = ds[var].pint.to("1")  # convert % to 1
    if var == "salinity":
        ds[var] = ds[var].pint.to("1")  # convert g/kg to 1
    if var == "chlorophyll_a":
        ds[var] = ds[var].pint.to("kg/m^3")  # convert 10^-6kg/m^3 to kg/m^3
    if var == "ship_speed":
        ds[var] = ds[var].pint.to("m/s")  # convert to m/s
    ds = ds.pint.dequantify()
    return ds


def var_ordering(ship):
    if ship == "Merian":
        ordering = [
            "t_air",
            "p_air",
            "wspd",
            "wdir",
            "rh",
            "lwr",
            "swr",
            "sea_floor_depth",
            "sst_2m",
            "sst_7m",
            "conductivity",
            "salinity",
            "chlorophyll_a",
            "current_speed",
            "current_dir",
            "wave_period",
            "wave_length",
            "wave_height",
            "wave_dir",
            "ship_speed",
            "ship_heading",
            "ship_heave",
            "ship_heave_std",
            "ship_pitch",
            "ship_pitch_std",
            "ship_roll",
            "ship_roll_std",
        ]
    elif ship == "Sonne":
        ordering = [
            "t_air",
            "p_air",
            "wspd",
            "wdir",
            "rh",
            "lwr",
            "swr",
            "sea_floor_depth",
            "sst_2m",
            "sst_4m",
            "conductivity",
            "salinity",
            "chlorophyll",
            "ship_speed",
            "ship_heading",
            "ship_heave",
            "ship_heave_std",
            "ship_pitch",
            "ship_pitch_std",
            "ship_roll",
            "ship_roll_std",
        ]
    elif ship == "Meteor":
        ordering = [
            "lwr",
            "swr",
            "sea_floor_depth",
            "ship_speed",
            "ship_heading",
            "ship_heave",
            "ship_heave_std",
            "ship_pitch",
            "ship_pitch_std",
            "ship_roll",
            "ship_roll_std",
        ]
    elif ship == "Test":
        ordering = ["t_air", "p_air"]
    return ordering


def run(args):
    fn = args.inputfile
    outname = args.outputfile
    ship_name = args.ship
    attrs_dict_name = args.attributes
    with open(attrs_dict_name, "r") as stream:
        rename_dict = yaml.safe_load(stream)
    varname_dict = {
        key: rename_dict[key]["varname"] for key in list(rename_dict.keys())
    }
    varname_swap = {v: k for k, v in varname_dict.items()}
    std_dict = {
        "ship_heave": {
            "varname": "ship_heave_std",
            "long_name": "standard deviation of ship heave",
            "units": "m",
            "instrument": "Seapath (DShip)",
            "cell_method": "minutely standard deviation",
        },
        "ship_roll": {
            "varname": "ship_roll_std",
            "long_name": "standard deviation of ship roll",
            "units": "degree",
            "instrument": "Seapath (DShip)",
            "cell_method": "minutely standard deviation",
        },
        "ship_pitch": {
            "varname": "ship_pitch_std",
            "long_name": "standard deviation of ship pitch",
            "units": "degree",
            "instrument": "Seapath (DShip)",
            "cell_method": "minutely standard deviation",
        },
    }

    def get_attrs(ds, var):
        attrs = rename_dict[varname_swap[var]]
        del attrs["varname"]
        ds[var].attrs = attrs
        return ds

    sw_dir = lambda ew_str: -1 if re.match(".*[SW].*", ew_str) else 1

    def format_loc(loc_string):
        """Turn string of location into float"""

        if type(loc_string) == str:
            deg_loc, min_loc, dir_loc = re.split("°|'", loc_string)
            return (int(deg_loc) + float(min_loc) / 60.0) * sw_dir(dir_loc)
        else:
            return np.nan

    dship_data = pd.read_csv(fn, delimiter=";", encoding="ISO-8859-1", low_memory=False)
    units_dict = dict(dship_data.loc[1])
    dship_data = dship_data.drop(axis=0, index=[0, 1])
    dship_data = dship_data.rename(columns=varname_dict)
    if ship_name == "Merian":
        dship_data.lat, dship_data.lon = dship_data.lat.apply(
            format_loc
        ), dship_data.lon.apply(format_loc)
    dship_data.where(dship_data != "NODATA", np.nan, inplace=True)
    dship_data.where(dship_data != "np.nan", np.nan, inplace=True)
    secs_since_1970 = dship_data["seconds since 1970"].values
    timestamps = np.datetime64(
        "1970-01-01T00:00:00"
    ) + secs_since_1970 * np.timedelta64(1, "s")
    timestamps = timestamps.astype("datetime64[ns]")
    var_list_full = list(dship_data.keys())
    for var in var_list_full:
        if var not in list(varname_dict.values()):
            dship_data = dship_data.drop(var, axis=1)
    dship_data = dship_data.astype("float", errors="ignore")
    dship_data["time"] = timestamps
    dship_sec = dship_data.set_index("time").to_xarray()

    if "salinity" in varname_swap:
        if units_dict[varname_swap["salinity"]] == "PSU":
            units_dict[varname_swap["salinity"]] = "g/kg"
    if "ship_speed" in varname_swap:
        if units_dict[varname_swap["ship_speed"]] == "kn":
            units_dict[varname_swap["ship_speed"]] = "knot"

    for v in dship_sec.keys():
        dship_sec[v].attrs["units"] = units_dict[varname_swap[v]]

    for v in varname_swap.keys():
        dship_sec = remove_fill_vals(dship_sec, ship_name, v)

    dship_sec = dship_sec.pint.quantify(unit_registry=units)
    for v in varname_swap.keys():
        minutely_data = min_means(dship_sec, v)
        minutely_data = fix(minutely_data, ship_name, v)
        minutely_data = get_attrs(minutely_data, v)
        if v == list(varname_swap.keys())[0]:
            dship = minutely_data
        else:
            dship = xr.merge([dship, minutely_data])

    for v in list(["ship_heave", "ship_roll", "ship_pitch"]):
        if ship_name == "Test":
            break
        stdname = std_dict[v]["varname"]
        prepare_data = remove_fill_vals(dship_sec, ship_name, v)
        minutely_data = min_std(prepare_data, v)
        minutely_data = minutely_data.rename({v: stdname})
        attrs = std_dict[v]
        del attrs["varname"]
        minutely_data[stdname].attrs = attrs
        dship = xr.merge([dship, minutely_data])

    dship = dship.set_coords(["lat", "lon"])
    dship = dship[var_ordering(ship_name)]
    dship = dship.pint.dequantify()
    dship.to_netcdf(outname)
