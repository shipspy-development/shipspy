import numpy as np
import xarray as xr
import pandas as pd
import yaml


def configure_sections_parser(parser):
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
        "-s",
        "--sectionfile",
        metavar="SECTION_FILE",
        help="txt file with campaign and section information",
        default=None,
        required=True,
    )

    parser.add_argument(
        "-t",
        "--time_dimensioname",
        metavar="start_time",
        help="Name of dimension along which section is assigned",
        default="time",
        required=False,
    )

    parser.add_argument(
        "-a",
        "--attributes",
        metavar="GLOBAL_ATTRIBUTES",
        help="Dictionary with global attributes (yaml)",
        required=False,
    )

    parser.add_argument(
        "-v",
        "--verbose",
        metavar="DEBUG",
        help="Set the level of verbosity [DEBUG, INFO," " WARNING, ERROR]",
        required=False,
        default="INFO",
    )

    parser.set_defaults(func=run)


def run(args):
    infilename = args.inputfile
    ds = xr.open_dataset(infilename)
    outfilename = args.outputfile
    secfilename = args.sectionfile
    dim_name = args.time_dimensioname

    complete_period = pd.read_csv(secfilename, names=["dates"])
    sct_header = np.argwhere(
        complete_period.dates.values == "section_number;section_cuts"
    ).flatten()[0]
    cutout_header = np.argwhere(
        complete_period.dates.values == "cutout_start;cutout_end"
    ).flatten()[0]
    cutouts = pd.read_csv(secfilename, header=cutout_header, delimiter=";")[
        : sct_header - cutout_header - 1
    ]
    for cut_start, cut_end in zip(
        cutouts.cutout_start.values, cutouts.cutout_end.values
    ):
        ds_before = ds.sel({dim_name: slice(None, cut_start)})
        ds_after = ds.sel({dim_name: slice(cut_end, None)})
        ds = xr.concat([ds_before, ds_after], dim_name)

    sections = pd.read_csv(secfilename, header=sct_header, delimiter=";")
    sec_cuts = sections["section_cuts"].values.astype("datetime64[ns]")
    sec_nums = sections["section_number"].values
    section_coords = np.zeros(
        len(
            ds.sel({dim_name: slice(None, sec_cuts[0] - np.timedelta64(1, "s"))})[
                dim_name
            ]
        )
    )
    for i, n in zip(np.arange(len(sec_cuts[:-1])), sec_nums[:-1]):
        sec = n * np.ones(
            len(
                ds.sel(
                    {
                        dim_name: slice(
                            sec_cuts[i], sec_cuts[i + 1] - np.timedelta64(1, "s")
                        )
                    }
                )[dim_name]
            )
        )
        section_coords = np.append(section_coords, sec)
    sec = sec_nums[-1] * np.ones(
        len(ds.sel({dim_name: slice(np.datetime64(sec_cuts[-1]), None)})[dim_name])
    )
    section_coords = np.append(section_coords, sec)
    section_coords = section_coords.astype("int")
    ds = ds.assign_coords(section=(dim_name, section_coords))
    ds.section.attrs["long_name"] = "data section"

    campaign_start = np.datetime64(complete_period["dates"].values[1])
    campaign_end = np.datetime64(complete_period["dates"].values[3])
    ds = ds.sel({dim_name: slice(campaign_start, campaign_end)})

    if args.attributes is not None:
        global_attrs_name = args.attributes
        with open(global_attrs_name, "r") as stream:
            global_attrs = yaml.safe_load(stream)
        ds.attrs = global_attrs

    ds.to_netcdf(
        outfilename,
        #        encoding={
        #            dim_name: {"dtype": "<i4", "units": f"seconds since {campaign_start}"}
        #        },
    )
