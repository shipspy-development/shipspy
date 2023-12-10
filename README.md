# shipspy: ship campaign data processing and standardisation

With shipspy, data from ship campaigns can be processed and standardised. A number of instruments are included.

## Setup

Install shipspy (coming soon):

```
pip install shipspy
```
## Using shipspy

The following processing options are available.

### DShip

With the dship subcommand, data from the ship integrated system can be processed. For an exemplary order to download underway data, see [1]. The following steps are done:
* Fixing time stamps 
* Averaging to minutely time stamps
* Removing unphysical values
* Renaming variables, changing units, and adding attributes
* Converting to netCDF

To process dship data run
```
shipspy dship -i <input file> -o <output file> -a <attribute dictionary> -s <ship>
```
with
* `input file`: downloaded data from BSH or GEOMAR (or directly from the ship) as txt or dat file with unix time stamp (seconds since 1970-01-01)
* `output file`: file name of netCDF output file
* `attribute dictionary`: yaml dictionary with variable names and attribute. For example dictionaries for the different research vessels see `rv_information` 
* `ship`: name of the ship. Options are `Merian`, `Sonne`, and `Meteor`

### Renaming

The rename command can be used to 
* fix the time stamps
* rename variables and add attributes
* convert units to SI
* remove unphysical values
* Convert and save the file as netCDF

To rename a data file run
```
shipspy rename -i <input file> -o <output file> -a <attribute dictionary> -d <instrument>
```
with 
* `input file`: input file, file format depends on the instrument
* `output file`: file name of netCDF output file
* `attribute dictionary`: yaml dictionary with variable names and attribute. For examples see [[1]](https://github.com/LauraKoehler/arc_processing)
* `instrument`: instrument name. Options are `calitoo`, `ceilometer`, `ctd`, `hatpro`, `microtops`, `radiosondes`, `test`, `uav`

### Sections

The section command adds a new coordinate and specifies the time period of the campaign. To use it run
```
shipspy section -i <input file> -o <output file> -s <section file> -t <time dimension name> -a <global attribute dictionary>
```
with
* `input file`: input netCDF file
* `output file`: file name of netCDF output file
* `section file`: txt file specifying the campaign dates (start, end, break) and the sections in which it should be divided. For an example, see [[1]](https://github.com/LauraKoehler/arc_processing).
* `time dimension name`: name of the time dimension. Default is `time` but sometimes it can be something like `start_time`.
* `global attribute dictionary`: yaml file with global attributes if wanted. For example files see [[1]](https://github.com/LauraKoehler/arc_processing).

The [repository [1]](https://github.com/LauraKoehler/arc_processing) with the settings for the ARC and additional scripts can serve as a template.

# References
[1] Köhler, L. (2023). ARC: Processing of atmospheric and oceanographic measurements (Version v1.0.0) [Computer software]. https://github.com/LauraKoehler/arc_processing
