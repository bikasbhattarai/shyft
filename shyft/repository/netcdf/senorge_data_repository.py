# This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
# See file COPYING for more details **/
# questions: sven.decker@geo.uio.no
import re
import os
import numpy as np
from netCDF4 import Dataset
from shyft import api
from .. import interfaces
from .time_conversion import convert_netcdf_time
from .utils import _limit_2D, _slice_var_2D, _numpy_to_geo_ts_vec, _make_time_slice, _get_files, dummy_var

class SeNorgeDataRepositoryError(Exception):
    pass


class SeNorgeDataRepository(interfaces.GeoTsRepository):
    """
    Repository for geo located timeseries given as WRF(*) data in
    netCDF(3) files.
    NetCDF dataset assumptions:
        * Dimensions:
           Time = UNLIMITED ; // (1 currently)
           DateStrLen = 19 ;
           west_east = 73 ;
           south_north = 60 ;
           bottom_top = 29 ;
           bottom_top_stag = 30 ;
           soil_layers_stag = 4 ;
           west_east_stag = 74 ;
           south_north_stag = 61 ;
        * Variables:
          TODO: A lot.  We really want to list them here?
    (*) WRF model output is from:
        http://www2.mmm.ucar.edu/wrf/users/docs/user_guide_V3/users_guide_chap5.htm
    """

    def __init__(self, epsg, directory, filename=None, elevation_file=None, padding=5000., allow_subset=False):
        """
        Construct the netCDF4 dataset reader for data from WRF NWP model,
        and initialize data retrieval.
        Parameters
        ----------
        epsg: string
            Unique coordinate system id for result coordinates.
            Currently "32632" and "32633" are supported.
        directory: string
            Path to directory holding one or possibly more WRF data files.
            os.path.isdir(directory) should be true, or exception is raised.
        filename: string, optional
            Name of netcdf file in directory that contains spatially
            distributed input data. Can be a regex pattern as well
        padding: float, optional
            padding in meters
        allow_subset: bool
            Allow extraction of a subset of the given source fields
            instead of raising exception.
        """
        directory = os.path.expandvars(directory)
        self._directory = directory
        if filename is None:
            filename = "seNorge2_PREC1d_grid_2015"
#            filename = "wrfout_d03_(\d{4})-(\d{2})"
        self._filename = filename
        self.allow_subset = allow_subset
        if not os.path.isdir(directory):
            raise SeNorgeDataRepositoryError("No such directory '{}'".format(directory))

        if elevation_file is not None:
            self.elevation_file = os.path.join(self._directory, elevation_file)
            if not os.path.isfile(self.elevation_file):
                raise SeNorgeDataRepositoryError(
                    "Elevation file '{}' not found".format(self.elevation_file))
        else:
            self.elevation_file = "/home/sven/workspace/shyft-data/repository/senorge_data_repository/seNorge2_dem_UTM33_comp.nc"

        self.shyft_cs = "+init=EPSG:{}".format(epsg)
        self._padding = padding

        # Field names and mappings
        self.senorge_shyft_map = {
            "mean_temperature": "temperature",
            "precipitation_amount": "precipitation",
            "radiation": "radiation",
            "wind_speed": "wind_speed"
            }

        self.var_units = {"mean_temperature": ['C'],
                          "precipitation_amount": ['kg/m^2', 'Mg/m^2', 'm', 'mm'],
                          "x_wind_10m": ['m/s'],
                          "y_wind_10m": ['m/s'],
                          "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time": ['W s/m^2'],
                          'dew_point_temperature_2m': ['K'],
                          'surface_air_pressure': ['Pa'],
                          'sea_level_pressure': ['Pa'],
                          }

        # Fields that need an additional timeslice because the measure average values
        # self._shift_fields = ("PREC_ACC_NC", "SWDOWN")
        self._shift_fields = ()

    """
            self.wrf_shyft_map = {
                "T2": "temperature",
                "HGT": "z",
                "PREC_ACC_NC": "precipitation",
                "U10": "x_wind",
                "V10": "y_wind",
                "SWDOWN": "radiation",
                "Q2": "mixing_ratio",
                "PSFC": "pressure"}
    """

    def get_timeseries(self, input_source_types, utc_period, geo_location_criteria=None):
        """
        see shyft.repository.interfaces.GeoTsRepository
        """

        filename = os.path.join(self._directory, self._filename)
        if not os.path.isfile(filename):
            if re.compile(self._filename).groups > 0:  # check if it is a filename-pattern
                filename = _get_files(self._directory, self._filename, utc_period.start, SeNorgeDataRepositoryError)
            else:
                raise SeNorgeDataRepositoryError("File '{}' not found".format(filename))

        if any(a in input_source_types for a in ('radiation', 'wind_speed', 'relative_hum')):
            return dummy_var(input_source_types, utc_period, geo_location_criteria)
        else:
            with Dataset(filename) as dataset:
                return self._get_data_from_dataset(dataset, input_source_types, utc_period, geo_location_criteria)



    """
    def _calculate_rel_hum(self, T2, PSFC, Q2):
        # constants
        EZERO = 6.112
        ESLCON1 = 17.67
        ESLCON2 = 29.65
        CELKEL = 273.15
        RD = 287.
        RV = 461.6
        EPS = 0.622

        # calculation
        RH = np.empty_like(T2)
        es = EZERO * np.exp(ESLCON1 * (T2 - CELKEL) / (T2 - ESLCON2))
        qvs = EPS * es / (0.01 * PSFC - (1.0 - EPS) * es)
        RH = Q2 / qvs
        RH[RH > 1.0] = 1.0
        RH[RH < 0.0] = 0.0
        return RH
    """

    def _check_and_get_coord_vars(self, dataset, var_types):
        cs = []
        coord_names = []
        for k, v in self.senorge_shyft_map.items():
            if v in var_types and k in dataset.variables:
                cs.append(dataset.variables[k].getncattr('grid_mapping'))
                coord_names.append([d for d in dataset.variables[k].dimensions if d in ['time', 'X', 'Y', 'latitude', 'longitude']])
        if not all(elem == cs[0] for elem in cs):
            SeNorgeDataRepositoryError('Requested vars have different coord_sys. Do index extraction per var.')
        if not all(elem == coord_names[0] for elem in coord_names):
            SeNorgeDataRepositoryError('Requested vars have different coords. Do index extraction per var.')
        time = dataset.variables.get("time", None)
        if not time:
            raise SeNorgeDataRepositoryError("Time variable not found in dataset.")
        time = convert_netcdf_time(time.units, time)


        if 'Y' in coord_names[0]:
            x = dataset.variables.get("X", None)
            y = dataset.variables.get("Y", None)
        elif 'latitude' in coord_names[0]:
            x = dataset.variables.get("longitude", None)
            y = dataset.variables.get("latitude", None)
        else:
            SeNorgeDataRepositoryError('No recognized coordinate dimension names found.')

        if not all([x, y]):
            raise SeNorgeDataRepositoryError("Spatial Coordinate variables not found in dataset.")
        if 'Y' in coord_names[0]:
            if not all([var.units in ['km', 'm', 'meters'] for var in [x, y]]) and x.units == y.units:
                raise SeNorgeDataRepositoryError("The unit for x and y coordinates should be either m or km.")
        else:
            if not (y.units == 'degrees_north' and x.units == 'degrees_east'):
                raise SeNorgeDataRepositoryError("The unit for latitude and longitude coordinates should be "
                                               "'degrees_north' and 'degrees_east' repectively.")
        coord_conv = 1.
        if y.units == 'km':
            coord_conv = 1000.

        data_cs = dataset.variables.get(cs[0], None)
        if data_cs is None:
            raise SeNorgeDataRepositoryError("No coordinate system information in dataset.")
        return time, x, y, data_cs, coord_conv


    def _get_data_from_dataset(self, dataset, input_source_types, utc_period, geo_location_criteria, ensemble_member=None):


         #x_var = dataset.variables.get("X", None)
        #y_var = dataset.variables.get("Y", None)
        #time = dataset.variables.get("time", None)
        #if not all([x_var, y_var, time]):
        #    raise SeNorgeDataRepositoryError("Something is wrong with the dataset."
        #                                 " x/y coords or time not found.")

        #time = convert_netcdf_time(time.units, time)
        #data_cs_proj4 = "+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
        #if data_cs_proj4 is None:
        #    raise SeNorgeDataRepositoryError("No coordinate system information in dataset.")


        #copied from met_netcdf repo
        # Check for presence and consistency of coordinate variables
        time, x_var, y_var, data_cs, coord_conv = self._check_and_get_coord_vars(dataset, input_source_types)


        # Check units of meteorological variables
#        unit_ok = {k: dataset.variables[k].units in self.var_units[k]
#                   for k in dataset.variables.keys() if self.senorge_shyft_map.get(k, None) in input_source_types}
#        if not all(unit_ok.values()):
#            raise SeNorgeDataRepositoryError("The following variables have wrong unit: {}.".format(
#                ', '.join([k for k, v in unit_ok.items() if not v])))
        # Make spatial slice
        x, y, (x_inds, y_inds), (x_slice, y_slice) = _limit_2D(x_var[:] * coord_conv, y_var[:] * coord_conv,
                                                               data_cs.proj4, self.shyft_cs, geo_location_criteria,
                                                               self._padding, SeNorgeDataRepositoryError)


        # Make temporal slilce
        time_slice, issubset = _make_time_slice(time, utc_period, SeNorgeDataRepositoryError)

        # from wrf_repo
#        x, y, (x_inds, y_inds), (x_slice, y_slice) = _limit_2D(x_var[:], y_var[:], data_cs_proj4, self.shyft_cs, geo_location_criteria, self._padding, SeNorgeDataRepositoryError, clip_in_data_cs=False)
        # end from wrf repo


        raw_data = {}
        for k in dataset.variables.keys():

            if k in self.senorge_shyft_map.keys():
                var = self.senorge_shyft_map.get(k, None)
                if var in input_source_types:
                # if k in self._shift_fields and issubset:  # Add one to time slice
                #     data_time_slice = slice(time_slice.start, time_slice.stop + 1)  # to supply the extra value that is needed for accumulated variables
                # else:
                #     data_time_slice = time_slice
                    data = dataset.variables[k]
                    pure_arr = _slice_var_2D(data, x_var.name, y_var.name, x_slice, y_slice, x_inds, y_inds, SeNorgeDataRepositoryError,
                                         # slices={'time': data_time_slice, 'ensemble_member': ensemble_member})
                                         slices={'time': time_slice, 'ensemble_member': ensemble_member})
                    raw_data[self.senorge_shyft_map[k]] = pure_arr, k, data.units

        if self.elevation_file is not None:
            _x, _y, z = self._read_elevation_file(self.elevation_file, x_var.name, y_var.name,
                                                  geo_location_criteria)
            assert np.linalg.norm(x - _x) < 1.0e-10  # x/y coordinates should match
            assert np.linalg.norm(y - _y) < 1.0e-10
        elif any([nm in dataset.variables.keys() for nm in ['altitude', 'surface_geopotential']]):
            var_nm = ['altitude', 'surface_geopotential'][
                [nm in dataset.variables.keys() for nm in ['altitude', 'surface_geopotential']].index(True)]
            z_data = dataset.variables[var_nm]
            z = _slice_var_2D(z_data, x_var.name, y_var.name, x_slice, y_slice, x_inds, y_inds,
                              SeNorgeDataRepositoryError)
            if var_nm == 'surface_geopotential':
                z /= self._G
        else:
            raise SeNorgeDataRepositoryError("No elevations found in dataset"
                                               ", and no elevation file given.")

        # Make sure requested fields are valid, and that dataset contains the requested data.

        if not self.allow_subset and not (set(raw_data.keys()).issuperset(input_source_types)):
            raise SeNorgeDataRepositoryError("Could not find all data fields")

        if ensemble_member is None and 'ensemble_member' in data.dimensions:
            dims_flat = [d for d in data.dimensions if d in ['time', 'ensemble_member', x_var.name]]
            ens_dim_idx = dims_flat.index('ensemble_member')
            ens_slice = len(dims_flat) * [slice(None)]
            returned_data = []
            for i in range(dataset.dimensions['ensemble_member'].size):
                ens_slice[ens_dim_idx] = i
                #print([(k, raw_data[k][0].shape) for k in raw_data])
                ensemble_raw = {k: (raw_data[k][0][ens_slice], raw_data[k][1], raw_data[k][2]) for k in raw_data.keys()}
                #print([(k,ensemble_raw[k][0].shape) for k in ensemble_raw])
                returned_data.append(_numpy_to_geo_ts_vec(self._transform_raw(ensemble_raw, time[time_slice]),#, issubset=issubset)
                                                          x, y, z, SeNorgeDataRepositoryError))
        else:
            returned_data = _numpy_to_geo_ts_vec(self._transform_raw(raw_data, time[time_slice]),#, issubset=issubset),
                                                 x, y, z, SeNorgeDataRepositoryError)
        return returned_data

        #end of copy

#        time_slice, issubset = _make_time_slice(time, utc_period, SeNorgeDataRepositoryError)
#        x, y, (x_inds, y_inds), (x_slice, y_slice) = _limit_2D(x_var[:], y_var[:], data_cs_proj4, self.shyft_cs, geo_location_criteria, self._padding, SeNorgeDataRepositoryError, clip_in_data_cs=False)

#        raw_data = {}

#        for k in dataset.variables.keys():
#            if self.senorge_shyft_map.get(k, None) in input_source_types:
#                if k in self._shift_fields and issubset:  # Add one to time slice
#                    data_time_slice = slice(time_slice.start, time_slice.stop + 1)
#                else:
#                    data_time_slice = time_slice
#                data = dataset.variables[k]
#                pure_arr = _slice_var_2D(data, x_var.dimensions[2], y_var.dimensions[1], x_slice, y_slice, x_inds, y_inds, SeNorgeDataRepositoryError,
#                                         slices={'Time': data_time_slice, 'ensemble_member': ensemble_member}
#                )
#                raw_data[self.senorge_shyft_map[k]] = pure_arr, k


        # Make sure requested fields are valid, and that dataset contains the requested data.
#        if not self.allow_subset and not (set(raw_data.keys()).issuperset(input_source_types)):
#            raise SeNorgeDataRepositoryError("Could not find all data fields")

        #extracted_data = self._transform_raw(raw_data, time[time_slice], issubset=issubset)
        #return _numpy_to_geo_ts_vec(extracted_data, x, y, SeNorgeDataRepositoryError)


    ######
    # from metnetcdf repo

    def _read_elevation_file(self, filename, x_var_name, y_var_name, geo_location_criteria):
        with Dataset(filename) as dataset:
            elev = dataset.variables["elevation"]
            if "elevation" not in dataset.variables.keys():
                raise interfaces.InterfaceError(
                    "File '{}' does not contain altitudes".format(filename))
            x, y, (x_inds, y_inds), (x_slice, y_slice) = _limit_2D(dataset.variables['easting'][:], dataset.variables['northing'][:],
                                                                  elev.projection, self.shyft_cs, geo_location_criteria,
                                                                   self._padding, SeNorgeDataRepositoryError)
            x_var_name = 'easting'
            y_var_name = 'northing'
            z = _slice_var_2D(elev, x_var_name, y_var_name, x_slice, y_slice, x_inds, y_inds, SeNorgeDataRepositoryError)
            return x, y, z


    def _transform_raw(self, data, time, issubset=False):
        """
        We need full time if deaccumulating
        """

        def noop_time(t):
            t0 = int(t[0])
            t1 = int(t[1])
            return api.TimeAxis(t0, t1 - t0, len(t))

        def dacc_time(t):
            t0 = int(t[0])
            t1 = int(t[1])
            return noop_time(t) if issubset else api.TimeAxis(t0, t1 - t0, len(t) - 1)

        def noop_space(x):
            return x

        def air_temp_conv(T):
            return T - 273.16  # definition says -273.15, but regression test says -273.16..

        def prec_conv(p):
            # return p[1:]
            return p

        # def prec_acc_conv(p):
        #    return np.clip(p[1:] - p[:-1], 0.0, 1000.0)

        def rad_conv(r):
            # dr = r[1:] - r[:-1]
            # return np.clip(dr/(time[1] - time[0]), 0.0, 5000.0)
            return r
        """
        convert_map = {"wind_speed": lambda x, t: (noop_space(x), noop_time(t)),
                       "relative_humidity_2m": lambda x, t: (noop_space(x), noop_time(t)),
                       "T2": lambda x, t: (air_temp_conv(x), noop_time(t)),
                       "SWDOWN": lambda x, t: (rad_conv(x), noop_time(t)),
                       "PREC_ACC_NC": lambda x, t: (prec_conv(x), noop_time(t))}
        # "precipitation_amount_acc": lambda x, t: (prec_acc_conv(x), dacc_time(t))}
        """

        convert_map = {"mean_temperature": lambda x, t, u: (air_temp_conv(x), noop_time(t)),
                       "precipitation_amount": lambda x, t, u: (prec_conv(x), noop_time(t))}

        res = {}
        for k, (v, ak, unit) in data.items():
            res[k] = convert_map[ak](v, time, unit)
        return res


