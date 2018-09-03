from __future__ import print_function
import unittest
from os import path

from shyft import shyftdata_dir
from shyft import api
from shyft.repository.netcdf.met_netcdf_data_repository import MetNetcdfDataRepository
from shyft.repository.netcdf.met_netcdf_data_repository import MetNetcdfDataRepositoryError
from shapely.geometry import box
import netCDF4
import numpy as np


class MetNetcdfDataRepositoryTestCase(unittest.TestCase):

    def test_get_timeseries(self):
        """
        Simple regression test of arome data respository.
        """
        EPSG, bbox, bpoly = self.arome_epsg_bbox

        # Period start
        n_hours = 30
        t0 = api.YMDhms(2015, 8, 24, 0)
        date_str = "{}{:02}{:02}_{:02}".format(t0.year, t0.month, t0.day, t0.hour)
        utc = api.Calendar()  # No offset gives Utc
        period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_hours))

        base_dir = path.join(shyftdata_dir, "repository", "arome_data_repository")
        f1 = "arome_metcoop_red_default2_5km_{}_diff_time_unit.nc".format(date_str)
        f2 = "arome_metcoop_red_test2_5km_{}.nc".format(date_str)

        ar1 = MetNetcdfDataRepository(EPSG, base_dir, filename=f1)
        ar2 = MetNetcdfDataRepository(EPSG, base_dir, filename=f2, elevation_file=f1)
        ar1_data_names = ("temperature", "wind_speed", "precipitation", "relative_humidity")
        ar2_data_names = ("radiation",)
        sources = ar1.get_timeseries(ar1_data_names, period, geo_location_criteria=bpoly)
        self.assertTrue(len(sources) > 0)
        sources2 = ar2.get_timeseries(ar2_data_names, period, geo_location_criteria=bpoly)

        self.assertTrue(set(sources) == set(ar1_data_names))
        self.assertTrue(set(sources2) == set(ar2_data_names))
        self.assertTrue(sources["temperature"][0].ts.size() == n_hours + 1)
        r0 = sources2["radiation"][0].ts
        p0 = sources["precipitation"][0].ts
        temp0 = sources["temperature"][0].ts
        self.assertTrue(r0.size() == n_hours)
        self.assertTrue(p0.size() == n_hours)
        self.assertTrue(r0.time(0) == temp0.time(0))
        self.assertTrue(p0.time(0) == temp0.time(0))
        self.assertTrue(r0.time_axis.total_period().end == temp0.time(temp0.size() - 1))
        self.assertTrue(p0.time_axis.total_period().end == temp0.time(temp0.size() - 1))
        self.assertTrue(p0.time(0), period.start)

    @property
    def arome_epsg_bbox(self):
        """A slice of test-data located in shyft-data repository/arome."""
        EPSG = 32632
        x0 = 436100.0  # lower left
        y0 = 6823000.0  # lower right
        nx = 74
        ny = 24
        dx = 1000.0
        dy = 1000.0
        return EPSG, ([x0, x0 + nx*dx, x0 + nx*dx, x0], [y0, y0, y0 + ny*dy, y0 + ny*dy]), box(x0, y0, x0 + dx*nx, y0 + dy*ny)

    def test_get_forecast(self):
        # Period start
        n_hours = 65
        utc = api.Calendar()  # No offset gives Utc
        t0 = utc.time(2015, 8, 24, 6)
        period1 = api.UtcPeriod(t0, t0 + api.deltahours(n_hours))
        period2 = api.UtcPeriod(t0 + api.deltahours(6), t0 + api.deltahours(6) + api.deltahours(n_hours))
        t_c1 = t0 + api.deltahours(1)
        t_c2 = t0 + api.deltahours(7)

        base_dir = path.join(shyftdata_dir, "repository", "arome_data_repository")
        # pattern = "arome_metcoop*default2_5km_*.nc"
        pattern = "arome_metcoop_red_default2_5km_(\d{4})(\d{2})(\d{2})[T_](\d{2})Z?.nc$"
        EPSG, bbox, bpoly = self.arome_epsg_bbox

        repos = MetNetcdfDataRepository(EPSG, base_dir, filename=pattern)
        data_names = ("temperature", "wind_speed", "precipitation", "relative_humidity")
        tc1_sources = repos.get_forecast(data_names, period1, t_c1, geo_location_criteria=bpoly)
        tc2_sources = repos.get_forecast(data_names, period2, t_c2, geo_location_criteria=bpoly)

        self.assertTrue(len(tc1_sources) == len(tc2_sources))
        self.assertTrue(set(tc1_sources) == set(data_names))
        self.assertTrue(tc1_sources["temperature"][0].ts.size() == n_hours + 1)

        tc1_precip = tc1_sources["precipitation"][0].ts
        tc2_precip = tc2_sources["precipitation"][0].ts

        self.assertEqual(tc1_precip.size(), n_hours)
        self.assertTrue(tc1_precip.time(0) != tc2_precip.time(0))

    def test_get_ensemble(self):
        EPSG = 32633
        upper_left_x = 436100.0
        upper_left_y = 7417800.0
        nx = 74
        ny = 94
        dx = 1000.0
        dy = 1000.0
        # Period start
        n_hours = 30
        utc = api.Calendar()  # No offset gives Utc
        t0 = utc.time(2015, 7, 26)
        period = api.UtcPeriod(t0, t0 + api.deltahours(n_hours))
        t_c = t0 + api.deltahours(1)

        base_dir = path.join(shyftdata_dir, "netcdf", "arome")
        # pattern = "fc*.nc"
        pattern = "fc_(\d{4})(\d{2})(\d{2})[T_](\d{2})Z?.nc$"
        bpoly = box(upper_left_x, upper_left_y - ny*dy, upper_left_x + nx*dx, upper_left_y)
        try:
            repos = MetNetcdfDataRepository(EPSG, base_dir, filename=pattern)
            data_names = ("temperature", "wind_speed", "relative_humidity")
            ensemble = repos.get_forecast_ensemble(data_names, period, t_c, geo_location_criteria=bpoly)
            self.assertTrue(isinstance(ensemble, list))
            self.assertEqual(len(ensemble), 10)
        except MetNetcdfDataRepositoryError as adre:
            self.skipTest("(test inconclusive- missing arome-data {0})".format(adre))

    def test_wrong_directory(self):
        with self.assertRaises(MetNetcdfDataRepositoryError) as context:
            MetNetcdfDataRepository(32632, "Foobar", filename="")
        self.assertEqual("No such directory 'Foobar'", context.exception.args[0])

    def test_wrong_elevation_file(self):
        with self.assertRaises(MetNetcdfDataRepositoryError) as context:
            MetNetcdfDataRepository(32632, shyftdata_dir, filename="", elevation_file="plain_wrong.nc")
        self.assertTrue(all(x in context.exception.args[0] for x in ["Elevation file",
                                                                     "not found"]))

    def test_wrong_file(self):
        with self.assertRaises(MetNetcdfDataRepositoryError) as context:
            utc = api.Calendar()  # No offset gives Utc
            t0 = api.YMDhms(2015, 12, 25, 18)
            period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(30))
            ar1 = MetNetcdfDataRepository(32632, shyftdata_dir, filename="plain_wrong.nc")
            ar1.get_timeseries(("temperature",), period, geo_location_criteria=None)
        self.assertTrue(all(x in context.exception.args[0] for x in ["File", "not found"]))

    def test_wrong_forecast(self):
        with self.assertRaises(MetNetcdfDataRepositoryError) as context:
            utc = api.Calendar()  # No offset gives Utc
            t0 = api.YMDhms(2015, 12, 25, 18)
            period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(30))
            ar1 = MetNetcdfDataRepository(32632, shyftdata_dir, filename="plain_wrong_(\d{4})(\d{2})(\d{2})[T_](\d{2})Z?.nc")
            ar1.get_forecast(("temperature",), period, utc.time(t0), geo_location_criteria=None)
        self.assertTrue(all(x in context.exception.args[0] for x in
                            ["No matches found for file_pattern = ", "and t_c = "]))

    # TODO: The test below verifies that error is raised when no point is found within polygon envelop. Adde a test for the case where there are points within the envelop but now within the polygon itself.
    def test_no_point_inside_polygon_bounds(self):
        EPSG, bbox, bpoly = self.arome_epsg_bbox
        bounds = bpoly.bounds
        bpoly = box(bounds[0], 6010000.0, bounds[2], 6035000.0)
        # Period start
        year = 2015
        month = 8
        day = 24
        hour = 6
        n_hours = 30
        date_str = "{}{:02}{:02}_{:02}".format(year, month, day, hour)
        utc = api.Calendar()  # No offset gives Utc
        t0 = api.YMDhms(year, month, day, hour)
        period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_hours))

        base_dir = path.join(shyftdata_dir, "repository", "arome_data_repository")
        filename = "arome_metcoop_red_default2_5km_{}.nc".format(date_str)
        reader = MetNetcdfDataRepository(EPSG, base_dir, filename=filename, padding=0.0)
        data_names = ("temperature", "wind_speed", "precipitation", "relative_humidity")
        with self.assertRaises(MetNetcdfDataRepositoryError) as context:
            reader.get_timeseries(data_names, period, geo_location_criteria=bpoly)
        self.assertEqual("No points in dataset which are within the bounding box of the geo_location_criteria polygon.",
                         context.exception.args[0])

    def test_geo_location_criteria_is_None(self):
        EPSG, _, _ = self.arome_epsg_bbox
        # Period start
        year = 2015
        month = 8
        day = 24
        hour = 6
        n_hours = 30
        date_str = "{}{:02}{:02}_{:02}".format(year, month, day, hour)
        utc = api.Calendar()  # No offset gives Utc
        t0 = api.YMDhms(year, month, day, hour)
        period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_hours))

        base_dir = path.join(shyftdata_dir, "repository", "arome_data_repository")
        filename = "arome_metcoop_red_default2_5km_{}.nc".format(date_str)
        reader = MetNetcdfDataRepository(EPSG, base_dir, filename=filename)
        data_names = ("temperature", "wind_speed", "precipitation", "relative_humidity")
        with netCDF4.Dataset(path.join(base_dir, filename)) as ds:
            nb_pts_in_file = ds.dimensions['x'].size*ds.dimensions['y'].size
        srcs = reader.get_timeseries(data_names, period, None)
        self.assertEqual(len(srcs['temperature']), nb_pts_in_file)

    def test_tiny_bbox(self):
        EPSG, _, _ = self.arome_epsg_bbox

        x = 432425.910493  # x coord of one pt in test file
        y = 6819847.92879  # y coord of one pt in test file
        dxy = 1000.  # should be less than the grid resolution (2500 m) to enclose only one point
        bpoly = box(x - dxy, y - dxy, x + dxy, y + dxy)  # a polygon containing only tht above point

        # Period start
        year = 2015
        month = 8
        day = 24
        hour = 6
        n_hours = 30
        date_str = "{}{:02}{:02}_{:02}".format(year, month, day, hour)
        utc = api.Calendar()  # No offset gives Utc
        t0 = api.YMDhms(year, month, day, hour)
        period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_hours))

        base_dir = path.join(shyftdata_dir, "repository", "arome_data_repository")
        filename = "arome_metcoop_red_default2_5km_{}.nc".format(date_str)
        reader = MetNetcdfDataRepository(EPSG, base_dir, filename=filename, padding=0.0)
        data_names = ("temperature", "wind_speed", "precipitation", "relative_humidity")
        try:
            tss = reader.get_timeseries(data_names, period, geo_location_criteria=bpoly)
        except MetNetcdfDataRepositoryError as err:
            self.fail("reader.get_timeseries raised MetNetcdfDataRepositoryError('{}') "
                      "unexpectedly.".format(err.args[0]))
        self.assertEqual(len(tss['temperature']), 1)

    def test_subsets(self):
        EPSG, bbox, bpoly = self.arome_epsg_bbox
        # Period start
        year = 2015
        month = 8
        day = 24
        hour = 6
        n_hours = 30
        date_str = "{}{:02}{:02}_{:02}".format(year, month, day, hour)
        utc = api.Calendar()  # No offset gives Utc
        t0 = api.YMDhms(year, month, day, hour)
        period = api.UtcPeriod(utc.time(t0), utc.time(t0) + api.deltahours(n_hours))

        base_dir = path.join(shyftdata_dir, "repository", "arome_data_repository")
        filename = "arome_metcoop_red_default2_5km_{}.nc".format(date_str)

        data_names = ("temperature", "wind_speed", "precipitation", "relative_humidity", "radiation")
        allow_subset = False
        reader = MetNetcdfDataRepository(EPSG, base_dir, filename=filename,
                                         allow_subset=allow_subset)
        with self.assertRaises(MetNetcdfDataRepositoryError) as context:
            reader.get_timeseries(data_names, period, None)
        self.assertEqual("Could not find all data fields", context.exception.args[0])
        allow_subset = True
        reader = MetNetcdfDataRepository(EPSG, base_dir, filename=filename,
                                         allow_subset=allow_subset)
        try:
            sources = reader.get_timeseries(data_names, period, geo_location_criteria=bpoly)
        except MetNetcdfDataRepositoryError as e:
            self.fail("MetNetcdfDataRepository.get_timeseries(data_names, period, None) "
                      "raised MetNetcdfDataRepositoryError unexpectedly.")
        self.assertEqual(len(sources), len(data_names) - 1)

    def test_utc_period_is_None(self):
        EPSG, bbox, bpoly = self.arome_epsg_bbox
        # Period start
        utc = api.Calendar()  # No offset gives Utc
        t0 = utc.time(2015, 8, 24, 6)
        t0_ymdhs = utc.calendar_units(t0)
        date_str = "{}{:02}{:02}_{:02}".format(*[getattr(t0_ymdhs, k) for k in ['year', 'month', 'day', 'hour']])
        period = None

        base_dir = path.join(shyftdata_dir, "repository", "arome_data_repository")
        filename = "arome_metcoop_red_default2_5km_{}.nc".format(date_str)
        reader = MetNetcdfDataRepository(EPSG, base_dir, filename=filename)
        src_name = "temperature"
        var_name_in_file = [k for k, v in reader._arome_shyft_map.items() if v == src_name][0]
        with netCDF4.Dataset(path.join(base_dir, filename)) as ds:
            var = ds.variables[var_name_in_file]
            nb_timesteps = var.shape[var.dimensions.index('time')]
        srcs = reader.get_timeseries((src_name,), period, geo_location_criteria=bpoly)
        self.assertEqual(srcs[src_name][0].ts.size(), nb_timesteps)

    def test_transform_functions_fixed_interval(self):
        """
        test the _transform_raw function.
        """
        EPSG, bbox, bpoly = self.arome_epsg_bbox

        # Period start
        t0 = api.YMDhms(2015, 8, 24, 0)
        date_str = "{}{:02}{:02}_{:02}".format(t0.year, t0.month, t0.day, t0.hour)
        utc = api.Calendar()  # No offset gives Utc

        base_dir = path.join(shyftdata_dir, "repository", "arome_data_repository")
        f1 = "arome_metcoop_red_default2_5km_{}_diff_time_unit.nc".format(date_str)
        ar1 = MetNetcdfDataRepository(EPSG, base_dir, filename=f1)
        np_raw_array = np.array(
                    [  # 0  # 1 #  2 #  3
                        [1.0, 2.0, 3.0, 4.0],
                        [1.1, 2.1, 3.1, 4.1],
                        [1.2, 2.2, 3.2, 4.2],
                        [1.4, 2.5, 3.6, 4.7]
                    ], dtype=np.float64
                )
        raw_values = {
            'wind_speed': (np_raw_array, 'wind_speed', 'm/s'),
            'rel_hum': (np_raw_array, 'relative_humidity_2m', '?'),
            'temperature': (273.15 + np_raw_array, 'air_temperature_2m', 'K'),
            'radiation': (3600.0*np_raw_array, 'integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time', 'W s/m2'),
            'prepitation_acc': (np_raw_array, 'precipitation_amount_acc', 'Mg/m^2'),
            'prepitation': (np_raw_array, 'precipitation_amount', 'mm')
        }
        raw_time = np.array([0, 3600, 7200, 10800], dtype=np.int64)

        rd = ar1._transform_raw(raw_values, raw_time)
        ta3 = api.TimeAxis(api.time(0), api.time(3600), 3)
        ta4 = api.TimeAxis(api.time(0), api.time(3600), 4)
        e_precip_acc = np.array(
                    [  # 0  # 1 #  2 #  3
                        [100.0, 100.0, 100.0, 100.0],
                        [100.0, 100.0, 100.0, 100.0],
                        [200.0, 300.0, 400.0, 500.0],
                    ], dtype=np.float64
        )
        e_precip = np.array(
                    [  # 0  # 1 #  2 #  3
                        [1.1, 2.1, 3.1, 4.1],
                        [1.2, 2.2, 3.2, 4.2],
                        [1.4, 2.5, 3.6, 4.7]
                    ], dtype=np.float64
                )
        e_rad = np.array(
            [  # 0  # 1 #  2 #  3
                [0.1, 0.1, 0.1, 0.1],
                [0.1, 0.1, 0.1, 0.1],
                [0.2, 0.3, 0.4, 0.5],
            ], dtype=np.float64
        )
        e = {
            'wind_speed': (np_raw_array, ta4),
            'rel_hum': (np_raw_array, ta4),
            'temperature': (np_raw_array, ta4),
            'radiation': (e_rad, ta3),
            'prepitation_acc': (e_precip_acc, ta3),
            'prepitation': (e_precip, ta3)
        }

        self.assertIsNotNone(rd)
        for k, r in rd.items():
            self.assertTrue(k in e)
            self.assertEqual(r[1], e[k][1], "expect correct time-axis")
            self.assertTrue(np.allclose(r[0], e[k][0]), "expect exact correct values")

    def test_transform_functions_variable_interval(self):
        """
        test the _transform_raw function.
        """
        EPSG, bbox, bpoly = self.arome_epsg_bbox

        # Period start
        n_hours = 30
        t0 = api.YMDhms(2015, 8, 24, 0)
        date_str = "{}{:02}{:02}_{:02}".format(t0.year, t0.month, t0.day, t0.hour)
        utc = api.Calendar()  # No offset gives Utc

        base_dir = path.join(shyftdata_dir, "repository", "arome_data_repository")
        f1 = "arome_metcoop_red_default2_5km_{}_diff_time_unit.nc".format(date_str)
        ar1 = MetNetcdfDataRepository(EPSG, base_dir, filename=f1)
        np_raw_array = np.array(
                    [  # 0  # 1 #  2 #  3
                        [1.0, 2.0, 3.0, 4.0],
                        [1.1, 2.1, 3.1, 4.1],
                        [1.2, 2.2, 3.2, 4.2],
                        [1.4, 2.5, 3.6, 4.7]
                    ], dtype=np.float64
                )
        raw_values = {
            'wind_speed': (np_raw_array, 'wind_speed', 'm/s'),
            'rel_hum': (np_raw_array, 'relative_humidity_2m', '?'),
            'temperature': (273.15 + np_raw_array, 'air_temperature_2m', 'K'),
            'radiation': (3600.0*np_raw_array, 'integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time', 'W s/m2'),
            'prepitation_acc': (np_raw_array, 'precipitation_amount_acc', 'Mg/m^2'),
            'prepitation': (np_raw_array, 'precipitation_amount', 'mm')
        }
        raw_time = np.array([0, 3600, 7200, 7200+2*3600], dtype=np.int64)  # last step is 2 hours!

        rd = ar1._transform_raw(raw_values, raw_time)
        ta3 = api.TimeAxis(api.UtcTimeVector(raw_time[:-1]), api.time(int(raw_time[-1])))
        ta4 = api.TimeAxis(api.UtcTimeVector(raw_time), api.time(int(raw_time[-1]+2*3600)))  # assume last step is also 2 hours
        e_precip_acc = np.array(
                    [  # 0  # 1 #  2 #  3
                        [100.0, 100.0, 100.0, 100.0],
                        [100.0, 100.0, 100.0, 100.0],
                        [100.0, 150.0, 200.0, 250.0],
                    ], dtype=np.float64
        )
        e_precip = np.array(
                    [  # 0  # 1 #  2 #  3
                        [1.1, 2.1, 3.1, 4.1],
                        [1.2, 2.2, 3.2, 4.2],
                        [1.4, 2.5, 3.6, 4.7]
                    ], dtype=np.float64
                )
        e_rad = np.array(
            [  # 0  # 1 #  2 #  3
                [0.1, 0.1, 0.1, 0.1],
                [0.1, 0.1, 0.1, 0.1],
                [0.1, 0.15, 0.2, 0.25],
            ], dtype=np.float64
        )
        e = {
            'wind_speed': (np_raw_array, ta4),
            'rel_hum': (np_raw_array, ta4),
            'temperature': (np_raw_array, ta4),
            'radiation': (e_rad, ta3),
            'prepitation_acc': (e_precip_acc, ta3),
            'prepitation': (e_precip, ta3)
        }

        self.assertIsNotNone(rd)
        for k, r in rd.items():
            self.assertTrue(k in e)
            self.assertEqual(r[1], e[k][1], "expect correct time-axis")
            self.assertTrue(np.allclose(r[0], e[k][0]), "expect exact correct values")


if __name__ == "__main__":
    unittest.main()
