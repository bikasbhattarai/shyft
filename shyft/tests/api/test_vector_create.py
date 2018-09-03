import unittest
import math
import numpy as np
from shyft.api import Calendar
from shyft.api import TsVector
from shyft.api import TimeSeries
from shyft.api import TimeAxis
from shyft.api import point_interpretation_policy as ts_point_fx
from shyft.api import deltahours
from shyft.api import GeoPointVector
from shyft.api import GeoPoint
from shyft.api import TemperatureSourceVector, TemperatureSource
from shyft.api import PrecipitationSourceVector
from shyft.api import RelHumSourceVector
from shyft.api import WindSpeedSourceVector
from shyft.api import RadiationSourceVector
from shyft.api import create_ts_vector_from_np_array
from shyft.api import create_temperature_source_vector_from_np_array
from shyft.api import create_precipitation_source_vector_from_np_array
from shyft.api import create_wind_speed_source_vector_from_np_array
from shyft.api import create_rel_hum_source_vector_from_np_array
from shyft.api import create_radiation_source_vector_from_np_array


class VectorCreate(unittest.TestCase):

    def test_create_basic(self):
        a = np.array([[1.1, 1.2, 1.3], [2.1, 2.2, 2.3]], dtype=np.float64)
        ta = TimeAxis(deltahours(0), deltahours(1), 3)
        tsv = create_ts_vector_from_np_array(ta, a, ts_point_fx.POINT_AVERAGE_VALUE)
        self.assertIsNotNone(tsv)
        for i in range(2):
            self.assertTrue(np.allclose(tsv[i].values.to_numpy(), a[i]))
            self.assertTrue(ta == tsv[i].time_axis)
            self.assertEqual(tsv[i].point_interpretation(), ts_point_fx.POINT_AVERAGE_VALUE)
        # create missmatch throws
        b = np.array([[], []], dtype=np.float64)
        try:
            create_ts_vector_from_np_array(ta, b, ts_point_fx.POINT_AVERAGE_VALUE)
            self.assertTrue(False, "Should throw for missmatch time-axis")
        except RuntimeError as e:
            pass
        # create empty ts works
        tb = TimeAxis(0, 0, 0)
        r = create_ts_vector_from_np_array(tb, b, ts_point_fx.POINT_AVERAGE_VALUE)
        self.assertEqual(len(r), 2)
        for ts in r:
            self.assertFalse(ts)
        # create empty returns empty
        c = np.empty(shape=(0, 0), dtype=np.float64)
        z = create_ts_vector_from_np_array(tb, c, ts_point_fx.POINT_AVERAGE_VALUE)
        self.assertEqual(len(z), 0)

    def test_create_xx_source_vector(self):
        # arrange the setup
        a = np.array([[1.1, 1.2, 1.3], [2.1, 2.2, 2.3]], dtype=np.float64)
        ta = TimeAxis(deltahours(0), deltahours(1), 3)
        gpv = GeoPointVector()
        gpv[:] = [GeoPoint(1, 2, 3), GeoPoint(4, 5, 6)]
        cfs = [(create_precipitation_source_vector_from_np_array, PrecipitationSourceVector),
               (create_temperature_source_vector_from_np_array, TemperatureSourceVector),
               (create_radiation_source_vector_from_np_array, RadiationSourceVector),
               (create_rel_hum_source_vector_from_np_array, RelHumSourceVector),
               (create_radiation_source_vector_from_np_array, RadiationSourceVector)]
        # test all creation types:
        for cf in cfs:
            r = cf[0](ta, gpv, a, ts_point_fx.POINT_AVERAGE_VALUE)  # act here
            self.assertTrue(isinstance(r, cf[1]))  # then the asserts
            self.assertEqual(len(r), len(gpv))
            for i in range(len(gpv)):
                self.assertEqual(r[i].mid_point(), gpv[i])
                self.assertTrue(np.allclose(r[i].ts.values.to_numpy(), a[i]))
                self.assertEqual(r[i].ts.point_interpretation(), ts_point_fx.POINT_AVERAGE_VALUE)

    def test_create_xx_vector_from_list(self):
        """ verify we can construct from list of xxSource"""
        ts = TemperatureSource(
            GeoPoint(1.0, 2.0, 3.0),
            TimeSeries(TimeAxis(0, 3600, 10), fill_value=1.0, point_fx=ts_point_fx.POINT_AVERAGE_VALUE)
        )
        tsv = TemperatureSourceVector([ts])
        self.assertEqual(len(tsv), 1)
        self.assertEqual(len(TemperatureSourceVector([])), 0)
        self.assertEqual(len(TemperatureSourceVector([ts, ts])), 2)

    def test_create_tsvector_from_ts_list(self):
        ts_list = [TimeSeries(TimeAxis(0, 3600, 10), fill_value=float(i), point_fx=ts_point_fx.POINT_AVERAGE_VALUE) for i in range(3)]
        tsv = TsVector(ts_list)
        assert tsv
        assert len(tsv) == 3
        assert tsv[0].value(0) == 0.0
        assert tsv[1].value(0) == 1.0
        assert tsv[2].value(0) == 2.0
