﻿from shyft import api
import unittest
import datetime as dt
import math


class Calendar(unittest.TestCase):
    """Verify and illustrate the Calendar & time from the api core, using
    pyunit. Note that the Calendar not yet support local/dst semantics (but
    plan to do so) Nevertheless, keeping it here allow users of api-Core to
    use start practicing time/calendar perimeter.
    """

    def setUp(self):
        self.utc = api.Calendar()  # A utc calendar
        self.std = api.Calendar(3600)  # UTC+01

    def tearDown(self):
        self.utc = None
        self.std = None

    def test_time_zone_region_id_list(self):
        self.assertGreater(len(api.Calendar.region_id_list()), 1)

    def test_create_calendar_from_region_id(self):
        osl = api.Calendar("Europe/Oslo")
        self.assertIsNotNone(osl)
        self.assertEqual(osl.tz_info.name(), "Europe/Oslo")
        self.assertEqual(osl.tz_info.base_offset(), 3600)
        t = osl.time(2015, 6, 1)
        self.assertTrue(osl.tz_info.is_dst(t))
        self.assertTrue(osl.tz_info.utc_offset(t), 7200)

    def test_calendar_add_and_diff_units(self):
        osl = api.Calendar("Europe/Oslo")
        t0 = osl.time(2016, 6, 1, 12, 0, 0)
        t1 = osl.add(t0, api.Calendar.DAY, 7)
        t2 = osl.add(t1, api.Calendar.WEEK, -1)
        self.assertEqual(t0, t2)
        self.assertEqual(7, osl.diff_units(t0, t1, api.Calendar.DAY))
        self.assertEqual(1, osl.diff_units(t0, t1, api.Calendar.WEEK))
        self.assertEqual(0, osl.diff_units(t0, t1, api.Calendar.MONTH))
        self.assertEqual(7*24, osl.diff_units(t0, t1, api.deltahours(1)))

    def test_calendar_add_during_dst(self):
        osl = api.Calendar("Europe/Oslo")
        t0 = osl.time(2016, 3, 27)  # dst change during spring
        t1 = osl.add(t0, api.Calendar.DAY, 1)
        t2 = osl.add(t1, api.Calendar.DAY, -1)
        self.assertEqual(t0, t2)
        self.assertEqual("2016-03-28T00:00:00+02", osl.to_string(t1))
        self.assertEqual(1, osl.diff_units(t0, t1, api.Calendar.DAY))
        self.assertEqual(23, osl.diff_units(t0, t1, api.Calendar.HOUR))
        t0 = osl.time(2016, 10, 30)
        t1 = osl.add(t0, api.Calendar.WEEK, 1)
        t2 = osl.add(t1, api.Calendar.WEEK, -1)
        self.assertEqual(t0, t2)
        self.assertEqual("2016-11-06T00:00:00+01", osl.to_string(t1))
        self.assertEqual(168 + 1, osl.diff_units(t0, t1, api.Calendar.HOUR))

    def test_calendar_add_3h_during_dst(self):
        osl = api.Calendar("Europe/Oslo")
        t0 = osl.time(2016, 3, 27)  # dst change during spring
        t1 = osl.add(t0, api.Calendar.DAY, 1)
        dt3h = api.deltahours(3)
        d3h = osl.diff_units(t0, t1, dt3h)
        self.assertEqual(8, d3h)

    def test_trim_day(self):
        t = api.utctime_now()
        td = self.std.trim(t, api.Calendar.DAY)
        c = self.std.calendar_units(td)
        a = self.std.calendar_units(t)
        self.assertEqual(c.second, 0, 'incorrect seconds should be 0')
        self.assertEqual(c.minute, 0, 'trim day should set minutes to 0')
        self.assertEqual(c.hour, 0, 'trim day should set hours to 0')
        self.assertEqual(a.year, c.year, 'trim day Should leave same year')
        self.assertEqual(a.month, c.month, 'trim day  Should leave month')
        self.assertEqual(a.day, c.day, 'should leave same day')

    def test_quarter(self):
        t = self.std.time(2017, 2, 28, 1, 2, 3)
        tt = self.std.trim(t, api.Calendar.QUARTER)
        self.assertEqual(tt, self.std.time(2017, 1, 1))
        self.assertEqual(1, self.std.quarter(t))

    def test_conversion_roundtrip(self):
        c1 = api.YMDhms(1960, 1, 2, 3, 4, 5)
        t1 = self.std.time(c1)
        c2 = self.std.calendar_units(t1)
        cw = self.std.calendar_week_units(t1)
        tw2 = self.std.time_from_week(1959, 53, 6, 3, 4, 5)
        tw1 = self.std.time(cw)
        self.assertEqual(tw2, t1)
        self.assertEqual(tw1, t1)
        self.assertEqual(cw.iso_year, 1959)
        self.assertEqual(cw.iso_week, 53)
        self.assertEqual(cw.week_day, 6)
        self.assertEqual(cw.hour, 3)
        self.assertEqual(cw.minute, 4)
        self.assertEqual(cw.second, 5)
        self.assertEqual(c1.year, c2.year, 'roundtrip should keep year')
        self.assertEqual(c1.month, c2.month)
        self.assertEqual(c1.day, c2.day)
        self.assertEqual(c1.hour, c2.hour)
        self.assertEqual(c1.second, c2.second)

    def test_utctime_now(self):
        a = api.utctime_now()
        x = dt.datetime.utcnow()
        b = self.utc.time(x.year, x.month, x.day,
                          x.hour, x.minute, x.second)
        self.assertLess(abs(a - b), 2, 'Should be less than 2 seconds')

    def test_utc_time_to_string(self):
        t = self.std.time(2000, 1, 2, 3, 4, 5, 6)
        s = self.std.to_string(t)
        self.assertEqual(s, "2000-01-02T03:04:05.000006+01")
        t = self.std.time(1960, 1, 2, 3, 4, 5, 6)
        s = self.std.to_string(t)
        self.assertEqual(s, "1960-01-02T03:04:05.000006+01")

    def test_UtcPeriod_to_string(self):
        c1 = api.YMDhms(2000, 1, 2, 3, 4, 5)
        t = self.utc.time(c1)
        p = api.UtcPeriod(t, t + api.deltahours(1))
        s = p.to_string()
        self.assertEqual(s, "[2000-01-02T03:04:05Z,2000-01-02T04:04:05Z>")
        s2 = self.std.to_string(p)
        self.assertEqual(s2, "[2000-01-02T04:04:05+01,2000-01-02T05:04:05+01>")

    def test_UtcPeriod_str(self):
        c1 = api.YMDhms(2000, 1, 2, 3, 4, 5)
        t = self.utc.time(c1)
        p = api.UtcPeriod(t, t + api.deltahours(1))
        s = str(p)
        self.assertEqual(s, "[2000-01-02T03:04:05Z,2000-01-02T04:04:05Z>")
        s = repr(p)
        self.assertEqual(s, "[2000-01-02T03:04:05Z,2000-01-02T04:04:05Z>")

    def test_UtcPeriod_methods(self):
        p0 = api.UtcPeriod()
        px = api.UtcPeriod(self.utc.time(2010), self.utc.time(2000))
        py = api.UtcPeriod(self.utc.time(2010), self.utc.time(2010))
        p1 = api.UtcPeriod(self.utc.time(2015), self.utc.time(2016))
        p2 = api.UtcPeriod(self.utc.time(2016), self.utc.time(2017))
        p3 = api.UtcPeriod(p1.start, p2.end)
        self.assertFalse(p0.valid())
        self.assertFalse(px.valid())
        self.assertTrue(py.valid())
        self.assertTrue(p1.contains(p1.start))
        self.assertFalse(p1.contains(p1.end))
        self.assertFalse(p1.overlaps(p2))
        self.assertFalse(p2.overlaps(p1))
        self.assertTrue(p3.overlaps(p1))
        self.assertTrue(p3.overlaps(p2))

    def test_UtcPeriod_trim(self):
        utc = api.Calendar()
        t0 = utc.time(2018, 1, 1, 0, 0, 0)
        t1 = utc.time(2018, 2, 1, 0, 0, 0)
        t2 = utc.time(2018, 3, 1, 0, 0, 0)
        t3 = utc.time(2018, 4, 1, 0, 0, 0)
        p_0_1 = api.UtcPeriod(t0, t1)
        p_0_3 = api.UtcPeriod(t0, t3)
        assert api.UtcPeriod(t0, t1).trim(utc, utc.MONTH, api.trim_policy.TRIM_IN) == p_0_1
        assert api.UtcPeriod(t0, t1).trim(utc, 3600*24*30, api.trim_policy.TRIM_IN) == p_0_1
        assert api.UtcPeriod(t0, t1 + 1).trim(utc, utc.MONTH) == p_0_1
        assert api.UtcPeriod(t1 - 1, t2 + 1).trim(utc, utc.MONTH, api.trim_policy.TRIM_OUT) == p_0_3

        # diff_units
        assert api.UtcPeriod(t0, t3).diff_units(utc, utc.MONTH) == 3
        # trim null period, should trow
        p_null = api.UtcPeriod()
        try:
            p_null.trim(utc, utc.DAY)
            did_raise = False
        except RuntimeError:
            did_raise = True
        assert did_raise

        assert p_null.diff_units(utc, utc.DAY) == 0
        assert api.UtcPeriod(t3, t0).diff_units(utc, utc.MONTH) == -3

    def test_UtcPeriod_intersection(self):
        utc = api.Calendar()
        t0 = utc.time(2018, 1, 1, 0, 0, 0)
        t1 = utc.time(2018, 1, 2, 0, 0, 0)
        t2 = utc.time(2018, 1, 3, 0, 0, 0)
        t3 = utc.time(2018, 1, 4, 0, 0, 0)

        p0 = api.UtcPeriod(t0, t2)
        p1 = api.UtcPeriod(t1, t3)
        p2 = api.UtcPeriod(t0, t1)
        p3 = api.UtcPeriod(t2, t3)

        assert api.UtcPeriod.intersection(p0, p1) == api.UtcPeriod(t1, t2)
        assert not api.UtcPeriod.intersection(p0, p3).valid()
        assert not api.UtcPeriod.intersection(p2, p3).valid()
        assert api.intersection(p0, p1) == api.UtcPeriod(t1, t2)

    def test_swig_python_time(self):
        """
        This particular test is here to point out a platform specific bug detected
        on windows.

        """
        c1 = api.YMDhms(1969, 12, 31, 23, 0, 0)
        t = self.utc.time(c1)  # at this point, the value returned from c++ is correct, but on its
        # way through swig layer to python it goes via int32 and then to int64, proper signhandling
        #
        t_str = self.utc.to_string(t)  # this one just to show it's still working as it should internally
        self.assertEqual(t_str, "1969-12-31T23:00:00Z")
        self.assertEqual(t, api.deltahours(-1))

    def test_time_construct(self):
        self.assertAlmostEqual(api.time().seconds, 0.0, msg='default value should be 0.0')
        self.assertAlmostEqual(api.time(3123456).seconds, 3123456, msg='should be constructible from integer type')
        self.assertAlmostEqual(api.time(3.123456).seconds, 3.123456, msg='should be constructible from float type')
        self.assertAlmostEqual(api.time('1970-01-01T00:00:23Z').seconds, 23.0, msg='should be constructible from iso 8601 string')

    def test_time_compare(self):
        self.assertTrue(api.time(123) == api.time(123))
        self.assertTrue(api.time(123) == 123)
        self.assertTrue(api.time(123) != api.time(123.2))
        self.assertTrue(api.time(123) != 123.4)
        self.assertTrue(api.time(123) <= api.time(123.2))
        self.assertTrue(api.time(123) <= api.time(123))
        self.assertTrue(api.time(1234) >= 134)
        self.assertTrue(api.time(1234) >= 1234)
        self.assertTrue(api.time(123) < api.time(123.2))
        self.assertTrue(api.time(1234) > 134)

    def test_time_math(self):
        t = api.time
        self.assertAlmostEqual(t(2) + t(2), 4)
        self.assertAlmostEqual(t(2) + 2, 4)
        self.assertAlmostEqual(t(2) + 2.2, 4.2)

        self.assertAlmostEqual(t(2) - t(2), 0.0)
        self.assertAlmostEqual(t(2) - 2, 0)
        self.assertAlmostEqual(t(2) - 2.2, -0.2)

        self.assertAlmostEqual(t(2)*t(3), 6.0)
        self.assertAlmostEqual(t(2)*0.1, 0.2)

        self.assertAlmostEqual(t(2)/t(3), 0.666667)
        self.assertAlmostEqual(t(2)/0.1, 20.0)

        self.assertAlmostEqual(t(2)//t(3), 0)
        self.assertAlmostEqual(t(10)//3, 3)

        self.assertAlmostEqual(abs(t(-3)), 3)
        self.assertAlmostEqual(abs(t(3.2)), 3.2)

        self.assertAlmostEqual(t(10)%3, 1.0)

    def test_time_floor(self):
        self.assertAlmostEqual(math.floor(api.time(3.2)), 3.0)
        self.assertAlmostEqual(math.floor(api.time(-3.2)), -4.0)

    def test_time_round(self):
        self.assertAlmostEqual(round(api.time(3.2)), 3.0)
        self.assertAlmostEqual(round(api.time(-3.7)), -4.0)

    def test_time_cast(self):
        self.assertAlmostEqual(int(api.time(10.23)), 10)
        self.assertAlmostEqual(float(api.time(10.23)), 10.23)
        self.assertAlmostEqual(api.time(1.23).seconds, 1.23)

    def test_time_large_number(self):
        a = 1000.0*api.time(1534832966.984426)
        self.assertIsNotNone(a)
        b = api.time(a)/1000.0
        self.assertIsNotNone(b)
        sb = str(b)
        self.assertTrue(len(sb) > 0)
        pass

    def test_time_hash(self):
        t0 = api.time('2018-01-01T01:02:03Z')
        t1 = api.time('2018-01-01T01:02:03Z')
        h0 = hash(t0)
        h1 = hash(t1)
        self.assertEqual(h0, h1)
        d = {t0: 'A'}
        self.assertTrue(t0 in d)
        self.assertTrue(t1 in d)


if __name__ == "__main__":
    unittest.main()
