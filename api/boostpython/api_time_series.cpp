#include "boostpython_pch.h"

#include "core/utctime_utilities.h"
#include "core/time_axis.h"
#include "core/time_series.h"
#include "api/api.h"
#include "api/time_series.h"

namespace expose {
    using namespace shyft;
    using namespace shyft::core;
    using namespace boost::python;
    using namespace std;

    namespace py = boost::python;

    shyft::api::ats_vector quantile_map_forecast_5(vector<shyft::api::ats_vector> const & forecast_set, vector<double> const& set_weights, shyft::api::ats_vector const& historical_data, shyft::api::gta_t const&time_axis, utctime interpolation_start ) {
        return shyft::api::quantile_map_forecast(forecast_set, set_weights, historical_data, time_axis, interpolation_start);
    }
    shyft::api::ats_vector quantile_map_forecast_6(vector<shyft::api::ats_vector> const & forecast_set, vector<double> const& set_weights, shyft::api::ats_vector const& historical_data, shyft::api::gta_t const&time_axis, utctime interpolation_start, utctime interpolation_end) {
        return shyft::api::quantile_map_forecast(forecast_set, set_weights, historical_data, time_axis, interpolation_start, interpolation_end);
    }
    shyft::api::ats_vector quantile_map_forecast_7(vector<shyft::api::ats_vector> const & forecast_set, vector<double> const& set_weights, shyft::api::ats_vector const& historical_data, shyft::api::gta_t const&time_axis, utctime interpolation_start, utctime interpolation_end, bool interpolated_quantiles) {
        return shyft::api::quantile_map_forecast(forecast_set, set_weights, historical_data, time_axis, interpolation_start, interpolation_end, interpolated_quantiles);
    }


    static void expose_ats_vector() {
        using namespace shyft::api;
        typedef ats_vector(ats_vector::*m_double)(double)const;
        typedef ats_vector(ats_vector::*m_ts)(apoint_ts const&)const;
        typedef ats_vector(ats_vector::*m_tsv) (ats_vector const&)const;

        class_<ats_vector>("TsVector",
                doc_intro("A vector of time-series that supports ts-math operations.")
                doc_intro("")
                doc_intro("Like:")
                doc_intro("  number bin_op ts_vector -> ts_vector")
                doc_intro("  ts_vector bin_op ts_vector -> ts_vector ")
                doc_intro("  ts bin_op ts_vector -> ts_vector")
                doc_intro("  where bin_op is any of (*,/,+,-)")
                doc_intro("")
                doc_intro("In addition, .average(..),.integral(..),.accumulate(..),.time_shift(..), .percentiles(..)")
                doc_intro("  is also supported")
                doc_intro("")
                doc_intro("All operation returns a *new* ts-vector, containing the resulting expressions")
            )
            .def(vector_indexing_suite<ats_vector>())
            .def(init<ats_vector const&>(args("clone_me")))
            .def("values_at",&ats_vector::values_at_time,args("t"),
                 doc_intro("Computes the value at specified time t for all time-series")
                 doc_parameters()
                 doc_parameter("t","int","seconds since epoch 1970 UTC")
            )
            .def("percentiles",&ats_vector::percentiles,args("time_axis","percentiles"),
                doc_intro("Calculate the percentiles, NIST R7, excel,R definition, of the timeseries")
                doc_intro("over the specified time-axis.")
                doc_intro("The time-series point_fx interpretation is used when performing")
                doc_intro("the true-average over the time_axis periods.")
                doc_parameters()
                doc_parameter("percentiles","IntVector","A list of numbers,[ 0, 25,50,-1,75,100] will return 6 time-series,\n -1 -> arithmetic average\n -1000 -> min extreme value\n +1000 max extreme value")
                doc_parameter("time_axis","TimeAxis","The time-axis used when applying true-average to the time-series")
                doc_returns("calculated_percentiles","TsVector","Time-series list with evaluated percentile results, same length as input")
            )
            .def("percentiles",&ats_vector::percentiles_f,args("time_axis","percentiles"),
                doc_intro("Calculate the percentiles, NIST R7, excel,R definition, of the timeseries")
                doc_intro("over the specified time-axis.")
                doc_intro("The time-series point_fx interpretation is used when performing")
                doc_intro("the true-average over the time_axis periods.")
                doc_parameters()
                doc_parameter("percentiles","IntVector","A list of numbers,[ 0, 25,50,-1,75,100] will return 6 time-series,\n -1 -> arithmetic average\n -1000 -> min extreme value\n +1000 max extreme value")
                doc_parameter("time_axis","TimeAxisFixedDeltaT","The time-axis used when applying true-average to the time-series")
                doc_returns("calculated_percentiles","TsVector","Time-series list with evaluated percentile results, same length as input")
            )
            .def("slice",&ats_vector::slice,args("indexes"),
                 doc_intro("returns a slice of self, specified by indexes")
                 doc_parameters()
                 doc_parameter("indexes","IntVector","the indicies to pick out from self, if indexes is empty, then all is returned")
                 doc_returns("slice","TsVector","a new TsVector, with content according to indexes specified")
            )
            .def("abs", &ats_vector::abs,
                doc_intro("create a new ts-vector, with all members equal to abs(self")
                doc_returns("tsv", "TsVector", "a new TsVector expression, that will provide the abs-values of self.values")
            )

            .def("average", &ats_vector::average, args("ta"),
                doc_intro("create a new vector of ts that is the true average of self")
                doc_intro("over the specified time-axis ta.")
                doc_parameters()
                doc_parameter("ta","TimeAxis","time-axis that specifies the periods where true-average is applied")
                doc_returns("tsv","TsVector","a new time-series expression, that will provide the true-average when requested")
                doc_notes()
                doc_note("the self point interpretation policy is used when calculating the true average")
			)
            .def("integral", &ats_vector::integral, args("ta"),
                doc_intro("create a new vector of ts that is the true integral of self")
                doc_intro("over the specified time-axis ta.")
                doc_intro(" defined as integral of the non-nan part of each time-axis interval")
                doc_parameters()
                doc_parameter("ta", "TimeAxis", "time-axis that specifies the periods where true-integral is applied")
                doc_returns("tsv", "TsVector", "a new time-series expression, that will provide the true-integral when requested")
                doc_notes()
                doc_note("the self point interpretation policy is used when calculating the true average")
            )
            .def("accumulate", &ats_vector::accumulate, args("ta"),
                doc_intro("create a new  vector of ts where each i'th value is the ")
                doc_intro("    integral f(t) *dt, from t0..ti,")
                doc_intro("given the specified time-axis ta")
                doc_parameters()
                doc_parameter("ta","TimeAxis","time-axis that specifies the periods where accumulated integral is applied")
                doc_returns("tsv","TsVector","a new time-series expression, that will provide the accumulated values when requested")
                doc_notes()
                doc_note("the self point interpretation policy is used when calculating the accumulated values")
            )
			.def("time_shift", &ats_vector::time_shift,args("delta_t"),
				doc_intro("create a new vector of ts that is a the time-shift'ed  version of self")
				doc_parameters()
                doc_parameter("delta_t","int","number of seconds to time-shift, positive values moves forward")
				doc_returns("tsv","TsVector",	"a new time-series, that appears as time-shifted version of self")
			)
            .def("extend", &ats_vector::extend_ts, (py::arg("ts"), py::arg("split_policy") = extend_ts_split_policy::EPS_LHS_LAST, py::arg("fill_policy") = extend_ts_fill_policy::EPF_NAN, py::arg("split_at") = utctime(0), py::arg("fill_value") = shyft::nan),
                doc_intro("create a new ats_vector where all time-series are extended by ts")
                doc_parameters()
                doc_parameter("ts", "TimeSeries", "time-series to extend each time-series in self with")
                doc_parameter("split_policy", "extend_ts_split_policy", "policy determining where to split between self and ts")
                doc_parameter("fill_policy", "extend_ts_fill_policy", "policy determining how to fill any gap between self and ts")
                doc_parameter("split_at", "utctime", "time at which to split if split_policy == EPS_VALUE")
                doc_parameter("fill_value", "float", "value to fill any gap with if fill_policy == EPF_FILL")
                doc_returns("new_ts_vec" ,"TsVector", "a new time-series vector where all time-series in self have been extended by ts")
            )
            .def("extend", &ats_vector::extend_vec, (py::arg("ts"), py::arg("split_policy") = extend_ts_split_policy::EPS_LHS_LAST, py::arg("fill_policy") = extend_ts_fill_policy::EPF_NAN, py::arg("split_at") = utctime(0), py::arg("fill_value") = shyft::nan),
                doc_intro("create a new ats_vector where all ts' are extended by the matching ts from ts_vec")
                doc_parameters()
                doc_parameter("ts_vec", "TsVector", "time-series vector to extend time-series in self with")
                doc_parameter("split_policy", "extend_ts_split_policy", "policy determining where to split between self and ts")
                doc_parameter("fill_policy", "extend_ts_fill_policy", "policy determining how to fill any gap between self and ts")
                doc_parameter("split_at", "utctime", "time at which to split if split_policy == EPS_VALUE")
                doc_parameter("fill_value", "float", "value to fill any gap with if fill_policy == EPF_FILL")
                doc_returns("new_ts_vec" ,"TsVector", "a new time-series vector where all time-series in self have been extended by the corresponding time-series in ts_vec")
            )
            .def("min",(m_double)&ats_vector::min,args("number"),"returns min of vector and a number")
            .def("min", (m_ts)&ats_vector::min, args("ts"), "returns min of ts-vector and a ts")
            .def("min", (m_tsv)&ats_vector::min, args("tsv"), "returns min of ts-vector and another ts-vector")
            .def("max", (m_double)&ats_vector::max, args("number"), "returns max of vector and a number")
            .def("max", (m_ts)&ats_vector::max, args("ts"), "returns max of ts-vector and a ts")
            .def("max", (m_tsv)&ats_vector::max, args("tsv"), "returns max of ts-vector and another ts-vector")
            .def("forecast_merge",&ats_vector::forecast_merge,args("lead_time","fc_interval"),
                 doc_intro("merge the forecasts in this vector into a time-series that is constructed")
                 doc_intro("taking a slice of length fc_interval starting lead_time into each of the forecasts")
                 doc_intro("of this time-series vector.")
                 doc_intro("The content of the vector should be ordered in forecast-time, each entry at least")
                 doc_intro("fc_interval separated from the previous.")
                 doc_intro("If there is missing forecasts (larger than fc_interval between two forecasts) this is")
                 doc_intro("automagically repaired using extended slices from the existing forecasts")
                 doc_parameters()
                 doc_parameter("lead_time","int","start slice number of seconds from t0 of each forecast")
                 doc_parameter("fc_interval","int","length of each slice in seconds, and thus also gives the forecast-interval separation")
                 doc_returns("merged time-series","TimeSeries","A merged forecast time-series")
                 )
             .def("nash_sutcliffe",&ats_vector::nash_sutcliffe,args("observation_ts","lead_time","delta_t","n"),
                doc_intro("Computes the nash-sutcliffe (wiki nash-sutcliffe) criteria between the")
                doc_intro("observation_ts over the slice of each time-series in the vector.")
                doc_intro("The slice for each ts is specified by the lead_time, delta_t and n")
                doc_intro("parameters. The function is provided to ease evaluation of forecast performance")
                doc_intro("for different lead-time periods into each forecast.")
                doc_intro("The returned value range is 1.0 for perfect match -oo for no match, or nan if observations is constant or data missing")
                doc_parameters()
                doc_parameter("observation_ts","TimeSeries","the observation time-series")
                doc_parameter("lead_time","int","number of seconds lead-time offset from each ts .time(0)")
                doc_parameter("delta_t","int","delta-time seconds to average as basis for n.s. simulation and observation values")
                doc_parameter("n","int","number of time-steps of length delta_t to slice out of each forecast/simulation ts")
                doc_returns("nash-sutcliffe value","double","the nash-sutcliffe criteria evaluated over all time-series in the TsVector for the specified lead-time, delta_t and number of elements")
                  doc_notes()
                  doc_see_also("nash_sutcliffe_goal_function")
             )
             .def("average_slice",&ats_vector::average_slice,args("lead_time","delta_t","n"),
                doc_intro("Returns a ts-vector with the average time-series of the specified slice")
                doc_intro("The slice for each ts is specified by the lead_time, delta_t and n")
                doc_intro("parameters. ")
                doc_parameters()
                doc_parameter("lead_time","int","number of seconds lead-time offset from each ts .time(0)")
                doc_parameter("delta_t","int","delta-time seconds to average as basis for n.s. simulation and observation values")
                doc_parameter("n","int","number of time-steps of length delta_t to slice out of each forecast/simulation ts")
                doc_returns("ts_vector_sliced","TsVector","a ts-vector with average ts of each slice specified.")
                  doc_notes()
                  doc_see_also("nash_sutcliffe,forecast_merge")
             )
            // defining vector math-operations goes here
            .def(-self)
            .def(self*double())
            .def(double()*self)
            .def(self*self)
            .def(shyft::api::apoint_ts()*self)
            .def(self*shyft::api::apoint_ts())

            .def(self/double())
            .def(double()/self)
            .def(self/self)
            .def(shyft::api::apoint_ts()/self)
            .def(self/shyft::api::apoint_ts())

            .def(self+double())
            .def(double()+self)
            .def(self+self)
            .def(shyft::api::apoint_ts()+self)
            .def(self+shyft::api::apoint_ts())

            .def(self-double())
            .def(double()-self)
            .def(self-self)
            .def(shyft::api::apoint_ts()-self)
            .def(self-shyft::api::apoint_ts())
            ;
            // expose min-max functions:
            typedef ats_vector(*f_atsv_double)(ats_vector const &, double b);
            typedef ats_vector(*f_double_atsv)(double b, ats_vector const &a);
            typedef ats_vector(*f_atsv_ats)(ats_vector const &, apoint_ts const& );
            typedef ats_vector(*f_ats_atsv)(apoint_ts const &b, ats_vector const& a);
            typedef ats_vector(*f_atsv_atsv)(ats_vector const &b, ats_vector const &a);

            def("min", (f_ats_atsv)min, args("ts", "ts_vector"), "return minimum of ts and ts_vector");
            def("min", (f_atsv_ats)min, args("ts_vector", "ts"), "return minimum of ts_vector and ts");
            def("min", (f_atsv_double)min, args("ts_vector", "number"), "return minimum of ts_vector and number");
            def("min", (f_double_atsv)min, args("number","ts_vector"), "return minimum of number and ts_vector");
            def("min", (f_atsv_atsv)min, args("a", "b"), "return minimum of ts_vectors a and b (requires equal size!)");
            def("max", (f_ats_atsv)max, args("ts", "ts_vector"), "return max of ts and ts_vector");
            def("max", (f_atsv_ats)max, args("ts_vector", "ts"), "return max of ts_vector and ts");
            def("max", (f_atsv_double)max, args("ts_vector", "number"), "return max of ts_vector and number");
            def("max", (f_double_atsv)max, args("number", "ts_vector"), "return max of number and ts_vector");
            def("max", (f_atsv_atsv)max, args("a", "b"), "return max of ts_vectors a and b (requires equal size!)");

            // we also need a vector of ats_vector for quantile_map_forecast function
            typedef std::vector<ats_vector> TsVectorSet;
            class_<TsVectorSet>("TsVectorSet",
                doc_intro("A set of TsVector")
                doc_intro("")
                doc_see_also("quantile_map_forecast,TsVector")
            )
            .def(vector_indexing_suite<TsVectorSet>())
            .def(init<TsVectorSet const&>(args("clone_me")))
            ;
            const char* qm_doc =
                doc_intro("Computes the quantile-mapped forecast from the supplied input.")
                doc_intro(" TBD:detailed description with references ")
                doc_parameters()
                doc_parameter("forecast_sets", "TsVectorSet", "forecast sets, each of them a TsVector with n forecasts (might differ in size and length)")
                doc_parameter("set_weights", "DoubleVector", "a weight for each of the forecast set in forecast-sets,correlated by same index)")
                doc_parameter("historical_data", "TsVector", "historical time-series that should cover the requested time-axis")
                doc_parameter("time_axis", "TimeAxis", "the time-axis that the resulting time-series are mapped into")
                doc_parameter("interpolation_start", "int", "time where the historical to forecast interpolation should start, 1970 utc seconds since epoch")
                doc_parameter("interpolation_end", "int", "time where the interpolation should end, if no_utctime, use end of forecast-set")
                doc_parameter("interpolated_quantiles", "bool", "whether the quantile values should be interpolated or assigned the values lower than or equal to the current quantile")
                doc_returns("qm_forecast", "TsVector", "quantile mapped forecast with the requested time-axis")
                ;

            def("quantile_map_forecast",quantile_map_forecast_5,
                args("forecast_sets","set_weights","historical_data","time_axis","interpolation_start"),
                qm_doc
                );
            def("quantile_map_forecast", quantile_map_forecast_6,
                args("forecast_sets", "set_weights", "historical_data", "time_axis", "interpolation_start", "interpolation_end"),
                qm_doc
            );
            def("quantile_map_forecast", quantile_map_forecast_7,
                args("forecast_sets", "set_weights", "historical_data", "time_axis", "interpolation_start", "interpolation_end", "interpolated_quantiles"),
                qm_doc
            );

    }

    #define DEF_STD_TS_STUFF() \
            .def("point_interpretation",&pts_t::point_interpretation,"returns the point interpretation policy")\
            .def("set_point_interpretation",&pts_t::set_point_interpretation,args("policy"),"set new policy")\
            .def("value",&pts_t::value,args("i"),"returns the value at the i'th time point")\
            .def("time",&pts_t::time,args("i"),"returns the time at the i'th point")\
            .def("get",&pts_t::get,args("t"),"returns the i'th Point ")\
            .def("set",&pts_t::set,args("i","v"),"set the i'th value")\
            .def("fill",&pts_t::fill,args("v"),"all values with v")\
            .def("scale_by",&pts_t::scale_by,args("v"),"scale all values by specified factor")\
            .def("size",&pts_t::size,"returns number of points")\
            .def("index_of",&pts_t::index_of,args("t"),"return the index of the intervall that contains t, or npos if not found")\
            .def("total_period",&pts_t::total_period,"returns the total period covered by the time-axis of this time-series")\
            .def("__call__",&pts_t::operator(),args("t"),"return the f(t) value for the time-series")


    template <class TA>
    static void point_ts(const char *ts_type_name,const char *doc) {
        typedef time_series::point_ts<TA> pts_t;
        class_<pts_t,bases<>,shared_ptr<pts_t>,boost::noncopyable>(ts_type_name, doc)
            .def(init<const TA&,const vector<double>&,optional<time_series::ts_point_fx>>(args("ta","v","policy"),"constructs a new timeseries from timeaxis and points"))
            .def(init<const TA&,double,optional<time_series::ts_point_fx>>(args("ta","fill_value","policy"),"constructs a new timeseries from timeaxis and fill-value"))
            DEF_STD_TS_STUFF()
            .def_readonly("v",&pts_t::v,"the point vector<double>, same as .values, kept around for backward compatibility")
			.def("get_time_axis", &pts_t::time_axis, "returns the time-axis", return_internal_reference<>()) // have to use func plus init.py fixup due to boost py policy
            ;
    }


    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(point_ts_overloads     ,shyft::api::TsFactory::create_point_ts,4,5);
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(time_point_ts_overloads,shyft::api::TsFactory::create_time_point_ts,3,4);

    static void TsFactory() {
        class_<shyft::api::TsFactory>("TsFactory","TsFactory is used to create point time-series that exposes the ITimeSeriesOfPoint interface, using the internal ts-implementations")
            .def("create_point_ts",&shyft::api::TsFactory::create_point_ts,point_ts_overloads())//args("n","tStart","dt","values","interpretation"),"returns a new fixed interval ts from specified arguments")
            .def("create_time_point_ts",&shyft::api::TsFactory::create_time_point_ts,time_point_ts_overloads())//args("period","times","values","interpretation"),"return a point ts from specified arguments")
            ;
    }


    static void expose_apoint_ts() {
        using namespace shyft::api;

        typedef shyft::api::apoint_ts pts_t;
        typedef pts_t (pts_t::*self_dbl_t)(double) const;
        typedef pts_t (pts_t::*self_ts_t)(const pts_t &)const;

        self_dbl_t min_double_f=&pts_t::min;
        self_ts_t  min_ts_f =&pts_t::min;
        self_dbl_t max_double_f=&pts_t::max;
        self_ts_t  max_ts_f =&pts_t::max;
        typedef shyft::api::ts_bind_info TsBindInfo;
        class_<TsBindInfo>("TsBindInfo",
            "TsBindInfo gives information about the time-series and it's binding\n"
            "represented by encoded string reference\n"
            "Given that you have a concrete ts,\n"
            "you can bind that the bind_info.ts\n"
            "using bind_info.ts.bind()\n"
            "see also Timeseries.find_ts_bind_info() and Timeseries.bind()\n"
            )
            .def_readwrite("id", &shyft::api::ts_bind_info::reference, "a unique id/url that identifies a time-series in a ts-database/file-store/service")
            .def_readwrite("ts", &shyft::api::ts_bind_info::ts,"the ts, provides .bind(another_ts) to set the concrete values")
            ;

        typedef vector<TsBindInfo> TsBindInfoVector;
        class_<TsBindInfoVector>("TsBindInfoVector", "A vector of TsBindInfo\nsee also TsBindInfo")
            .def(vector_indexing_suite<TsBindInfoVector>())
            ;


		class_<shyft::api::apoint_ts>("TimeSeries",
                doc_intro("A time-series providing mathematical and statistical operations and functionality.")
                doc_intro("")
                doc_intro("A time-series can be an expression, or a concrete point time-series.")
                doc_intro("All time-series do have a time-axis, values, and a point fx policy.")
                doc_intro("")
                doc_intro("The time-series can provide a value for all the intervals, and the point_fx policy")
                doc_intro("defines how the values should be interpreted:")
                doc_intro("POINT_INSTANT_VALUE:")
                doc_intro("    the point value is valid at the start of the period, linear between points")
                doc_intro("    -extend flat from last point to +oo, nan before first value")
                doc_intro("    typical for state-variables, like water-level, temperature measured at 12:00 etc.")
                doc_intro("POINT_AVERAGE_VALUE:")
                doc_intro("    the point represents an average or constant value over the period")
                doc_intro("    typical for model-input and results, precipitation mm/h, discharge m^3/s")
                doc_intro("")
                doc_intro("Example:")
                doc_intro("import numpy as np")
                doc_intro("from shyft.api import Calendar,deltahours,TimeAxis,TimeSeries,POINT_AVERAGE_VALUE as fx_avg,DoubleVector as dv")
                doc_intro("")
                doc_intro("utc = Calendar()  # ensure easy consistent explicit handling of calendar and time")
                doc_intro("ta = TimeAxis(utc.time(2016, 9, 1, 8, 0, 0), deltahours(1), 10)  # create a time-axis to use")
                doc_intro("a = TimeSeries(ta, dv.from_numpy(np.linspace(0, 10, num=len(ta))), fx_avg)")
                doc_intro("b = TimeSeries(ta, dv.from_numpy(np.linspace(0,  1, num=len(ta))), fx_avg)")
                doc_intro("c = a + b*3.0  # c is now an expression, time-axis is the overlap of a and b, lazy evaluation")
                doc_intro("c_values = c.values.to_numpy()  # compute and extract the values, as numpy array")
                doc_intro("")
                doc_intro("The TimeSeries functionality includes ")
                doc_intro(" resampling:average,accumulate,time_shift")
                doc_intro(" statistics: min/max,correlation by nash-sutcliffe, kling-gupta")
                doc_intro(" filtering: convolution,average")
                doc_intro(" partitioning and percentiles ")
                doc_intro("Please check notebooks, examples and api-tests for usage.")
                doc_see_also("TimeAxis,DoubleVector,Calendar,point_interpretation_policy")

            )
			.def(init<const time_axis::generic_dt&, double, optional<time_series::ts_point_fx> >(args("ta", "fill_value", "point_fx"), "construct a timeseries with timeaxis ta and specified fill-value, default point_fx=POINT_INSTANT_VALUE"))
			.def(init<const time_axis::generic_dt&, const std::vector<double>&, optional<time_series::ts_point_fx> >(args("ta", "values", "point_fx"), "construct a timeseries timeaxis ta and corresponding values, default point_fx=POINT_INSTANT_VALUE"))

			.def(init<const time_axis::fixed_dt&, double, optional<time_series::ts_point_fx> >(args("ta", "fill_value", "point_fx"), "construct a timeseries with timeaxis ta and specified fill-value, default point_fx=POINT_INSTANT_VALUE"))
			.def(init<const time_axis::fixed_dt&, const std::vector<double>&, optional<time_series::ts_point_fx> >(args("ta", "values", "point_fx"), "construct a timeseries timeaxis ta and corresponding values, default point_fx=POINT_INSTANT_VALUE"))

			.def(init<const time_axis::point_dt&, double, optional<time_series::ts_point_fx> >(args("ta", "fill_value", "point_fx"), "construct a timeseries with timeaxis ta and specified fill-value, default point_fx=POINT_INSTANT_VALUE"))
			.def(init<const time_axis::point_dt&, const std::vector<double>&, optional<time_series::ts_point_fx> >(args("ta", "values", "point_fx"), "construct a timeseries timeaxis ta and corresponding values, default point_fx=POINT_INSTANT_VALUE"))
            .def(init<const shyft::api::rts_t &>(args("core_result_ts"),"construct a timeseries from a shyft core time-series, to allow full ts-functionality in python"))

			.def(init<const shyft::api::apoint_ts&>(args("clone"), "creates a shallow copy of clone"))

			.def(init<const vector<double>&, utctimespan, const time_axis::generic_dt&>(args("pattern", "dt", "ta"), "construct a timeseries given a equally spaced dt pattern and a timeaxis ta"))
			.def(init<const vector<double>&, utctimespan,utctime, const time_axis::generic_dt&>(args("pattern", "dt","t0", "ta"), "construct a timeseries given a equally spaced dt pattern, starting at t0, and a timeaxis ta"))
            .def(init<std::string>(args("ts_id"),
                doc_intro("constructs a bind-able ts,")
                doc_intro("providing a symbolic possibly unique id that at a later time")
                doc_intro("can be bound, using the .bind(ts) method to concrete values")
                doc_intro("if the ts is used as ts, like size(),.value(),time() before it")
                doc_intro("is bound, then a runtime-exception is raised")
                doc_parameters()
                doc_parameter("ts_id","str","url-like identifier for the time-series,notice that shyft://<container>/<path> is for shyft-internal store")
                )
            )
            .def(init<std::string,const apoint_ts&>(args("ts_id","bts"),
                doc_intro("constructs a ready bound ts,")
                doc_intro("providing a symbolic possibly unique id that at a later time")
                doc_intro("can be used to correlate with back-end store\n")
                doc_parameters()
                doc_parameter("ts_id","str","url-type of id, notice that shyft://<container>/<path> is for shyft-internal store")
                doc_parameter("bts","TimeSeries","A concrete time-series, with point_fx policy, time_axis and values")
                )
            )
            .def("ts_id",&apoint_ts::id,
                doc_intro("returns ts_id of symbolic ts, or empty string if not symbolic ts")
                doc_returns("ts_id","str","url-like ts_id as passed to constructor or empty if the ts is not a ts with ts_id")
                doc_see_also("TimeSeries('url://like/id'),TimeSeries('url://like/id',ts_with_values)")
            )
			DEF_STD_TS_STUFF()
			// expose time_axis sih: would like to use property, but no return value policy, so we use get_ + fixup in init.py
			.def("get_time_axis", &shyft::api::apoint_ts::time_axis, "returns the time-axis", return_internal_reference<>())
			.add_property("values", &shyft::api::apoint_ts::values, "return the values (possibly calculated on the fly)")
			// operators
			.def(self * self)
			.def(double() * self)
			.def(self * double())

			.def(self + self)
			.def(double() + self)
			.def(self + double())

			.def(self / self)
			.def(double() / self)
			.def(self / double())

			.def(self - self)
			.def(double() - self)
			.def(self - double())

			.def(-self)
            .def(operator!(self))
            .def("abs", &shyft::api::apoint_ts::abs,
                doc_intro("create a new ts, abs(self")
                doc_returns("ts", "TimeSeries", "a new time-series expression, that will provide the abs-values of self.values")
            )
			.def("average", &shyft::api::apoint_ts::average, args("ta"),
                doc_intro("create a new ts that is the true average of self")
                doc_intro("over the specified time-axis ta.")
                doc_parameters()
                doc_parameter("ta","TimeAxis","time-axis that specifies the periods where true-average is applied")
                doc_returns("ts","TimeSeries","a new time-series expression, that will provide the true-average when requested")
                doc_notes()
                doc_note("the self point interpretation policy is used when calculating the true average")
			)
            .def("integral", &shyft::api::apoint_ts::integral, args("ta"),
                doc_intro("create a new ts that is the true integral of self")
                doc_intro("over the specified time-axis ta.")
                doc_intro(" defined as integral of the non-nan part of each time-axis interval")
                doc_parameters()
                doc_parameter("ta", "TimeAxis", "time-axis that specifies the periods where true-integral is applied")
                doc_returns("ts", "TimeSeries", "a new time-series expression, that will provide the true-integral when requested")
                doc_notes()
                doc_note("the self point interpretation policy is used when calculating the true average")
            )
            .def("accumulate", &shyft::api::apoint_ts::accumulate, args("ta"),
                doc_intro("create a new ts where each i'th value is the ")
                doc_intro("    integral f(t) *dt, from t0..ti,")
                doc_intro("given the specified time-axis ta")
                doc_parameters()
                doc_parameter("ta","TimeAxis","time-axis that specifies the periods where accumulated integral is applied")
                doc_returns("ts","TimeSeries","a new time-series expression, that will provide the accumulated values when requested")
                doc_notes()
                doc_note("the self point interpretation policy is used when calculating the accumulated values")
            )
			.def("time_shift", &shyft::api::apoint_ts::time_shift,args("delta_t"),
				doc_intro("create a new ts that is a the time-shift'ed  version of self")
				doc_parameters()
                doc_parameter("delta_t","int","number of seconds to time-shift, positive values moves forward")
				doc_returns("ts","TimeSeries",	"a new time-series, that appears as time-shifted version of self")
			)
            .def("convolve_w", &shyft::api::apoint_ts::convolve_w, args("weights", "policy"),
                doc_intro("create a new ts that is the convolved ts with the given weights list")
                doc_parameters()
                doc_parameter("weights","DoubleVector","the weights profile, use DoubleVector.from_numpy(...) to create these.\n"
                                "\t it's the callers responsibility to ensure the sum of weights are 1.0\n")
                doc_parameter("policy","convolve_policy","(.USE_FIRST|USE_ZERO|USE_NAN)\n"
                "\t Specifies how to handle initial weight.size()-1 values\n")
                doc_returns("ts","TimeSeries","a new time-series that is evaluated on request to the convolution of self")
                doc_see_also("ConvolvePolicy")
            )
			.def("rating_curve", &shyft::api::apoint_ts::rating_curve, py::arg("rc_param"),
				doc_intro("Create a new TimeSeries that is computed using a RatingCurveParameter instance.")
				doc_intro("")
				doc_intro("Example:")
				doc_intro("")
				doc_intro(">>> import numpy as np")
				doc_intro(">>> from shyft.api import (")
				doc_intro("...     utctime_now, deltaminutes,")
				doc_intro("...     TimeAxis, TimeSeries,")
				doc_intro("...     RatingCurveFunction, RatingCurveParameters")
				doc_intro("... )")
				doc_intro(">>>")
				doc_intro(">>> # parameters")
				doc_intro(">>> t0 = utctime_now()")
				doc_intro(">>> dt = deltaminutes(30)")
				doc_intro(">>> n = 48*2")
				doc_intro(">>>")
				doc_intro(">>> # make rating function, each with two segments")
				doc_intro(">>> rcf_1 = RatingCurveFunction()")
				doc_intro(">>> rcf_1.add_segment(0, 2, 0, 1)    # add segment from level 0, computing f(h) = 2*(h - 0)**1")
				doc_intro(">>> rcf_1.add_segment(5.3, 1, 1, 1.4)  # add segment from level 5.3, computing f(h) = 1.3*(h - 1)**1.4")
				doc_intro(">>> rcf_2 = RatingCurveFunction()")
				doc_intro(">>> rcf_2.add_segment(0, 1, 1, 1)    # add segment from level 0, computing f(h) = 1*(h - 1)**1")
				doc_intro(">>> rcf_2.add_segment(8.0, 0.5, 0, 2)  # add segment from level 8.0, computing f(h) = 0.5*(h - 0)**2")
				doc_intro(">>>")
				doc_intro(">>> # add rating curves to a parameter pack")
				doc_intro(">>> rcp = RatingCurveParameters()")
				doc_intro(">>> rcp.add_curve(t0, rcf_1)  # rcf_1 is active from t0")
				doc_intro(">>> rcp.add_curve(t0+dt*n//2, rcf_2)  # rcf_2 takes over from t0 + dt*n/2")
				doc_intro(">>>")
				doc_intro(">>> # create a time-axis/-series")
				doc_intro(">>> ta = TimeAxis(t0, dt, n)")
				doc_intro(">>> ts = TimeSeries(ta, np.linspace(0, 12, n))")
				doc_intro(">>> rc_ts = ts.rating_curve(rcp)  # create a new time series computed using the rating curve functions")
				doc_intro(">>>")
				doc_parameters()
				doc_parameter("rc_param", "RatingCurveParameter", "RatingCurveParameter instance.")
				doc_returns("rcts", "TimeSeries", "A new TimeSeries computed using self and rc_param.")
			)
            .def("extend", &shyft::api::apoint_ts::extend, (py::arg("ts"), py::arg("split_policy") = extend_ts_split_policy::EPS_LHS_LAST, py::arg("fill_policy") = extend_ts_fill_policy::EPF_NAN, py::arg("split_at") = utctime(0), py::arg("fill_value") = shyft::nan),
                doc_intro("create a new time-series that is self extended with ts")
                doc_parameters()
                doc_parameter("ts", "TimeSeries", "time-series to extend self with, only values after both the start of self, and split_at is used")
                doc_parameter("split_policy", "extend_split_policy", "policy determining where to split between self and ts")
                doc_parameter("fill_policy", "extend_fill_policy", "policy determining how to fill any gap between self and ts")
                doc_parameter("split_at", "utctime", "time at which to split if split_policy == EPS_VALUE")
                doc_parameter("fill_value", "float", "value to fill any gap with if fill_policy == EPF_FILL")
                doc_returns("extended_ts" ,"TimeSeries", "a new time-series that is the extension of self with ts")
            )
            .def("min",min_double_f,args("number"),"create a new ts that contains the min of self and number for each time-step")
            .def("min",min_ts_f,args("ts_other"),"create a new ts that contains the min of self and ts_other")
            .def("max",max_double_f,args("number"),"create a new ts that contains the max of self and number for each time-step")
            .def("max",max_ts_f,args("ts_other"),"create a new ts that contains the max of self and ts_other")
            //.def("max",max_stat_ts_ts_f,args("ts_a","ts_b"),"create a new ts that is the max(ts_a,ts_b)").staticmethod("max")
            //.def("min",min_stat_ts_ts_f,args("ts_a","ts_b"),"create a new ts that is the max(ts_a,ts_b)").staticmethod("min")
			.def("partition_by",&shyft::api::apoint_ts::partition_by,
                args("calendar","t", "partition_interval", "n_partitions","common_t0"),
				doc_intro("convert ts to a list of n_partitions partition-ts.")
				doc_intro("each partition covers partition_interval, starting from utctime t")
				doc_parameters()
				doc_parameter("cal","Calendar","The calendar to use, typically utc")
				doc_parameter("t","utctime","specifies where to pick the first partition")
				doc_parameter("partition_interval","utctimespan","the length of each partition, Calendar.YEAR,Calendar.DAY etc.")
				doc_parameter("n_partitions","int","number of partitions")
				doc_parameter("common_t0","utctime","specifies the time to correlate all the partitions")
				doc_returns("ts-partitions","TsVector","with length n_partitions, each ts is time-shifted and averaged expressions")
                doc_see_also("time_shift,average,TsVector")
				)
            .def("bind",&shyft::api::apoint_ts::bind,args("bts"),
                doc_intro("given that this ts,self, is a bind-able ts (aref_ts)")
                doc_intro("and that bts is a concrete point TimeSeries, make")
                doc_intro("a *copy* of bts and use it as representation")
                doc_intro("for the values of this ts")
                doc_parameters()
                doc_parameter("bts","TimeSeries","a concrete point ts, with time-axis and values")
                doc_notes()
                doc_note("raises runtime_error if any of preconditions is not true")
                doc_see_also("find_ts_bind_info,TimeSeries('a-ref-string')")
            )
            .def("bind_done",&shyft::api::apoint_ts::do_bind,
                 doc_intro("after bind operations on unbound time-series of an expression is done, call bind_done()")
                 doc_intro("to prepare the expression for use")
                 doc_notes()
                 doc_note("Usually this is done automatically by the dtss framework, but if not using dtss")
                 doc_note("this function is needed *after* the symbolic ts's are bound")
                 doc_see_also(".bind(), .find_ts_bind_info(), needs_bind()")
            )
            .def("needs_bind",&shyft::api::apoint_ts::needs_bind,
                 doc_intro("returns true if there are any unbound time-series in the expression")
                 doc_intro("this time-series represent")
                 doc_see_also(".find_ts_bind_info(),bind() and bind_done()")

            )
            .def("find_ts_bind_info",&shyft::api::apoint_ts::find_ts_bind_info,
                doc_intro("recursive search through the expression that this ts represents,\n")
                doc_intro("and return a list of TsBindInfo that can be used to\n")
                doc_intro("inspect and possibly 'bind' to ts-values \ref bind.\n")
                doc_returns("bind_info","TsBindInfoVector","A list of BindInfo where each entry contains a symbolic-ref and a ts that needs binding")
                doc_see_also("bind() method")

            )
            .def("serialize",&shyft::api::apoint_ts::serialize_to_bytes,
                "convert ts (expression) into a binary blob\n"
            )
            .def("deserialize",&shyft::api::apoint_ts::deserialize_from_bytes,args("blob"),
               "convert a blob, as returned by .serialize() into a Timeseries"
            ).staticmethod("deserialize")

        ;
        typedef shyft::api::apoint_ts (*avg_func_t)(const shyft::api::apoint_ts&,const shyft::time_axis::generic_dt&);
        typedef shyft::api::apoint_ts(*int_func_t)(const shyft::api::apoint_ts&, const shyft::time_axis::generic_dt&);
        avg_func_t avg=shyft::api::average;
        int_func_t intfnc = shyft::api::integral;
		avg_func_t acc = shyft::api::accumulate;
        def("average",avg,args("ts","time_axis"),"creates a true average time-series of ts for intervals as specified by time_axis");
        def("integral", intfnc, args("ts", "time_axis"), "creates a true integral time-series of ts for intervals as specified by time_axis");
		def("accumulate", acc, args("ts", "time_axis"), "create a new ts that is the integral f(t) *dt, t0..ti, the specified time-axis");
        //def("max",shyft::api::max,(boost::python::arg("ts_a"),boost::python::arg("ts_b")),"creates a new time-series that is the max of the supplied ts_a and ts_b");

        typedef shyft::api::apoint_ts (*ts_op_ts_t)(const shyft::api::apoint_ts&a, const shyft::api::apoint_ts&b);
        typedef shyft::api::apoint_ts (*double_op_ts_t)(double, const shyft::api::apoint_ts&b);
        typedef shyft::api::apoint_ts (*ts_op_double_t)(const shyft::api::apoint_ts&a, double);

        ts_op_ts_t max_ts_ts         = shyft::api::max;
        double_op_ts_t max_double_ts = shyft::api::max;
        ts_op_double_t max_ts_double = shyft::api::max;
        def("max",max_ts_ts    ,args("ts_a","ts_b"),"returns a new ts as max(ts_a,ts_b)");
        def("max",max_double_ts,args("a"   ,"ts_b"),"returns a new ts as max(a,ts_b)");
        def("max",max_ts_double,args("ts_a","b"   ),"returns a new ts as max(ts_a,b)");

        ts_op_ts_t min_ts_ts         = shyft::api::min;
        double_op_ts_t min_double_ts = shyft::api::min;
        ts_op_double_t min_ts_double = shyft::api::min;
        def("min",min_ts_ts    ,args("ts_a","ts_b"),"returns a new ts as min(ts_a,ts_b)");
        def("min",min_double_ts,args("a"   ,"ts_b"),"returns a new ts as min(a,ts_b)");
        def("min",min_ts_double,args("ts_a","b"   ),"returns a new ts as min(ts_a,b)");

        def("time_shift", shyft::api::time_shift,args("timeseries","delta_t"),
            "returns a delta_t time-shifted time-series\n"
            " the values are the same as the original,\n"
            " but the time_axis equals the original + delta_t\n");

        def("create_glacier_melt_ts_m3s", shyft::api::create_glacier_melt_ts_m3s, args("temperature", "sca_m2", "glacier_area_m2", "dtf"),
            doc_intro("create a ts that provide the glacier-melt algorithm based on the inputs")
            doc_parameters()
            doc_parameter("temperature", "TimeSeries", "a temperature time-series, unit [deg.Celcius]")
            doc_parameter("sca_m2", "TimeSeries", "a snow covered area (sca) time-series, unit [m2]")
            doc_parameter("glacier_area_m2", "float", "the glacier area, unit[m2]")
            doc_parameter("dtf","float","degree timestep factor [mm/day/deg.C]; lit. values for Norway: 5.5 - 6.4 in Hock, R. (2003), J. Hydrol., 282, 104-115")
            doc_returns("glacier_melt","TimeSeries","an expression computing the glacier melt based on the inputs")
        );
		/* local scope */ {

			typedef shyft::time_axis::fixed_dt ta_t;
			typedef shyft::time_series::average_accessor<pts_t, ta_t> AverageAccessorTs;
			class_<AverageAccessorTs>("AverageAccessorTs", "Accessor to get out true average for the time-axis intervals for a point time-series", no_init)
				.def(init<const pts_t&, const ta_t&>(args("ts", "ta"), "construct accessor from ts and time-axis ta"))
				.def(init<shared_ptr<pts_t>, const ta_t&>(args("ts", "ta"), "constructor from ref ts and time-axis ta"))
				.def("value", &AverageAccessorTs::value, args("i"), "returns the i'th true average value")
				.def("size", &AverageAccessorTs::size, "returns number of intervals in the time-axis for this accessor")
				;
		}
    }

	static void expose_rating_curve_classes() {

		// overloads for rating_curve_segment::flow
		double (shyft::core::rating_curve_segment::*rcs_flow_1)(double) const = &shyft::core::rating_curve_segment::flow;
		std::vector<double> (shyft::core::rating_curve_segment::*rcs_flow_2)(const std::vector<double> &, std::size_t, std::size_t) const = &shyft::core::rating_curve_segment::flow;

		class_<shyft::core::rating_curve_segment>("RatingCurveSegment",
				doc_intro("Represent a single rating-curve equation.")
				doc_intro("")
				doc_intro("The rating curve function is `a*(h - b)^c` where `a`, `b`, and `c` are parameters")
				doc_intro("for the segment and `h` is the water level to compute flow for. Additionally there")
				doc_intro("is a `lower` parameter for the least water level the segment is valid for. Seen")
				doc_intro("separatly a segment is considered valid for any level greater than `lower`.")
				doc_intro("")
				doc_intro("The function segments are gathered into `RatingCurveFunction`s to represent a")
				doc_intro("set of different rating functions for different levels.")
				doc_see_also("RatingCurveFunction, RatingCurveParameters")
			)
			.def_readonly("lower", &shyft::core::rating_curve_segment::lower,
						  "Least valid water level. Not mutable after constructing a segment.")
			.def_readwrite("a", &shyft::core::rating_curve_segment::a, "Parameter a")
			.def_readwrite("b", &shyft::core::rating_curve_segment::b, "Parameter b")
			.def_readwrite("c", &shyft::core::rating_curve_segment::c, "Parameter c")
			.def(init<double, double, double, double>(args("lower", "a", "b", "c"), "Defines a new RatingCurveSegment with the specified parameters"))
			.def("valid", &shyft::core::rating_curve_segment::valid, (py::args("level")),
					doc_intro("Check if a water level is valid for the curve segment")
					doc_parameter("level", "float", "water level")
					doc_returns("valid", "bool", "True if level is greater or equal to lower")
				)
            //NOTE: For some reason boost 1.65 needs this def *before* the other simpler def, otherwise it fails finding the simple one
            .def("flow", rcs_flow_2, (py::arg("levels"), py::arg("i0") = 0u, py::arg("iN") = std::numeric_limits<std::size_t>::max()),
                doc_intro("Compute the flow for a range of water levels")
                doc_parameters()
                doc_parameter("levels", "DoubleVector", "Vector of water levels")
                doc_parameter("i0", "int", "first index to use from levels, defaults to 0")
                doc_parameter("iN", "int", "first index _not_ to use from levels, defaults to std::size_t maximum.")
                doc_returns("flow", "DoubleVector", "Vector of flow values.")
            )
            .def("flow", rcs_flow_1, (py::arg("level")),
					doc_intro("Compute the flow for the given water level.")
					doc_notes()
					doc_note("There is _no_ check to see if level is valid. It's up to the user to call")
					doc_note("with a correct level.")
					doc_parameters()
					doc_parameter("level", "float", "water level")
					doc_returns("flow", "double", "the flow for the given water level")
				)
			.def("__str__", &shyft::core::rating_curve_segment::operator std::string, "Stringify the segment.")
			;

		// overloads for rating_curve_function::flow
		double (shyft::core::rating_curve_function::*rcf_flow_val)(double) const = &shyft::core::rating_curve_function::flow;
		std::vector<double> (shyft::core::rating_curve_function::*rcf_flow_vec)(const std::vector<double> & ) const = &shyft::core::rating_curve_function::flow;
		// overloads for rating_curve_function::add_segment
		void (shyft::core::rating_curve_function::*rcf_add_args)(double, double, double, double) = &shyft::core::rating_curve_function::add_segment;
		void (shyft::core::rating_curve_function::*rcf_add_obj)(const rating_curve_segment & ) = &shyft::core::rating_curve_function::add_segment;

		class_<shyft::core::rating_curve_function>("RatingCurveFunction",
				doc_intro("Combine multiple RatingCurveSegments into a rating function.")
				doc_intro("")
				doc_intro("RatingCurveFunction aggregates multiple RatingCurveSegments and routes.")
				doc_intro("computation calls to the correct segment based on the water level to compute for.")
				doc_see_also("RatingCurveSegment, RatingCurveParameters")
			)
			.def(init<>("Defines a new empty rating curve function."))
			.def("size", &shyft::core::rating_curve_function::size, "Get the number of RatingCurveSegments composing the function.")
			.def("add_segment", rcf_add_args, (py::arg("lower"), py::arg("a"), py::arg("b"), py::arg("c")),
					doc_intro("Add a new curve segment with the given parameters.")
					doc_see_also("RatingCurveSegment")
				)
			.def("add_segment", rcf_add_obj, py::arg("segment"),
					doc_intro("Add a new curve segment as a copy of an exting.")
					doc_see_also("RatingCurveSegment")
				)
            // ref. note above regarding the order of overloaded member functions
            .def("flow", rcf_flow_vec, py::arg("levels"),
                doc_intro("Compute flow for a range of water levels.")
                doc_parameters()
                doc_parameter("levels", "DoubleVector", "Range of water levels to compute flow for.")
            )
            .def("flow", rcf_flow_val, py::arg("level"),
					doc_intro("Compute flow for the given level.")
					doc_parameters()
					doc_parameter("level", "float", "Water level to compute flow for.")
				)
			.def("__iter__", py::range(&shyft::core::rating_curve_function::cbegin,
									   &shyft::core::rating_curve_function::cend),
				 "Constant iterator. Invalidated on calls to .add_segment")
			.def("__str__", &shyft::core::rating_curve_function::operator std::string, "Stringify the function.")
			;

		// overloads for rating_curve_function::flow
		double (shyft::core::rating_curve_parameters::*rcp_flow_val)(utctime, double) const = &shyft::core::rating_curve_parameters::flow;
		std::vector<double> (shyft::core::rating_curve_parameters::*rcp_flow_ts)(const shyft::api::apoint_ts & ) const = &shyft::core::rating_curve_parameters::flow<shyft::api::apoint_ts>;
		// overloads for rating_curve_function::add_segment
		void (shyft::core::rating_curve_parameters::*rcp_add_obj)(utctime, const rating_curve_function & ) = &shyft::core::rating_curve_parameters::add_curve;

		class_<shyft::core::rating_curve_parameters>("RatingCurveParameters",
				doc_intro("Parameter pack controlling rating level computations.")
				doc_intro("")
				doc_intro("A parameter pack encapsulates multiple RatingCurveFunction's with time-points.")
				doc_intro("When used with a TimeSeries representing level values it maps computations for")
				doc_intro("each level value onto the correct RatingCurveFunction, which again maps onto the")
				doc_intro("correct RatingCurveSegment for the level value.")
				doc_see_also("RatingCurveSegment, RatingCurveFunction, TimeSeries.rating_curve")
			)
			.def(init<>("Defines a empty RatingCurveParameter instance"))
			.def("add_curve", rcp_add_obj, (py::arg("t"), py::arg("curve")),
					doc_intro("Add a curve to the parameter pack.")
					doc_parameters()
					doc_parameter("t", "RatingCurveFunction", "First time-point the curve is valid for.")
					doc_parameter("curve", "RatingCurveFunction", "RatingCurveFunction to add at t.")
				)
			.def("flow", rcp_flow_val, (py::arg("t"), py::arg("level")),
					doc_intro("Compute the flow at a specific time point.")
					doc_parameters()
					doc_parameter("t", "utctime", "Time-point of the level value.")
					doc_parameter("level", "float", "Level value at t.")
					doc_returns("flow", "float", "Flow correcponding to input level at t, `nan` if level is less than the least water level of the first segment or before the time of the first rating curve function.")
				)
			.def("flow", rcp_flow_ts, py::arg("ts"),
					doc_intro("Compute the flow at a specific time point.")
					doc_parameters()
					doc_parameter("ts", "TimeSeries", "Time series of level values.")
					doc_returns("flow", "DoubleVector", "Flow correcponding to the input levels of the time-series, `nan` where the level is less than the least water level of the first segment and for time-points before the first rating curve function.")
				)
			.def("__iter__", py::range(&shyft::core::rating_curve_parameters::cbegin,
									   &shyft::core::rating_curve_parameters::cend),
				 "Constant iterator. Invalidated on calls to .add_curve")
			.def("__str__", &shyft::core::rating_curve_parameters::operator std::string, "Stringify the parameters.")
			;
	}

	static void expose_correlation_functions() {
		const char * kg_doc =
			doc_intro("Computes the kling-gupta KGEs correlation for the two time-series over the specified time_axis")
			doc_parameters()
			doc_parameter("observed_ts","TimeSeries","the observed time-series")
			doc_parameter("model_ts","TimeSeries","the time-series that is the model simulated / calculated ts")
			doc_parameter("time_axis","TimeAxis","the time-axis that is used for the computation")
			doc_parameter("s_r","float","the kling gupta scale r factor(weight the correlation of goal function)")
			doc_parameter("s_a","float","the kling gupta scale a factor(weight the relative average of the goal function)")
			doc_parameter("s_b","float","the kling gupta scale b factor(weight the relative standard deviation of the goal function)")
            doc_returns("KGEs","float","The  KGEs= 1-EDs that have a maximum at 1.0");
		def("kling_gupta", shyft::api::kling_gupta, args("observation_ts", "model_ts", "time_axis", "s_r", "s_a", "s_b"),
			kg_doc
		);

		const char *ns_doc =
            doc_intro("Computes the Nash-Sutcliffe model effiency coefficient (n.s) ")
			doc_intro("for the two time-series over the specified time_axis\n")
			doc_intro("Ref:  http://en.wikipedia.org/wiki/Nash%E2%80%93Sutcliffe_model_efficiency_coefficient \n")
            doc_parameters()
			doc_parameter("observed_ts","TimeSeries","the observed time-series")
			doc_parameter("model_ts","TimeSeries","the time-series that is the model simulated / calculated ts")
			doc_parameter("time_axis","TimeAxis","the time-axis that is used for the computation")
			doc_returns("ns","float","The  n.s performance, that have a maximum at 1.0");

		def("nash_sutcliffe", shyft::api::nash_sutcliffe, args("observation_ts", "model_ts", "time_axis"),
			ns_doc
		);


	}
	static void expose_periodic_ts() {
		const char *docstr =
			doc_intro("Create a Timeseries by repeating the pattern-specification")
			doc_parameters()
			doc_parameter("pattern","DoubleVector","the value-pattern as a sequence of values")
            doc_parameter("dt","int","number of seconds between the pattern values, e.g. deltahours(3)")
            doc_parameter("t0","utctime","specifies the start-time of the pattern")
            doc_parameter("ta","TimeAxis","the time-axis for which the pattern is repeated\n\te.g. your pattern might be 8 3h values,and you could supply\n\ta time-axis 'ta' at hourly resolution")
			;
		def("create_periodic_pattern_ts", shyft::api::create_periodic_pattern_ts, args("pattern","dt","t0","ta"), docstr);

	}
    void timeseries() {
        enum_<time_series::ts_point_fx>("point_interpretation_policy")
            .value("POINT_INSTANT_VALUE",time_series::POINT_INSTANT_VALUE)
            .value("POINT_AVERAGE_VALUE",time_series::POINT_AVERAGE_VALUE)
            .export_values()
            ;
        enum_<time_series::statistics_property>("statistics_property")
            .value("AVERAGE",time_series::statistics_property::AVERAGE)
            .value("MIN_EXTREME",time_series::statistics_property::MIN_EXTREME)
            .value("MAX_EXTREME",time_series::statistics_property::MAX_EXTREME)
            ;

        enum_<time_series::api::extend_ts_fill_policy>(
            "extend_fill_policy",
            "Ref TimeSeries.extend function, this policy determines how to represent values in a gap\n"
            "EPF_NAN : use nan values in the gap\n"
            "EPF_LAST: use the last value before the gap\n"
            "EPF_FILL: use a supplied value in the gap\n"
            )
            .value("FILL_NAN",   time_series::api::extend_ts_fill_policy::EPF_NAN)
            .value("USE_LAST",   time_series::api::extend_ts_fill_policy::EPF_LAST)
            .value("FILL_VALUE", time_series::api::extend_ts_fill_policy::EPF_FILL)
            ;
        enum_<time_series::api::extend_ts_split_policy>(
            "extend_split_policy",
            "Ref TimeSeries.extend function, this policy determines where to split/shift from one ts to the other\n"
            "EPS_LHS_LAST : use nan values in the gap\n"
            "EPS_RHS_FIRST: use the last value before the gap\n"
            "EPS_VALUE    : use a supplied value in the gap\n"
            )
            .value("LHS_LAST",  time_series::api::extend_ts_split_policy::EPS_LHS_LAST)
            .value("RHS_FIRST", time_series::api::extend_ts_split_policy::EPS_RHS_FIRST)
            .value("AT_VALUE",  time_series::api::extend_ts_split_policy::EPS_VALUE)
            ;

        enum_<time_series::convolve_policy>(
            "convolve_policy",
            "Ref Timeseries.convolve_w function, this policy determinte how to handle initial conditions\n"
            "USE_FIRST: value(0) is used for all values before value(0), 'mass preserving'\n"
            "USE_ZERO : fill in zero for all values before value(0):shape preserving\n"
            "USE_NAN  : nan filled in for the first length-1 values of the filter\n"
            )
            .value("USE_FIRST", time_series::convolve_policy::USE_FIRST)
            .value("USE_ZERO", time_series::convolve_policy::USE_ZERO)
            .value("USE_NAN", time_series::convolve_policy::USE_NAN)
            .export_values()
            ;
        class_<time_series::point> ("Point", "A timeseries point specifying utctime t and value v")
            .def(init<utctime,double>(args("t","v")))
            .def_readwrite("t",&time_series::point::t)
            .def_readwrite("v",&time_series::point::v)
            ;
        point_ts<time_axis::fixed_dt>("TsFixed","A time-series with a fixed delta t time-axis, used by the Shyft core,see also TimeSeries for end-user ts");
        point_ts<time_axis::point_dt>("TsPoint","A time-series with a variable delta time-axis, used by the Shyft core,see also TimeSeries for end-user ts");
        TsFactory();
		expose_rating_curve_classes();
        expose_apoint_ts();
		expose_periodic_ts();
		expose_correlation_functions();
		expose_ats_vector();
    }
}
