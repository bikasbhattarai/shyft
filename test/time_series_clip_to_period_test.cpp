#include "test_pch.h"
#define _USE_MATH_DEFINES
#include "core/time_series_dd.h"

using shyft::time_series::dd::apoint_ts;
using shyft::time_series::dd::ats_vector;
using time_axis=shyft::time_series::dd::gta_t;
using shyft::core::utctime;
using shyft::core::calendar;
using shyft::core::deltahours;
using shyft::time_series::ts_point_fx;
using std::vector;

namespace shyft::time_series::dd {
    
    /** @brief clip a concrete time-series to period
     *
     * @note not valid for expression time-series (yet).
     * 
     */
    apoint_ts clip_to_period(apoint_ts const& ts, utcperiod p) {
        //TODO: implement test and clip algorithm here.
        
        return ts;
    }
    
    /** @brief clip all time-series in a tsvector to specified clip_to_period
     *
     */
    ats_vector clip_to_period(const ats_vector& tsv, utcperiod p) {
        //TODO: construct a new ts-vector,
        // then copy, append clipped time-series  into the new ts-vector
        return tsv;
    }
}
TEST_SUITE("time_series") {
    TEST_CASE("ts_clip_to_period") {
        //TODO: verify that clip_to_period(ts,p) works.
    }
    TEST_CASE("tsv_clip_to_period") {
        /// TODO:
        /// Arrange:
        /// construct a tsv with several types of concrete time-time_series (time_axis fixed, calendar, point)
        /// different lengths
        ///
        /// act: call clip_to_period
        /// assert: verify that all are within rules and limits as set by period and time-resolution of the time-series.
        
    }   
}
