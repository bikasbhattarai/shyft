/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <future>

#include "utctime_utilities.h"
#include "time_series.h"
#include "time_series_average.h"

namespace shyft {
    namespace qm {
        using namespace std;
        using shyft::core::to_seconds;

        /**\brief quantile_index generates a pr. time-step index for order by value
        * \tparam tsa_t time-series accessor type that fits to tsv_t::value_type and ta_t, thread-safe fast access to the value aspect
        *               of the time-series for each period in the specified time-axis
        *               requirement to the type is:
        *               constructor tsa_t(ts,ta) including just enough stuff for trivial copy-ct, move-ct etc.
        *                tsa_t.value(size_t i) -> double gives the i'th time-axis period value of ts (cache/lookup if needed)
        *
        * \tparam tsv_t time-series vector type
        * \tparam ta_t time-axis type, require .size() (the tsa_t hides other methods required)
        * \param tsv time-series vector that keeps .size() time-series
        * \param ta time-axis
        * \return vector<vector<int>> that have dimension [ta.size()][tsv.size()] holding the indexes of the by-value in time-step ordered list
        *
        */
        template <class tsa_t, class tsv_t, class ta_t>
        vector<vector<int>> quantile_index(tsv_t const &tsv, ta_t const &ta) {
            vector < tsa_t> tsa; tsa.reserve(tsv.size()); // create ts-accessors to enable fast thread-safe access to ts value aspect
            for (const auto& ts : tsv) tsa.emplace_back(ts, ta, shyft::time_series::extension_policy::USE_NAN);

            vector<vector<int>> qi(ta.size()); // result vector, qi[i] -> a vector of index-order tsv[i].value(t)
            //vector<int> pi(tsv.size()); for (size_t i = 0; i<tsv.size(); ++i) pi[i] = i;//initial order 0..n-1, and we re-use it in the loop
            for (size_t t = 0; t<ta.size(); ++t) {
				vector<int> pi;pi.reserve(pi.size()); // get rid of nan-members before sort
				for (size_t i = 0; i<tsa.size(); ++i) {
					if (isfinite(tsa[i].value(t))) {
						pi.emplace_back(i);
					}
				}
                sort(begin(pi), end(pi), [&tsa, t](int a, int b)->int {return tsa[a].value(t)<tsa[b].value(t); });
                // Filter away the nans, which signal that we have time series that don't extend this far
                qi[t] = pi; // it's a copy assignment, ok, since sort above will start out with previous time-step order
            }
            return qi;
        }


        /** \brief compute quantile-values based on weigh_value ordered (wvo) items
        *
        *
        * \note The requirement to the wvo_accessor
        *  (1) sum weights == 1.0
        *  (2) weight(i), value(i) are already ordered by value asc
        *   are the callers responsibility
        *   no checks/asserts are done
        *
        *  CONSIDER: we might want to rework the output of this algorithm
        *        so that it streams/put data directly into the
        *        target-time-series (instead of building a vector
        *        and then spread that across time-series)
        *
        *
        * \tparam WVO Weight Value Ordered accessor, see also wvo_accessor
        *  .size() -> number of weight,value items
        *  .value(i)-> double the i'th value
        *  .weight(i)-> double the i'th weight
        * \param n_q number of quantiles,>1, evenly distributed over 0..1.0
        * \param items the weighted items, weight,value, ordered asc value by i'th index,
        *        assuming that sum weights = 1.0
        * \return vector<double> of size n_q, where the i'th value is the i'th computed quantile
        */
        template <class WVO> // to-consider: maybe just iterators to weight-value pairs ?
        vector<double> compute_weighted_quantiles(size_t n_q, WVO const& wv) {
            const double q_step = 1.0 / (n_q - 1); // iteration is 0-based..
            vector<double> q; q.reserve(n_q);
            // tip: 'j' denotes the 'source' quantile material, 'i' denotes the wanted quantiles
            size_t j = 0;// the j'th item in wv
            double q_j = 0.0;// since the 'weights' are normalized, each weight represent the length of quantile-segment the value is valid for
            double w_j = wv.weight(j);// We try to have only one access for each indexed element
            double v_j = wv.value(j); // (we *could* trust the compiler optimizer here)
            for (size_t i = 0; i<n_q; ++i) { // ensure we are precise on number of quantiles we provide
                double q_i = q_step*i;  // multiplication or addition, almost same performance(?)
                while (q_j + w_j  < q_i) { //  as long as current q_j+w_j are LESS than target q_i
                    q_j += w_j;            // climb along the 'j' quantile axis until we reach wanted level
                    if (j + 1 < wv.size()) {  // only get next element if there are more left
                        w_j = wv.weight(++j);
                        v_j = wv.value(j);
                    }
                }
                q.emplace_back(v_j);// note: we *could* do weighted interpolation q_j < q_i <q_j+1, but we do simple rank here
            }
            return q;
        }


        /** \brief compute interpolated quantile-values based on weigh_value ordered (wvo) items
        *
        *
        * \note The requirement to the wvo_accessor
        *  (1) sum weights == 1.0
        *  (2) weight(i), value(i) are already ordered by value asc
        *   are the callers responsibility
        *   no checks/asserts are done
        *
        *  CONSIDER: we might want to rework the output of this algorithm
        *        so that it streams/put data directly into the
        *        target-time-series (instead of building a vector
        *        and then spread that across time-series)
        *
        *  The difference between this function and compute_weighted_quantiles
        *  is that this function assigns a value to each desired quantile that
        *  is the interpolated value between the two closest input quantiles,
        *  while compute_weighted_quantiles assigns a value to each desired
        *  quantile that is just the value of the input quantile that is less
        *  than or equal to the desired quantile. (This function treats the edge
        *  cases similarly - i.e. whenever the desired quantile is lower than
        *  half of the first input quantile or greater than the last half of
        *  the last input quantile, the value assigned is the first and last
        *  input quantile values, respectively, and no interpolation is
        *  performed).
        *
        * \tparam WVO Weight Value Ordered accessor, see also wvo_accessor
        *  .size() -> number of weight,value items
        *  .value(i)-> double the i'th value
        *  .weight(i)-> double the i'th weight
        * \param n_q number of quantiles,>1, evenly distributed over 0..1.0
        * \param items the weighted items, weight,value, ordered asc value by i'th index,
        *        assuming that sum weights = 1.0
        * \return vector<double> of size n_q, where the i'th value is the i'th computed quantile
        */
        template <class WVO>
        vector<double> compute_interp_weighted_quantiles(size_t n_q, WVO const& wv) {
            const double q_step = 1.0 / (n_q - 1); // iteration is 0-based..
            vector<double> q; q.reserve(n_q);
            size_t j = -1;
            double q_j = 0.0;
            double w_j = 0;
            double width = 0;
            for (size_t i = 0; i< n_q; ++i) {
                double q_i = q_step * i;
                while (q_j <= q_i && j + 1 != wv.size()) {
                    width = 0.5 * w_j;
                    w_j = wv.weight(++j);
                    width += 0.5 * w_j;
                    q_j += width;
                }
                // Calculate interpolation
                double v_j;

                if (j == 0 || (j + 1 == wv.size() && q_i >= q_j)) {
                    // Flat values beyond interpolation range
                    v_j = wv.value(j);
                } else {
                    double interp_start = wv.value(j-1);
                    double interp_end = wv.value(j);
                    double inclination = (interp_end - interp_start) / width;
                    double offset = q_j - width;
                    v_j = interp_start + inclination * (q_i - offset);
                }
                q.emplace_back(v_j);
            }
            return q;
        }

        /**\brief a Weight Value Ordered collection for ts-vector type
        *
        * This type is a stack-context only, light weight wrapper,
        *  for ts-vector type, to be feed into the
        * computed_weighted_quantiles algorithm.
        *
        * The idea here is to give a class that is useful for
        * anything that resembles typical shyft time-series or time-series accessor
        * Note that we hide the mechanism for how 'things' are sorted
        * keeping the indirect-indexing array an internal detail away from the
        * compute_weighted_quantile function.
        *
        * \see computed_weighted_quantiles
        */
        template <class tsa_t>
        struct wvo_accessor {
            vector<vector<int>>const& ordered_ix; ///<index ordering of the i'th timestep
            vector<double>const & w;///< weights of each ts - assumed to be same for all timesteps
            vector<tsa_t>const& tsv;///< access to .value(t), where t is size_t
            size_t t_ix = 0;///< current time step index
            vector<double> w_sum;///< sum of weights per timestep so we can return normalized .weight(i). Note that it must be per timestep even if we have the same weights for all timestep - this is because ordered_idx can contain varying length-vectors of indices, since certain time series might be invalid at certain times.
                                 //-------------------
            wvo_accessor(vector<vector<int>> const&ordered_ix, vector<double> const&w, vector<tsa_t> const &tsv)
                :ordered_ix(ordered_ix), w(w), tsv(tsv)
            {
                for (size_t t = 0; t<ordered_ix.size(); ++t) {
                    double sum = 0.0;
                    for (size_t i = 0; i<ordered_ix[t].size(); ++i) {
                        sum += w[ordered_ix[t][i]];
                    }
                    w_sum.emplace_back(sum);
                }
            }

            //-- here goes the template contract requirements
            size_t size() const { return ordered_ix[t_ix].size(); } //<< number of forecasts available at time-point t_ix
            double weight(size_t i) const { return w[ordered_ix[t_ix][i]] / w_sum[t_ix]; }
            double value(size_t i) const { return tsv[ordered_ix[t_ix][i]].value(t_ix); }
        };

        /** compute the utcperiod that covers all ts in tsv
         *
         * Given a tsv, find the largest period that covers all time-axis.
         * \tparam tsv_t any kind of iterable that returns object that supports
         *               .time_axis(), that provides .total_period().
         * \param tsv - the ts-vector with forecasts
         * \return utcperiod, the period that covers all time_axis().total_period()
         */
        template<class tsv_t>
        inline core::utcperiod compute_cover_period(const tsv_t&tsv) {
            core::utcperiod fc_period;
            for (const auto&ts : tsv) {
                auto p = ts.time_axis().total_period();
                if (!fc_period.valid())
                    fc_period = p;
                else {
                    fc_period.start = std::min(fc_period.start, p.start);
                    fc_period.end = std::max(fc_period.end, p.end);
                }
            }
            return fc_period;
        }

        /** given interpolation start,end, and max time t_end, return empty period, or ready clipped period */
        inline core::utcperiod compute_interpolation_period(core::utcperiod fc_period,core::utctime interpolation_start,core::utctime interpolation_end,core::utctime t_end) {
            core::utcperiod interpolation_period;
            if (core::is_valid(interpolation_start)) {
                interpolation_period = core::utcperiod(interpolation_start, t_end);
                if (core::is_valid(interpolation_end))
                    interpolation_period.end = interpolation_end;
                // now clip i.start <= fc_end
                // and clip i.end <= fc_end
                interpolation_period.start = std::min(fc_period.end, interpolation_period.start);
                interpolation_period.end = std::min(fc_period.end, interpolation_period.end);
            }
            return interpolation_period;
        }

        /**\brief The main quantile mapping function, which, using quantile
        * calculations, maps the values of the weighted 'forecast' time
        * series vectors onto the 'prior' time series. This mapping is done
        * for each time step in the specified time axis.
        * \tparam tsa_t The time-series accessor type.
        * \tparam tsv_t The time-series vector type.
        * \tparam ta_t The time-axis type.
        * \param pri_tsv The time-series vector containing the 'prior'
        *      observations. The values in this vector will be more or less
        *      overwritten by the values in fc_tsv that, using the quantile
        *      mapping, correspond to the values in pri_tsv.
        * \param fc_tsv The time-series vector containing the 'forecast'.
        *      The values in this vector will overwrite the corresponding ones
        *      (in the quantile mapping sense) in pri_tsv.
        * \param pri_idx_v The vector containing the indices that will order
        *      the observations in pri_tsv according to their value, so that
        *      if the current time axis index is t, then pri_idx_v[t] will
        *      give a vector sorting pri_tsv for that timestep.
        * \param fc_idx_v Same as pri_idx_v, but applies to fc_tsv.
        * \param fc_weights Contains the weights to apply to each
        *      value in fc_tsv. Assumed to be the same for all timesteps.
        * \param time_axis The time axis over which to perform the mapping.
        * \param interpolation_start The start of the period over which we
        *      want to interpolate between the forecast and prior. We assume
        *      that the end of the interpolation period coincides with the
        *      last point of time_axis. For all times before
        *      interpolation_start, the corresponding value in pri_tsv will
        *      only contain the corresponding quantile-mapped value in
        *      fc_tsv, while for times after interpolation_start, the value
        *      in pri_tsv will contain a linearly interpolated mix of the
        *      value already in pri_tsv and the corresponding quantile-mapped
        *      value in fc_tsv. If this parameter is set to
        *      core::no_utctime, interpolation will not take place.
        * \param interpolation_end As interpolation_start, but denotes the end
        *      of the interpolation period.
        * \param interpolated_quantiles Whether to interpolate between forecast
        *      quantile values or use the quantile values that are less than or
        *      equal to the current quantile (default) when mapping the forecast
        *      observations to the priors. Note that this interpolation is
        *      something else than what is referred to wrt the
        *      interpolation_start and interpolation_end.
        * \return tsv_t Containing the quantile mapped values at the times
        *      indicated by time_axis, and interpolated after
        *      interpolation_start.
        **/
        template <class tsa_t,class rtsv_t, class tsv_t, class ta_t>
        rtsv_t quantile_mapping(tsv_t const &pri_tsv, tsv_t const &fc_tsv,
                vector<vector<int>> const &pri_idx_v,
                vector<vector<int>> const &fc_idx_v,
                vector<double> const &fc_weights,
                ta_t const &time_axis,
                core::utctime const &interpolation_start,
                core::utctime const interpolation_end = core::no_utctime,
                bool interpolated_quantiles = false) {

            vector<tsa_t> pri_accessor_vec;
            pri_accessor_vec.reserve(pri_tsv.size());
            for (const auto& ts : pri_tsv)
                pri_accessor_vec.emplace_back(ts, time_axis);

            auto fc_period=compute_cover_period(fc_tsv); // compute the maximum forecast period
            auto interpolation_period=compute_interpolation_period(fc_period,interpolation_start,interpolation_end,time_axis.time(time_axis.size() - 1));

            vector<tsa_t> fc_accessor_vec;
            fc_accessor_vec.reserve(fc_tsv.size());
            for (const auto& ts : fc_tsv)
                fc_accessor_vec.emplace_back(ts, time_axis);

            auto wvo_fc = wvo_accessor<tsa_t>(fc_idx_v, fc_weights, fc_accessor_vec);

            rtsv_t output;
            output.reserve(pri_tsv.size());
            for (size_t i = 0; i<pri_tsv.size(); ++i) {
                output.emplace_back(time_axis, nan, time_series::POINT_AVERAGE_VALUE);
            }
            for (size_t t = 0; t<time_axis.size(); ++t) {
                wvo_fc.t_ix = t;
                size_t num_pri_cases = pri_idx_v[t].size();
                if (wvo_fc.size() > 0 && (!core::is_valid(interpolation_period.end) || time_axis.time(t)<interpolation_period.end)) {
                    vector<double> quantile_vals;
                    if (interpolated_quantiles) {
                        quantile_vals = compute_interp_weighted_quantiles(num_pri_cases, wvo_fc);
                    } else {
                        quantile_vals = compute_weighted_quantiles(num_pri_cases, wvo_fc);
                    }
                    if ( (interpolation_period.contains(time_axis.time(t)) ||
                            interpolation_period.end == time_axis.time(t))) {
                        core::utctime start = interpolation_period.start;
                        core::utctime end = interpolation_period.end;
                        double interp_weight = to_seconds(time_axis.time(t) - start)/to_seconds(end - start);
                        for (size_t i = 0; i < num_pri_cases; ++i)
                            output[pri_idx_v[t][i]].set(t,(1.0 - interp_weight)*quantile_vals[i] + interp_weight*pri_accessor_vec[pri_idx_v[t][i]].value(t));
                    } else {
                        for (size_t i = 0; i < num_pri_cases; ++i)  output[pri_idx_v[t][i]].set(t, quantile_vals[i]);
                    }
                } else { // if no more forecast available, or after valid end, use the prior scenario value for the specified time-points
                    for (size_t i = 0; i < output.size(); ++i)  output[i].set(t, pri_accessor_vec[i].value(t));
                }
            }

            return output;
        }

        /** append historical_ts to the end qm ts,
         * inplace construct the resulting ts into the result vector.
         * 
         * \param historical_ts with values that should be appended to the qm, could be of any time-resolution.
         * \param qm the quantile mapped result, starts the result series, the time-axis should be aligned with the ta-parameter first part
         * \param ta the total time-axis for the resulting time-series, the qm ts should be aligned with the start of this one
         * \param ta_hist the time-axis  that should be appended
         * \param r the result vector where the result should be .emplace_back constructed into
         */
         template <class ta_t,class ts_t,class rts_t=ts_t>
         void merge_qm_result( const ts_t &historical_ts, const rts_t& qm,const ta_t &ta, const ta_t&ta_hist,vector<ts_t>&r) {
                vector<double> v(ta.size(),nan);
                auto qm_v{qm.values()};
                auto b=begin(qm_v);
				const auto p_fx = historical_ts.point_interpretation();
                std::copy(b,b+qm.size(),begin(v));
				vector<double> hv = p_fx==time_series::POINT_AVERAGE_VALUE?
					time_series::accumulate_stair_case(ta_hist, historical_ts, true):
					time_series::accumulate_linear(ta_hist, historical_ts, true);
				std::copy(begin(hv), end(hv), begin(v) + qm.size());
                r.emplace_back(ta, move(v), time_series::POINT_AVERAGE_VALUE); // because the operation is true average
            };

        /** \brief the quantile_map_forecast applies quantile_mapping to weighted forecast_set and historical data
        *
        *
        * \see quantile_mapping function
        *
        * \tparam tsa_t time-series accessor type to extract average data from the ts-vectors elements
        * \tparam tsv_t time-series vector type, ala std::vector<ts_t>
        * \tparam ta_t time-axis type, usually a fixed time-step time-axis
        *
        * \param forecast_set a set of forecasts with n-elements
        * \param set_weights the weights for each of the n-elements in the forecast_set
        * \param historical_data historical time-series that should cover the forecast period
        * \param time_axis the time-axis that we would like the resulting time-series to be mapped into
        * \param interpolation_start the time within the forecast period where the interpolation period should start
        * \param interpolation_end   the time within the forecast period where the interpolation period should end, default last fc
        * \param interpolated_quantiles Whether to interpolate between forecast
        *      quantile values or use the quantile values that are less than or
        *      equal to the current quantile (default) when mapping the forecast
        *      observations to the priors. Note that this interpolation is
        *      something else than what is referred to wrt the
        *      interpolation_start and interpolation_end.

        * \return a time-series vector with the resulting quantile-mapped time-series
        */
        template <class tsa_t, class tsv_t, class ta_t>
        tsv_t quantile_map_forecast(vector<tsv_t> const &forecast_sets,
            vector<double> const &set_weights, tsv_t const &historical_data,
            ta_t const &time_axis,
            core::utctime interpolation_start,core::utctime interpolation_end=core::no_utctime,
            bool interpolated_quantiles=false) {
            tsv_t forecasts_unpacked;
            vector<double> weights_unpacked;
            for (size_t i = 0; i<forecast_sets.size(); ++i) {
                forecasts_unpacked.reserve(forecasts_unpacked.size() + forecast_sets[i].size());
                weights_unpacked.reserve(weights_unpacked.size() + forecast_sets[i].size());
                for (size_t j = 0; j<forecast_sets[i].size(); ++j) {
                    forecasts_unpacked.emplace_back(forecast_sets[i][j]);
                    weights_unpacked.emplace_back(set_weights[i]);
                }
            }
            auto fc_period=compute_cover_period(forecasts_unpacked); // compute the maximum forecast period
            auto ta_last_step = time_axis.time(time_axis.size() - 1);
            auto interpolation_period=compute_interpolation_period(fc_period,interpolation_start,interpolation_end,ta_last_step);
            auto qm_end = interpolation_period.valid()?interpolation_period.end:fc_period.end;
            auto qm_ix_end = time_axis.index_of(qm_end);
            auto qm_n_steps = qm_ix_end==std::string::npos?time_axis.size():qm_ix_end+1;
            auto qm_time_axis = time_axis.slice(0,qm_n_steps); // +1, ->include time-step we end into
            auto historical_indices_handle = async(launch::async, quantile_index<tsa_t, tsv_t, ta_t>, historical_data, qm_time_axis);
            auto forecast_indices = quantile_index<tsa_t>(forecasts_unpacked, qm_time_axis);
            using rtsv_t=vector<time_series::point_ts<ta_t>>;
            auto qm_tsv=quantile_mapping<tsa_t,rtsv_t>(
                            historical_data,
                            forecasts_unpacked,
                            historical_indices_handle.get(),
                            forecast_indices,
                            weights_unpacked,
                            qm_time_axis,
                            interpolation_start,
                            interpolation_end,
                            interpolated_quantiles
                    );
            if(qm_time_axis==time_axis) {// early exit, zero copy if possible
                tsv_t r;r.reserve(qm_tsv.size());
                for(auto &qm:qm_tsv) {
                    r.emplace_back(qm.ta,move(qm.v),time_series::POINT_AVERAGE_VALUE);
                }
                return r;
            }

            tsv_t r;
            r.reserve(historical_data.size());
			auto ta_ext = time_axis.slice(qm_n_steps, time_axis.size() - qm_n_steps);// the time-axis extension with historical data
            for (size_t i = 0; i<historical_data.size(); ++i)  merge_qm_result(historical_data[i],qm_tsv[i],time_axis,ta_ext,r);
            return r;
        }
    }
}
