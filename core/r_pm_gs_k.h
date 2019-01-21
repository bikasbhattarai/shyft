/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#pragma once

#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "core_serialization.h"

#include "radiation.h"
#include "penman_monteith.h"
#include "kirchner.h"
#include "gamma_snow.h"
#include "actual_evapotranspiration.h"
#include "precipitation_correction.h"
#include "glacier_melt.h"
#include "unit_conversion.h"
#include "routing.h"
#include "mstack_param.h"

namespace shyft {
  namespace core {
    namespace r_pm_gs_k {
        using namespace std;

        /** \brief Simple parameter struct for the RPMGSK method stack
         *
         * This struct contains the parameters to the methods used in the RPMGSK assembly.
         *
         * \tparma RParameter Radiation parameter type that implements the interface:
         *    - RParameter.albedo const --> double, land albedo parameter
         *    - RParameter.turbidity const --> double, turbidity parameter
         * \tparam PMParameter PenmanMonteith parameter type that implements the interface:
         *    - PMParameter.lai const --> double, leaf area index parameter in PenmanMonteith.
         * \tparam GSState GammaSnow parameter type that implements the parameter interface for GammaSnow.
         * \tparam KState Kirchner parameter type that implements the parameter interface for Kirchner.
         * \sa GammaSnowParameter \sa Kirchner \sa RPMGSK  \sa PenmanMonteith
         */
        struct parameter {
            typedef radiation::parameter rad_parameter_t;
            typedef penman_monteith::parameter pm_parameter_t;
            typedef gamma_snow::parameter gs_parameter_t;
            typedef actual_evapotranspiration::parameter ae_parameter_t;
            typedef kirchner::parameter kirchner_parameter_t;
            typedef precipitation_correction::parameter precipitation_correction_parameter_t;
            typedef glacier_melt::parameter glacier_melt_parameter_t;
            typedef routing::uhg_parameter routing_parameter_t;
            typedef mstack_parameter mstack_parameter_t;

            parameter(const rad_parameter_t& rad,
                      const pm_parameter_t& pm,
                      const gs_parameter_t& gs,
                      const ae_parameter_t& ae,
                      const kirchner_parameter_t& k,
                      const precipitation_correction_parameter_t& p_corr,
                      glacier_melt_parameter_t gm=glacier_melt_parameter_t(),
                      routing_parameter_t routing=routing_parameter_t(),
                      mstack_parameter_t msp=mstack_parameter_t()
                     )
             : rad(rad), pm(pm), gs(gs), ae(ae), kirchner(k), p_corr(p_corr) ,gm(gm),routing(routing),msp{msp}{ /*Do nothing */ }
			parameter()=default;
			parameter(const parameter&)=default;
			parameter(parameter&&)=default;
			parameter& operator=(const parameter &c)=default;
			parameter& operator=(parameter&&c)=default;


            rad_parameter_t rad;
            pm_parameter_t pm;
            gs_parameter_t gs;
            ae_parameter_t ae;
            kirchner_parameter_t  kirchner;
            precipitation_correction_parameter_t p_corr;
            glacier_melt_parameter_t gm;
            routing_parameter_t routing;
            mstack_parameter_t msp; ///< method stack parameter(s)
            ///<calibration support, needs vector interface to params, size is the total count
            size_t size() const { return 32; }
            ///<calibration support, need to set values from ordered vector
            void set(const vector<double>& p) {
                if (p.size() != size())
                    throw runtime_error("RPMGSK Parameter Accessor: .set size missmatch");
                int i = 0;
                kirchner.c1 = p[i++];
                kirchner.c2 = p[i++];
                kirchner.c3 = p[i++];
                ae.ae_scale_factor = p[i++];
                gs.tx = p[i++];
                gs.wind_scale = p[i++];
                gs.max_water = p[i++];
                gs.wind_const = p[i++];
                gs.fast_albedo_decay_rate = p[i++];
                gs.slow_albedo_decay_rate = p[i++];
                gs.surface_magnitude = p[i++];
                gs.max_albedo = p[i++];
                gs.min_albedo = p[i++];
                gs.snowfall_reset_depth = p[i++];
                gs.snow_cv = p[i++];
                gs.glacier_albedo = p[i++];
                p_corr.scale_factor = p[i++];
                gs.snow_cv_forest_factor=p[i++];
                gs.snow_cv_altitude_factor=p[i++];
                rad.albedo = p[i++];
                rad.turbidity = p[i++];
                pm.lai = p[i++];
				gs.initial_bare_ground_fraction = p[i++];
				gs.winter_end_day_of_year = size_t(p[i++]);
				gs.calculate_iso_pot_energy = p[i++] != 0.0 ? true : false;
                gm.dtf = p[i++];
                routing.velocity = p[i++];
                routing.alpha = p[i++];
                routing.beta  = p[i++];
                gs.n_winter_days= p[i++];
                gm.direct_response = p[i++];
                msp.reservoir_direct_response_fraction=p[i++];
            }

            ///< calibration support, get the value of i'th parameter
            double get(size_t i) const {
                switch (i) {
                    case  0:return kirchner.c1;
                    case  1:return kirchner.c2;
                    case  2:return kirchner.c3;
                    case  3:return ae.ae_scale_factor;
                    case  4:return gs.tx;
                    case  5:return gs.wind_scale;
                    case  6:return gs.max_water;
                    case  7:return gs.wind_const;
                    case  8:return gs.fast_albedo_decay_rate;
                    case  9:return gs.slow_albedo_decay_rate;
                    case 10:return gs.surface_magnitude;
                    case 11:return gs.max_albedo;
                    case 12:return gs.min_albedo;
                    case 13:return gs.snowfall_reset_depth;
                    case 14:return gs.snow_cv;
                    case 15:return gs.glacier_albedo;
                    case 16:return p_corr.scale_factor;
                    case 17:return gs.snow_cv_forest_factor;
                    case 18:return gs.snow_cv_altitude_factor;
					case 19:return rad.albedo;
					case 20:return rad.turbidity;
					case 21:return pm.lai;
					case 22:return gs.initial_bare_ground_fraction;
					case 23:return (double)gs.winter_end_day_of_year;
					case 24:return gs.calculate_iso_pot_energy ? 1.0 : 0.0;
                    case 25:return gm.dtf;
                    case 26:return routing.velocity;
                    case 27:return routing.alpha;
                    case 28:return routing.beta;
                    case 29:return double(gs.n_winter_days);
                    case 30:return gm.direct_response;
                    case 31:return msp.reservoir_direct_response_fraction;

                default:
                    throw runtime_error("PTGSK Parameter Accessor:.get(i) Out of range.");
                }
                return 0;
            }

            ///< calibration and python support, get the i'th parameter name
            string get_name(size_t i) const {
                static const char *names[] = {
                    "kirchner.c1",
                    "kirchner.c2",
                    "kirchner.c3",
                    "ae.ae_scale_factor",
                    "gs.tx",
                    "gs.wind_scale",
                    "gs.max_water",
                    "gs.wind_const",
                    "gs.fast_albedo_decay_rate",
                    "gs.slow_albedo_decay_rate",
                    "gs.surface_magnitude",
                    "gs.max_albedo",
                    "gs.min_albedo",
                    "gs.snowfall_reset_depth",
                    "gs.snow_cv",
                    "gs.glacier_albedo",
                    "p_corr.scale_factor",
                    "gs.snow_cv_forest_factor",
                    "gs.snow_cv_altitude_factor",
                    "rad.albedo",
                    "rad.turbidity",
                    "pm.lai",
					"gs.initial_bare_ground_fraction",
					"gs.winter_end_day_of_year",
					"gs.calculate_iso_pot_energy",
                    "gm.dtf",
                    "routing.velocity",
                    "routing.alpha",
                    "routing.beta",
                    "gs.n_winter_days",
                    "gm.direct_response",
                    "msp.reservoir_direct_response_fraction"
                };
                if (i >= size())
                    throw runtime_error("RPMGSK Parameter Accessor:.get_name(i) Out of range.");
                return names[i];
            }
        };

        /** \brief Simple state struct for the RPMGSK method stack
         *
         * This struct contains the states of the methods used in the RPMGSK assembly.
         *
         * \tparam GSState GammaSnow state type that implements the state interface for GammaSnow.
         * \tparam KState Kirchner state type that implements the state interface for Kirchner.
         * \sa GammaSnowState \sa Kirchner \sa RPMGSK \sa PenmanMonteith \sa GammaSnow
         */
        struct state {
            typedef gamma_snow::state gs_state_t;
            typedef kirchner::state kirchner_state_t;
            state() {}
            state(const gs_state_t& gs, const kirchner_state_t& k) : gs(gs), kirchner(k) {}
            gs_state_t gs;
            kirchner_state_t kirchner;
            bool operator==(const state& x) const {return gs==x.gs && kirchner==x.kirchner;}

            state scale_snow(const double& snow_storage_area_fraction) const {
                state c{*this};
                c.gs.temp_swe *= snow_storage_area_fraction;
                c.gs.lwc *=snow_storage_area_fraction;// is this reasonable?
                return c;
            }

            /**adjust kirchner q with the  specified scale-factor
            * to support the process of tuning output of a step
            * to a specified observed/wanted average
            */
            void adjust_q(double scale_factor) {kirchner.adjust_q(scale_factor);}
            x_serialize_decl();
        };


        /** \brief Simple response struct for the RPMGSK method stack
         *
         * This struct contains the responses of the methods used in the PRPMGSK assembly.
         */
        struct response {
            // Model responses
            typedef radiation::response  rad_response_t;
            typedef penman_monteith::response pm_response_t;
            typedef gamma_snow::response gs_response_t;
            typedef actual_evapotranspiration::response  ae_response_t;
            typedef kirchner::response kirchner_response_t;
            rad_response_t rad;
            pm_response_t pm;
            gs_response_t gs;
            ae_response_t ae;
            kirchner_response_t kirchner;
            double gm_melt_m3s;
            // Stack response
            double total_discharge;
            double charge_m3s;
            // scale snow parts relative snow_storage_area
            response scale_snow(const double& snow_storage_area_fraction) const {
                response c{*this};
                c.gs.storage *= snow_storage_area_fraction;
                c.gs.outflow *=snow_storage_area_fraction;
                // are there others that we should also scale, sca, is a still meaningful, unscaled ?
                return c;
            }
        };

        /** \brief Calculation Model using assembly of Radiatiom, PenmanMonteith, GammaSnow and Kirchner
         *
         * This model first uses Radiation for calculating net radiation (predicted or translated depending
         * on availability of the radiation data) and soil heat flux,
         * than PenmanMonteith for calculating the reference
         * evapotranspiration based on radiation result and time series data for temperature
         * and relative humidity. Then it uses the GammaSnow method
         * to calculate the snow/ice adjusted runoff using time series data for
         * precipitation and wind speed in addition to the time series used in
         * the PriestleyTaylor method. The PriestleyTaylor potential evaporation is
         * is used to calculate the actual evapotranspiration that is passed on to the
         * last step, Kirchner.
         * Kirchner is run with the gamma snow output
         * and actual evaporation response from the two methods above to
         * calculate the discharge.
         *
         * TODO: This method should construct an internal time stepping rpmgsk struct.
         * This stepper should be used as an iterator in time integration loops. This
         * paves the way for inter-cell communication (e.g. snow transport) without
         * touching this simple interface.
         *
         * \tparam TS Time series type that implements:
         *    - TS::source_accessor_type --> Type of accessor used to retrieve time series data.
         *    - TS.accessor(const T& time_axis) const --> TS::source_accessor_type. Accessor object
         *      for this time series.
         *
         * \tparam TS::source_accessor_type Time series source accessor type that implements:
         *    - TS::source_accessor_type(const TS& time_series, const T& time_axis) --> Construct
         *      accessor for the given time series and time axis.
         *    - TS::source_accessor_type.value(size_t i) --> double, -value of the time series for
         *      period i in the time axis.
         * \tparam T Time axis type that implements:
         *    - T.size() const --> Number of periods in the time axis.
         *    - T(size_t i) const --> shyft::core::utcperiod, - Start and end as shyft::core::utctime
         *      of the i'th period.
         * \tparam S State type that implements:
         *    - S::gs_state_type --> State type for the GammaSnow method.
         *    - S::kirchner_state_type --> State type for the Kirchner method.
         *    - S.gs --> S::gs_state_type, - State variables for the GammaSnow method
         *    - S.kirchner --> S::kirchner_state_type, - State variables for the Kirchner method
         * \tparam R Response type that implements:
         *    - R::gs_response_type --> Response type for the GammaSnow routine.
         *    - R.gs --> R::gs_response_type, -Response object passed to the GammaSnow routine.
         * \tparam P Parameter type that implements:
         *    - P::rad_parameter_type --> Parameter type for the Radiation method
         *    - P::pm_parameter_type --> Parameter type for the PenmanMonteith method
         *    - P::gs_parameter_type --> Parameter type for the GammaSnow method.
         *    - P::ae_parameter_type --> Parameter type for the ActualEvapotranspiration method.
         *    - P::kirchner_parameter_type --> Parameter type for the Kirchner method.
         *    - P.rad --> P::rad_parameter_type --> Parameters for Radiation method.
         *    - P.pm --> P::pm_parameter_type --> Parameters for the PenmanMonteith method.
         *    - P.gs --> P::gs_parameter_type --> Parameters for thge GammaSnow method.
         *    - P.ae --> P::ae_parameter_type --> Parameters for thge ActualEvapotranspiration method.
         *    - P.kirchner --> P::kirchner_parameter_type --> Parameters for the Kirchner method.
         * \tparam SC State collector type that implements:
         *    - SC.collect(utctime t, const S& state) --> Possibly save some states at time t.
         * \tparam RC Response collector type that implements:
         *    - RC.collect(utctime t, const R& response) --> Possibly save some responses at time t.
         */
        template<template <typename, typename> class A, class R, class T_TS, class P_TS, class WS_TS, class RH_TS, class RAD_TS, class T,
        class S, class GCD, class P, class SC, class RC>
        void run_r_pm_gs_k(const GCD& geo_cell_data,
            const P& parameter,
            const T& time_axis, int start_step, int  n_steps,
            const T_TS& temp,
            const P_TS& prec,
            const WS_TS& wind_speed,
            const RH_TS& rel_hum,
            const RAD_TS& rad,
            S& state,
            SC& state_collector,
            RC& response_collector
            ) {
            // Access time series input data through accessors of template A (typically a direct accessor).
            using temp_accessor_t = A<T_TS, T>;
            using prec_accessor_t = A<P_TS, T>;
            using wind_speed_accessor_t = A<WS_TS, T>;
            using rel_hum_accessor_t = A<RH_TS, T>;
            using rad_accessor_t = A<RAD_TS, T>; //TODO: think about how to evaluate the availability of radiation data

            auto temp_accessor = temp_accessor_t(temp, time_axis);
            auto prec_accessor = prec_accessor_t(prec, time_axis);
            auto wind_speed_accessor = wind_speed_accessor_t(wind_speed, time_axis);
            auto rel_hum_accessor = rel_hum_accessor_t(rel_hum, time_axis);
            auto rad_accessor = rad_accessor_t(rad, time_axis);

            // Initialize the method stack
            precipitation_correction::calculator p_corr(parameter.p_corr.scale_factor);
            R response;
            radiation::calculator < typename P::rad_parameter_t, typename R::rad_response_t> radc(parameter.rad);
            penman_monteith::calculator< typename P::pm_parameter_t, typename R::pm_response_t> pm(parameter.pm);
            gamma_snow::calculator<typename P::gs_parameter_t, typename S::gs_state_t, typename R::gs_response_t> gs;
            kirchner::calculator<kirchner::trapezoidal_average, typename P::kirchner_parameter_t> kirchner(parameter.kirchner);
            //
            // Get the initial states

            const double forest_fraction=geo_cell_data.land_type_fractions_info().forest();
            const double glacier_fraction = geo_cell_data.land_type_fractions_info().glacier();
            const double gm_direct = parameter.gm.direct_response; //glacier melt directly out of cell
            const double gm_routed = 1-gm_direct; // glacier melt routed through kirchner
            const double snow_storage_fraction = geo_cell_data.land_type_fractions_info().snow_storage();// on this part, snow builds up, and melts.-> season time-response
            const double kirchner_routed_prec =  geo_cell_data.land_type_fractions_info().reservoir()*(1.0-parameter.msp.reservoir_direct_response_fraction) + geo_cell_data.land_type_fractions_info().lake();
            const double direct_response_fraction = glacier_fraction*gm_direct + geo_cell_data.land_type_fractions_info().reservoir()*parameter.msp.reservoir_direct_response_fraction;// only direct response on reservoirs
            const double kirchner_fraction = 1 - direct_response_fraction;
            const double cell_area_m2 = geo_cell_data.area();
            const double glacier_area_m2 = geo_cell_data.area()*glacier_fraction;
            const double altitude= geo_cell_data.mid_point().z;
            // Step through times in axis
            size_t i_begin = n_steps > 0 ? start_step : 0;
            size_t i_end = n_steps > 0 ? start_step + n_steps : time_axis.size();
            std::cout<<i_end<<std::endl;
            for (size_t i = i_begin ; i < i_end ; ++i) {
                utcperiod period = time_axis.period(i);
                utctime t = time_axis.time(i);
                double temp = temp_accessor.value(i);
                double rad = rad_accessor.value(i);
                double rel_hum = rel_hum_accessor.value(i);
                double prec = p_corr.calc(prec_accessor.value(i));
                double ws = wind_speed_accessor.value(i);
                std::cout<<"temp: " <<temp<<std::endl;
                std::cout<<"rh: "<<rel_hum<<std::endl;
                std::cout<<"rad: "<<rad<<std::endl;
                state_collector.collect(i, state.scale_snow(snow_storage_fraction));///< \note collect the state at the beginning of each period (the end state is saved anyway)

                gs.step(state.gs, response.gs, period.start, period.timespan(), parameter.gs,
                        temp, rad, prec, wind_speed_accessor.value(i), rel_hum,forest_fraction,altitude);
                response.gm_melt_m3s = glacier_melt::step(parameter.gm.dtf, temp, cell_area_m2*response.gs.sca, glacier_area_m2);
//                radc.net_radiation(response.rad, geo_cell_data.mid_point().x, t, 0.0, 0.0, temp, rel_hum, altitude);// radiation model
//                radc.net_radiation(response.rad, 40.4, t, 0.0, 0.0, temp, rel_hum, altitude, rad/0.0864);// radiation model TODO: check translation, it doesn't work now
                radc.net_radiation(response.rad, 40.4, t, 0.0, 0.0, temp, rel_hum, altitude);// radiation model
                std::cout<<"sw radiation:" <<response.rad.sw_radiation<<std::endl;
                std::cout<<"lw radiation:" <<response.rad.lw_radiation<<std::endl;
                std::cout<<"net radiation:" <<response.rad.net_radiation<<std::endl;
                pm.reference_evapotranspiration_asce(response.pm, response.rad.net_radiation*0.0036,temp,rel_hum,altitude,ws); // TODO: move height of measurements into parameters
                std::cout<<"et_ref: "<<response.pm.et_ref<<std::endl;
                std::cout<<"----------------------------------------"<<std::endl;
                //response.pt.pot_evapotranspiration = pt.potential_evapotranspiration(temp, rad, rel_hum)*to_seconds(calendar::HOUR); //mm/s -> mm/h
                response.ae.ae = actual_evapotranspiration::calculate_step(
                                  state.kirchner.q,
                                  response.pm.et_ref,
                                  parameter.ae.ae_scale_factor,
                                  std::max(response.gs.sca,glacier_fraction), // a evap only on non-snow/non-glac area
                                  period.timespan()
                                );
                double gm_mmh= shyft::m3s_to_mmh(response.gm_melt_m3s, cell_area_m2);
                kirchner.step(period.start, period.end, state.kirchner.q, response.kirchner.q_avg, response.gs.outflow*snow_storage_fraction + prec*kirchner_routed_prec + gm_routed*gm_mmh, response.ae.ae); // all units mm/h over 'same' area

                response.total_discharge =
                      std::max(0.0,prec - response.ae.ae)*direct_response_fraction // when it rains, remove ae. from direct response
                    + gm_direct*gm_mmh  // glacier melt direct response
                    + response.kirchner.q_avg*kirchner_fraction;
                response.charge_m3s =
                    + shyft::mmh_to_m3s(prec, cell_area_m2)
                    - shyft::mmh_to_m3s(response.ae.ae, cell_area_m2)
                    + response.gm_melt_m3s
                    - shyft::mmh_to_m3s(response.total_discharge, cell_area_m2);
                // Possibly save the calculated values using the collector callbacks.
                response_collector.collect(i, response.scale_snow(snow_storage_fraction));///< \note collect the response valid for the i'th period (current state is now at the end of period)
                if(i+1==i_end)
                    state_collector.collect(i+1, state.scale_snow(snow_storage_fraction));///< \note last iteration,collect the  final state as well.
            }
            response_collector.set_end_response(response.scale_snow(snow_storage_fraction));
        }
    } // pt_gs_k
  } // core
} // shyft
  //-- serialization support shyft
x_serialize_export_key(shyft::core::r_pm_gs_k::state);
