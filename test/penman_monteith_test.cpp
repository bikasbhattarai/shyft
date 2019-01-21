#include "test_pch.h"
#include "core/utctime_utilities.h"
#include <vector>
#include <chrono>
#include <boost/math/constants/constants.hpp>
#include "core/radiation.h"
#include <cmath>
#include <random>
#include <tuple>
#include "core/hydro_functions.h"
#include "core/penman_monteith.h"


namespace shyft::test {

    class trapezoidal_average {
    private:
        double area = 0.0;
        double f_a = 0.0;; // Left hand side of next integration subinterval
        double t_start = 0.0; // Start of integration period
        double t_a = 0.0; // Left hand side time of next integration subinterval
    public:
        explicit trapezoidal_average() {}

        /** \brief initialize must be called to reset states before being used during ode integration.
         */
        void initialize(double f0, double t_start) {
            this->f_a = f0;
            this->t_start = t_start;
            t_a = t_start;
            area = 0.0;
        }

        /** \brief Add contribution to average using a simple trapezoidal rule
         *
         * See: http://en.wikipedia.org/wiki/Numerical_integration
         */
        void add(double f, double t) {
            area += 0.5*(f_a + f)*(t - t_a);
            f_a = f;
            t_a = t;
        }

        double result() const { return area/(t_a - t_start); }
    };

}


TEST_SUITE("penman_monteith") {
    using namespace shyft::core;
//    using shyft::core::radiation::surface_normal;
    using shyft::core::calendar;
    using shyft::core::utctime;
    using namespace shyft::test;
    // test basics: creation, etc



    TEST_CASE("hourly_full"){

        //===========================//
        // getting radiation  //
        radiation::parameter rad_p;
        radiation::response rad_r;
        rad_p.albedo = 0.2;
        rad_p.turbidity = 1.0;
        radiation::calculator<radiation::parameter,radiation::response> rad(rad_p);
        calendar utc_cal;
        // Greeley, Colorado weather station
        double lat = 40.41;
        double elevation = 1462.4;
        double ht = 1.68;
        double hws = 3.0;
        utctime t;
        // checking for horizontal surface Eugene, OR, p.64, fig.1b
        arma::vec surface_normal({0.0,0.0,1.0});
        double slope = 0.0;
        double aspect = 0.0;
        utctime ta;

        ta = utc_cal.time(2000, 06, 2, 00, 00, 0, 0);
        //rad.psw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);




        //evapotranspiraiton PM
        penman_monteith::parameter pm_p(2.0,3.0,1.68);
        penman_monteith::response pm_r;
        penman_monteith::calculator<penman_monteith::parameter,penman_monteith::response> pm_calculator(pm_p);
        trapezoidal_average av_et_ref;
        av_et_ref.initialize(pm_r.et_ref, 0.0);
        // ref.: ASCE=EWRI Appendix C: Example Calculation of ET
        double temperature[23] = {16.5, 15.4, 15.5, 13.5, 13.2, 16.2, 20.0, 22.9, 26.4, 28.2, 29.8, 30.9, 31.8, 32.5, 32.9, 32.4, 30.2, 30.6, 28.3, 25.9, 23.9, 20.1, 19.9};
        double vap_pressure[23] = {1.26, 1.34, 1.31, 1.26, 1.24, 1.31, 1.36, 1.39, 1.25, 1.17, 1.03, 1.02, 0.98, 0.87, 0.86, 0.93, 1.14, 1.27, 1.27, 1.17, 1.20, 1.10, 1.05};
        double windspeed[23] = {0.5, 1.0, 0.68, 0.69, 0.29, 1.24, 1.28, 0.88, 0.72, 1.52, 1.97, 2.07, 2.76, 2.990, 3.10, 2.77, 3.41, 2.78, 2.95, 3.27, 2.86, 2.7, 2.0};
        double svp[23];
        double rhumidity[23];
        for (int i=0;i<23;i++){
            svp[i] = hydro_functions::svp(temperature[i]);
            //std::cout<<"svp: "<<svp[i]<<std::endl; // matched
            rhumidity[i] = svp[i]*100/vap_pressure[i];
        }

        rad.net_radiation(rad_r, lat, ta, slope, aspect, temperature[0], rhumidity[0], elevation);
        for (int h = 1; h < 24; ++h) {
            t = utc_cal.time(2000, 06, 2, h, 00, 0, 0); // June
            rad.net_radiation(rad_r, lat, t, slope,aspect, temperature[h-1], rhumidity[h-1], elevation);
            std::cout<<"lw: "<<rad_r.lw_radiation<<std::endl;
            std::cout<<"sw: "<<rad_r.sw_radiation<<std::endl;
            std::cout<<"net: "<<rad_r.net_radiation<<std::endl;
            std::cout<<"-----------------"<<std::endl;
            // TODO: move time coefficient to radiation model, think how to choose it, should I have time_axis everywhere?
            pm_calculator.reference_evapotranspiration_asce(pm_r,rad_r.net_radiation*0.0036,temperature[h-1], rhumidity[h-1],elevation, windspeed[h-1]);
            std::cout<<"et_ref: "<<pm_r.et_ref<<std::endl;
            std::cout<<"-----------------"<<std::endl;
            av_et_ref.add(pm_r.et_ref, h);

        }
        std::cout << "et_ref_soh_av: " << av_et_ref.result() << std::endl;
        FAST_CHECK_EQ(av_et_ref.result(), doctest::Approx(0.4).epsilon(0.005));

    }
    TEST_CASE("hourly_standartized"){

        //===========================//
        // getting radiation  //
        radiation::parameter rad_p;
        radiation::response rad_r;
        rad_p.albedo = 0.23;
        rad_p.turbidity = 1.0;
        radiation::calculator<radiation::parameter,radiation::response> rad(rad_p);
        calendar utc_cal;
        // Greeley, Colorado weather station
        double lat = 40.41;
        double elevation = 1462.4;
        double ht = 1.68;
        double hws = 3.0;
        utctime t;
        // checking for horizontal surface Eugene, OR, p.64, fig.1b
        arma::vec surface_normal({0.0,0.0,1.0});
        double slope = 0.0;
        double aspect = 0.0;
        utctime ta;

        ta = utc_cal.time(2000, 06, 2, 00, 00, 0, 0);
        //rad.psw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);




        //evapotranspiraiton PM
        penman_monteith::parameter pm_p;
        penman_monteith::response pm_r;
        penman_monteith::calculator<penman_monteith::parameter,penman_monteith::response> pm_calculator(pm_p);
        trapezoidal_average av_et_ref;
        av_et_ref.initialize(pm_r.et_ref, 0.0);
        // ref.: ASCE=EWRI Appendix C: Example Calculation of ET
        double temperature[23] = {16.5, 15.4, 15.5, 13.5, 13.2, 16.2, 20.0, 22.9, 26.4, 28.2, 29.8, 30.9, 31.8, 32.5, 32.9, 32.4, 30.2, 30.6, 28.3, 25.9, 23.9, 20.1, 19.9};
        double vap_pressure[23] = {1.26, 1.34, 1.31, 1.26, 1.24, 1.31, 1.36, 1.39, 1.25, 1.17, 1.03, 1.02, 0.98, 0.87, 0.86, 0.93, 1.14, 1.27, 1.27, 1.17, 1.20, 1.10, 1.05};
        double windspeed[23] = {0.5, 1.0, 0.68, 0.69, 0.29, 1.24, 1.28, 0.88, 0.72, 1.52, 1.97, 2.07, 2.76, 2.990, 3.10, 2.77, 3.41, 2.78, 2.95, 3.27, 2.86, 2.7, 2.0};
        double radiation_m[23] = {0.0, 0.0, 0.0, 0.0, 0.03, 0.46, 1.09, 1.74, 2.34, 2.84, 3.25, 3.21, 3.34, 2.96, 2.25, 1.35, 0.88, 0.79, 0.27, 0.03, 0.0, 0.0};
        double svp[23];
        double rhumidity[23];
        for (int i=0;i<23;i++){
            svp[i] = hydro_functions::svp(temperature[i]);
            std::cout<<"svp: "<<svp[i]<<std::endl; // matched
            rhumidity[i] = svp[i]*100/vap_pressure[i];
        }

//        rad.net_radiation_asce_ewri(rad_r, lat, ta, slope, aspect, temperature[0], rhumidity[0], elevation, radiation_m[0]);
        rad.net_radiation(rad_r, lat, ta, slope, aspect, temperature[0], rhumidity[0], elevation, radiation_m[0]/0.0036);
        std::cout<<"========== standartized ===========" <<std::endl;
        for (int h = 1; h < 24; ++h) {
            t = utc_cal.time(2000, 06, 2, h, 00, 0, 0); // June
//            rad.net_radiation_asce_ewri(rad_r, lat, t, slope,aspect, temperature[h-1], rhumidity[h-1], elevation, radiation_m[h-1]);
            rad.net_radiation(rad_r, lat, t, slope,aspect, temperature[h-1], rhumidity[h-1], elevation, radiation_m[h-1]/0.0036);
            //std::cout<<"-----------------"<<std::endl;
            // TODO: move time coefficient to radiation model, think how to choose it, should I have time_axis everywhere?
            pm_calculator.reference_evapotranspiration_st(pm_r,rad_r.net_radiation*0.0036,temperature[h-1], rhumidity[h-1],elevation, hws, windspeed[h-1],ht);
            std::cout<<"lw: "<<rad_r.lw_radiation<<" | "<<"sw: "<<rad_r.sw_radiation<<" | "<<"net: "<<rad_r.net_radiation<<" || "<<"et_ref: "<<pm_r.et_ref<<std::endl;
            std::cout<<"--------------------------------------------------------------------------"<<std::endl;
            av_et_ref.add(pm_r.et_ref, h);

        }
        std::cout << "et_ref_soh_av: " << av_et_ref.result() << std::endl;
        FAST_CHECK_EQ(av_et_ref.result(), doctest::Approx(0.4).epsilon(0.005));

    }

}
