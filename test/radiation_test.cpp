#include "test_pch.h"
#include "core/utctime_utilities.h"
#include <armadillo>
#include <vector>
#include <chrono>
#include <boost/math/constants/constants.hpp>
#include "core/radiation.h"
#include <cmath>
#include <random>
#include <tuple>

//namespace shyft::core::radiation{
//    /** \tparam C is a cell, like shyft-cell, */
//    template <class C>
//    vector<arma::vec> surface_normal( const vector<C>& cells){
//        vector<arma::vec> r;
//        for(const auto&c:cells) {
//            double x=c.geo.mid_point().x;
//            r.push_back(arma::vec({1.0,1.0,1.0}));
//        }
//        return r;
//    }
//}
namespace rasputin {
    using namespace std;
    using Point = std::array<double, 3>;
    using PointList = std::vector<Point>;
    using VectorList = PointList;
    using Vector = Point;
    using Face = std::array<int, 3>;
    using FaceList = std::vector<Face>;

    VectorList surface_normals(const PointList &pts, const FaceList &faces) {
        VectorList result;
        result.reserve(faces.size());
        for (const auto face: faces) {
            const auto p0 = pts[face[0]];
            const auto p1 = pts[face[1]];
            const auto p2 = pts[face[2]];
            const arma::vec::fixed<3> v0 = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
            const arma::vec::fixed<3> v1 = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
            const arma::vec::fixed<3> n = arma::cross(v0/arma::norm(v0), v1/arma::norm(v1));
            result.emplace_back(n.at(2) >= 0.0 ? Vector{n.at(0), n.at(1), n.at(2)} : Vector{-n.at(0), -n.at(1), -n.at(2)});
        }
        return result;
    }
    std::vector<double> slopes(const PointList &pts, const FaceList &faces){
        VectorList normals = surface_normals(pts,faces);
        std::vector<double> result;
        result.reserve(normals.size());
        for (const auto &normal: normals) {
            result.push_back(atan2(pow(pow(normal.at(0),2)+pow(normal.at(1),2),0.5),normal.at(2)));
        }
        return result;
    }
    std::vector<double> aspects(const PointList &pts, const FaceList &faces){
        VectorList normals = surface_normals(pts,faces);
        std::vector<double> result;
        result.reserve(normals.size());
        for (const auto &normal: normals) {
            result.push_back(atan2(normal.at(1),normal.at(0)));
        }
        return result;
    }

    double elevation{0.0};
}

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

        double result() const { /*std::cout<<" times "<<(t_a-t_start)<<std::endl; */return area/(t_a - t_start); }
    };

}


TEST_SUITE("radiation") {
    using shyft::core::radiation::parameter;
    using shyft::core::radiation::response;
    using shyft::core::radiation::calculator;
//    using shyft::core::radiation::surface_normal;
    using shyft::core::calendar;
    using shyft::core::utctime;
    using shyft::test::trapezoidal_average;
    bool verbose= getenv("SHYFT_VERBOSE")!=nullptr;
//    bool verbose = true;
    // test basics: creation, etc



    TEST_CASE("geometry"){

        std::array<double,3> point1{{0.0,0.0,0.0}};
        std::array<double,3> point2{{1.0,0.0,0.0}};
        std::array<double,3> point3{{0.0,1.0,0.0}};
        std::array<int,3> face{{0,1,2}};
        std::vector<std::array<double,3>> points;
        points.push_back(point1);
        points.push_back(point2);
        points.push_back(point3);

        std::vector<std::array<int,3>> faces;
        faces.push_back(face);

        std::vector<std::array<double,3>> normals;
        normals = rasputin::surface_normals(points,faces);
        if (verbose){
            for (const auto n:normals)
                std::cout<<n.at(0)<<"; "<<n.at(1) << "; "<<n.at(2)<<std::endl;

            std::vector<double> slopes;

            slopes = rasputin::slopes(points, faces);
            for (const auto s:slopes)
                std::cout<<s<<std::endl;
        }


    }

    TEST_CASE("check_solar_radiation_horizontal_inst"){
        parameter p;
        response r;
        p.albedo = 0.2;
        p.turbidity = 1.0;
        calculator<parameter,response> rad(p);
        calendar utc_cal;
        double lat = 44.0;
        utctime t;
        // checking for horizontal surface Eugene, OR, p.64, fig.1b
        double slope = 0.0;
        double aspect = 0.0;
        utctime ta;
        trapezoidal_average av_rahor;
        trapezoidal_average av_ra;
        trapezoidal_average av_rs;
        std::uniform_real_distribution<double> ur(150.0, 390.0);
        std::default_random_engine gen;
//
        SUBCASE("June_translated") {
//            std::cout << "========= June translated ========" << std::endl;
            ta = utc_cal.time(2002, 06, 21, 00, 00, 0, 0);
            //rad.psw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);
            rad.net_radiation(r, lat, ta, slope, aspect, 20.0, 50.0, 150.0, ur(gen));
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rs.initialize(r.sw_radiation, 0.0);
            for (int h = 0; h < 24; ++h) {
                t = utc_cal.time(2002, 06, 21, h, 00, 0, 0); // June
                rad.net_radiation(r, lat, t, slope,aspect, 20.0, 50.0, 150.0, ur(gen));
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rs.add(r.sw_radiation, h);
            }

//            std::cout << "ra: " << av_ra.result() << std::endl;
//            std::cout << "rs: " << av_rs.result() << std::endl;
            FAST_CHECK_EQ(av_ra.result(), doctest::Approx(500.0).epsilon(0.05));
            FAST_CHECK_EQ(av_rahor.result(), doctest::Approx(av_ra.result()).epsilon(0.05));
            //FAST_CHECK_EQ(av_rs.result(), doctest::Approx(370.0).epsilon(0.05));

        }
        SUBCASE("June") {
//
            ta = utc_cal.time(2002, 06, 21, 00, 00, 0, 0);
            //rad.psw_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);
            rad.net_radiation(r, lat, ta, slope,aspect, 20.0, 50.0, 150.0);
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rs.initialize(r.sw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(2002, 06, 21, h, 00, 0, 0); // June
                //rad.psw_radiation(r, lat, t, slope,aspect, 20.0, 50.0, 150.0);
                rad.net_radiation(r, lat, t, slope, aspect, 20.0, 50.0, 150.0);
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rs.add(r.sw_radiation, h);
            }

            if (verbose){
std::cout << "========= Horizontal =======" << std::endl;
std::cout << "========= June ========" << std::endl;
                std::cout << "ra: " << av_ra.result()<<std::endl;
                std::cout << "rs: " << av_rs.result()<<std::endl;
            }
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(500.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rahor.result(), doctest::Approx(av_ra.result()).epsilon(0.05));
                    FAST_CHECK_EQ(av_rs.result(), doctest::Approx(370.0).epsilon(0.05));

        }
        SUBCASE("January") {
//
            ta = utc_cal.time(2002, 01, 1, 00, 00, 0, 0);
            rad.net_radiation(r, lat, ta, slope,aspect, 20.0, 50.0, 150.0);
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rs.initialize(r.sw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(2002, 01, 1, h, 00, 0, 0); // January
                rad.net_radiation(r, lat, t, slope,aspect, 20.0, 50.0, 150.0);
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rs.add(r.sw_radiation, h);
            }

            if (verbose){
std::cout << "========= January =======" << std::endl;
            std::cout << "ra: " << av_ra.result()<<std::endl;
            std::cout << "rs: " << av_rs.result()<<std::endl;
            }
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(130.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rahor.result(), doctest::Approx(av_ra.result()).epsilon(0.05));
                    FAST_CHECK_EQ(av_rs.result(), doctest::Approx(80.0).epsilon(0.05));
        }
        SUBCASE("December") {
//
            ta = utc_cal.time(2002, 12, 21, 00, 00, 0, 0);
            rad.net_radiation(r, lat, ta, slope,aspect, 20.0, 50.0, 150.0);
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rs.initialize(r.sw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(2002, 12, 21, h, 00, 0, 0); // January
                rad.net_radiation(r, lat, t, slope,aspect, 20.0, 50.0, 150.0);
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rs.add(r.sw_radiation, h);
            }

            if (verbose){
std::cout << "========= December =======" << std::endl;
            std::cout << "ra: " << av_ra.result()<<std::endl;
            std::cout << "rs: " << av_rs.result()<<std::endl;
            }
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(130.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rahor.result(), doctest::Approx(av_ra.result()).epsilon(0.05));
                    FAST_CHECK_EQ(av_rs.result(), doctest::Approx(80.0).epsilon(0.05));
        }

    }
    TEST_CASE("check_solar_radiation_horizontal_step"){
    parameter p;
    response r;
    p.albedo = 0.2;
    p.turbidity = 1.0;
    calculator<parameter,response> rad(p);
    calendar utc_cal;
    double lat = 44.0;
    utctime t1,t2;
    // checking for horizontal surface Eugene, OR, p.64, fig.1b
    double slope = 0.0;
    double aspect = 0.0;
    utctime ta1,ta2;
    std::uniform_real_distribution<double> ur(150.0, 390.0);
    std::default_random_engine gen;
    //
    SUBCASE("June_translated") {
        //
        double rastep = 0.0;
        double rsostep = 0.0;
        for (int h = 1; h < 24; ++h) {
            t1 = utc_cal.time(2002, 06, 21, h-1, 00, 0, 0); // June
            t2 = utc_cal.time(2002, 06, 21, h, 00, 0, 0); // June
            rad.net_radiation_step(r, lat, t1,t2, slope, aspect, 20.0, 50.0, 150.0,ur(gen));
            rastep+=r.ra;
            rsostep+=r.sw_radiation;
        }

        if (verbose){
std::cout << "========= Horizontal =======" << std::endl;
std::cout << "========= June translated ========" << std::endl;
        std::cout << "ra: " << rastep<<std::endl;
        std::cout << "rs: " << rsostep<<std::endl;
        }


    }
    SUBCASE("June") {
//
        double rastep = 0.0;
        double rsostep = 0.0;
        for (int h = 1; h < 24; ++h) {
            t1 = utc_cal.time(2002, 06, 21, h-1, 00, 0, 0); // June
            t2 = utc_cal.time(2002, 06, 21, h, 00, 0, 0); // June
            rad.net_radiation_step(r, lat, t1,t2, slope, aspect, 20.0, 50.0, 150.0);
            rastep+=r.ra;
            rsostep+=r.sw_radiation;
        }

        if (verbose){
std::cout << "========= June step ========" << std::endl;
        std::cout << "ra: " << rastep<<std::endl;
        std::cout << "rs: " << rsostep<<std::endl;
        }
        FAST_CHECK_EQ(rastep, doctest::Approx(500.0).epsilon(0.05));
        FAST_CHECK_EQ(rsostep, doctest::Approx(370.0).epsilon(0.05));


    }
    SUBCASE("January") {
//
        double rastep = 0.0;
        double rsostep = 0.0;
        for (int h = 1; h < 24; ++h) {
            t1 = utc_cal.time(2002, 01, 1, h-1, 00, 0, 0); // June
            t2 = utc_cal.time(2002, 01, 1, h, 00, 0, 0); // June
            rad.net_radiation_step(r, lat, t1,t2, slope, aspect, 20.0, 50.0, 150.0);
            rastep+=r.ra;
            rsostep+=r.sw_radiation;
        }
        if (verbose){
std::cout << "========= January =======" << std::endl;
        std::cout << "ra: " << rastep<<std::endl;
        std::cout << "rs: " << rsostep<<std::endl;
        }
        FAST_CHECK_EQ(rastep, doctest::Approx(130.0).epsilon(0.05));
        FAST_CHECK_EQ(rsostep, doctest::Approx(75.0).epsilon(0.05));
    }
    SUBCASE("December") {
//
        double rastep = 0.0;
        double rsostep = 0.0;
        for (int h = 1; h < 24; ++h) {
            t1 = utc_cal.time(2002, 12, 21, h-1, 00, 0, 0); // June
            t2 = utc_cal.time(2002, 12, 21, h, 00, 0, 0); // June
            rad.net_radiation_step(r, lat, t1,t2, slope, aspect, 20.0, 50.0, 150.0);
            rastep+=r.ra;
            rsostep+=r.sw_radiation;
        }
        if (verbose){
std::cout << "========= December =======" << std::endl;
        std::cout << "ra: " << rastep<<std::endl;
        std::cout << "rs: " << rsostep<<std::endl;
        }
        FAST_CHECK_EQ(rastep, doctest::Approx(130.0).epsilon(0.05));
        FAST_CHECK_EQ(rsostep, doctest::Approx(75.0).epsilon(0.05));
    }

}


   TEST_CASE("check_solar_radiation_slope_45s"){
        parameter p;
        response r;
        p.albedo = 0.2;
        p.turbidity = 1.0;
        calculator<parameter,response> rad(p);
        calendar utc_cal;
        double lat = 44.0;
        utctime t;
        // checking for horizontal surface Eugene, OR, p.64, fig.1d
        // 24h  average radiation
        double slope = 45;//*shyft::core::radiation::pi/180; // 45 S
       // double proj = sin(slope);
        double aspect = 0.0;//*shyft::core::radiation::pi/180;// facing south
        //arma::vec surface_normal({proj*cos(aspect),proj*sin(aspect),cos(slope)});
        utctime ta;
        trapezoidal_average av_rahor;
        trapezoidal_average av_ra;
        trapezoidal_average av_rs;
        std::uniform_real_distribution<double> ur(100.0, 390.0);
        std::default_random_engine gen;
//
        SUBCASE("June_translated") {
//
            ta = utc_cal.time(2002, 06, 21, 00, 00, 0, 0);
            rad.net_radiation(r, lat, ta, slope,aspect, 20.0, 50.0, 150.0);
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rs.initialize(r.sw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(2002, 06, 21, h, 00, 0, 0); // June
                rad.net_radiation(r, lat, t, slope,aspect, 20.0, 50.0, 150.0);
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rs.add(r.sw_radiation, h);
            }

            if (verbose){
std::cout << "========= Slope 45S =======" << std::endl;
std::cout << "========= June translated ========" << std::endl;
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
            std::cout << "sun_rise: " << rad.sun_rise() << std::endl;
            std::cout << "sun_set: " << rad.sun_set() << std::endl;
            }
            FAST_CHECK_EQ(av_ra.result(), doctest::Approx(390.0).epsilon(0.05));
            //FAST_CHECK_EQ(av_rs.result(), doctest::Approx(310.0).epsilon(0.05));
        }
        SUBCASE("June") {
//
            ta = utc_cal.time(2002, 06, 21, 00, 00, 0, 0);
//            rad.net_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);
            rad.net_radiation(r, lat, ta, slope,aspect, 20.0, 50.0, 150.0);
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rs.initialize(r.sw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(2002, 06, 21, h, 00, 0, 0); // June
                rad.net_radiation(r, lat, t, slope,aspect, 20.0, 50.0, 150.0);
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rs.add(r.sw_radiation, h);
            }
            if (verbose){
std::cout << "========= June ========" << std::endl;
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
            std::cout << "sun_rise: " << rad.sun_rise() << std::endl;
            std::cout << "sun_set: " << rad.sun_set() << std::endl;
            }
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(390.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rs.result(), doctest::Approx(310.0).epsilon(0.05));
        }
        SUBCASE("January") {
//
            ta = utc_cal.time(2002, 01, 1, 00, 00, 0, 0);
//            rad.net_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);
            rad.net_radiation(r, lat, ta, slope,aspect, 20.0, 50.0, 150.0);
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rs.initialize(r.sw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(2002, 01, 1, h, 00, 0, 0); // June
                rad.net_radiation(r, lat, t, slope,aspect, 20.0, 50.0, 150.0);
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rs.add(r.sw_radiation, h);
            }
            if (verbose){
std::cout << "========= January ========" << std::endl;
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
            std::cout << "sun_rise: " << rad.sun_rise() << std::endl;
            std::cout << "sun_set: " << rad.sun_set() << std::endl;
            }
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(390.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rs.result(), doctest::Approx(270.0).epsilon(0.05));
        }
        SUBCASE("December") {
//
            ta = utc_cal.time(2002, 12, 12, 00, 00, 0, 0);
//            rad.net_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);
            rad.net_radiation(r, lat, ta, slope,aspect, 20.0, 50.0, 150.0);
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rs.initialize(r.sw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(2002, 12, 12, h, 00, 0, 0); // June
                rad.net_radiation(r, lat, t, slope,aspect, 20.0, 50.0, 150.0);
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rs.add(r.sw_radiation, h);
            }
            if (verbose){
std::cout << "========= December ========" << std::endl;
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
            }
                    FAST_CHECK_EQ(av_ra.result(), doctest::Approx(390.0).epsilon(0.05));
                    FAST_CHECK_EQ(av_rs.result(), doctest::Approx(270.0).epsilon(0.05));
        }

    }
    TEST_CASE("check_solar_radiation_45s_step"){
        parameter p;
        response r;
        p.albedo = 0.2;
        p.turbidity = 1.0;
        calculator<parameter,response> rad(p);
        calendar utc_cal;
        double lat = 44.0;
        utctime t1,t2;
        // checking for horizontal surface Eugene, OR, p.64, fig.1b
        double slope = 45.0;
        double aspect = 0.0;
        utctime ta1,ta2;
        std::uniform_real_distribution<double> ur(150.0, 390.0);
        std::default_random_engine gen;

        SUBCASE("June") {
//
            double rastep = 0.0;
            double rsostep = 0.0;
            for (int h = 1; h < 24; ++h) {
            t1 = utc_cal.time(2002, 06, 21, h-1, 00, 0, 0); // June
            t2 = utc_cal.time(2002, 06, 21, h, 00, 0, 0); // June
            rad.net_radiation_step(r, lat, t1,t2, slope, aspect, 20.0, 50.0, 150.0);
            rastep+=r.ra;
            rsostep+=r.sw_radiation;
            }

            if (verbose){
            std::cout << "========= June step ========" << std::endl;
            std::cout << "rastep: " << rastep << std::endl;
            std::cout << "rsostep: " << rsostep << std::endl;
            }
            FAST_CHECK_EQ(rastep, doctest::Approx(390.0).epsilon(0.05));
            FAST_CHECK_EQ(rsostep, doctest::Approx(310.0).epsilon(0.05));


        }
        SUBCASE("January") {
//
            double rastep = 0.0;
            double rsostep = 0.0;
            for (int h = 1; h < 24; ++h) {
            t1 = utc_cal.time(2002, 01, 1, h-1, 00, 0, 0); // June
            t2 = utc_cal.time(2002, 01, 1, h, 00, 0, 0); // June
            rad.net_radiation_step(r, lat, t1,t2, slope, aspect, 20.0, 50.0, 150.0);
            rastep+=r.ra;
            rsostep+=r.sw_radiation;
            }
            if (verbose){
            std::cout << "========= January =======" << std::endl;
            std::cout << "rastep: " << rastep << std::endl;
            std::cout << "rsostep: " << rsostep << std::endl;
            }
            FAST_CHECK_EQ(rastep, doctest::Approx(370.0).epsilon(0.05));
            FAST_CHECK_EQ(rsostep, doctest::Approx(180.0).epsilon(0.05));
        }
        SUBCASE("December") {

            double rastep = 0.0;
            double rsostep = 0.0;
            for (int h = 1; h < 24; ++h) {
            t1 = utc_cal.time(2002, 12, 21, h-1, 00, 0, 0); // June
            t2 = utc_cal.time(2002, 12, 21, h, 00, 0, 0); // June
            rad.net_radiation_step(r, lat, t1,t2, slope, aspect, 20.0, 50.0, 150.0);
            rastep+=r.ra;
            rsostep+=r.sw_radiation;
            }
            if (verbose){
            std::cout << "========= December =======" << std::endl;
            std::cout << "rastep: " << rastep << std::endl;
            std::cout << "rsostep: " << rsostep << std::endl;
            }
            FAST_CHECK_EQ(rastep, doctest::Approx(370.0).epsilon(0.05));
            FAST_CHECK_EQ(rsostep, doctest::Approx(180.0).epsilon(0.05));
        }

    }

    TEST_CASE("check_solar_radiation_slope_90s"){
        parameter p;
        response r;
        p.albedo = 0.05;
        p.turbidity = 1.0;
        calculator<parameter,response> rad(p);
        calendar utc_cal;
        double lat = 44.0;
        utctime t;
        // checking for horizontal surface Eugene, OR, p.64, fig.1d
        // 24h  average radiation
        double slope = 90;//*shyft::core::radiation::pi/180; // 45 S
        // double proj = sin(slope);
        double aspect = 0.0;//*shyft::core::radiation::pi/180;// facing south
        //arma::vec surface_normal({proj*cos(aspect),proj*sin(aspect),cos(slope)});
        utctime ta;
        trapezoidal_average av_rahor;
        trapezoidal_average av_ra;
        trapezoidal_average av_rs;
        std::uniform_real_distribution<double> ur(100.0, 390.0);
        std::default_random_engine gen;

        SUBCASE("June_translated") {

            ta = utc_cal.time(2002, 06, 21, 00, 00, 0, 0);
            rad.net_radiation(r, lat, ta, slope,aspect, 20.0, 50.0, 150.0);
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rs.initialize(r.sw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(2002, 06, 21, h, 00, 0, 0); // June
                rad.net_radiation(r, lat, t, slope,aspect, 20.0, 50.0, 150.0);
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rs.add(r.sw_radiation, h);
            }
            if (verbose){
            std::cout << "========= Slope 90S =======" << std::endl;
            std::cout << "========= June translated ========" << std::endl;
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
            std::cout << "sun_rise: " << rad.sun_rise() << std::endl;
            std::cout << "sun_set: " << rad.sun_set() << std::endl;
            }
            FAST_CHECK_EQ(av_ra.result(), doctest::Approx(110.0).epsilon(0.05));
        //FAST_CHECK_EQ(av_rs.result(), doctest::Approx(310.0).epsilon(0.05));
        }
        SUBCASE("June") {

            ta = utc_cal.time(2002, 06, 21, 00, 00, 0, 0);
            //            rad.net_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);
            rad.net_radiation(r, lat, ta, slope,aspect, 20.0, 50.0, 150.0);
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rs.initialize(r.sw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(2002, 06, 21, h, 00, 0, 0); // June
                rad.net_radiation(r, lat, t, slope,aspect, 20.0, 50.0, 150.0);
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rs.add(r.sw_radiation, h);
            }
            if (verbose){
            std::cout << "========= June ========" << std::endl;

            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
            std::cout << "sun_rise: " << rad.sun_rise() << std::endl;
            std::cout << "sun_set: " << rad.sun_set() << std::endl;
            }
            FAST_CHECK_EQ(av_ra.result(), doctest::Approx(110.0).epsilon(0.05));
            FAST_CHECK_EQ(av_rs.result(), doctest::Approx(90.0).epsilon(0.05));
        }
        SUBCASE("January") {

            ta = utc_cal.time(2002, 01, 1, 00, 00, 0, 0);
            //            rad.net_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);
            rad.net_radiation(r, lat, ta, slope,aspect, 20.0, 50.0, 150.0);
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rs.initialize(r.sw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(2002, 01, 1, h, 00, 0, 0); // June
                rad.net_radiation(r, lat, t, slope,aspect, 20.0, 50.0, 150.0);
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rs.add(r.sw_radiation, h);
            }
            if (verbose){
            std::cout << "========= January ========" << std::endl;
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
            std::cout << "sun_rise: " << rad.sun_rise() << std::endl;
            std::cout << "sun_set: " << rad.sun_set() << std::endl;
            }
            FAST_CHECK_EQ(av_ra.result(), doctest::Approx(410.0).epsilon(0.05));
            FAST_CHECK_EQ(av_rs.result(), doctest::Approx(300.0).epsilon(0.05));
        }
        SUBCASE("December") {
            ;
            ta = utc_cal.time(2002, 12, 12, 00, 00, 0, 0);
            //            rad.net_radiation(r, lat, ta, surface_normal, 20.0, 50.0, 150.0);
            rad.net_radiation(r, lat, ta, slope,aspect, 20.0, 50.0, 150.0);
            av_rahor.initialize(rad.ra_radiation_hor(), 0.0);
            av_ra.initialize(rad.ra_radiation(), 0.0);
            av_rs.initialize(r.sw_radiation, 0.0);
            for (int h = 1; h < 24; ++h) {
                t = utc_cal.time(2002, 12, 12, h, 00, 0, 0); // June
                rad.net_radiation(r, lat, t, slope,aspect, 20.0, 50.0, 150.0);
                av_rahor.add(rad.ra_radiation_hor(), h);
                av_ra.add(rad.ra_radiation(), h);
                av_rs.add(r.sw_radiation, h);
            }
            if (verbose){
            std::cout << "========= December ========" << std::endl;
            std::cout << "rahor: " << av_rahor.result() << std::endl;
            std::cout << "ra: " << av_ra.result() << std::endl;
            std::cout << "rs: " << av_rs.result() << std::endl;
            }
            FAST_CHECK_EQ(av_ra.result(), doctest::Approx(410.0).epsilon(0.05));
            FAST_CHECK_EQ(av_rs.result(), doctest::Approx(300.0).epsilon(0.05));
        }

    }

    TEST_CASE("check_solar_radiation_90s_step"){
        parameter p;
        response r;
        p.albedo = 0.05;
        p.turbidity = 1.0;
        calculator<parameter,response> rad(p);
        calendar utc_cal;
        double lat = 44.0;
        utctime t1,t2;
        // checking for horizontal surface Eugene, OR, p.64, fig.1b
        double slope = 90.0;
        double aspect = 0.0;
        utctime ta1,ta2;
        std::uniform_real_distribution<double> ur(150.0, 390.0);
        std::default_random_engine gen;

        SUBCASE("June") {

            double rastep = 0.0;
            double rsostep = 0.0;
            for (int h = 1; h < 24; ++h) {
            t1 = utc_cal.time(2002, 06, 21, h-1, 00, 0, 0); // June
            t2 = utc_cal.time(2002, 06, 21, h, 00, 0, 0); // June
            rad.net_radiation_step(r, lat, t1,t2, slope, aspect, 20.0, 50.0, 150.0);
            rastep+=r.ra;
            rsostep+=r.sw_radiation;
            }
            if (verbose){
            std::cout << "========= slope 90s =======" << std::endl;
            std::cout << "========= June step ========" << std::endl;
            std::cout << "rastep: " << rastep << std::endl;
            std::cout << "rsostep: " << rsostep << std::endl;
            }
            FAST_CHECK_EQ(rastep, doctest::Approx(110.0).epsilon(0.05));
            FAST_CHECK_EQ(rsostep, doctest::Approx(100.0).epsilon(0.05));


        }
        SUBCASE("January") {

            double rastep = 0.0;
            double rsostep = 0.0;
            for (int h = 1; h < 24; ++h) {
            t1 = utc_cal.time(2002, 01, 1, h-1, 00, 0, 0); // June
            t2 = utc_cal.time(2002, 01, 1, h, 00, 0, 0); // June
            rad.net_radiation_step(r, lat, t1,t2, slope, aspect, 20.0, 50.0, 150.0);
            rastep+=r.ra;
            rsostep+=r.sw_radiation;
            }
            if (verbose){
            std::cout << "========= January =======" << std::endl;
            std::cout << "rastep: " << rastep << std::endl;
            std::cout << "rsostep: " << rsostep << std::endl;
            }
            FAST_CHECK_EQ(rastep, doctest::Approx(390.0).epsilon(0.05));
            FAST_CHECK_EQ(rsostep, doctest::Approx(180.0).epsilon(0.05));
        }
        SUBCASE("December") {

            double rastep = 0.0;
            double rsostep = 0.0;
            for (int h = 1; h < 24; ++h) {
            t1 = utc_cal.time(2002, 12, 21, h-1, 00, 0, 0); // June
            t2 = utc_cal.time(2002, 12, 21, h, 00, 0, 0); // June
            rad.net_radiation_step(r, lat, t1,t2, slope, aspect, 20.0, 50.0, 150.0);
            rastep+=r.ra;
            rsostep+=r.sw_radiation;
            }
            if (verbose){
            std::cout << "========= December =======" << std::endl;
            std::cout << "rastep: " << rastep << std::endl;
            std::cout << "rsostep: " << rsostep << std::endl;
            }
            FAST_CHECK_EQ(rastep, doctest::Approx(390.0).epsilon(0.05));
            FAST_CHECK_EQ(rsostep, doctest::Approx(180.0).epsilon(0.05));
        }

    }
    TEST_CASE("check_solar_radiation_90n_step"){
        parameter p;
        response r;
        p.albedo = 0.05;
        p.turbidity = 1.0;
        calculator<parameter,response> rad(p);
        calendar utc_cal;
        double lat = 44.0;
        utctime t1,t2;
        // checking for horizontal surface Eugene, OR, p.64, fig.1b
        double slope = 90.0;
        double aspect = 180.0;
        utctime ta1,ta2;
        std::uniform_real_distribution<double> ur(150.0, 390.0);
        std::default_random_engine gen;

        SUBCASE("June") {

            double rastep = 0.0;
            double rsostep = 0.0;
            for (int h = 1; h < 24; ++h) {
                t1 = utc_cal.time(2002, 06, 21, h-1, 00, 0, 0); // June
                t2 = utc_cal.time(2002, 06, 21, h, 00, 0, 0); // June
                rad.net_radiation_step(r, lat, t1,t2, slope, aspect, 20.0, 50.0, 150.0);
                rastep+=r.ra;
                rsostep+=r.sw_radiation;
            }
            if (verbose){
            std::cout << "========= slope 90n =======" << std::endl;
            std::cout << "========= June step ========" << std::endl;
            std::cout << "rastep: " << rastep << std::endl;
            std::cout << "rsostep: " << rsostep << std::endl;
            }
            FAST_CHECK_EQ(rastep, doctest::Approx(110.0).epsilon(0.05));
            FAST_CHECK_EQ(rsostep, doctest::Approx(50.0).epsilon(0.05));


        }
        SUBCASE("January") {

            double rastep = 0.0;
            double rsostep = 0.0;
            for (int h = 1; h < 24; ++h) {
                t1 = utc_cal.time(2002, 01, 1, h-1, 00, 0, 0); // June
                t2 = utc_cal.time(2002, 01, 1, h, 00, 0, 0); // June
                rad.net_radiation_step(r, lat, t1,t2, slope, aspect, 20.0, 50.0, 150.0);
                rastep+=r.ra;
                rsostep+=r.sw_radiation;
            }
            if (verbose){
            std::cout << "========= January =======" << std::endl;
            std::cout << "rastep: " << rastep << std::endl;
            std::cout << "rsostep: " << rsostep << std::endl;
            }
            FAST_CHECK_EQ(rastep, doctest::Approx(0.0).epsilon(0.05));
            FAST_CHECK_EQ(rsostep, doctest::Approx(5.0).epsilon(0.05));
        }
        SUBCASE("December") {

            double rastep = 0.0;
            double rsostep = 0.0;
            for (int h = 1; h < 24; ++h) {
                t1 = utc_cal.time(2002, 12, 21, h-1, 00, 0, 0); // June
                t2 = utc_cal.time(2002, 12, 21, h, 00, 0, 0); // June
                rad.net_radiation_step(r, lat, t1,t2, slope, aspect, 20.0, 50.0, 150.0);
                rastep+=r.ra;
                rsostep+=r.sw_radiation;
            }
            if (verbose){
            std::cout << "========= December =======" << std::endl;
            std::cout << "rastep: " << rastep << std::endl;
            std::cout << "rsostep: " << rsostep << std::endl;
            }
            FAST_CHECK_EQ(rastep, doctest::Approx(0.0).epsilon(0.05));
            FAST_CHECK_EQ(rsostep, doctest::Approx(5.0).epsilon(0.05));
        }

    }
//
}
