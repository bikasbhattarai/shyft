/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
///	Copyright 2012 Statkraft Energi A/S
///
///	This file is part of Shyft.
///
///	Shyft is free software: you can redistribute it and/or modify it under the terms of
/// the GNU Lesser General Public License as published by the Free Software Foundation,
/// either version 3 of the License, or (at your option) any later version.
///
///	Shyft is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
/// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
/// PURPOSE. See the GNU Lesser General Public License for more details.
///
///	You should have received a copy of the GNU Lesser General Public License along with
/// Shyft, usually located under the Shyft root directory in two files named COPYING.txt
/// and COPYING_LESSER.txt.	If not, see <http://www.gnu.org/licenses/>.
///
/// Adapted from early enki method programmed by Kolbjørn Engeland and Sjur Kolberg
///
#pragma once

#include "core/utctime_utilities.h"
#include "core/geo_cell_data.h"
#include <armadillo>
#include <vector>

#include <chrono>
#include <boost/math/constants/constants.hpp>
#include "core/hydro_functions.h"
//#include <proj.h>
namespace shyft {

    namespace core {
        // TODO use rasputin namespace

        namespace radiation {
            using namespace std;
            const double pi = boost::math::constants::pi<double>();
            using namespace hydro_functions;

            struct parameter {
                double albedo = 0.1; // average albedo of the surrounding ground surface:0.15-0.25 -- grass, 0.10-0.15 -- coniferous forest, 0.15 - 0.25 -- deciduous forest, 0.04-0.08 -- open water, 0.15-0.35 -- bare soil
                double turbidity = 1.0; // 1.0 - clean, 0.5 -- dusty
                parameter(double albedo = 0.1, double turbidity = 1.0) : albedo(albedo), turbidity(turbidity) {}
            };

            //struct state {}; // No state variables for this method

            struct response {
                double sw_radiation = 0.0; // translated  solar radiation on a sloping surface based on measured horizontal radiation [W/m^2]
                double lw_radiation = 0.0; // long-wave radiation [W/m^2]
                double net_radiation = 0.0; // net radiation [W/m^2]
                double ra = 0.0; // temporary output for extensive testing on python side
                double rah = 0.0;
                double omega1 = 0.0;
                double omega2 = 0.0;
            };

            template<class P, class R>
            struct calculator {
                P param;

                explicit calculator(const P &p) : param(p) {}

                double latitude() const {
                    return phi_ * rad2deg;
                }// latitude, [deg] should be available from cell?/// TODO: add PROJ4 for conversion from cartesian to wgs84

                double ra_radiation() const { return ra_; } // extraterrestrial solar radiation for inclined surface[W/m2]
                double ra_radiation_hor() { return rahor_; } // extraterrestrial solar radiation for horizontal surfaces

                double sun_rise() const { return omega1_24_ * rad2deg; }

                double sun_set() const { return omega2_24_ * rad2deg; }

                /**\brief computes instantaneous net radiation, [W/m2]*/
                void net_radiation(R &response, double latitude, utctime t, double slope=0.0, double aspect = 0.0,
                                    double temperature = 0.0, double rhumidity = 40.0, double elevation = 0.0,
                                    double rsm = 0.0){

                    response.sw_radiation = tsw_radiation(latitude, t, slope, aspect, temperature, rhumidity, elevation,rsm);
                    response.lw_radiation = lw_radiation(temperature, rhumidity);
                    response.net_radiation = response.sw_radiation+response.lw_radiation;
                    response.ra = ra_radiation();
                    response.rah = ra_radiation_hor();
                    response.omega1 = sun_rise();
                    response.omega2 = sun_set();
//                    std::cout<<"calendar time: " <<(utc.calendar_units(t).hour + utc.calendar_units(t).minute / 60.0)<<std::endl;

                }
                /**\brief computes instantaneous net radiation, [W/m2]*/
                void net_radiation_step(R &response, double latitude, utctime t1, utctime t2, double slope=0.0, double aspect = 0.0,
                                   double temperature = 0.0, double rhumidity = 40.0, double elevation = 0.0,
                                   double rsm = 0.0){

                    response.sw_radiation = psw_radiation_step(latitude, t1,t2, slope, aspect, temperature, rhumidity, elevation);
                    response.lw_radiation = lw_radiation(temperature, rhumidity);
                    response.net_radiation = response.sw_radiation+response.lw_radiation;
                    response.ra = ra_radiation();
                    response.rah = ra_radiation_hor();
                    response.omega1 = sun_rise();
                    response.omega2 = sun_set();
//                    std::cout<<"calendar time: " <<(utc.calendar_units(t).hour + utc.calendar_units(t).minute / 60.0)<<std::endl;

                }

                /**\brief net radiation asce-ewri*/
                void net_radiation_asce_ewri(R &response, double latitude, utctime t, double slope=0.0, double aspect = 0.0,
                                             double temperature = 0.0, double rhumidity = 40.0, double elevation = 0.0,
                                             double rsm = 0.0){


                    double phi = latitude;
                    double longz = 105.0;
                    double longitudem = 104.78;
                    double doy = utc.day_of_year(t);
                    double delta = 0.409*sin(2*pi/365*doy-1.39);
                    double oms = acos(-tan(phi)*tan(delta));
                    double b = 2*pi/364*(doy - 81);
                    double t1 = 1;
                    double Sc = 0.1645*sin(2*b) - 0.1255*cos(b) - 0.025*sin(b);
                    double hour = utc.calendar_units(t).hour + utc.calendar_units(t).minute / 60.0;
                    double om = pi/12.0*((hour+0.06667*(longz-longitudem)+Sc)-12.0);
                    double om1 = om - pi*t1/24;
                    double om2 = om + pi*t1/24;
                    double Gsc = 4.92;
                    double dr = 1+0.033*cos(2*pi/365*doy);
                    if (om1 < -oms)
                        om1 = -oms;
                    if (om2 <-oms)
                        om2 = -oms;
                    if (om1>oms)
                        om1 = oms;
                    if (om2>oms)
                        om2  = oms;
                    if (om1>om2)
                        om1 = om2;

                    double ra = 12.0/pi*Gsc*dr*(om2-om1*sin(phi)*sin(delta)+cos(phi)*cos(delta)*(sin(om2)-sin(om1)));
                    if  ((om<-oms) or (om>oms))
                        ra = 0.0;
                    double rso = (0.75+0.00002*elevation)*ra;
                    double avp = actual_vp(temperature,rhumidity);
                    double fcd = std::max(0.05,std::min(1.0,1.35*std::max(0.3,std::min(1.0,rsm/rso))-0.35)); // eq.45
//                    response.sw_radiation = rsm*(1 - param.albedo);
                    response.sw_radiation = psw_radiation(latitude,t,slope,aspect,temperature,rhumidity,elevation)*0.0036;
                    response.lw_radiation = 2.042*pow(10,-10)*fcd*(0.34-0.14*sqrt(avp))*pow(temperature+273.16,4);
                    response.net_radiation = response.sw_radiation-response.lw_radiation;
                }


            private:
                double delta_;
                double omega_;
                double phi_;
                double slope_;
                double aspect_;

                double ra_ = 0.0; // extraterrestrial solar radiation for inclined surface[W/m2]
                double rahor_ = 0.0; // extraterrestrial solar radiation for horizontal surfaces
                double omega1_;
                double omega2_;
                void update_omega12(double omega1, double omega2){
                    omega1_ = omega1;
                    omega2_ = omega2;
                }

                calendar utc;
                double doy_; // day of the yearI
                double lt_; // local time
                double costt_, costthor_;
                double a_, b_, c_,g_,h_;
                double f1_, f2_,f3_,f4_,f5_;
                double omega1_24_, omega2_24_, omega1_24b_, omega2_24b_; //integration limits, actual periods of sun
                double fb_;

                /** \brief computes necessary geometric parameters
                 * \param omega, [rad] -- hour angle
                 * \return cos(theta) -- theta -- angle of incidence             * */
                double costt(double omega) {
                    return -a_ + b_ * cos(omega) + c_ * sin(omega); // eq.14
                }

                /** \brief computes necessary geometric parameters
                 * \param phi, [rad] -- latitude
                 * \param s, [rad] -- slope angle
                 * \param gamma, [rad] -- aspect angle
                 * \return cos(theta) -- theta -- angle of incidence             * */
                void compute_abc(double delta, double phi, double s = 0.0, double gamma = 0.0) {
                    a_ = sin(delta) * cos(phi) * sin(s) * cos(gamma) - sin(delta) * sin(phi) * cos(s); // eq.11a
                    b_ = cos(delta) * cos(phi) * cos(s) + cos(delta) * sin(phi) * sin(s) * cos(gamma);//eq 11b
                    c_ = cos(delta) * sin(s) * sin(gamma);
                    g_ = sin(delta) * sin(phi);
                    h_ = cos(delta) * cos(phi);
                }

                /** \brief computes necessary geometric parameters
                 * \param omega1, [rad] -- tstart
                 * \param omega2, [rad] -- tend
                 * \return f1,f2,f3,f4,f5             * */
                void compute_fs(double omega1, double omega2) {
                    f1_ = sin(omega2)-sin(omega1); // eq.11a
                    f2_ = cos(omega2)-cos(omega1);//eq 11b
                    f3_ = omega2-omega1;
                    f4_ =sin(2*omega2)-sin(2*omega1);
                    f5_ = sin(omega2)*sin(omega2)-sin(omega1)*sin(omega1);

                }
                /** \brief computes necessary geometric parameters
                 * \param phi, [rad] -- latitude
                 * \param s, [rad] -- slope angle
                 * \param gamma, [rad] -- aspect angle
                 * \return cos(theta) -- theta -- angle of incidence             * */
                double compute_beta_step() {
                    double nominator = (b_*g_-a_*h_)*f1_-c_*g_*f2_+(0.5*b_*h_-a_*g_)*f3_+0.25*b_*h_*f4_+0.5*c_*h_*f5_;
                    double denominator = b_*f1_-c_*f2_-a_*f3_>0?b_*f1_-c_*f2_-a_*f3_:0.0001;
                    return nominator/denominator;
                }

                /**\brief compute sun rise and sun set values
                 * \param delta,[rad] - earrh declination angle
                 * \param phi,[rad] -- latitude
                 * \param slope,[rad]
                 * \param aspect,[rad]
                 * calculates  local variables omega1_24_, omega2_24_, omega1_24b_, omega2_24b_*/
                void compute_sun_rise_set(double delta, double phi, double slope, double aspect) {
                    ///TODO get info about hemisphere from data, don't see anything in geopoint, but it is inside netcdf file -- add geopoint_with_crs to interface
                    double omega_s; // omega_s -- time of potential horizontal sunset, -omega_s --time of horizontal sunrize
                    // this solar noon and sunrise taken from ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.575
                    if (abs(phi - delta) >= pi / 2) { omega_s = pi; }
                    if (abs(phi - delta) < pi / 2 and (abs(phi + delta) >= pi / 2)) { omega_s = 0; }
                    if (abs(phi - delta) < pi / 2 and (abs(phi + delta) < pi / 2)) {
                        omega_s = acos(-tan(delta) * tan(phi));
                    } // for horizontal surface

                    compute_abc(delta, phi, slope, aspect);
                    double costt_sunset = costt(omega_s);
                    double costt_sunrise = costt(-omega_s);

                    /// TODO: verify this weird sun rise and sunset procedure

                    // Lower integration limit
                    double bbcc = b_ * b_ + c_ * c_ > 0.0 ? b_ * b_ + c_ * c_ : 0.0001;
                    double sqrt_bca = bbcc - a_ * a_ > 0.0 ? bbcc - a_ * a_ : 0.0; // the authors suggest here 0.0001???
                    double sin_omega1 = std::min(1.0, std::max(-1.0, (a_ * c_ - b_ * std::pow(sqrt_bca, 0.5)) / bbcc));//eq.(13a)
                    double omega1 = asin(sin_omega1);
                    omega1_24_ = omega1;
                    double costt_omega1 = costt(omega1);
                    if ((costt_sunrise <= costt_omega1) and (costt_omega1 < 0.001)) {omega1_24_ = omega1;}
                    else {
                        omega1 = -pi - omega1;
                        if (costt(omega1) > 0.001) {omega1_24_ = -omega_s;}
                        else {
                            if (omega1 <= -omega_s) {omega1_24_ = -omega_s;}
                            else {omega1_24_ = omega1;}
                        }
                    }
                    omega1_24_ = std::max(-omega_s,omega1_24_);

                    // Upper integration limit
                    double sin_omega2 = std::min(1.0, std::max(-1.0, (a_ * c_ + b_ * std::pow(sqrt_bca, 0.5)) / bbcc));//eq.(13b)
                    double omega2 = asin(sin_omega2);
                    omega2_24_ = omega2;
                    double costt_omega2 = costt(omega2);
                    if (costt_sunset <= costt_omega2 and costt_omega2 < 0.001) {omega2_24_ = omega2;}
                    else {
                        omega2 = pi - omega2;
                        if (costt(omega2) > 0.001) {omega2_24_ = omega_s;}
                        else {
                            if (omega2 >= omega_s) {omega2_24_ = omega_s;}
                            else {omega2_24_ = omega2;}
                        }
                    }
                    omega2_24_ = std::min(omega_s,omega2_24_);

                    if (omega1_24_>omega2_24_){omega1_24_=omega2_24_;}// slope is always shaded

                    // two periods of direct beam radiation (eq.7)
                    if (sin(slope) > sin(phi) * cos(delta) + cos(phi) * sin(delta)) {
                        double sinA = std::min(1.0, std::max(-1.0, (a_ * c_ + b_ * std::pow(sqrt_bca, 0.5)) / bbcc));
                        double A = asin(sinA);
                        double sinB = std::min(1.0, std::max(-1.0, (a_ * c_ + b_ * std::pow(sqrt_bca, 0.5)) / bbcc));
                        double B = std::asin(sinB);
                        omega2_24b_ = std::min(A, B);
                        omega1_24b_ = std::max(A, B);
                        compute_abc(delta_, phi_, slope_, aspect_);
                        double costt_omega2_24b = costt(omega2_24b_);
                        if (costt_omega2_24b < -0.001 or costt_omega2_24b > 0.001) { omega2_24b_ = -pi - omega2_24b_; }
                        double costt_omega1_24b = costt(omega1_24b_);
                        if ((costt_omega1_24b < -0.001 or costt_omega1_24b > 0.001)) { omega1_24b_ = pi - omega1_24b_; }
                        if ((omega2_24b_ > omega1_24_) or (omega1_24b_ < omega2_24_)) {
                            omega2_24b_ = omega1_24_;
                            omega1_24b_ = omega1_24_;
                        } // single period of sun
                    }
                }

                /**\brief computes solar hour angle from local time
                 * ref.: https://en.wikipedia.org/wiki/Equation_of_time
                 * ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574, no EOT correction provided
                 * \param lt -- local time (LT) [h]
                 * \param longitute = (0 for UTC time)
                 * \param tz = 0 -- time zone [h], difference of LT from UTC
                 * we use UTC, so longitude = 0.0, tz = 0.0
                 * \return HRA, [rad] -- Earth hour angle*/
                double hour_angle(double lt, double longitude = 0.0, double tz = 0.0) {
                    double LSTM = 15 * tz; // local standard time meridian
//                double B = 360/365*(doy)*(-81);
//                double EoT = 9.87*sin(2*B) - 7.3*cos(B) - 1.5*sin(B);// equation of time
                    double M = 2 * pi / 365.2596358 * doy_ * pi / 180;
                    double EoT = -7.659 * sin(M) +
                                 9.863 * sin(2 * M + 3.5932);//https://en.wikipedia.org/wiki/Equation_of_time*/
                    double TC = 4 * (longitude - LSTM) +
                                EoT; //time correction factor including EoT correction for Earth eccentricity
                    double LST = lt - TC / 60; // local solar time
                    return 15 * (lt - 12) *
                           deg2rad; ///TODO: find right EOT and data for validation, so it should return: 15*(LST-12)
                }

                /**\brief computes earth declination angle
                 * // ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.574, eq.(D.5)
                 * \param doy -- day of the year
                 * \return declination, [rad] */
                double compute_earth_declination(double doy) {
                    double G = 2 * pi / 365.0 * (doy - 1);
                    return 0.006918 - 0.399912 * cos(G) + 0.070257 * sin(G) - 0.006758 * cos(2 * G) +
                           0.000907 * sin(2 * G) - 0.002697 * cos(3 * G) + 0.00148 * sin(3 * G);
                }

                /**\brief computes extraterrestrial solar radiation
                 * \param cos(theta)
                 * \param doy -- day of the year
                 * return ra, [W/m^2]*/
                double compute_ra(double cos_theta, double doy) {
                    return gsc * cos_theta * (1 + 0.0033 * cos(doy * 2 * pi / 365)); // eq.(1)
                }

                double compute_ra_step(double cos_theta, double doy) {
                    return gsc * cos_theta * (1 + 0.0033 * cos(doy * 2 * pi / 365))/pi/2; // eq.(1)
                }

                double fi() { return 0.75 + 0.25 * cos(slope_) - 0.5 / pi * slope_;/*eq.(32)*/}

                double fia(double kb, double kd) {
                    return (1 - kb) *
                           (1 + pow(kb / (abs(kb + kd) > 0.0 ? (kb + kd) : 0.0001*(kb>=0?1:-1)), 0.5) * pow(sin(slope_ / 2), 3)) *
                           fi() + fb_ * kb; /*eq.(33)*/
                }

                /**\brief compute 24h parameter*/
                double costt_step(double omega1, double omega2) {
                    return -a_*(omega2-omega1) + b_ * (sin(omega2)-sin(omega1))+c_*(cos(omega2)-cos(omega1));
                }


                /** \brief computes instantaneous predicted short-wave clear-sky radiation (direct, diffuse, reflected) for inclined surfaces
                 * ref.: Allen, R. G.; Trezza, R. & Tasumi, M. Analytical integrated functions for daily solar radiation on slopes Agricultural and Forest Meteorology, 2006, 139, 55-73
                 * \param latitude, [deg]
                 * \param utctime,
                 * \param slope, [deg]
                 * \param aspect, [deg]
                 * \param temperature, [degC]
                 * \param rhumidity, [percent]
                 * \param elevation
                 * \return */
                double psw_radiation(double latitude, utctime t, double slope=0.0, double aspect = 0.0,
                                   double temperature = 0.0, double rhumidity = 40.0, double elevation = 0.0) {
                    doy_ = utc.day_of_year(t);
                    lt_ = utc.calendar_units(t).hour + utc.calendar_units(t).minute / 60.0;
                    delta_ = compute_earth_declination(doy_);
                    omega_ = hour_angle(lt_); // earth hour angle

                    slope_ = slope*pi/180.0;
                    aspect_ = aspect*pi/180;
                    phi_ = latitude * pi / 180;

                    compute_sun_rise_set(delta_, phi_, 0.0, 0.0); // for horizontal surface
                    compute_abc(delta_, phi_, 0.0, 0.0);
                    costthor_ = costt(omega_);

                    if (omega_ >= omega1_24_ and omega_ <= omega2_24_) {
                        rahor_ = std::max(0.0, compute_ra(costthor_, doy_)); // eq.(1) with cos(theta)hor
                    } else {
                        rahor_ = 0.0;
                    };


                    compute_sun_rise_set(delta_, phi_, slope, aspect);
                    compute_abc(delta_, phi_, slope, aspect);
                    costt_ = costt(omega_); // eq.(14)

                    if (omega_ >= omega1_24_ and omega_ <= omega2_24_) {
                        ra_ = std::max(0.0, compute_ra(costt_, doy_)); // eq.(1)
                    } else {
                        ra_ = 0.0;
                    };


                    double W; //equivalent depth of precipitable water in the atmosphere[mm]
                    double eatm = atm_pressure(
                            elevation); // [kPa] atmospheric pressure as a function of elevation ///TODO: get elevation from cell.midpoint().z
                    double ea = actual_vp(temperature, rhumidity); //[kPa] actual vapor pressure
                    W = 0.14 * ea * eatm + 2.1; // eq.(18)

                    double Kbo;
                    double sin_beta, sin_betahor;
                    sin_betahor = costthor_; // eq.(20) equal to (4), cos_tthor = sin_betahor /// TODO: check if cos_tt = sin_beta is true for inclined surface
                    sin_beta = costt_;
                    // clearness index for direct beam radiation

                    Kbo = std::min(1.0, std::max(-0.4, 0.98 * exp(-0.00146 * eatm / param.turbidity / sin_beta -
                                                        0.075 * pow((W / sin_beta), 0.4)))); // eq.(17)
                    double Kbohor = std::min(1.0, std::max(-0.4, 0.98 * exp(-0.00146 * eatm / param.turbidity / sin_betahor -
                                                                  0.075 * pow((W / sin_betahor), 0.4)))); // eq.(17)

                    double Kdo; // transmissivity of diffuse radiation, eq.(19)a,b,c
                    if (Kbo >= 0.15) { Kdo = 0.35 - 0.36 * Kbo; }
                    else if (Kbo < 0.15 and Kbo > 0.065) { Kdo = 0.18 + 0.82 * Kbo; }
                    else { Kdo = 0.10 + 2.08 * Kbo; }

                    double Kdohor;
                    if (Kbohor >= 0.15) { Kdohor = 0.35 - 0.36 * Kbohor; }
                    else if (Kbohor < 0.15 and Kbohor > 0.065) { Kdohor = 0.18 + 0.82 * Kbohor; }
                    else { Kdohor = 0.10 + 2.08 * Kbohor; }

                    double fi_ = fi();//eq.(32)

                    // fb_ = min(5.0,Kbo/Kbohor*ra_/(rahor_>0.0?rahor_:max(0.00001,ra_)));//eq.(34)
                    fb_ = ra_ / (rahor_ > 0.0 ? rahor_ : std::max(0.00001, ra_));//eq.(34)

                    double fia_ = fia(Kbohor, Kdohor); //eq.(33)


                    double dir_radiation = Kbo * ra_;
                    double dif_radiation = fia_ * Kdo * rahor_;
                    double ref_radiation = param.albedo * (1 - fi_) * (Kbo + Kdo) * rahor_;

                    return dir_radiation + dif_radiation + ref_radiation; // predicted clear sky solar radiation for inclined surface [W/m2]
                }

                /**\brief computes extraterrestrial radiation*/
                double compute_ra(double step,double delta, double phi, double s, double gamma, double omega1, double omega2, double doy){

                    // NEXT Calculations is done to define the integration omegalim properly, based on the step
                    compute_sun_rise_set(delta, phi, s, gamma);

                    // 24-h step: integrate over full period
                    if (step>=23) {
                        omega1 = std::max(omega1, omega1_24_);
                        omega2 = std::min(omega2, omega2_24_);
                    }

                    // 1-h or 3-h step, if not the whle period has sun
                    if (omega2>omega1_24_){
                        if ((omega1<omega1_24_) and (omega2<=omega2_24_)){
                            omega1 = omega1_24_;
                        }
                    }
                    if (omega1<omega2_24_) {
                        if ((omega2 > omega2_24_) and (omega1 >= omega1_24_)) {
                            omega2 = omega2_24_;
                        }
                    }
                    update_omega12(omega1,omega2);

                    double costt = costt_step(omega1,omega2);//-a_*(omega2-omega1) + b_ * (sin(omega2)-sin(omega1))+c_*(cos(omega2)-cos(omega1));

                    if (omega1 >= omega1_24_ and omega2 <= omega2_24_) {
                        return std::max(0.0, gsc * costt * (1 + 0.0033 * cos(doy * 2 * pi / 365))/pi/2); // eq.(1) with cos(theta)hor
                    } else {
                        return 0.0;
                    };
                }
                /** \brief computes instantaneous predicted short-wave clear-sky radiation (direct, diffuse, reflected) for inclined surfaces
                * ref.: Allen, R. G.; Trezza, R. & Tasumi, M. Analytical integrated functions for daily solar radiation on slopes Agricultural and Forest Meteorology, 2006, 139, 55-73
                * \param latitude, [deg]
                * \param utctime tstart,
                * \param ucttime tend,
                * \param slope, [deg]
                * \param aspect, [deg]
                * \param temperature, [degC]
                * \param rhumidity, [percent]
                * \param elevation
                * \return */
                double psw_radiation_step(double latitude, utctime t1, utctime t2, double slope=0.0, double aspect = 0.0,
                                     double temperature = 0.0, double rhumidity = 40.0, double elevation = 0.0) {
                    double doy = utc.day_of_year(t1);
                    double lt1 = utc.calendar_units(t1).hour + utc.calendar_units(t1).minute / 60.0;
                    double lt2 = utc.calendar_units(t2).hour + utc.calendar_units(t2).minute / 60.0;
                    double delta = compute_earth_declination(doy);
                    double omega1 = hour_angle(lt1); // earth hour angle
                    double omega2 = hour_angle(lt2); // earth hour angle

                    slope_= slope*pi/180.0;
                    aspect_ = aspect*pi/180;
                    double phi = latitude * pi / 180;


                    compute_abc(delta,phi,0.0,0.0);
                    rahor_ = compute_ra(lt2-lt1,delta,phi,0.0,0.0,omega1,omega2,doy);
                    compute_fs(omega1_,omega2_);
                    double sin_betahor = std::min(1.0,std::max(-1.0,compute_beta_step())); // eq.(22),

                    compute_abc(delta,phi,slope_,aspect_);
                    ra_ = compute_ra(lt2-lt1,delta,phi,slope_,aspect_,omega1,omega2,doy);
                    compute_fs(omega1_,omega2_);
                    double sin_beta = std::min(1.0,std::max(-1.0,compute_beta_step())); // eq.(22),


                    double W; //equivalent depth of precipitable water in the atmosphere[mm]
                    double eatm = atm_pressure(elevation); // [kPa] atmospheric pressure as a function of elevation
                    double ea = actual_vp(temperature, rhumidity); //[kPa] actual vapor pressure
                    W = 0.14 * ea * eatm + 2.1; // eq.(18)

                    double Kbo = std::min(1.0, std::max(-0.4, 0.98 * exp(-0.00146 * eatm / param.turbidity / sin_beta -
                                                                  0.075 * pow((W / sin_beta), 0.4)))); // eq.(17)
                    double Kbohor = std::min(1.0, std::max(-0.4, 0.98 * exp(-0.00146 * eatm / param.turbidity / sin_betahor -
                                                                            0.075 * pow((W / sin_betahor), 0.4)))); // eq.(17)

                    double Kdo; // transmissivity of diffuse radiation, eq.(19)a,b,c
                    if (Kbo >= 0.15) { Kdo = 0.35 - 0.36 * Kbo; }
                    else if (Kbo < 0.15 and Kbo > 0.065) { Kdo = 0.18 + 0.82 * Kbo; }
                    else { Kdo = 0.10 + 2.08 * Kbo; }

                    double Kdohor;
                    if (Kbohor >= 0.15) { Kdohor = 0.35 - 0.36 * Kbohor; }
                    else if (Kbohor < 0.15 and Kbohor > 0.065) { Kdohor = 0.18 + 0.82 * Kbohor; }
                    else { Kdohor = 0.10 + 2.08 * Kbohor; }

                    double fi_ = fi();//eq.(32)

                    // fb_ = min(5.0,Kbo/Kbohor*ra_/(rahor_>0.0?rahor_:max(0.00001,ra_)));//eq.(34)
                    fb_ = ra_ / (rahor_ > 0.0 ? rahor_ : std::max(0.00001, ra_));//eq.(34)
                    double fia_ = fia(Kbohor, Kdohor); //eq.(33)

                    double dir_radiation = Kbo * ra_;
                    double dif_radiation = fia_ * Kdo * rahor_;
                    double ref_radiation = param.albedo * (1 - fi_) * (Kbo + Kdo) * rahor_;
                    return dir_radiation + dif_radiation + ref_radiation; // predicted clear sky solar radiation for inclined surface [W/m2]
                }
                /**\brief translates measured solar radiation from horizontal surfaces to slopes
                 * ref.: Allen, R. G.; Trezza, R. & Tasumi, M. Analytical integrated functions for daily solar radiation on slopes Agricultural and Forest Meteorology, 2006, 139, 55-73
                 * \param latitude, [deg]
                 * \param utctime,
                 * \param slope, [deg]
                 * \param aspect, [deg]
                 * \param temperature, [degC]
                 * \param rhumidity, [percent]
                 * \param elevation, [m]
                 * \param rsm,[W/m^2] -- measured solar radiation
                 * \return */
                double tsw_radiation(double latitude, utctime t, double slope = 0.0, double aspect = 0.0,
                                   double temperature = 0.0, double rhumidity = 40.0, double elevation = 0.0,
                                   double rsm = 0.0) {
                    // first calculate all predicted values
                    double tsw_rad = psw_radiation(latitude, t, slope, aspect, temperature, rhumidity, elevation);
                    double tauswhor = rsm > 0.0 ? rsm / (rahor_ > 0.0 ? rahor_ : rsm)
                                                : 1.0; //? not sure if we use a theoretical rahor here
                    double KBhor;
                    if (tauswhor >= 0.42) { KBhor = 1.56 * tauswhor - 0.55; }
                    else if (tauswhor > 0.175 and tauswhor < 0.42) {
                        KBhor = 0.022 - 0.280 * tauswhor + 0.828 * tauswhor * tauswhor + 0.765 * pow(tauswhor, 3);
                    }
                    else { KBhor = 0.016 * tauswhor; }
                    double KDhor = tauswhor - KBhor;
                    if (rsm>0.0)
                        tsw_rad = rsm * (fb_ * KBhor / tauswhor + fia(KBhor, KDhor) * KDhor / tauswhor +
                                                    param.albedo * (1 - fi()));
                    return tsw_rad;
                }
                /**\brief clear-sky longwave raditiation
                 * ref.: Lawrence Dingman Physical Hydrology, Third Edition, 2015, p.261
                 * \param temperature, [degC] -- air temperature
                 * \param rhumidity, [persent] -- relative humidity
                 * \param ss_temp, [K] -- surface temperature
                 * response.lw_radiation W/m^2 */
                /// TODO https://www.hydrol-earth-syst-sci.net/17/1331/2013/hess-17-1331-2013-supplement.pdf
                // TODO discuss the option to have different formulations here.
                double lw_radiation(double temperature, double rhumidity){
                    double Lin = 2.7*actual_vp(temperature,rhumidity)+0.245*(temperature+273.15)-45.14;
                    double ss_temp = std::min(temperature+273.15-2.5,273.16);
                    double epsilon_ss = 0.95;//water TODO: as parameter
                    double Lout = epsilon_ss*sigma*pow(ss_temp,4)+(1-epsilon_ss)*Lin;
                    return (Lin-Lout)*MJm2d2Wm2;
                }
            };


        }

    }
}