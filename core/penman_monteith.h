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

#include <chrono>
#include <boost/math/constants/constants.hpp>
#include "core/hydro_functions.h"
namespace shyft::core{
   // namespace evapotranspiration{
        /**\brief Evapotranspiration models
         * Penman-Monteith, ref.: ASCE-EWRI The ASCE Standardized Reference Evapotranspiration Equation Environmental and
         * Water Resources Institute (EWRI) of the American Society of Civil Engineers Task Com- mittee on Standardization of
         * Reference Evapotranspiration Calculation, ASCE, Washington, DC, Environmental and Water Resources Institute (EWRI) of
         * the American Society of Civil Engineers Task Com- mittee on Standardization of Reference Evapotranspiration Calculation,
         * ASCE, Washington, DC, 2005*/
        namespace penman_monteith {
            using namespace std;
            using namespace shyft::core::hydro_functions;
            /**ASCE-EWRI hourly based model
             * it is stated taht this model can be used with the smaller periods as well, but we need to account for fractional hour
             * TODO: discuss how to manipulate between models and timesteps*/
            struct parameter{
                double lai = 2.0; //leaf area index, defaulted value taken from MORECS-metoff 1981 model for grass
                double hveg = 0.15; // vegetation height, [m] (grass from MORECS)
                double rl = 50.0; // effective stomatal resistance, [s/m] (calculated for grass from MORECS data)
                parameter(double hveg=0.12, double lai = 2.0,double rl=50.0) : hveg(hveg), lai(lai), rl(rl) {}
            };
            struct response {
                double et_ref = 0.0; // reference evapotranspiration [mm/h]
            };

            template<class P, class R>
            struct calculator {
                P param;

                explicit calculator(const P &p) : param(p) {}

                // manipulating with time-steps
//                 utctimespan delta = time_axis.delta();
//                std::cout<< to_seconds(delta)  <<"  ===== "<<std::endl;

                /**\brief computes ASCE PENMAN-MONTEITH reference evapotranspiation, [mm/h] based on eq. B.1 from ASCE-EWRI
                 * \param response -- updating response structure
                 * \param net_radiation, [MJ/m2/h] -- net radiation
                 * \param temperature [degC] -- air temperature
                 * \param rhumidity [-] -- relative humidity
                 * \param elevation [m]
                 * \param height_ws [m] -- wind speed measurement height
                 * \param windspeed [m/s at height_ws]
                 * \param height_t [m] -- height of humidity,temperature measuerements*/
                void reference_evapotranspiration_asce(R& response, double net_radiation, double temperature, double rhumidity, double elevation=0.1, double height_ws=2.0, double windspeed=0.1, double height_t=2.0){

                    double pressure = atm_pressure(elevation, 293.15);// recommended value for T0 during growing season is 20 gradC, see eq. B.8
                    double lambda_rho = vaporization_latent_heat(temperature)*rho_w;
//                     std::cout<<"lambda_rho: "<<lambda_rho<<std::endl;
                    double delta = svp_slope(temperature);
//                     std::cout<<"delta: "<<delta<<std::endl;
                    double avp = actual_vp(temperature, rhumidity);
//                     std::cout<<"avp: "<<avp<<std::endl;
                    double G = soil_heat_flux(net_radiation);
//                     std::cout<<"G: "<<G<<std::endl;
                    double rho = density_air(pressure, temperature, avp);
//                     std::cout<<"rho: "<<rho<<std::endl;
                    double ra = resistance_aerodynamic(height_ws,windspeed);
//                     std::cout<<"ra: "<<ra<<std::endl;
                    double sat_vp = svp(temperature);
//                     std::cout<<"sat_vp: "<<sat_vp<<std::endl;
                    double nominator = delta * (net_radiation + soil_heat_flux(net_radiation)) +
                                       ktime*density_air(pressure, temperature, avp)*cp/resistance_aerodynamic(height_ws,windspeed, height_t)*(svp(temperature)-avp);
//                     std::cout<<"nominator: "<<nominator<<std::endl;
                    double denominator = delta + gamma_pm(pressure, temperature)*(1+resistance_surface(param.hveg)/resistance_aerodynamic(height_ws,windspeed));
//                     std::cout<<"denominator: "<<denominator<<std::endl;
                    response.et_ref = nominator/denominator/lambda_rho;
                    return;
                };
                /**\brief computes ASCE PENMAN-MONTEITH reference evapotranspiation, [mm/h] Standartised eq.1 from ASCE-EWRI
                * \param response -- updating response structure
                * \param net_radiation, [MJ/m2/h] -- net radiation
                * \param temperature [degC] -- air temperature
                * \param rhumidity [-] -- relative humidity
                * \param elevation [m]
                * \param height_ws [m] -- wind speed measurement height
                * \param windspeed [m/s at height_ws]
                * \param height_t [m] -- height of humidity,temperature measuerements*/
                void reference_evapotranspiration_st(R& response, double net_radiation, double temperature, double rhumidity, double elevation=0.1, double height_ws=2.0, double windspeed=0.1, double height_t=2.0){

//                     double pressure = atm_pressure(elevation, 293.15);// recommended value for T0 during growing season is 20 gradC, see eq. B.8
                    double pressure = 101.3*pow((293.0-0.0065*elevation)/293,5.26); // eq.34
//                     double gamma = gamma_pm(pressure, temperature);
                    double gamma = 0.000665*pressure;
//                     double lambda_rho = vaporization_latent_heat(temperature)*rho_w;
//                     std::cout<<"lambda_rho: "<<lambda_rho<<std::endl;
//                     double delta = svp_slope(temperature);
                    double delta = 2593.0*exp(17.27*temperature/(temperature+237.3))/(temperature+237.3)/(temperature+237.3);
                    std::cout<<"delta: "<<delta<<std::endl;
                    double avp = actual_vp(temperature, rhumidity);
                    std::cout<<"avp: "<<avp<<std::endl;
//                     double G = soil_heat_flux(net_radiation);

                    double G = 0.1*net_radiation;
                    if (net_radiation<0.0)
                        G= 0.5*net_radiation;


                    std::cout<<"G: "<<G<<std::endl;
                    double rho = density_air(pressure, temperature, avp);
//                     std::cout<<"rho: "<<rho<<std::endl;
                    double ra = resistance_aerodynamic(height_ws,windspeed);
//                     std::cout<<"ra: "<<ra<<std::endl;
                    double sat_vp = svp(temperature);
                    std::cout<<"sat_vp: "<<sat_vp<<std::endl;
                    double Cn = 37.0;
                    double Cd = 0.24;
                    if (net_radiation<=0.0){
                        Cd = 0.96;
                    }

                    double nominator = 0.408* delta * (net_radiation + soil_heat_flux(net_radiation)) +
                                       gamma*Cn*ws_adjustment(height_ws,windspeed)*(sat_vp-avp)/(temperature+273);
//                     std::cout<<"nominator: "<<nominator<<std::endl;
                    double denominator = delta + gamma*(1+Cd*ws_adjustment(height_ws,windspeed));
//                     std::cout<<"denominator: "<<denominator<<std::endl;
                    response.et_ref = max(0.01,nominator/denominator);
                    return;
                };

            private:
                // ============================================ //
                // constants for asce-ewri standard equation
                /**Values for Cn and Cd in Eq. 1 ref.: ASCE-EWRI p.21*/
                const double Cni[2][3] = {{900.0,37.0,37.0}, {1600.0,66.0,66.0}};// Cn short reference: daily -- 900; hourly day -- 37 hourly night -- 37
                // CN tall reference: daily 1600, hourly day 66; houry night
                const double Cdi[2][3] = {{0.34, 0.24, 0.96}, {0.38, 0.25, 1.7}}; // Cd short reference: daily -- 0.34; hourly 0.24; hourly night 0.96
                // Cd tall reference: dai;y 0.38, houry day 0.25, hiurly night 1.7
                /**ASCE Penman-Monteith Terms Standardized for application of the RET equation*/
                const double zveg[2] = {0.12,0.5}; // vegetation height short and tall reference
                const double zh[2] = {2.0,2.0}; //[m], height of air temperature and humidity measurements 1.5-2.5 m
                const double zw[2] = {2.0,2.0};//[m], height corresponding to wind speed
                const double zo[2] = {0.08,0.08}; //[m], zero plane displacement height
                const double lambda[2] = {2.45,2.45}; //[MJ/kg], latent heat of vaorization
                const double rs_dailly[2] = {70.0, 45.0}; //[s/m], surface resistance, daily, short and tall reference
                const double rs_h_day[2] = {50.0,30.0}; //[s/m], surface resistance, hourly during daytime, short and tall
                const double rs_h_night[2] = {200.0,200.0};//[s/m], surface resistance, hourly during nighttime, short and tall
                // predicting daytime Rnet > 0; nighttime Rnet <= 0.0
                // ============================================ //
                // constants for asce penman-monteith method
                const double rho_w = 1;//[Mg/m^3]
                const double cp = 1.013*0.001; // specific heat of moist air, [MJ/kg/gradC]
                const double ktime = 3600; // unit conversion for ET in mm/h
                //double h = 0.1; // mean height of vegetation --> parameter
                double d =  0.67 * param.hveg; // zero displacement height, [m]
                double zom = 0.123*param.hveg; // roughness length governing momentum transfer, [m]
                double zoh =  0.0123*param.hveg;// roughness height for transfer of heat and vapor, [m]
                double kappa = 0.41; // von Karman's constant

                /**\brief aerodynamic resistance, [s/m], eq. B.2
                 * \param hws -- wind speed measurements height, [m]
                 * \param ws -- wind speed at height hws, [m/s]
                 * \param ht [m] -- height of humidity and or temperature measurements*/
                double resistance_aerodynamic(double hws, double ws, double ht=2.0){
                    return log((hws - d)/zom)*log((ht-d)/zoh)/kappa/kappa/min(ws,0.01);
                }

                /**\brief active leaf-area index (LAI), eq.B.4
                 * \param hveg -- vegetation height, [m]*/
                double lai_a(double hveg){
                    double lai = 24.0*hveg; // for clipped grass
                    lai = 5.5 + 1.5*log(hveg);// for alfalfa /// TODO find condition to switch
                    return 0.5*lai;
                }
                /**\brief bulk surface resistance, [s/m], eq. B.3
                 * \param zveg -- vegetation height,[m]*/
                double resistance_surface(double hveg=0.1){ /// TODO discuss the incorporation of MORECS model here, is it worth to account for seasons?
                    return param.rl/lai_a(hveg);
                }

                /// TODO: move all to the hydro_functions???

                /**\brief Latent heat of vaporization, eq. B.7; ASCE-EWRI
                 * \param temperature [degC]*/
                double vaporization_latent_heat(double temperature = 20.0){
                    return 2.501 - (2.361*0.001)*temperature;//with default value gives 2.45;
                }
//                 /**\brief atmospheric pressure, ASCE-EWRI, eq.B.8
//                  * */
//                 double pressure_air(double elevation){
//                     const double p0 = 101.325;// pressure at sea level, [kPa]
//                     const double T0 = 293.15; // temperature at sea level for growing season, [K]
//                     const double alpha1 = 0.0065; // lapse rate, [K/m]
//                     const double z0 = 0; // sea level, [m]
//                     const double Rgas = 287; // universal gas constant, [J/kg/K]
//                     const double g = 9.807; // gravitational acceleration, [m/s2]
//                     return p0*pow((T0-alpha1*(elevation-z0))/T0,g/alpha1/Rgas);
//                 }

                /**\brief atmospheric air density, ASCE-EWRI, eq. B.10
                 * \param pressure [kPa]
                 * \param temperature [degC]
                 * \param actual vapor pressure, [kPa]*/
                double density_air(double pressure, double temperature, double ea){
                    double Tk = 273.16+temperature;
                    double Tkv = Tk/(1-0.378*ea/pressure); // B.11
                    return 3.486*pressure/Tkv; // B.10
                }

                /**\brief psychrometric constant, [kPa/gradC]
                 * \param atmospheric pressure, [kPa]
                 * \param air temperature, [degC]*/
                double gamma_pm(double pressure, double temperature=20.0){

                    double lambda = vaporization_latent_heat(temperature);
                    const double epsilon = 0.622; // ratio of molecular weight of water vapor/dry air
                    return cp*pressure/epsilon/lambda;
                }

                /**\brief soil heat flux density (G) for hourly periods, eq.B.13
                 * [net_radiation] should be either MJ/m2/h or Wt/m2*/
                double soil_heat_flux(double net_radiation){
                    double kg = 0.4;
                    if (net_radiation <= 0.0)
                        kg = 2.0;
                    return kg*exp(-0.5*lai_a(param.hveg))*net_radiation;
                }

                /**\brief wind speed adjustment for measurement height, eq.B.14*/
                double ws_adjustment(double height_measure, double ws_measure){
                    return ws_measure*log((2-d)/zom)/log((height_measure-d)/zom);
                }
            };


        }
   // } //namespace evapotranspiration

}