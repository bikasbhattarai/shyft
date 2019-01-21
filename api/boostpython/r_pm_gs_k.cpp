/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#include "boostpython_pch.h"

#include "core/utctime_utilities.h"
#include "core/radiation.h"
#include "core/penman_monteith.h"
#include "core/actual_evapotranspiration.h"
#include "core/precipitation_correction.h"
#include "core/gamma_snow.h"
#include "core/kirchner.h"
#include "core/r_pm_gs_k.h"
#include "api/api.h"
#include "core/r_pm_gs_k_cell_model.h"
#include "core/region_model.h"
#include "core/model_calibration.h"
#include "expose_statistics.h"
#include "expose.h"

static char const* version() {
   return "v1.0";
}

namespace expose {
    namespace r_pm_gs_k {
        using namespace boost::python;
        namespace py = boost::python;
        using namespace shyft::core;
        using namespace shyft::core::r_pm_gs_k;
        using std::string;
        using std::vector;

        typedef vector<state> RPMGSKStateVector;

        static void
        parameter_state_response() {

            class_<parameter,bases<>,std::shared_ptr<parameter>>("RPMGSKParameter",
                              "Contains the parameters to the methods used in the RPMGSK assembly\n"
                              "radiation, penman_monteith, gamma_snow,actual_evapotranspiration,precipitation_correction,kirchner\n"
                )
                .def(init<radiation::parameter,penman_monteith::parameter,gamma_snow::parameter,actual_evapotranspiration::parameter,kirchner::parameter,precipitation_correction::parameter,py::optional<glacier_melt::parameter,routing::uhg_parameter,mstack_parameter>>(args("rad","pm","gs","ae","k","p_corr","gm","routing","msp"),"create object with specified parameters"))
                .def(init<const parameter&>(args("p"),"clone a parameter"))
                .def_readwrite("rad",&parameter::rad,"radiation parameter")
                .def_readwrite("pm",&parameter::pm,"penman_monteith parameter")
                .def_readwrite("gs",&parameter::gs,"gamma-snow parameter")
                .def_readwrite("gm", &parameter::gm, "glacier melt parameter")
				.def_readwrite("ae",&parameter::ae,"actual evapotranspiration parameter")
                .def_readwrite("kirchner",&parameter::kirchner,"kirchner parameter")
                .def_readwrite("p_corr",&parameter::p_corr,"precipitation correction parameter")
                .def_readwrite("routing",&parameter::routing,"routing cell-to-river catchment specific parameters")
                .def_readwrite("msp",&parameter::msp,"contains the method stack parameters")
                .def("size",&parameter::size,"returns total number of calibration parameters")
                .def("set",&parameter::set,args("p"),"set parameters from vector/list of float, ordered as by get_name(i)")
                .def("get",&parameter::get,args("i"),"return the value of the i'th parameter, name given by .get_name(i)")
                .def("get_name",&parameter::get_name,args("i"),"returns the i'th parameter name, see also .get()/.set() and .size()")
                ;

            typedef std::map<int,parameter> RPMGSKParameterMap;
            class_<RPMGSKParameterMap>("RPMGSKParameterMap","dict (int,parameter)  where the int is the catchment_id")
                .def(map_indexing_suite<RPMGSKParameterMap>())
            ;

            class_<state>("RPMGSKState")
                .def(init<gamma_snow::state,kirchner::state>(args("gs","k"),"initializes state with gamma-snow gs and kirchner k"))
                .def_readwrite("gs",&state::gs,"gamma-snow state")
                .def_readwrite("kirchner",&state::kirchner,"kirchner state")
                ;


            class_<RPMGSKStateVector,bases<>,std::shared_ptr<RPMGSKStateVector> >("RPMGSKStateVector")
                .def(vector_indexing_suite<RPMGSKStateVector>())
                ;


            class_<response>("RPMGSKResponse","This struct contains the responses of the methods used in the RPMGSK assembly")
                .def_readwrite("rad",&response::rad,"radiation response")
                .def_readwrite("pm",&response::pm,"penman_monteith response")
                .def_readwrite("gs",&response::gs,"gamma-snnow response")
                .def_readwrite("gm_melt_m3s", &response::gm_melt_m3s, "glacier melt response[m3s]")
                .def_readwrite("ae",&response::ae,"actual evapotranspiration response")
                .def_readwrite("kirchner",&response::kirchner,"kirchner response")
                .def_readwrite("total_discharge",&response::total_discharge,"total stack response")
                ;
        }

        static void
        collectors() {
            typedef shyft::core::r_pm_gs_k::all_response_collector RPMGSKAllCollector;
            class_<RPMGSKAllCollector>("RPMGSKAllCollector", "collect all cell response from a run")
                .def_readonly("destination_area",&RPMGSKAllCollector::destination_area,"a copy of cell area [m2]")
                .def_readonly("avg_discharge",&RPMGSKAllCollector::avg_discharge,"Kirchner Discharge given in [m^3/s] for the timestep")
                .def_readonly("snow_sca",&RPMGSKAllCollector::snow_sca," gamma snow covered area fraction, sca.. 0..1 - at the end of timestep (state)")
                .def_readonly("snow_swe",&RPMGSKAllCollector::snow_swe,"gamma snow swe, [mm] over the cell sca.. area, - at the end of timestep")
                .def_readonly("snow_outflow",&RPMGSKAllCollector::snow_outflow," gamma snow output [m^3/s] for the timestep")
                .def_readonly("glacier_melt", &RPMGSKAllCollector::glacier_melt, " glacier melt (outflow) [m3/s] for the timestep")
                .def_readonly("ae_output",&RPMGSKAllCollector::ae_output,"actual evap mm/h")
                .def_readonly("pe_output",&RPMGSKAllCollector::pe_output,"pot evap mm/h")
                .def_readonly("end_reponse",&RPMGSKAllCollector::end_reponse,"end_response, at the end of collected")
                .def_readonly("avg_charge",&RPMGSKAllCollector::charge_m3s,"average charge in [m^3/s]")
            ;

            typedef shyft::core::r_pm_gs_k::discharge_collector RPMGSKDischargeCollector;
            class_<RPMGSKDischargeCollector>("RPMGSKDischargeCollector", "collect all cell response from a run")
                .def_readonly("cell_area",&RPMGSKDischargeCollector::cell_area,"a copy of cell area [m2]")
                .def_readonly("avg_discharge",&RPMGSKDischargeCollector::avg_discharge,"Kirchner Discharge given in [m^3/s] for the timestep")
                .def_readonly("snow_sca",&RPMGSKDischargeCollector::snow_sca," gamma snow covered area fraction, sca.. 0..1 - at the end of timestep (state)")
                .def_readonly("snow_swe",&RPMGSKDischargeCollector::snow_swe,"gamma snow swe, [mm] over the cell sca.. area, - at the end of timestep")
                .def_readonly("end_reponse",&RPMGSKDischargeCollector::end_response,"end_response, at the end of collected")
                .def_readwrite("collect_snow",&RPMGSKDischargeCollector::collect_snow,"controls collection of snow routine")
                .def_readonly("avg_charge", &RPMGSKDischargeCollector::charge_m3s, "average charge in [m^3/s]")
                ;
            typedef shyft::core::r_pm_gs_k::null_collector RPMGSKNullCollector;
            class_<RPMGSKNullCollector>("RPMGSKNullCollector","collector that does not collect anything, useful during calibration to minimize memory&maximize speed")
                ;

            typedef shyft::core::r_pm_gs_k::state_collector RPMGSKStateCollector;
            class_<RPMGSKStateCollector>("RPMGSKStateCollector","collects state, if collect_state flag is set to true")
                .def_readwrite("collect_state",&RPMGSKStateCollector::collect_state,"if true, collect state, otherwise ignore (and the state of time-series are undefined/zero)")
                .def_readonly("kirchner_discharge",&RPMGSKStateCollector::kirchner_discharge,"Kirchner state instant Discharge given in m^3/s")
                .def_readonly("gs_albedo",&RPMGSKStateCollector::gs_albedo,"")
                .def_readonly("gs_lwc",&RPMGSKStateCollector::gs_lwc,"")
                .def_readonly("gs_surface_heat",&RPMGSKStateCollector::gs_surface_heat,"")
                .def_readonly("gs_alpha",&RPMGSKStateCollector::gs_alpha,"")
                .def_readonly("gs_sdc_melt_mean",&RPMGSKStateCollector::gs_sdc_melt_mean,"")
                .def_readonly("gs_acc_melt",&RPMGSKStateCollector::gs_acc_melt,"")
                .def_readonly("gs_iso_pot_energy",&RPMGSKStateCollector::gs_iso_pot_energy,"")
                .def_readonly("gs_temp_swe",&RPMGSKStateCollector::gs_temp_swe,"")
            ;

        }

        static void
        cells() {
              typedef shyft::core::cell<parameter, environment_t, state, state_collector, all_response_collector> RPMGSKCellAll;
              typedef shyft::core::cell<parameter, environment_t, state, null_collector, discharge_collector> RPMGSKCellOpt;
              expose::cell<RPMGSKCellAll>("RPMGSKCellAll","tbd: RPMGSKCellAll doc");
              expose::cell<RPMGSKCellOpt>("RPMGSKCellOpt","tbd: RPMGSKCellOpt doc");
              expose::statistics::gamma_snow<RPMGSKCellAll>("RPMGSKCell");//it only gives meaning to expose the *All collect cell-type
              expose::statistics::actual_evapotranspiration<RPMGSKCellAll>("RPMGSKCell");
              expose::statistics::penman_monteith<RPMGSKCellAll>("RPMGSKCell");
              //expose::statistics::radiation<RPMGSKCellAll>("RPMGSKCell");
              expose::statistics::kirchner<RPMGSKCellAll>("RPMGSKCell");
              expose::cell_state_etc<RPMGSKCellAll>("RPMGSK");// just one expose of state

        }

        static void
        models() {
            typedef shyft::core::region_model<r_pm_gs_k::cell_discharge_response_t, shyft::api::a_region_environment> RPMGSKOptModel;
            typedef shyft::core::region_model<r_pm_gs_k::cell_complete_response_t, shyft::api::a_region_environment> RPMGSKModel;
            expose::model<RPMGSKModel>("RPMGSKModel","RPMGSK");
            expose::model<RPMGSKOptModel>("RPMGSKOptModel","RPMGSK");
            def_clone_to_similar_model<RPMGSKModel, RPMGSKOptModel>("create_opt_model_clone");
            def_clone_to_similar_model<RPMGSKOptModel,RPMGSKModel>("create_full_model_clone");
        }


        static void
        model_calibrator() {
            expose::model_calibrator<shyft::core::region_model<r_pm_gs_k::cell_discharge_response_t,shyft::api::a_region_environment>>("RPMGSKOptimizer");
        }
    }
}


BOOST_PYTHON_MODULE(_r_pm_gs_k)
{

    boost::python::scope().attr("__doc__")="Shyft python api for the r_pm_gs_k model";
    boost::python::def("version", version);
	boost::python::docstring_options doc_options(true, true, false);// all except c++ signatures
    expose::r_pm_gs_k::parameter_state_response();
    expose::r_pm_gs_k::cells();
    expose::r_pm_gs_k::models();
    expose::r_pm_gs_k::collectors();
    expose::r_pm_gs_k::model_calibrator();
}
