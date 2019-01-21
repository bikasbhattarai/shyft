/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#pragma once

namespace expose {
    namespace statistics {
        typedef shyft::time_series::dd::apoint_ts rts_;
        typedef std::vector<double> vd_;
        typedef const std::vector<int64_t>& cids_;
        typedef size_t ix_;
        using shyft::core::stat_scope;
        using namespace boost::python;
        namespace py=boost::python;
        
        template<class cell>
        static void kirchner(const char *cell_name) {
            char state_name[200];sprintf(state_name,"%sKirchnerStateStatistics",cell_name);
            typedef typename shyft::api::kirchner_cell_state_statistics<cell>    sc_stat;

            rts_ (sc_stat::*discharge_ts)(cids_, stat_scope) const = &sc_stat::discharge;
            vd_  (sc_stat::*discharge_vd)(cids_, ix_, stat_scope) const =&sc_stat::discharge;
            class_<sc_stat>(state_name,"Kirchner response statistics",no_init)
                .def(init<std::shared_ptr<std::vector<cell>> >(args("cells"),"construct Kirchner cell response statistics object"))
                .def("discharge",discharge_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("discharge",discharge_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("discharge_value",&sc_stat::discharge_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum discharge[m3/s]  for cells matching catchments_ids at the i'th timestep")
            ;
        }

		template<class cell>
		static void hbv_soil(const char *cell_name) {
			char state_name[200]; sprintf(state_name, "%sHbvSoilStateStatistics", cell_name);
			typedef typename shyft::api::hbv_soil_cell_state_statistics<cell>    sc_stat;

			rts_(sc_stat::*discharge_ts)(cids_, stat_scope) const = &sc_stat::discharge;
			vd_(sc_stat::*discharge_vd)(cids_, ix_, stat_scope) const = &sc_stat::discharge;
			class_<sc_stat>(state_name, "HbvSoil response statistics", no_init)
				.def(init<std::shared_ptr<std::vector<cell>> >(args("cells"), "construct Kirchner cell response statistics object"))
				.def("discharge", discharge_ts, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
				.def("discharge", discharge_vd, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("discharge_value", &sc_stat::discharge_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum discharge[m3/s]  for cells matching catchments_ids at the i'th timestep")
				;
		}

		template<class cell>
		static void hbv_tank(const char *cell_name) {
			char state_name[200]; sprintf(state_name, "%sHbvTankStateStatistics", cell_name);
			typedef typename shyft::api::hbv_tank_cell_state_statistics<cell>    sc_stat;

			rts_(sc_stat::*discharge_ts)(cids_, stat_scope) const = &sc_stat::discharge;
			vd_(sc_stat::*discharge_vd)(cids_, ix_, stat_scope) const = &sc_stat::discharge;
			class_<sc_stat>(state_name, "HbvSoil response statistics", no_init)
				.def(init<std::shared_ptr<std::vector<cell>> >(args("cells"), "construct Kirchner cell response statistics object"))
				.def("discharge", discharge_ts, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
				.def("discharge", discharge_vd, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("discharge_value", &sc_stat::discharge_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum discharge[m3/s]  for cells matching catchments_ids at the i'th timestep")
				;
		}

        template <class cell>
        static void priestley_taylor(const char *cell_name) {
            char response_name[200];sprintf(response_name,"%sPriestleyTaylorResponseStatistics",cell_name);
            typedef typename shyft::api::priestley_taylor_cell_response_statistics<cell> rc_stat;

            rts_ (rc_stat::*output_ts)(cids_, stat_scope) const = &rc_stat::output;
            vd_  (rc_stat::*output_vd)(cids_, ix_, stat_scope) const =&rc_stat::output;
            class_<rc_stat>(response_name,"PriestleyTaylor response statistics",no_init)
                .def(init<std::shared_ptr<std::vector<cell>> >(args("cells"),"construct PriestleyTaylor cell response statistics object"))
                .def("output",output_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("output",output_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("output_value", &rc_stat::output_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns for cells matching catchments_ids at the i'th timestep")
            ;
        }

        template <class cell>
        static void penman_monteith(const char *cell_name) {
            char response_name[200];sprintf(response_name,"%sPenmanMonteithResponseStatistics",cell_name);
            typedef typename shyft::api::penman_monteith_cell_response_statistics<cell> rc_stat;

            rts_ (rc_stat::*output_ts)(cids_, stat_scope) const = &rc_stat::output;
            vd_  (rc_stat::*output_vd)(cids_, ix_, stat_scope) const =&rc_stat::output;
            class_<rc_stat>(response_name,"PenmanMonteith response statistics",no_init)
                    .def(init<std::shared_ptr<std::vector<cell>> >(args("cells"),"construct PenmanMonteith cell response statistics object"))
                    .def("output",output_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                    .def("output",output_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
                    .def("output_value", &rc_stat::output_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns for cells matching catchments_ids at the i'th timestep")
                    ;
        }


        template <class cell>
        static void actual_evapotranspiration(const char *cell_name) {
            char response_name[200];sprintf(response_name,"%sActualEvapotranspirationResponseStatistics",cell_name);
            typedef typename shyft::api::actual_evapotranspiration_cell_response_statistics<cell> rc_stat;

            rts_ (rc_stat::*output_ts)(cids_, stat_scope) const = &rc_stat::output;
            vd_  (rc_stat::*output_vd)(cids_, ix_, stat_scope) const =&rc_stat::output;
            rts_ (rc_stat::*pot_ratio_ts)(cids_, stat_scope) const = &rc_stat::pot_ratio;
            vd_  (rc_stat::*pot_ratio_vd)(cids_, ix_, stat_scope) const =&rc_stat::pot_ratio;
            class_<rc_stat>(response_name,"ActualEvapotranspiration response statistics",no_init)
                .def(init<std::shared_ptr<std::vector<cell>> >(args("cells"),"construct ActualEvapotranspiration cell response statistics object"))
                .def("output",output_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("output",output_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("output_value", &rc_stat::output_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns for cells matching catchments_ids at the i'th timestep")
                .def("pot_ratio",pot_ratio_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns the avg ratio (1-exp(-water_level*3/scale_factor)) for catcment_ids")
                .def("pot_ratio",pot_ratio_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns the ratio the ratio (1-exp(-water_level*3/scale_factor)) for cells matching catchments_ids at the i'th timestep")
				.def("pot_ratio_value", &rc_stat::pot_ratio_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns the ratio avg (1-exp(-water_level*3/scale_factor)) value for cells matching catchments_ids at the i'th timestep")
				;
        }
		template <class cell>
		static void hbv_actual_evapotranspiration(const char *cell_name) {
			char response_name[200]; sprintf(response_name, "%sHbvActualEvapotranspirationResponseStatistics", cell_name);
			typedef typename shyft::api::hbv_actual_evapotranspiration_cell_response_statistics<cell> rc_stat;

			rts_(rc_stat::*output_ts)(cids_, stat_scope) const = &rc_stat::output;
			vd_(rc_stat::*output_vd)(cids_, ix_, stat_scope) const = &rc_stat::output;
			class_<rc_stat>(response_name, "HbvActualEvapotranspiration response statistics", no_init)
				.def(init<std::shared_ptr<std::vector<cell>> >(args("cells"), "construct HbvActualEvapotranspiration cell response statistics object"))
				.def("output", output_ts, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
				.def("output", output_vd, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("output_value", &rc_stat::output_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns for cells matching catchments_ids at the i'th timestep")
				;
		}


        template <class cell>
        static void gamma_snow(const char *cell_name) {
            char state_name[200];sprintf(state_name,"%sGammaSnowStateStatistics",cell_name);
            char response_name[200];sprintf(response_name,"%sGammaSnowResponseStatistics",cell_name);
            typedef typename shyft::api::gamma_snow_cell_state_statistics<cell>    sc_stat;
            typedef typename shyft::api::gamma_snow_cell_response_statistics<cell> rc_stat;

            rts_ (sc_stat::*albedo_ts)(cids_, stat_scope) const = &sc_stat::albedo;
            vd_  (sc_stat::*albedo_vd)(cids_, ix_, stat_scope) const =&sc_stat::albedo;

            rts_ (sc_stat::*lwc_ts)(cids_, stat_scope) const = &sc_stat::lwc;
            vd_  (sc_stat::*lwc_vd)(cids_, ix_, stat_scope) const =&sc_stat::lwc;

            rts_ (sc_stat::*surface_heat_ts)(cids_, stat_scope) const = &sc_stat::surface_heat;
            vd_  (sc_stat::*surface_heat_vd)(cids_, ix_, stat_scope) const =&sc_stat::surface_heat;

            rts_ (sc_stat::*alpha_ts)(cids_, stat_scope) const = &sc_stat::alpha;
            vd_  (sc_stat::*alpha_vd)(cids_, ix_, stat_scope) const =&sc_stat::alpha;

            rts_ (sc_stat::*sdc_melt_mean_ts)(cids_, stat_scope) const = &sc_stat::sdc_melt_mean;
            vd_  (sc_stat::*sdc_melt_mean_vd)(cids_, ix_, stat_scope) const =&sc_stat::sdc_melt_mean;

            rts_ (sc_stat::*acc_melt_ts)(cids_, stat_scope) const = &sc_stat::acc_melt;
            vd_  (sc_stat::*acc_melt_vd)(cids_, ix_, stat_scope) const =&sc_stat::acc_melt;

            rts_ (sc_stat::*iso_pot_energy_ts)(cids_, stat_scope) const = &sc_stat::iso_pot_energy;
            vd_  (sc_stat::*iso_pot_energy_vd)(cids_, ix_, stat_scope) const =&sc_stat::iso_pot_energy;

            rts_ (sc_stat::*temp_swe_ts)(cids_, stat_scope) const = &sc_stat::temp_swe;
            vd_  (sc_stat::*temp_swe_vd)(cids_, ix_, stat_scope) const =&sc_stat::temp_swe;

            class_<sc_stat>(state_name,"GammaSnow state statistics",no_init)
                .def(init<std::shared_ptr<std::vector<cell>> >(args("cells"),"construct GammaSnow cell state statistics object"))
                .def("albedo",albedo_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("albedo",albedo_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("albedo_value", &sc_stat::albedo_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("lwc",lwc_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("lwc",lwc_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("lwc_value", &sc_stat::lwc_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("surface_heat",surface_heat_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("surface_heat",surface_heat_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("surface_heat_value", &sc_stat::surface_heat_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("alpha",alpha_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("alpha",alpha_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("alpha_value", &sc_stat::alpha_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("sdc_melt_mean",sdc_melt_mean_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("sdc_melt_mean",sdc_melt_mean_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("sdc_melt_mean_value", &sc_stat::sdc_melt_mean_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("acc_melt",acc_melt_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("acc_melt",acc_melt_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("acc_melt_value", &sc_stat::acc_melt_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("iso_pot_energy",iso_pot_energy_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("iso_pot_energy",iso_pot_energy_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("iso_pot_energy_value", &sc_stat::iso_pot_energy_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("temp_swe",temp_swe_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("temp_swe",temp_swe_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("temp_swe_value", &sc_stat::temp_swe_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
            ;


            rts_ (rc_stat::*sca_ts)(cids_, stat_scope) const = &rc_stat::sca;
            vd_  (rc_stat::*sca_vd)(cids_, ix_, stat_scope) const =&rc_stat::sca;

            rts_ (rc_stat::*swe_ts)(cids_, stat_scope) const = &rc_stat::swe;
            vd_  (rc_stat::*swe_vd)(cids_, ix_, stat_scope) const =&rc_stat::swe;

            rts_ (rc_stat::*outflow_ts)(cids_, stat_scope) const = &rc_stat::outflow;
            vd_  (rc_stat::*outflow_vd)(cids_, ix_, stat_scope) const =&rc_stat::outflow;

            rts_ (rc_stat::*glacier_melt_ts)(cids_, stat_scope) const = &rc_stat::glacier_melt;
            vd_  (rc_stat::*glacier_melt_vd)(cids_, ix_, stat_scope) const = &rc_stat::glacier_melt;

            class_<rc_stat>(response_name,"GammaSnow response statistics",no_init)
                .def(init<std::shared_ptr<std::vector<cell>> >(args("cells"),"construct GammaSnow cell response statistics object"))
                .def("outflow",outflow_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("outflow",outflow_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("outflow_value", &rc_stat::outflow_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("swe",swe_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("swe",swe_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("swe_value", &rc_stat::swe_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("sca",sca_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("sca",sca_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("sca_value", &rc_stat::sca_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("glacier_melt", glacier_melt_ts, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids[m3/s]")
                .def("glacier_melt", glacier_melt_vd, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep [m3/s]")
                .def("glacier_melt_value", &rc_stat::glacier_melt_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep[m3/s]")

                ;

        }
        template <class cell>
        static void universal_snow(const char *cell_name) {
            char state_name[200];sprintf(state_name,"%sUniversalSnowStateStatistics",cell_name);
            char response_name[200];sprintf(response_name,"%sUniversalSnowResponseStatistics",cell_name);
            typedef typename shyft::api::universal_snow_cell_state_statistics<cell>    sc_stat;
            typedef typename shyft::api::universal_snow_cell_response_statistics<cell> rc_stat;

            rts_ (sc_stat::*albedo_ts)(cids_, stat_scope) const = &sc_stat::albedo;
            vd_  (sc_stat::*albedo_vd)(cids_, ix_, stat_scope) const =&sc_stat::albedo;

            rts_ (sc_stat::*lwc_ts)(cids_, stat_scope) const = &sc_stat::lwc;
            vd_  (sc_stat::*lwc_vd)(cids_, ix_, stat_scope) const =&sc_stat::lwc;

            rts_ (sc_stat::*surface_heat_ts)(cids_, stat_scope) const = &sc_stat::surface_heat;
            vd_  (sc_stat::*surface_heat_vd)(cids_, ix_, stat_scope) const =&sc_stat::surface_heat;

            rts_ (sc_stat::*alpha_ts)(cids_, stat_scope) const = &sc_stat::alpha;
            vd_  (sc_stat::*alpha_vd)(cids_, ix_, stat_scope) const =&sc_stat::alpha;

            rts_ (sc_stat::*sdc_melt_mean_ts)(cids_, stat_scope) const = &sc_stat::sdc_melt_mean;
            vd_  (sc_stat::*sdc_melt_mean_vd)(cids_, ix_, stat_scope) const =&sc_stat::sdc_melt_mean;

            rts_ (sc_stat::*acc_melt_ts)(cids_, stat_scope) const = &sc_stat::acc_melt;
            vd_  (sc_stat::*acc_melt_vd)(cids_, ix_, stat_scope) const =&sc_stat::acc_melt;

            rts_ (sc_stat::*iso_pot_energy_ts)(cids_, stat_scope) const = &sc_stat::iso_pot_energy;
            vd_  (sc_stat::*iso_pot_energy_vd)(cids_, ix_, stat_scope) const =&sc_stat::iso_pot_energy;

            rts_ (sc_stat::*temp_swe_ts)(cids_, stat_scope) const = &sc_stat::temp_swe;
            vd_  (sc_stat::*temp_swe_vd)(cids_, ix_, stat_scope) const =&sc_stat::temp_swe;

            class_<sc_stat>(state_name,"UniversalSnow state statistics",no_init)
                .def(init<std::shared_ptr<std::vector<cell>> >(args("cells"),"construct UniversalSnow cell state statistics object"))
                .def("albedo",albedo_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("albedo",albedo_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("albedo_value", &sc_stat::albedo_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("lwc",lwc_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("lwc",lwc_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("lwc_value", &sc_stat::lwc_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("surface_heat",surface_heat_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("surface_heat",surface_heat_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("surface_heat_value", &sc_stat::surface_heat_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("alpha",alpha_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("alpha",alpha_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("alpha_value", &sc_stat::alpha_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("sdc_melt_mean",sdc_melt_mean_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("sdc_melt_mean",sdc_melt_mean_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("sdc_melt_mean_value", &sc_stat::sdc_melt_mean_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("acc_melt",acc_melt_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("acc_melt",acc_melt_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("acc_melt_value", &sc_stat::acc_melt_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("iso_pot_energy",iso_pot_energy_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("iso_pot_energy",iso_pot_energy_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("iso_pot_energy_value", &sc_stat::iso_pot_energy_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("temp_swe",temp_swe_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("temp_swe",temp_swe_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("temp_swe_value", &sc_stat::temp_swe_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
            ;


            rts_ (rc_stat::*sca_ts)(cids_, stat_scope) const = &rc_stat::sca;
            vd_  (rc_stat::*sca_vd)(cids_, ix_, stat_scope) const =&rc_stat::sca;

            rts_ (rc_stat::*swe_ts)(cids_, stat_scope) const = &rc_stat::swe;
            vd_  (rc_stat::*swe_vd)(cids_, ix_, stat_scope) const =&rc_stat::swe;

            rts_ (rc_stat::*outflow_ts)(cids_, stat_scope) const = &rc_stat::outflow;
            vd_  (rc_stat::*outflow_vd)(cids_, ix_, stat_scope) const =&rc_stat::outflow;

            rts_ (rc_stat::*glacier_melt_ts)(cids_, stat_scope) const = &rc_stat::glacier_melt;
            vd_  (rc_stat::*glacier_melt_vd)(cids_, ix_, stat_scope) const = &rc_stat::glacier_melt;

            class_<rc_stat>(response_name,"UniversalSnow response statistics",no_init)
                .def(init<std::shared_ptr<std::vector<cell>> >(args("cells"),"construct UniversalSnow cell response statistics object"))
                .def("outflow",outflow_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("outflow",outflow_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("outflow_value", &rc_stat::outflow_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("swe",swe_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("swe",swe_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("swe_value", &rc_stat::swe_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("sca",sca_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("sca",sca_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("sca_value", &rc_stat::sca_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("glacier_melt", glacier_melt_ts, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids[m3/s]")
                .def("glacier_melt", glacier_melt_vd, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep [m3/s]")
                .def("glacier_melt_value", &rc_stat::glacier_melt_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep[m3/s]")

                ;

        }

        template <class cell>
        static void hbv_snow(const char *cell_name) {
            char state_name[200];sprintf(state_name,"%sHBVSnowStateStatistics",cell_name);
            char response_name[200];sprintf(response_name,"%sHBVSnowResponseStatistics",cell_name);
            typedef typename shyft::api::hbv_snow_cell_state_statistics<cell>    sc_stat;
            typedef typename shyft::api::hbv_snow_cell_response_statistics<cell> rc_stat;

            rts_ (sc_stat::*swe_ts)(cids_, stat_scope) const = &sc_stat::swe;
            vd_  (sc_stat::*swe_vd)(cids_, ix_, stat_scope) const =&sc_stat::swe;
            rts_ (sc_stat::*sca_ts)(cids_, stat_scope) const = &sc_stat::sca;
            vd_  (sc_stat::*sca_vd)(cids_, ix_, stat_scope) const =&sc_stat::sca;

            class_<sc_stat>(state_name,"HBVSnow state statistics",no_init)
                .def(init<std::shared_ptr<std::vector<cell>> >(args("cells"),"construct HBVSnow cell state statistics object"))
                .def("swe",swe_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("swe",swe_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("swe_value", &sc_stat::swe_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("sca",sca_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("sca",sca_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("sca_value", &sc_stat::sca_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
            ;


            rts_ (rc_stat::*outflow_ts)(cids_, stat_scope) const = &rc_stat::outflow;
            vd_  (rc_stat::*outflow_vd)(cids_, ix_, stat_scope) const =&rc_stat::outflow;
            rts_ (rc_stat::*glacier_melt_ts)(cids_, stat_scope) const = &rc_stat::glacier_melt;
            vd_  (rc_stat::*glacier_melt_vd)(cids_, ix_, stat_scope) const = &rc_stat::glacier_melt;

            class_<rc_stat>(response_name,"HBVSnow response statistics",no_init)
                .def(init<std::shared_ptr<std::vector<cell>> >(args("cells"),"construct HBVSnow cell response statistics object"))
                .def("outflow",outflow_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("outflow",outflow_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("outflow_value", &rc_stat::outflow_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("glacier_melt", glacier_melt_ts, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids[m3/s]")
                .def("glacier_melt", glacier_melt_vd, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep [m3/s]")
                .def("glacier_melt_value", &rc_stat::glacier_melt_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep[m3/s]")
                ;
        }
        

        template <class cell>
        static void hbv_physical_snow(const char *cell_name) {
            char state_name[200];sprintf(state_name,"%sHBVPhysicalSnowStateStatistics",cell_name);
            char response_name[200];sprintf(response_name,"%sHBVPhysicalSnowResponseStatistics",cell_name);
            typedef typename shyft::api::hbv_physical_snow_cell_state_statistics<cell>    sc_stat;
            typedef typename shyft::api::hbv_physical_snow_cell_response_statistics<cell> rc_stat;

            rts_ (sc_stat::*swe_ts)(cids_, stat_scope) const = &sc_stat::swe;
            vd_  (sc_stat::*swe_vd)(cids_, ix_, stat_scope) const =&sc_stat::swe;
            rts_ (sc_stat::*sca_ts)(cids_, stat_scope) const = &sc_stat::sca;
            vd_  (sc_stat::*sca_vd)(cids_, ix_, stat_scope) const =&sc_stat::sca;
            rts_ (sc_stat::*surface_heat_ts)(cids_, stat_scope) const = &sc_stat::surface_heat;
            vd_  (sc_stat::*surface_heat_vd)(cids_, ix_, stat_scope) const =&sc_stat::surface_heat;

            class_<sc_stat>(state_name,"HBVPhysicalSnow state statistics",no_init)
                .def(init<std::shared_ptr<std::vector<cell>> >(args("cells"),"construct HBVPhysicalSnow cell state statistics object"))
                .def("swe",swe_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("swe",swe_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("swe_value", &sc_stat::swe_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("sca",sca_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("sca",sca_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("sca_value", &sc_stat::sca_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("surface_heat",surface_heat_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("surface_heat",surface_heat_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("surface_heat_value", &sc_stat::sca_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")

            ;


            rts_ (rc_stat::*outflow_ts)(cids_, stat_scope) const = &rc_stat::outflow;
            vd_  (rc_stat::*outflow_vd)(cids_, ix_, stat_scope) const =&rc_stat::outflow;
            rts_ (rc_stat::*glacier_melt_ts)(cids_, stat_scope) const = &rc_stat::glacier_melt;
            vd_  (rc_stat::*glacier_melt_vd)(cids_, ix_, stat_scope) const = &rc_stat::glacier_melt;

            class_<rc_stat>(response_name,"HBVSnow response statistics",no_init)
                .def(init<std::shared_ptr<std::vector<cell>> >(args("cells"),"construct HBVSnow cell response statistics object"))
                .def("outflow",outflow_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("outflow",outflow_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("outflow_value", &rc_stat::outflow_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("glacier_melt", glacier_melt_ts, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids[m3/s]")
                .def("glacier_melt", glacier_melt_vd, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep [m3/s]")
                .def("glacier_melt_value", &rc_stat::glacier_melt_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep[m3/s]")
                ;
        }


        template <class cell>
        static void skaugen(const char *cell_name) {
            char state_name[200];sprintf(state_name,"%sSkaugenStateStatistics",cell_name);
            char response_name[200];sprintf(response_name,"%sSkaugenResponseStatistics",cell_name);
            typedef typename shyft::api::skaugen_cell_state_statistics<cell>    sc_stat;
            typedef typename shyft::api::skaugen_cell_response_statistics<cell> rc_stat;

            rts_ (sc_stat::*alpha_ts)(cids_, stat_scope) const = &sc_stat::alpha;
            vd_  (sc_stat::*alpha_vd)(cids_, ix_, stat_scope) const =&sc_stat::alpha;
            rts_ (sc_stat::*nu_ts)(cids_, stat_scope) const = &sc_stat::nu;
            vd_  (sc_stat::*nu_vd)(cids_, ix_, stat_scope) const =&sc_stat::nu;
            rts_ (sc_stat::*lwc_ts)(cids_, stat_scope) const = &sc_stat::lwc;
            vd_  (sc_stat::*lwc_vd)(cids_, ix_, stat_scope) const =&sc_stat::lwc;
            rts_ (sc_stat::*residual_ts)(cids_, stat_scope) const = &sc_stat::residual;
            vd_  (sc_stat::*residual_vd)(cids_, ix_, stat_scope) const =&sc_stat::residual;
            rts_ (sc_stat::*swe_ts)(cids_, stat_scope) const = &sc_stat::swe;
            vd_  (sc_stat::*swe_vd)(cids_, ix_, stat_scope) const =&sc_stat::swe;
            rts_ (sc_stat::*sca_ts)(cids_, stat_scope) const = &sc_stat::sca;
            vd_  (sc_stat::*sca_vd)(cids_, ix_, stat_scope) const =&sc_stat::sca;

            class_<sc_stat>(state_name,"Skaugen snow state statistics",no_init)
                .def(init<std::shared_ptr<std::vector<cell>> >(args("cells"),"construct Skaugen snow cell state statistics object"))
                .def("alpha",alpha_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("alpha",alpha_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("alpha_value", &sc_stat::alpha_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("nu",nu_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("nu",nu_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("nu_value",&sc_stat::nu_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("lwc",lwc_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("lwc",lwc_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("lwc_value", &sc_stat::lwc_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("residual",residual_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("residual",residual_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("residual_value", &sc_stat::residual_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("swe",swe_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("swe",swe_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("swe_value", &sc_stat::swe_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("sca",sca_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("sca",sca_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("sca_value", &sc_stat::sca_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
            ;


            rts_ (rc_stat::*outflow_ts)(cids_, stat_scope) const = &rc_stat::outflow;
            vd_  (rc_stat::*outflow_vd)(cids_, ix_, stat_scope) const =&rc_stat::outflow;
            rts_ (rc_stat::*total_stored_water_ts)(cids_, stat_scope) const = &rc_stat::total_stored_water;
            vd_  (rc_stat::*total_stored_water_vd)(cids_, ix_, stat_scope) const =&rc_stat::total_stored_water;
            rts_ (rc_stat::*glacier_melt_ts)(cids_, stat_scope) const = &rc_stat::glacier_melt;
            vd_  (rc_stat::*glacier_melt_vd)(cids_, ix_, stat_scope) const = &rc_stat::glacier_melt;

            class_<rc_stat>(response_name,"Skaugen snow response statistics",no_init)
                .def(init<std::shared_ptr<std::vector<cell>> >(args("cells"),"construct Skaugen snow cell response statistics object"))
                .def("outflow",outflow_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("outflow",outflow_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("outflow_value", &rc_stat::outflow_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("total_stored_water",total_stored_water_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("total_stored_water",total_stored_water_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("total_stored_water_value", &rc_stat::total_stored_water_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("glacier_melt", glacier_melt_ts, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids[m3/s]")
                .def("glacier_melt", glacier_melt_vd, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep [m3/s]")
                .def("glacier_melt_value", &rc_stat::glacier_melt_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep[m3/s]")
                ;
        }

        template <class cell>
        static void basic_cell(const char *cell_name) {
            char base_name[200];sprintf(base_name,"%sStatistics",cell_name);
            typedef typename shyft::api::basic_cell_statistics<cell> bc_stat;

            rts_ (bc_stat::*discharge_ts)(cids_,stat_scope) const = &bc_stat::discharge;
            vd_  (bc_stat::*discharge_vd)(cids_,ix_,stat_scope) const =&bc_stat::discharge;

            rts_(bc_stat::*charge_ts)(cids_,stat_scope) const = &bc_stat::charge;
            vd_(bc_stat::*charge_vd)(cids_, ix_,stat_scope) const = &bc_stat::charge;

            rts_ (bc_stat::*temperature_ts)(cids_,stat_scope) const = &bc_stat::temperature;
            vd_  (bc_stat::*temperature_vd)(cids_,ix_,stat_scope) const =&bc_stat::temperature;

            rts_ (bc_stat::*radiation_ts)(cids_,stat_scope) const = &bc_stat::radiation;
            vd_  (bc_stat::*radiation_vd)(cids_,ix_,stat_scope) const =&bc_stat::radiation;

            rts_ (bc_stat::*wind_speed_ts)(cids_,stat_scope) const = &bc_stat::wind_speed;
            vd_  (bc_stat::*wind_speed_vd)(cids_,ix_,stat_scope) const =&bc_stat::wind_speed;

            rts_ (bc_stat::*rel_hum_ts)(cids_,stat_scope) const = &bc_stat::rel_hum;
            vd_  (bc_stat::*rel_hum_vd)(cids_,ix_,stat_scope) const =&bc_stat::rel_hum;

            rts_ (bc_stat::*precipitation_ts)(cids_,stat_scope) const = &bc_stat::precipitation;
            vd_  (bc_stat::*precipitation_vd)(cids_,ix_,stat_scope) const =&bc_stat::precipitation;



            class_<bc_stat>(base_name,
                    doc_intro("This class provides statistics for group of cells, as specified")
                    doc_intro("by the list of catchment identifiers, or list of cell-indexes passed to the methods.")
                    doc_intro("It is provided both as a separate class, but is also provided")
                    doc_intro("automagically through the region_model.statistics property."),
                    no_init
                )
                .def(init<std::shared_ptr<std::vector<cell>> >((py::arg("self"),py::arg("cells")),"construct basic cell statistics object for the list of cells, usually from the region-model"))
                .def("discharge",discharge_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("discharge",discharge_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("discharge_value", &bc_stat::discharge_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("charge", charge_ts, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum charge[m^3/s] for catcment_ids")
                .def("charge", charge_vd, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns charge[m^3/s]  for cells matching catchments_ids at the i'th timestep")
                .def("charge_value", &bc_stat::charge_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns charge[m^3/s] for cells matching catchments_ids at the i'th timestep")
                .def("temperature",temperature_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("temperature",temperature_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("temperature_value", &bc_stat::temperature_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("precipitation",precipitation_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("precipitation",precipitation_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("precipitation_value", &bc_stat::precipitation_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("radiation",radiation_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("radiation",radiation_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("radiation_value", &bc_stat::radiation_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("wind_speed",wind_speed_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("wind_speed",wind_speed_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("wind_speed_value", &bc_stat::wind_speed_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
                .def("rel_hum",rel_hum_ts,(py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns sum  for catcment_ids")
                .def("rel_hum",rel_hum_vd,(py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix),"returns  for cells matching catchments_ids at the i'th timestep")
				.def("rel_hum_value", &bc_stat::rel_hum_value, (py::arg("self"),py::arg("indexes"),py::arg("i"),py::arg("ix_type")=stat_scope::catchment_ix), "returns  for cells matching catchments_ids at the i'th timestep")
				.def("total_area", &bc_stat::total_area, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns total area[m2] for cells matching catchments_ids")
				.def("forest_area", &bc_stat::forest_area, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns forest area[m2] for cells matching catchments_ids")
				.def("glacier_area", &bc_stat::glacier_area, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns glacier area[m2] for cells matching catchments_ids")
				.def("lake_area", &bc_stat::lake_area, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns lake area[m2] for cells matching catchments_ids")
				.def("reservoir_area", &bc_stat::reservoir_area, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns reservoir area[m2] for cells matching catchments_ids")
				.def("unspecified_area", &bc_stat::unspecified_area, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns unspecified area[m2] for cells matching catchments_ids")
                .def("snow_storage_area", &bc_stat::snow_storage_area, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns snow_storage area where snow can build up[m2], eg total_area - lake and reservoir")
				.def("elevation", &bc_stat::elevation, (py::arg("self"),py::arg("indexes"),py::arg("ix_type")=stat_scope::catchment_ix), "returns area-average elevation[m.a.s.l] for cells matching catchments_ids")
				;
        }
    }
}
