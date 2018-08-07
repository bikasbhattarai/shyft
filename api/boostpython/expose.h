/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#pragma once
#include "expose_statistics.h"
#include "api/api.h"
#include "api/api_state.h"
namespace expose {
    using namespace boost::python;
    namespace py=boost::python;
    using namespace std;

    template <class T>
    static vector<double> geo_cell_data_vector(shared_ptr<vector<T>> cell_vector) {
        vector<double> r; r.reserve(shyft::api::geo_cell_data_io::size()*cell_vector->size());//Assume approx 200 chars pr. cell
        for(const auto& cell:*cell_vector)
            shyft::api::geo_cell_data_io::push_to_vector(r,cell.geo);
        return r;
    }

    template <class T>
    static vector<T> create_from_geo_cell_data_vector(const vector<double>& s) {
        if(s.size()==0 || s.size()% shyft::api::geo_cell_data_io::size())
            throw invalid_argument("create_from_geo_cell_data_vector: size of vector of double must be multiple of 11");
        vector<T> r; r.reserve(s.size()/shyft::api::geo_cell_data_io::size());// assume this is ok size for now
        for(size_t i=0;i<s.size();i+=shyft::api::geo_cell_data_io::size()) {
            T cell;
            cell.geo=shyft::api::geo_cell_data_io::from_raw_vector(s.data()+i);
            r.push_back(cell);
        }
        return r;
    }
    template <class S,class S_ID_V>
    static shared_ptr<vector<S>> extract_cell_state_vector(const S_ID_V&v) {
        auto r = make_shared<vector<S>>();
        r->reserve(v->size());
        for(const auto&sid:*v) {
            r->push_back(sid.state);
        }
        return r;
    }

    template<class C>
    static void cell_state_etc(const char *stack_name) {
        typedef typename C::state_t cstate_t;
        typedef typename shyft::api::cell_state_with_id<cstate_t> CellState;
        char cs_name[200];sprintf(cs_name, "%sStateWithId", stack_name);
        class_<CellState>(cs_name, "Keep the cell id and cell state")
            .def_readwrite("id", &CellState::id, "the cell identifier for the state")
            .def_readwrite("state", &CellState::state, "the cell state")
            .def("cell_state", &shyft::api::cell_state_id_of, args("geo_cell_data"), "create a cell state with id for the supplied cell.geo")
            .staticmethod("cell_state")
            ;
        char csv_name[200];sprintf(csv_name, "%sVector", cs_name);
        typedef std::vector<CellState> cell_state_with_id_vector;
        typedef std::shared_ptr<cell_state_with_id_vector> cell_state_with_id_vector_;
        class_<cell_state_with_id_vector, bases<>, cell_state_with_id_vector_  >(csv_name, "vector of cell state")
            .def(vector_indexing_suite<std::vector<CellState>>())
            ;
        def("extract_state_vector",extract_cell_state_vector<cstate_t,cell_state_with_id_vector_>,(py::arg("cell_state_id_vector")),
            doc_intro("Given a cell-state-with-id-vector, returns a pure state vector that can be inserted directly into region-model")
            doc_parameters()
            doc_parameter("cell_state_id_vector","xStateWithIdVector","a complete consistent with region-model vector, all states, as in cell-order")
            doc_returns("cell_state_vector","XStateVector","a vector with cell-id removed, order preserved")
        );
        def("serialize", shyft::api::serialize_to_bytes<CellState>, args("states"), "make a blob out of the states");
        def("deserialize", shyft::api::deserialize_from_bytes<CellState>, args("bytes", "states"), "from a blob, fill in states");
    }

    template <class C>
    static void cell_state_io(const char *cell_name) {

        char csh_name[200];sprintf(csh_name, "%sStateHandler", cell_name);
        typedef shyft::api::state_io_handler<C> CellStateHandler;
        class_<CellStateHandler>(csh_name, "Provides functionality to extract and restore state from cells")
            .def(init<std::shared_ptr<std::vector<C>> >(args("cells"),"construct a cell state handler for the supplied cells"))
            .def("extract_state", &CellStateHandler::extract_state,( py::arg("self"),py::arg("cids")),
                doc_intro("Extract cell state for the optionaly specified catchment ids, cids")
                doc_parameters()
                doc_parameter("cids","IntVector","list of catchment-id's, if empty, extract all")
                doc_returns("cell_states","CellStateIdVector","the state with identifier for the cells")
            )
            .def("apply_state", &CellStateHandler::apply_state,( py::arg("self"), py::arg("cell_id_state_vector"), py::arg("cids")),
                doc_intro("apply the supplied cell-identified state to the cells,")
                doc_intro("limited to the optionally supplied catchment id's")
                doc_intro("If no catchment-id's specified, it applies to all cells")
                doc_parameters()
                doc_parameter("cell_id_state_vector","","")
                doc_parameter("cids","IntVector","list of catchment-id's, if empty, apply all")
                doc_returns("not_applied_list","IntVector",
                            "a list of indices into cell_id_state_vector that did not match any cells\n"
                            "\t taken into account the optionally catchment-id specification\n")
            )
        ;


    }

    template <class T>
    static void cell(const char *cell_name,const char* cell_doc) {
      class_<T>(cell_name,cell_doc)
        .def_readwrite("geo",&T::geo,"geo_cell_data information for the cell, such as mid-point, forest-fraction and other cell-specific personalities.")
        .add_property("parameter",&T::get_parameter,&T::set_parameter,"reference to parameter for this cell, typically shared for a catchment")
        .def_readwrite("env_ts",&T::env_ts,"environment time-series as projected to the cell after the interpolation/preparation step")
        .def_readwrite("state",&T::state,"Current state of the cell")
        .def_readonly("sc",&T::sc,"state collector for the cell")
        .def_readonly("rc",&T::rc,"response collector for the cell")
        .def("set_parameter",&T::set_parameter,args("parameter"),"set the cell method stack parameters, typical operations at region_level, executed after the interpolation, before the run")
        .def("set_state_collection",&T::set_state_collection,args("on_or_off"),"collecting the state during run could be very useful to understand models")
        .def("set_snow_sca_swe_collection",&T::set_snow_sca_swe_collection,"collecting the snow sca and swe on for calibration scenario")
        .def("mid_point",&T::mid_point,"returns geo.mid_point()",return_internal_reference<>())
        .def("run",&T::run,(py::arg("self"),py::arg("time_axis"),py::arg("start_step"),py::arg("n_steps")),
             doc_intro("run the cell (given it's initialized)")
             doc_intro("before run, the caller must ensure the cell is ready to run, is initialized")
             doc_intro("after the run, the cell state, as well as resource collector/state-collector is updated")
             doc_parameters()
             doc_parameter("time_axis","TimeAxisFixedDeltaT","time-axis to run, should match the run-time-axis used for env_ts")
             doc_parameter("start_step","int","first interval, ref. time-axis to start run")
             doc_parameter("n_steps","int","number of time-steps to run")
             )
      ;
      char cv[200];sprintf(cv,"%sVector",cell_name);
      class_<vector<T>,bases<>,shared_ptr<vector<T>> > (cv,"vector of cells")
        .def(vector_indexing_suite<vector<T>>())
        .def("geo_cell_data_vector", geo_cell_data_vector<T>,
             "returns a persistable DoubleVector representation of of geo_cell_data for all cells.\n"
             "that object can in turn be used to construct a <Cell>Vector of any cell type\n"
             "using the <Cell>Vector.create_from_geo_cell_data_vector")
             .staticmethod("geo_cell_data_vector")
        .def("create_from_geo_cell_data_vector",create_from_geo_cell_data_vector<T>,
             "create a cell-vector filling in the geo_cell_data records as given by the DoubleVector.\n"
             "This function works together with the geo_cell_data_vector static method\n"
             "that provides a correctly formatted persistable vector\n"
             "Notice that the context and usage of these two functions is related\n"
             "to python orchestration and repository data-caching\n")
             .staticmethod("create_from_geo_cell_data_vector")

        ;
      expose::statistics::basic_cell<T>(cell_name);//common for all type of cells, so expose it here
      cell_state_io<T>(cell_name);
    }


    template <class M>
    static void model(const char *model_name,const char *model_doc) {
        char m_doc[5000];
        sprintf(m_doc,
            " %s , a region_model is the calculation model for a region, where we can have\n"
            "one or more catchments.\n"
            "The role of the region_model is to describe region, so that we can run the\n"
            "region computational model efficiently for a number of type of cells, interpolation and\n"
            "catchment level algorihtms.\n"
            "\n"
            "The region model keeps a list of cells, of specified type \n"
            "as well as parameters for the cells.\n"
            "The model also keeps state, such as region_env(forcing variables), time-axis and intial state\n"
            "- they are non-empty after initializing, and running the model\n"
                ,model_name);
        // NOTE: explicit expansion of the run_interpolate method is needed here, using this specific syntax
        auto run_interpolation_f= &M::run_interpolation;
        auto run_interpolation_f_g=&M::run_interpolation_g;
		auto interpolate_f = &M::interpolate;
        class_<M>(model_name,m_doc,no_init)
	     .def(init<const M&>((py::arg("self"),py::arg("other_model")),
            doc_intro("Create a copy of the other_model")
            doc_parameters()
            doc_parameter("other_model","RegionModel","region-model to copy")
          ))
         .def(init< const vector<shyft::core::geo_cell_data>&, const typename M::parameter_t& >( (py::arg("self"),py::arg("geo_data_vector"), py::arg("region_param")), 
            doc_intro("Creates a model from GeoCellDataVector and region model parameters")
            doc_parameters()
            doc_parameter("geo_data_vector","GeoCellDataVector","contains the geo-related characteristics for the cells")
            doc_parameter("region_param","Parameter","contains the parameters for all cells of this region model")
          ))
         .def(init< shared_ptr< vector<typename M::cell_t> >&, const typename M::parameter_t&, const map<int,typename M::parameter_t>& >( 
            (py::arg("self"),py::arg("cells"), py::arg("region_param"), py::arg("catchment_parameters")),
            doc_intro("Creates a model from cells and region model parameters, and specified catchment parameters")
            doc_intro("The cell-vector and catchment-id's should match those specified in the catchment_parameters mapping")
            doc_parameters()
            doc_parameter("cells","CellVector","contains the cells, each with geo-properties and type matching the region-model type")
            doc_parameter("region_param","Parameter","contains the parameters for cells that does not have catchment specific parameters")
            doc_parameter("catchment_parameters","ParameterMap","contains mapping (a kind of dict, where the key is catchment-id and value is parameters for cells matching catchment-id")
             
          ))
         .def_readonly("time_axis",&M::time_axis,"the time_axis (type TimeAxisFixedDeltaT) as set from run_interpolation, determines the time-axis for run")
		 .def_readwrite("interpolation_parameter",&M::ip_parameter,"the most recently used interpolation parameter as passed to run_interpolation or interpolate routine")
         .def_readwrite("initial_state",&M::initial_state,"empty or the the initial state as established on the first invokation of .set_states() or .run_cells()")
         .def_readwrite("ncore",&M::ncore,
                        "determines how many core to utilize during run_cell processing,\n"
                        "0(=default) means detect by hardware probe"
                        )
         .def_readwrite("region_env",&M::region_env,"empty or the region_env as passed to run_interpolation() or interpolate()")
         .def_readwrite("river_network",&M::river_network,
                        "river network that when enabled do the routing part of the region-model\n"
                        "See also RiverNetwork class for how to build a working river network\n"
                        "Then use the connect_catchment_to_river(cid,rid) method\n"
                        "to route cell discharge into the river-network\n")
         .def("has_routing",&M::has_routing,(py::arg("self")),"true if some cells routes to river-network")
         .def("river_output_flow_m3s",&M::river_output_flow_m3s,(py::arg("self"),py::arg("rid")),"returns the routed output flow of the specified river id (rid))")
         .def("river_upstream_inflow_m3s",&M::river_upstream_inflow_m3s, (py::arg("self"), py::arg("rid")),"returns the routed upstream inflow to the specified river id (rid))")
         .def("river_local_inflow_m3s",&M::river_local_inflow_m3s, (py::arg("self"), py::arg("rid")),"returns the routed local inflow from connected cells to the specified river id (rid))")
		 .def("connect_catchment_to_river",&M::connect_catchment_to_river, (py::arg("self"), py::arg("cid"),py::arg("rid")),
            doc_intro("Connect routing of all the cells in the specified catchment id to the specified river id")
            doc_intro("")
            doc_parameters()
            doc_parameter("cid","int","catchment identifier")
            doc_parameter("rid","int","river identifier, can be set to 0 to indicate disconnect from routing")
         )
         .def("number_of_catchments",&M::number_of_catchments, (py::arg("self")),"compute and return number of catchments using info in cells.geo.catchment_id()")
		 .def_readonly("catchment_ids",&M::catchment_ids,
             doc_intro("provides the list of catchment identifiers,'cids' within this model")
         )
         .def("extract_geo_cell_data",&M::extract_geo_cell_data,(py::arg("self")),
             "extracts the geo_cell_data and return it as GeoCellDataVector that can\n"
             "be passed into a the constructor of a new region-model (clone-operation)\n"
         )
         .def("initialize_cell_environment",&M::initialize_cell_environment,(py::arg("self"),py::arg("time_axis")),
                doc_intro("Initializes the cell enviroment (cell.env.ts* )")
                doc_intro("")
                doc_intro("The method initializes the cell environment, that keeps temperature, precipitation etc")
                doc_intro("that is local to the cell.The initial values of these time - series is set to zero.")
                doc_intro("The region-model time-axis is set to the supplied time-axis, so that")
                doc_intro("the any calculation steps will use the supplied time-axis.")
                doc_intro("This call is needed once prior to call to the .interpolate() or .run_cells() methods")
                doc_intro("")
                doc_intro("The call ensures that all cells.env ts are reset to zero, with a time-axis and")
                doc_intro(" value-vectors according to the supplied time-axis.")
                doc_intro(" Also note that the region-model.time_axis is set to the supplied time-axis.")
                doc_intro("")
                doc_parameters()
                doc_parameter("time_axis","TimeAxisFixedDeltaT","specifies the time-axis for the region-model, and thus the cells")
                doc_returns("nothing","","")
		 )
		 .def("initialize_cell_environment",&M::initialize_cell_environment_g,(py::arg("self"),py::arg("time_axis")),
                doc_intro("Initializes the cell enviroment (cell.env.ts* )")
                doc_intro("")
                doc_intro("The method initializes the cell environment, that keeps temperature, precipitation etc")
                doc_intro("that is local to the cell.The initial values of these time - series is set to zero.")
                doc_intro("The region-model time-axis is set to the supplied time-axis, so that")
                doc_intro("the any calculation steps will use the supplied time-axis.")
                doc_intro("This call is needed once prior to call to the .interpolate() or .run_cells() methods")
                doc_intro("")
                doc_intro("The call ensures that all cells.env ts are reset to zero, with a time-axis and")
                doc_intro(" value-vectors according to the supplied time-axis.")
                doc_intro(" Also note that the region-model.time_axis is set to the supplied time-axis.")
                doc_intro("")
                doc_parameters()
                doc_parameter("time_axis","TimeAxis","specifies the time-axis (fixed type) for the region-model, and thus the cells")
                doc_returns("nothing","","")
		 )

		 .def("interpolate", interpolate_f, (py::arg("self"),py::arg("interpolation_parameter"),py::arg("env"),py::arg("best_effort")=true),
                doc_intro("do interpolation interpolates region_environment temp,precip,rad.. point sources")
                doc_intro("to a value representative for the cell.mid_point().")
                doc_intro("")
                doc_intro("note: initialize_cell_environment should be called once prior to this function")
                doc_intro("")
                doc_intro("Only supplied vectors of temp, precip etc. are interpolated, thus")
                doc_intro("the user of the class can choose to put in place distributed series in stead.")
                doc_intro("")
                doc_parameters()
                doc_parameter("interpolation_parameter","InterpolationParameter","contains wanted parameters for the interpolation")
                doc_parameter("env","RegionEnvironment","contains the region environment with geo-localized time-series for P,T,R,W,Rh")
                doc_parameter("best_effort","bool","default=True, don't throw, just return True/False if problem, with best_effort, unfilled values is nan")
                doc_returns("success","bool","True if interpolation runs with no exceptions(btk,raises if to few neighbours)")
		 )
         .def("run_cells",&M::run_cells,(py::arg("self"),py::arg("use_ncore")=0,py::arg("start_step")=0,py::arg("n_steps")=0),
                doc_intro("run_cells calculations over specified time_axis,optionally with thread_cell_count, start_step and n_steps")
                doc_intro("require that initialize(time_axis) or run_interpolation is done first")
                doc_intro("If start_step and n_steps are specified, only the specified part of the time-axis is covered.")
                doc_intro("notice that in any case, the current model state is used as a starting point")
                doc_parameters()
                doc_parameter("use_ncore","int","number of worker threads, or cores to use, if 0 is passed, the the core-count is used to determine the count")
                doc_parameter("start_step","int","start_step in the time-axis to start at, default=0, meaning start at the beginning")
                doc_parameter("n_steps","int","number of steps to run in a partial run, default=0 indicating the complete time-axis is covered")
         )
         .def("run_interpolation",run_interpolation_f,(py::arg("self"),py::arg("interpolation_parameter"),py::arg("time_axis"),py::arg("env"),py::arg("best_effort")=true),
                doc_intro("run_interpolation interpolates region_environment temp,precip,rad.. point sources")
                doc_intro("to a value representative for the cell.mid_point().")
                doc_intro("")
                doc_intro("note: This function is equivalent to")
                doc_intro("    self.initialize_cell_environment(time_axis)")
                doc_intro("    self.interpolate(interpolation_parameter,env)")
                doc_parameters()
                doc_parameter("interpolation_parameter","InterpolationParameter","contains wanted parameters for the interpolation")
                doc_parameter("time_axis","TimeAxisFixedDeltaT","should be equal to the time-axis the region_model is prepared running for")
                doc_parameter("env","RegionEnvironment","contains the ref: region_environment type")
                doc_parameter("best_effort","bool","default=True, don't throw, just return True/False if problem, with best_effort, unfilled values is nan")
                doc_returns("success","bool","True if interpolation runs with no exceptions(btk,raises if to few neighbours)")
            )
         .def("run_interpolation",run_interpolation_f_g,(py::arg("self"),py::arg("interpolation_parameter"),py::arg("time_axis"),py::arg("env"),py::arg("best_effort")=true),
                doc_intro("run_interpolation interpolates region_environment temp,precip,rad.. point sources")
                doc_intro("to a value representative for the cell.mid_point().")
                doc_intro("")
                doc_intro("note: This function is equivalent to")
                doc_intro("    self.initialize_cell_environment(time_axis)")
                doc_intro("    self.interpolate(interpolation_parameter,env)")
                doc_parameters()
                doc_parameter("interpolation_parameter","InterpolationParameter","contains wanted parameters for the interpolation")
                doc_parameter("time_axis","TimeAxis","should be equal to the time-axis the region_model is prepared running for")
                doc_parameter("env","RegionEnvironment","contains the ref: region_environment type")
                doc_parameter("best_effort","bool","default=True, don't throw, just return True/False if problem, with best_effort, unfilled values is nan")
                doc_returns("success","bool","True if interpolation runs with no exceptions(btk,raises if to few neighbours)")
            )

         .def("set_region_parameter",&M::set_region_parameter,(py::arg("self"),py::arg("p")),
                    "set the region parameter, apply it to all cells \n"
                    "that do *not* have catchment specific parameters.\n")
         .def("get_region_parameter",&M::get_region_parameter,(py::arg("self")),"provide access to current region parameter-set",return_internal_reference<>())
         .def("set_catchment_parameter",&M::set_catchment_parameter,(py::arg("self"),py::arg("catchment_id"),py::arg("p")),
                    "creates/modifies a pr catchment override parameter\n"
                    "param catchment_id the 0 based catchment_id that correlates to the cells catchment_id\n"
                    "param a reference to the parameter that will be kept for those cells\n"
         )
         .def("remove_catchment_parameter",&M::remove_catchment_parameter,(py::arg("self"),py::arg("catchment_id")),"remove a catchment specific parameter override, if it exists.")
         .def("has_catchment_parameter",&M::has_catchment_parameter, (py::arg("self"), py::arg("catchment_id")),"returns true if there exist a specific parameter override for the specified 0-based catchment_id")

			 .def("get_catchment_parameter",&M::get_catchment_parameter, (py::arg("self"), py::arg("catchment_id")),
                    "return the parameter valid for specified catchment_id, or global parameter if not found.\n"
                    "note Be aware that if you change the returned parameter, it will affect the related cells.\n"
                    "param catchment_id 0 based catchment id as placed on each cell\n"
                    "returns reference to the real parameter structure for the catchment_id if exists,\n"
                    "otherwise the global parameters\n"
         ,return_internal_reference<>())
		.def("set_catchment_calculation_filter",&M::set_catchment_calculation_filter,(py::arg("self"),py::arg("catchment_id_list")),
                    "set/reset the catchment based calculation filter. This affects what get simulate/calculated during\n"
                    "the run command. Pass an empty list to reset/clear the filter (i.e. no filter).\n"
                    "\n"
                    "param catchment_id_list is a catchment id vector\n"
         )
		 .def("set_calculation_filter", &M::set_calculation_filter, (py::arg("self"), py::arg("catchment_id_list"), py::arg("river_id_list")),
                 "set/reset the catchment *and* river based calculation filter. This affects what get simulate/calculated during\n"
                 "the run command. Pass an empty list to reset/clear the filter (i.e. no filter).\n"
                 "\n"
                 "param catchment_id_list is a catchment id vector\n"
                "param river_id_list is a river id vector\n"
         )
         .def("is_calculated",&M::is_calculated,(py::arg("self"),py::arg("catchment_id")),"true if catchment id is calculated during runs, ref set_catchment_calculation_filter")
         .def("get_states",&M::get_states,(py::arg("self"),py::arg("end_states")),
                    "collects current state from all the cells\n"
                    "note that catchment filter can influence which states are calculated/updated.\n"
                    "param end_states a reference to the vector<state_t> that are filled with cell state, in order of appearance.\n"
        )

        .def("set_states",&M::set_states,(py::arg("self"), py::arg("states")),
                    "set current state for all the cells in the model.\n"
                    "states is a vector<state_t> of all states, must match size/order of cells.\n"
                    "note throws runtime-error if states.size is different from cells.size\n"
        )
        .def("revert_to_initial_state",&M::revert_to_initial_state,(py::arg("self")),
             "Given that the cell initial_states are established, these are \n"
             "copied back into the cells\n"
             "Note that the cell initial_states vector is established at the first call to \n"
             ".set_states() or run_cells()\n"
             )
        .def("set_state_collection",&M::set_state_collection,(py::arg("self"), py::arg("catchment_id"),py::arg("on_or_off")),
                    "enable state collection for specified or all cells\n"
                    "note that this only works if the underlying cell is configured to\n"
                    "do state collection. This is typically not the  case for\n"
                    "cell-types that are used during calibration/optimization\n"
        )
        .def("set_snow_sca_swe_collection",&M::set_snow_sca_swe_collection,(py::arg("self"), py::arg("catchment_id"),py::arg("on_or_off")),
                    "enable/disable collection of snow sca|sca for calibration purposes\n"
                    "param cachment_id to enable snow calibration for, -1 means turn on/off for all\n"
                    "param on_or_off true|or false.\n"
                    "note if the underlying cell do not support snow sca|swe collection, this \n"
        )
		.def("adjust_q",&M::adjust_q,(py::arg("self"),py::arg("q_scale"),py::arg("cids")),
            doc_intro("adjust the current state content q of ground storage by scale-factor")
            doc_intro("")
            doc_intro("Adjust the content of the ground storage, e.g. state.kirchner.q, or")
            doc_intro("hbv state.(tank|soil).(uz,lz|sm), by the specified scale factor.")
            doc_intro("The this function plays key role for adjusting the state to")
            doc_intro("achieve a specified/wanted average discharge flow output for the")
            doc_intro("model at the first time-step.")
            doc_parameters()
            doc_parameter("q_scale","float","the scale factor to apply to current storage state")
            doc_parameter("cids","IntVector","if empty, all cells are in scope, otherwise only cells that have specified catchment ids.")
        )

        .def("adjust_state_to_target_flow",&M::adjust_state_to_target_flow,(py::arg("self"),py::arg("wanted_flow_m3s"),py::arg("cids"),py::arg("start_step")=0,
            py::arg("scale_range")=10.0,py::arg("scale_eps")=1.0e-3,py::arg("max_iter")=300,py::arg("n_steps")=1
        ),
             doc_intro("state adjustment to achieve wanted/observed flow")
			 doc_intro("")
			 doc_intro("This function provides an easy and consistent way to adjust the")
			 doc_intro("state of the cells(kirchner, or hbv-tank-levels) so that the average output")
			 doc_intro("from next n_steps time-steps matches the wanted flow for the same period.")
			 doc_intro("")
			 doc_intro("This is quite complex, since the amount of adjustment needed is dependent of the")
			 doc_intro("cell-state, temperature/precipitation in time-step, glacier-melt, length of the time-step,")
			 doc_intro("and calibration factors sensitivity.")
			 doc_intro("")
			 doc_intro("The approach here is to use dlib::find_min_single_variable to solve")
			 doc_intro("the problem, instead of trying to reverse compute the needed state.")
			 doc_intro("")
			 doc_intro("This has several benefits, it deals with the full stack and state, and it can be made")
			 doc_intro("method stack independent.")
			 doc_intro("")
			 doc_intro("Notice that the model should be prepared for run prior to calling this function")
			 doc_intro("and that there should be a current model state that gives the starting point")
			 doc_intro("for the adjustment.")
			 doc_intro("Also note that when returning, the active state reflects the")
			 doc_intro("achieved flow returned, and that the current state  for the cells")
			 doc_intro ("belonging to the catchment-ids is modified as needed to provide this average-flow.")
			 doc_intro("The state when returning is set to the start of the i'th period specified")
			 doc_intro("to reach the desired flow.")
			 doc_intro("")
             doc_parameters()
			 doc_parameter("wanted_flow_m3s","float","the average flow first time-step we want to achieve")
			 doc_parameter("cids","IntVector"," catchments, represented by catchment-ids that should be adjusted")
             doc_parameter("start_step","int","what time-step number in the time-axis to use, default 0")
             doc_parameter("scale_range","float","optimizer boundaries is s_0/scale_range .. s_0*scale_range, s_0=wanted_flow_m3s/q_0 , default =10.0")
             doc_parameter("scale_eps","float","optimizer eps, stop criteria (ref. dlib), eps=s_0*scale_eps , default =1-e3")
             doc_parameter("max_iter","int","optimizer max evaluations before giving up to find optimal solution")
             doc_parameter("n_steps","int","number of time-steps in the time-axis to average the to the wanted_flow_m3s, default=1")
			 doc_returns("obtained flow in m3/s units.","FlowAdjustResult","note: this can deviate from wanted flow due to model and state constraints")
        )
        .def("get_cells",&M::get_cells, (py::arg("self")),"cells as shared_ptr<vector<cell_t>>")
        .def("size",&M::size,(py::arg("self")),"return number of cells")
        .add_property("cells",&M::get_cells,"cells of the model")
        .add_property("current_state",&M::current_state," a copy of the current model state")
        .def("is_cell_env_ts_ok",&M::is_cell_env_ts_ok,(py::arg("self")),
             doc_intro("Use this function after the interpolation step, before .run_cells(), to verify")
             doc_intro("that all cells selected for computation (calculation_filter), do have ")
             doc_intro("valid values.")
             doc_returns("all_ok","bool","return false if any nan is found, otherwise true")
         )
         ;

    }

    template <class F, class O>
    O clone_to_opt_impl(F const& f,bool with_catchment_params) {
        O o(f.extract_geo_cell_data(), f.get_region_parameter());
        o.time_axis = f.time_axis;
        o.ip_parameter = f.ip_parameter;
        o.region_env = f.region_env;
        o.initial_state = f.initial_state;
        o.river_network = f.river_network;
        auto fc = f.get_cells();
        auto oc = o.get_cells();
        for (size_t i = 0;i < f.size();++i) {
            (*oc)[i].env_ts = (*fc)[i].env_ts;
            (*oc)[i].state = (*fc)[i].state;
        }
        if(with_catchment_params) {
            auto cids=f.catchment_ids();
            for(const auto& cid:cids) {
                if(f.has_catchment_parameter(cid))
                    o.set_catchment_parameter(cid, f.get_catchment_parameter(cid));
            }
        }
        return o;
    }

    template <typename F, typename O>
    void def_clone_to_similar_model(const char *func_name) {
        auto pfi = &clone_to_opt_impl< F, O>;
        def(func_name, pfi, (py::arg("src_model"),py::arg("with_catchment_params")=false),
            doc_intro("Clone a model to a another similar type model, full to opt-model or vice-versa")
            doc_intro("The entire state except catchment-specific parameters, filter and result-series are cloned")
            doc_intro("The returned model is ready to run_cells(), state and interpolated enviroment is identical to the clone source")
            doc_parameters()
            doc_parameter("src_model","XXXX?Model","The model to be cloned, with state interpolation done, etc")
            doc_parameter("with_catchment_params","bool","default false, if true also copy catchment specific parameters")
            doc_returns("new_model","XXXX?Model","new_model ready to run_cells, or to put into the calibrator/optimizer")
        );
    }



    template<class RegionModel>
    static void
    model_calibrator(const char *optimizer_name) {

        typedef typename RegionModel::parameter_t parameter_t;
        typedef shyft::time_series::dd::apoint_ts pts_t;
        typedef shyft::core::model_calibration::optimizer<RegionModel, parameter_t, pts_t> Optimizer;
        typedef typename Optimizer::target_specification_t target_specification_t;

        // fix overloads mapping vs. vector& new parameter stuff
        std::vector<double>(Optimizer::*optimize_v)(const std::vector<double>&, size_t, double, double ) = &Optimizer::optimize;
        parameter_t(Optimizer::*optimize_p)(const parameter_t&, size_t, double, double) = &Optimizer::optimize;

        parameter_t(Optimizer::*optimize_global_p)(const parameter_t&, size_t, double, double) = &Optimizer::optimize_global;

        std::vector<double>(Optimizer::*optimize_dream_v)(const std::vector<double>&, size_t) = &Optimizer::optimize_dream;
        parameter_t(Optimizer::*optimize_dream_p)(const parameter_t&, size_t) = &Optimizer::optimize_dream;

        std::vector<double> (Optimizer::*optimize_sceua_v)(const std::vector<double>&,size_t,double,double)=&Optimizer::optimize_sceua;
        parameter_t(Optimizer::*optimize_sceua_p)(const parameter_t&, size_t, double, double) = &Optimizer::optimize_sceua;

        double (Optimizer::*calculate_goal_function_v)(const std::vector<double>&) = &Optimizer::calculate_goal_function;
        double (Optimizer::*calculate_goal_function_p)(const parameter_t&) = &Optimizer::calculate_goal_function;



        class_<Optimizer>(optimizer_name,
            doc_intro(
                "The optimizer for parameters for a region model\n"
                "It provides needed functionality to orchestrate a search for the optimal parameters so that the goal function\n"
                "specified by the target_specifications are minimized.\n"
                "The user can specify which parameters (model specific) to optimize, giving range min..max for each of the\n"
                "parameters. Only parameters with min != max are used, thus minimizing the parameter space.\n"
                "\n"
                "Target specification ref: TargetSpecificationVector allows a lot of flexibility when it comes to what\n"
                "goes into the goal-function.\n"
                "\n"
                "This class provides several goal-function search algorithms:\n"
                "    .optimize               min-bobyqa  a fast local optimizer, http://dlib.net/optimization.html#find_min_bobyqa\n"
                "    .optimize_global   a global optimizer, http://dlib.net/optimization.html#global_function_search\n"
                "    .optimize_sceua   a global optimizer,  https://www.sciencedirect.com/science/article/pii/0022169494900574\n"
                "    .optimize_dream  a global optimizer,\n"
                "                                                            Theory is found in: Vrugt, J. et al: Accelerating Markov Chain Monte Carlo\n"
                "                                                            simulations by Differential Evolution with Self-Adaptive Randomized Subspace\n"
                "                                                            Sampling. Int. J. of Nonlinear Sciences and Numerical Simulation 10(3) 2009.\n"
                "\n\n"
                "Each method searches for the optimum parameter-set, given the input-constraints and time-limit, max_iterations and accuracy(method dependent).\n"
                "Also note that after the optimization, you have a complete trace of the parameter-search with the corresponding goal-function value\n"
                "This enable you to analyze the search-function, and allows you to select other parameter-sets that based on \n"
                "hydrological criterias that is not captured in the goal-function specification\n"
            )
            ,no_init
        )
        .def(init< RegionModel&,
                   const std::vector<target_specification_t>&,
                   const std::vector<double>&,
                   const std::vector<double>& >(
                   (py::arg("self"), py::arg("model"), py::arg("targets"), py::arg("p_min"), py::arg("p_max")),
                    doc_intro(
                        "Construct an optimizer for the specified region model.\n"
                        "Set  p_min.param.x = p_max.param.x  to disable optimization for a parameter param.x\n"
                    )
                    doc_parameters()
                    doc_parameter("model","OptModel","the model to be optimized, the model should be initialized, interpolation/preparation  step done" )
                    doc_parameter("targets","TargetSpecificationVector","specifies how to calculate the goal-function")
                    doc_parameter("p_min","Parameter","minimum values for the parameters to be optimized")
                    doc_parameter("p_max","Parameter","maximum values for the parameters to be optimized")
                  )
        )
        .def(init<RegionModel&>(py::args("model"),
            doc_intro("Construct a parameter Optimizer for the supplied model\n"
            "Use method .set_target_specification(...) to provide the target specification,\n"
            "then invoke opt_param= o.optimize(p_starting_point..)\n"
            "to get back the optimized parameters for the supplied model and target-specification\n"
            )
            doc_parameters()
            doc_parameter("model","OptModel","the model to be optimized, the model should be initialized, interpolation/preparation  step done" )
            )
        )
        .def("set_target_specification",&Optimizer::set_target_specification,
             (py::arg("self"),py::arg("target_specification"),py::arg("parameter_lower_bound"),py::arg("parameter_upper_bound")),
            doc_intro(
                "Set the target specification, parameter lower and upper bound to be used during \n"
                "subsequent call to the .optimize() methods.\n"
                "Only parameters with lower_bound != upper_bound will be subject to optimization\n"
                "The object properties target_specification,lower and upper bound are updated and\n"
                "will reflect the current setting.\n"
            )
            doc_parameters()
            doc_parameter("target_specification","TargetSpecificationVector","the complete target specification composition of one or more criteria")
            doc_parameter("parameter_lower_bound", "Parameter", "the lower bounds of the parameters")
            doc_parameter("parameter_upper_bound", "Parameter", "the upper bounds of the parameters")
        )
        .def("establish_initial_state_from_model", &Optimizer::establish_initial_state_from_model,
             doc_intro(
                "Copies the Optimizer referenced region-model current state\n"
                "to a private store in the Optimizer object.\n"
                "This state is used to for restore prior to each run of the model during calibration\n"
                "notice that if you forget to call this method, it will be called automatically once you\n"
                "call one of the optimize methods.\n"
             )
        )
        .def("get_initial_state",&Optimizer::get_initial_state,(py::arg("self"),py::arg("i")),"returns a copy of the i'th cells initial state")
        .def("optimize",optimize_v,args("p","max_n_evaluations","tr_start","tr_stop"),
                "(deprecated)Call to optimize model, starting with p parameter set, using p_min..p_max as boundaries.\n"
                "where p is the full parameter vector.\n"
                "the p_min,p_max specified in constructor is used to reduce the parameterspace for the optimizer\n"
                "down to a minimum number to facilitate fast run.\n"
                "param p contains the starting point for the parameters\n"
                "param max_n_evaluations stop after n calls of the objective functions, i.e. simulations.\n"
                "param tr_start is the trust region start , default 0.1, ref bobyqa\n"
                "param tr_stop is the trust region stop, default 1e-5, ref bobyqa\n"
                "return the optimized parameter vector\n"
        )
        .def("optimize", optimize_p, (py::arg("self"), py::arg("p"),py::arg( "max_n_evaluations"),py::arg( "tr_start"), py::arg("tr_stop")),
            doc_intro(
                "Call to optimize model, using find_min_bobyqa,  starting with p parameters\n"
                "as the start point\n"
                "The current target specification, parameter lower and upper bound\n"
                "is taken into account\n"
            )
            doc_parameters()
            doc_parameter("p","Parameter","contains the starting point for the parameters")
            doc_parameter("max_n_evaluations","int","stop after n calls of the objective functions, i.e. simulations.")
            doc_parameter("tr_start","float","minbobyqa is the trust region start , default 0.1, ref bobyqa")
            doc_parameter("tr_stop","float"," is the trust region stop, default 1e-5, ref bobyqa")
            doc_returns("p_opt","Parameter","the the optimized parameters")
        )
        .def("optimize_global", optimize_global_p,(py::arg("self"),py::arg("p"),py::arg("max_n_evaluations"),py::arg("max_seconds"),py::arg("solver_eps")),
            doc_intro("Finds the global optimum parameters for the model.\n"
                "The current target specification, parameter lower and upper bound\n"
                "is taken into account\n"
                ".. refer to _dlib_global_search:\n"
                " http://dlib.net/optimization.html#global_function_search\n"
            )
            doc_parameters()
            doc_parameter("p","Parameter", "the potential starting point for the global search(currently not used by dlib impl)")
            doc_parameter("max_n_evaluations","int","stop after n calls of the objective functions, i.e. simulations.")
            doc_parameter("max_seconds","float", "stop search for for solution after specified time-limit")
            doc_parameter("solver_eps","float","search for minimum goal-function value at this accuracy, continue search for possibly other global minima when this accuracy is reached.")
            doc_returns("p_opt","Parameter","the optimal found minima given the inputs")
        )

        .def("optimize_dream",optimize_dream_v,(py::arg("self"), py::arg("p"),py::arg("max_n_evaluations")),
                doc_intro(
                    "(Deprecated)Call to optimize model, using DREAM alg., find p, using p_min..p_max as boundaries.\n"
                    "where p is the full parameter vector.\n"
                    "the p_min,p_max specified in constructor is used to reduce the parameterspace for the optimizer\n"
                    "down to a minimum number to facilitate fast run.\n"
                    "param p is used as start point (not really, DREAM use random, but we should be able to pass u and q....\n"
                    "param max_n_evaluations stop after n calls of the objective functions, i.e. simulations.\n"
                    "return the optimized parameter vector\n"
                )
        )
        .def("optimize_dream", optimize_dream_p, (py::arg("self"), py::arg("p"),py::arg("max_n_evaluations")),
            doc_intro(
                "Call to optimize model with the DREAM algorithm.\n"
                "The supplied p is ignored (DREAM selects starting point randomly)\n"
                "The current target specification, parameter lower and upper bound\n"
                "is taken into account\n"
            )
            doc_parameters()
            doc_parameter("p","Parameter", "the potential starting point for the global search(currently not used by dlib impl)")
            doc_parameter("max_n_evaluations","int","stop after n calls of the objective functions, i.e. simulations.")
            doc_returns("p_opt","Parameter","the optimal found minima given the inputs")
        )
        .def("optimize_sceua",optimize_sceua_v,(py::arg("self"), py::arg("p"),py::arg("max_n_evaluations"),py::arg("x_eps"),py::arg("y_eps")),
               doc_intro(
                    "(Deprecated)Call to optimize model, using SCE UA, using p as startpoint, find p, using p_min..p_max as boundaries.\n"
                    "where p is the full parameter vector.\n"
                    "the p_min,p_max specified in constructor is used to reduce the parameter-space for the optimizer\n"
                    "down to a minimum number to facilitate fast run.\n"
                    "param p is used as start point and is updated with the found optimal points\n"
                    "param max_n_evaluations stop after n calls of the objective functions, i.e. simulations.\n"
                    "param x_eps is stop condition when all changes in x's are within this range\n"
                    "param y_eps is stop condition, and search is stopped when goal function does not improve anymore within this range\n"
                    "return the optimized parameter vector\n"
               )
        )
        .def("optimize_sceua", optimize_sceua_p,(py::arg("self"), py::arg("p"),py::arg("max_n_evaluations"),py::arg("x_eps"),py::arg("y_eps")),
            doc_intro(
                "Call to optimize model using SCE UA algorithm, starting with p parameters\n"
                "as the start point\n"
                "The current target specification, parameter lower and upper bound\n"
                "is taken into account\n"
            )
            doc_parameters()
            doc_parameter("p","Parameter", "the potential starting point for the global search")
            doc_parameter("max_n_evaluations","int","stop after n calls of the objective functions, i.e. simulations.")
            doc_parameter("x_eps","float","is stop condition when all changes in x's are within this range")
            doc_parameter("y_eps","float","is stop condition, and search is stopped when goal function does not improve anymore within this range")
            doc_returns("p_opt","Parameter","the optimal found minima given the inputs")
        )

        .def("reset_states",&Optimizer::reset_states,(py::arg("self")),"reset the state of the model to the initial state before starting the run/optimize")
        .def("set_parameter_ranges",&Optimizer::set_parameter_ranges,(py::arg("self"),py::arg("p_min"), py::arg("p_max")),
             doc_intro(
                 "Set the parameter ranges for the optimization search.\n"
                " Set min=max=wanted parameter value for those not subject to change during optimization\n"
                " - changes/sets the parameter_lower_bound.. paramter_upper_bound as specified in constructor\n"
             )
             doc_parameters()
             doc_parameter("p_min", "Parameter", "the lower bounds of the parameters")
             doc_parameter("p_max", "Parameter", "the upper bounds of the parameters")
        )
        .def("set_verbose_level",&Optimizer::set_verbose_level,(py::arg("self"),py::arg("level")),"set verbose level on stdout during calibration,0 is silent,1 is more etc.")
        .def("calculate_goal_function",calculate_goal_function_v,(py::arg("self"),py::arg("full_vector_of_parameters")),
                doc_intro(
                    "(Deprecated)calculate the goal_function as used by minbobyqa,etc.,\n"
                    "using the full set of  parameters vectors (as passed to optimize())\n"
                    "and also ensures that the shyft state/cell/catchment result is consistent\n"
                    "with the passed parameters passed\n"
                    "param full_vector_of_parameters contains all parameters that will be applied to the run.\n"
                    "returns the goal-function, weigthed nash_sutcliffe|Kling-Gupta sum \n"
                )
        )
        .def("calculate_goal_function", calculate_goal_function_p, (py::arg("self"), py::arg("parameters")),
            doc_intro(
                "Calculate the goal_function as used by minbobyqa,etc.,\n"
                "using the supplied set of parameters\n"
                "and also ensures that the shyft state/cell/catchment result is consistent\n"
                "with the passed parameters passed\n"
                "param parameters contains all parameters that will be applied to the run.\n"
                "You can also use this function to build your own external supplied optimizer in python"
            )
            doc_parameters()
            doc_parameter("parameters", "Parameter", "the region model parameter to use when evaluating the goal-function")
            doc_returns("goal_function_value", "float", "the goal-function, weigthed nash_sutcliffe|Kling-Gupta sum etc. value ")
        )
        .def_readwrite("target_specification",&Optimizer::targets,
           doc_intro("The current target-specifications used during optimization")
        )
        .def_readwrite("parameter_lower_bound",&Optimizer::parameter_lower_bound,"the lower bound parameters\n")
        .def_readwrite("parameter_upper_bound",&Optimizer::parameter_upper_bound,"the upper bound parameters\n")
        .def("parameter_active",&Optimizer::active_parameter,(py::arg("self"), py::arg("i")),
            doc_intro("returns true if the i'th parameter is active, i.e. lower != upper bound\n")
            doc_parameters()
            doc_parameter("i","int","the index of the parameter")
            doc_returns("active","bool","True if the parameter abs(p[i].min -p[i].max)> zero_limit")
        )
        .def_readonly("trace_size",&Optimizer::trace_size,
            doc_intro("returns the size of the parameter-trace")
            doc_see_also("trace_goal_function_value,trace_parameter")
        )
        .def_readonly("trace_goal_function_values",&Optimizer::goal_fn_trace,
            doc_intro("the goal-function values in the order of searching for the minimum value")
            doc_intro("The trace_parameter(i) gives the corresponding i'th parameter")
            doc_see_also("trace_parameter,trace_value,trace_size")
        )
        .def("trace_goal_function_value",&Optimizer::trace_goal_fn,(py::arg("self"),py::args("i")),
            doc_intro("returns the i'th goal function value")
        )
        .def("trace_parameter",&Optimizer::trace_parameter,(py::arg("self"),py::args("i")),
            doc_intro("returns the i'th parameter tried, corresponding to the ")
            doc_intro("i'th trace_goal_function value")
            doc_see_also("trace_goal_function,trace_size")
        )
        ;

    }
}
