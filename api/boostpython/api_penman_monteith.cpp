/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#include "boostpython_pch.h"

#include "core/penman_monteith.h"
#include <armadillo>

namespace expose {

    void penman_monteith() {
        using namespace shyft::core::penman_monteith;
        using namespace boost::python;
        namespace py = boost::python;

        class_<parameter>("PenmanMonteithParameter")
            .def(init<double,double,double>(args("lai","height_ws","height_t"),"a new object with specified parameters"))
            .def_readwrite("lai",&parameter::lai,"typical value 2.0")
            .def_readwrite("height_ws", &parameter::height_ws,"typical value 2.0")
            .def_readwrite("height_t", &parameter::height_t,"typical value 2.0")
            .def_readwrite("height_veg", &parameter::height_t,"grass 0.15")
            ;

        class_<response>("PenmanMonteithResponse")
            .def_readwrite("et_ref",&response::et_ref)
            ;

        typedef calculator<parameter,response> PenmanMonteithCalculator;
        class_<PenmanMonteithCalculator>("PenmanMonteithCalculator",
                "Evapotranspiration model, Penman-Monteith equation, PM \n"
                "(ref.: ASCE-EWRI The ASCE Standardized Reference Evapotranspiration Equation Environmental and\n"
                "         * Water Resources Institute (EWRI) of the American Society of Civil Engineers Task Com- mittee on Standardization of\n"
                "         * Reference Evapotranspiration Calculation, ASCE, Washington, DC, Environmental and Water Resources Institute (EWRI) of\n"
                "         * the American Society of Civil Engineers Task Com- mittee on Standardization of Reference Evapotranspiration Calculation,\n"
                "         * ASCE, Washington, DC, 2005\n"
                "calculating reference evapotranspiration\n"
                "This function is plain and simple, taking albedo and turbidity\n"
                "into the constructor and provides 2 functions:\n"
                " reference_evapotranspiration_asce;\n"
                " reference_evapotranspiration_st;\n"
                "[mm/s] units.\n",no_init
            )
            .def(init<const parameter&>(args("param"),"create a calculator using supplied parameter"))
            .def("reference_evapotranspiration_asce",&PenmanMonteithCalculator::reference_evapotranspiration_asce,(py::arg("self"),py::arg("response"), py::arg("net_radiation"), py::arg("temperature"), py::arg("rhumidity"),  py::arg("elevation"),  py::arg("windspeed")),
                     doc_intro("calculates reference evapotranspiration, updating response")
            )
            .def("reference_evapotranspiration_st",&PenmanMonteithCalculator::reference_evapotranspiration_st,(py::arg("self"),py::arg("response"), py::arg("net_radiation"), py::arg("temperature"), py::arg("rhumidity"),  py::arg("elevation"),  py::arg("windspeed")),
                    doc_intro("calculates reference evapotranspiration, updating response")
                    )
            ;
    }
}

BOOST_PYTHON_MODULE(_penman_monteith)
{

    boost::python::scope().attr("__doc__")="Shyft python api for the penman-monteith evapotranspiration model";
    boost::python::docstring_options doc_options(true, true, false);// all except c++ signatures
    expose::penman_monteith();
}
