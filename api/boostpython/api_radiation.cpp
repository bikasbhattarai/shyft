/** This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
See file COPYING for more details **/
#include "boostpython_pch.h"

#include "core/radiation.h"
#include <armadillo>

namespace expose {

    void radiation() {
        using namespace shyft::core::radiation;
        using namespace boost::python;
        namespace py = boost::python;

        class_<parameter>("RadiationParameter")
            .def(init<double,double>(args("albedo","turbidity"),"a new object with specified parameters"))
            .def_readwrite("albedo",&parameter::albedo,"typical value 0.2")
            .def_readwrite("turbidity", &parameter::turbidity,"typical value 1.0")
            ;

        class_<response>("RadiationResponse")
            .def_readwrite("sw_radiation",&response::sw_radiation)
            .def_readwrite("lw_radiation",&response::lw_radiation)
            .def_readwrite("net_radiation",&response::net_radiation)
            .def_readwrite("ra",&response::ra)
            .def_readwrite("rah",&response::rah)
            .def_readwrite("omega1",&response::omega1)
            .def_readwrite("omega2",&response::omega2)
            ;

        typedef calculator<parameter,response> RadiationCalculator;
        class_<RadiationCalculator>("RadiationCalculator",
                "Radiation,R, (ref.: Allen, R. G.; Trezza, R. & Tasumi, M. Analytical integrated functions for daily solar radiation on slopes Agricultural and Forest Meteorology, 2006, 139, 55-73)\n"
                "primitive implementation for calculating predicted clear-sky short-wave solar radiation for inclined surfaces\n"
                "This function is plain and simple, taking albedo and turbidity\n"
                "into the constructor and provides 2 functions:\n"
                " psw_radiation calculates predicted solar radiation (if no measured data available);\n"
                " tsw_radiation translates measured horizontal radiation into sloped surface\n"
                "[mm/s] units.\n",no_init
            )
            .def(init<const parameter&>(args("param"),"create a calculator using supplied parameter"))
            .def("net_radiation",&RadiationCalculator::net_radiation,(py::arg("self"),py::arg("response"), py::arg("latitude"), py::arg("t"), py::arg("slope"), py::arg("aspect"), py::arg("temperature"), py::arg("rhumidity"), py::arg("elevation"),  py::arg("rsm")),
                     doc_intro("calculates net radiation, updating response")
            .def("net_radiation_step",&RadiationCalculator::net_radiation_step,(py::arg("self"),py::arg("response"), py::arg("latitude"), py::arg("t1"),py::arg("t2"), py::arg("slope"), py::arg("aspect"), py::arg("temperature"), py::arg("rhumidity"), py::arg("elevation")),
                     doc_intro("calculates net radiation, updating response")
            )
            ;
    }
}

BOOST_PYTHON_MODULE(_radiation)
{

    boost::python::scope().attr("__doc__")="Shyft python api for the radiation model";
    boost::python::docstring_options doc_options(true, true, false);// all except c++ signatures
    expose::radiation();
}
