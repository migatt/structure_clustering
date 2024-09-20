#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Atom.hpp"
#include "Machine.hpp"
#include "Result.hpp"
#include "Structure.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {

    py::class_<Machine>(m, "Machine")
        .def(py::init<>())
        .def("setOnlyConnectedGraphs", &Machine::setOnlyConnectedGraphs)
        .def("setCovalentRadius", &Machine::setCovalentRadius)
        .def("addPairDistance", &Machine::addPairDistance)
        .def("cluster", &Machine::cluster);

    py::class_<Structure>(m, "Structure")
        .def(py::init<int>())
        .def("addAtom", &Structure::addAtom)
        .def_property_readonly("numConnections", &Structure::getNumConnections);

    py::class_<Atom>(m, "Atom").def(
        py::init<const int, const double, const double, const double>());

    py::class_<Result>(m, "Result")
        .def("export", &Result::exportDat)
        .def_property_readonly("clusters", &Result::getClusters)
        .def_property_readonly("singles", &Result::getSingles)
        .def_property_readonly("structures", &Result::getStructures);

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
