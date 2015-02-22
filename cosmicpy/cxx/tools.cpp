// Copyright (c) 2014-2015, CosmicPy Developers
// Licensed under CeCILL 2.1 - see LICENSE.rst
#include "besselwindow.h"
using namespace boost::python;

BOOST_PYTHON_MODULE(tools)
{
	np::initialize();
	
	class_< BesselWindow >("BesselWindow", init<int64_t, int64_t, double, char* >())
	  .def("rgrid", &BesselWindow::rgrid)
	  .def("kgrid", &BesselWindow::kgrid)
	  .def("custom", &BesselWindow::custom)
	  .def("optimal", &BesselWindow::optimal)
	  .def("tabulated", &BesselWindow::tabulated)
	  .add_property("npoints", &BesselWindow::getNpoints, &BesselWindow::setNpoints)
	  .add_property("precision", &BesselWindow::getPrecision, &BesselWindow::setPrecision);
}