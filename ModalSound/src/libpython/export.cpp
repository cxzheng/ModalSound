#ifndef LIB_PYTHON_EXPORT_INC
#   define LIB_PYTHON_EXPORT_INC

#include <boost/python.hpp>

#include "geometry/Point3.hpp"
#include "modal/ModalShape.h"
#include "transfer/BoundaryIntegralTransfer.hpp"
#include "transfer/FMMBoundIntTransfer.hpp"

namespace bp = boost::python;


void export_tuple3()
{
    bp::class_< Tuple3d >("Tuple3d", bp::init<>())
        .def(bp::init<double, double, double>())
        .def("norm", &Tuple3d::norm)
        .def("norm_sqr", &Tuple3d::norm_sqr);

    bp::class_< Point3d >("Point3d", bp::init<>())
        .def(bp::init<double, double, double>())
        .def_readwrite("x", &Point3d::x)
        .def_readwrite("y", &Point3d::y)
        .def_readwrite("z", &Point3d::z)
        .def("distance", &Point3d::distance)
        .def("distance_sqr", &Point3d::distance_sqr)
        .def(bp::self += double())
        .def(bp::self *= double());
}

void export_modal_shape()
{
    // map the modal namespace to a sub-module
    // make "from sploosh.modal import <whatever>" work
    bp::object modalModule(bp::handle<>(bp::borrowed(PyImport_AddModule("sploosh.modal"))));

    // make "from sploosh import modal" work
    bp::scope().attr("modal") = modalModule;
    // set the current scope to the new sub-module
    bp::scope modal_scope = modalModule;

    bp::class_<ModalShape>("ModalShape", bp::init<const char*>())
        .add_property("num_modes", &ModalShape::num_modes)
        .def("eigen_mode", &ModalShape::eigen_mode);
}

void export_boundary_integral_transfer()
{
    // map the transfer namespace to a sub-module
    // make "from sploosh.transfer import <whatever>" work
    bp::object transferModule(bp::handle<>(bp::borrowed(PyImport_AddModule("sploosh.transfer"))));

    // make "from sploosh import transfer" work
    bp::scope().attr("transfer") = transferModule;
    // set the current scope to the new sub-module
    bp::scope transfer_scope = transferModule;

    typedef BoundaryIntegralTransfer<double> BITd;
    bp::class_< BITd >("BoundaryIntegralTransfer", bp::init<const char*, const char*>())
        .add_property("wave_number", &BITd::wave_number)
        .add_property("num_elements", &BITd::num_elements)
        .def("eval", &BITd::eval);

    typedef FMMBoundIntTransfer<double> FBITd;
    bp::class_< FBITd, bp::bases<BITd> >("FMMBoundIntTransfer", bp::init<const char*, const char*>())
        //.add_property("num_expan", &FBITd::expansion_num, &FBITd::set_expansion_num)
        .add_property("num_expan", &FBITd::expansion_num)
        .def("compute_moments", &FBITd::compute_moments)
        .def("store_moments", &FBITd::store_moments)
        .def("eval", &FBITd::eval);
}

//////////////////////////////////////////////////////////////////////////

BOOST_PYTHON_MODULE(sploosh)
{
    export_tuple3();
    export_boundary_integral_transfer();
    export_modal_shape();
}

#endif
