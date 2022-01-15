#include <string>

#include <Eigen/Eigen>

#include <cork/rawmesh/rawMesh.h>
#include <cork/mesh/mesh.h>
#include <cork/cork.h>

#include <pybind11/pybind11.h>

#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include <pybind11/complex.h>
#include <pybind11/eigen.h>

#include <tuple>

namespace py = pybind11;

namespace pycork {

typedef Eigen::Matrix<double,Eigen::Dynamic,3> EigenVecX3d;
typedef Eigen::Matrix<uint64_t,Eigen::Dynamic,3> EigenVecX3i;
typedef std::tuple<EigenVecX3d, EigenVecX3i> MeshTuple;

bool isSolid(const EigenVecX3d &verts,
             const EigenVecX3i &tris) {

    CorkMesh mesh;

    eigenToCorkMesh(verts, tris, &mesh);

    bool result = true;

    if(mesh.isSelfIntersecting())
        result = false;

    if(!mesh.isClosed())
        result = false;

    return result;

}


MeshTuple booleanUnion(const EigenVecX3d &vertsA,
                       const EigenVecX3i &trisA,
                       const EigenVecX3d &vertsB,
                       const EigenVecX3i &trisB) {

    CorkMesh meshA, meshB;

    eigenToCorkMesh(vertsA, trisA, &meshA);
    eigenToCorkMesh(vertsB, trisB, &meshB);

    MeshTuple meshOut;

    meshA.boolUnion(meshB);

    corkMesh2Eigen(meshA, std::get<0>(meshOut), std::get<1>(meshOut));

    return meshOut;
}

MeshTuple booleanDifference(const EigenVecX3d &vertsA,
                            const EigenVecX3i &trisA,
                            const EigenVecX3d &vertsB,
                            const EigenVecX3i &trisB) {

    CorkMesh meshA, meshB;

    eigenToCorkMesh(vertsA, trisA, &meshA);
    eigenToCorkMesh(vertsB, trisB, &meshB);

    MeshTuple meshOut;

    meshA.boolDiff(meshB);

    corkMesh2Eigen(meshA, std::get<0>(meshOut), std::get<1>(meshOut));

    return meshOut;
}

MeshTuple booleanIntersection(const EigenVecX3d &vertsA,
                              const EigenVecX3i &trisA,
                              const EigenVecX3d &vertsB,
                              const EigenVecX3i &trisB) {

    CorkMesh meshA, meshB;

    eigenToCorkMesh(vertsA, trisA, &meshA);
    eigenToCorkMesh(vertsB, trisB, &meshB);

    MeshTuple meshOut;

    meshA.boolIsct(meshB);

    corkMesh2Eigen(meshA, std::get<0>(meshOut), std::get<1>(meshOut));

    return meshOut;
}


MeshTuple booleanXor(const EigenVecX3d &vertsA,
                     const EigenVecX3i &trisA,
                     const EigenVecX3d &vertsB,
                     const EigenVecX3i &trisB) {

    CorkMesh meshA, meshB;

    eigenToCorkMesh(vertsA, trisA, &meshA);
    eigenToCorkMesh(vertsB, trisB, &meshB);

    MeshTuple meshOut;

    meshA.boolXor(meshB);

    corkMesh2Eigen(meshA, std::get<0>(meshOut), std::get<1>(meshOut));

    return meshOut;
}


MeshTuple resolveIntersection(const EigenVecX3d &vertsA,
                              const EigenVecX3i &trisA) {

    CorkMesh meshA;

    eigenToCorkMesh(vertsA, trisA, &meshA);

    MeshTuple meshOut;

    meshA.resolveIntersections();

    corkMesh2Eigen(meshA, std::get<0>(meshOut), std::get<1>(meshOut));

    return meshOut;
}

#if 0

void computeDifference(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out) {

    CorkMesh cmIn0, cmIn1;
    corkTriMesh2CorkMesh(in0, &cmIn0);
    corkTriMesh2CorkMesh(in1, &cmIn1);

    cmIn0.boolDiff(cmIn1);

    corkMesh2CorkTriMesh(&cmIn0, out);
}

void computeIntersection(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out) {

    CorkMesh cmIn0, cmIn1;
    corkTriMesh2CorkMesh(in0, &cmIn0);
    corkTriMesh2CorkMesh(in1, &cmIn1);

    cmIn0.boolIsct(cmIn1);

    corkMesh2CorkTriMesh(&cmIn0, out);
}

void computeSymmetricDifference(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out){

    CorkMesh cmIn0, cmIn1;

    corkTriMesh2CorkMesh(in0, &cmIn0);
    corkTriMesh2CorkMesh(in1, &cmIn1);

    cmIn0.boolXor(cmIn1);

    corkMesh2CorkTriMesh(&cmIn0, out);
}

void resolveIntersections(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out) {

    CorkMesh cmIn0, cmIn1;

    corkTriMesh2CorkMesh(in0, &cmIn0);
    corkTriMesh2CorkMesh(in1, &cmIn1);

    cmIn0.disjointUnion(cmIn1);
    cmIn0.resolveIntersections();

    corkMesh2CorkTriMesh(&cmIn0, out);
}
#endif


} // end of namespace

//PYBIND11_MAKE_OPAQUE(std::vector<slm::LayerGeometry::Ptr>)
//PYBIND11_MAKE_OPAQUE(std::vector<aerotech::Axis::Ptr>)

/*
struct CorkTriMesh
{
    uint    n_triangles;
    uint    n_vertices;
    uint    *triangles;
    float   *vertices;
};

*/


PYBIND11_MODULE(pycork, m) {

    m.doc() = R"pbdoc(
        Pycork Module
        -----------------------
        .. currentmodule:: pycork
        .. autosummary::
           :toctree: _generate

    )pbdoc";

    m.def("isSolid", &pycork::isSolid, "Determines if the mesh is manifold",
                     py::arg("vertices"), py::arg("tris"))
     .def("union", &pycork::booleanUnion, "Computes boolean Union between two meshes",
                    py::arg("vertsA"), py::arg("trisA"),
                    py::arg("vertsB"), py::arg("trisB"))
     .def("difference", &pycork::booleanDifference, "Computes boolean differnece between two meshes",
                   py::arg("vertsA"), py::arg("trisA"),
                   py::arg("vertsB"), py::arg("trisB"))
     .def("intersection", &pycork::booleanIntersection, "Computes boolean intersection between two meshes",
                   py::arg("vertsA"), py::arg("trisA"),
                   py::arg("vertsB"), py::arg("trisB"))
     .def("intersection", &pycork::booleanIntersection, "Computes boolean xor between two meshes",
                           py::arg("vertsA"), py::arg("trisA"),
                           py::arg("vertsB"), py::arg("trisB"))
     .def("resolveIntersection", &pycork::resolveIntersection, "Computes the intersection between two meshes",
                                  py::arg("vertsA"), py::arg("trisA"));


#if 0
    // This function will test whether or not a mesh is solid
     bool isSolid(CorkTriMesh mesh);

    // Boolean operations follow
    // result = A U B
     void computeUnion(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);

    // result = A - B
     void computeDifference(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);

    // result = A ^ B
     void computeIntersection(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);

    // result = A XOR B
     void computeSymmetricDifference(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);

    // Not a Boolean operation, but related:
    //  No portion of either surface is deleted.  However, the
    //  curve of intersection between the two surfaces is made explicit,
    //  such that the two surfaces are now connected.
     void resolveIntersections(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);

    py::enum_<aerotech::AxisId>(m, "AxisId")
    .value("X", AxisId::X)
    .value("Y", AxisId::Y)
    .value("Z", AxisId::Z)
    .value("A", AxisId::A)
    .value("B", AxisId::B)
    .value("C", AxisId::C)
    .export_values();


    py::enum_<aerotech::EdgeMode>(m, "EdgeMode")
    .value("None", EdgeMode::NONE)
    .value("Enter", EdgeMode::ENTER)
    .value("Exit", EdgeMode::EXIT)
    .value("Both", EdgeMode::BOTH)
    .export_values();

    py::enum_<aerotech::RampMode>(m, "RampMode")
    .value("Rate", RampMode::RATE)
    .value("Time", RampMode::TIME)
    .export_values();

    py::enum_<aerotech::RampType>(m, "RampType")
    .value("Linear", RampType::LINEAR)
    .value("SCurve", RampType::SCURVE)
    .value("Sine", RampType::SINE)
    .export_values();

    py::enum_<aerotech::MotionMode>(m, "MotionMode")
    .value("Absolute", MotionMode::MOTION_ABSOLUTE )
    .value("Incremental", MotionMode::MOTION_INCREMENTAL)
    .export_values();

    py::register_exception<A3200Exception>(m, "Exception");

    py::bind_vector<std::vector<aerotech::Axis::Ptr>>(m, "VectorAxis");

    py::implicitly_convertible<py::list, std::vector<aerotech::Axis::Ptr>>();

    py::class_<aerotech::PSO, std::shared_ptr<aerotech::PSO>>(m, "PSO")
        .def(py::init<aerotech::Axis::Ptr>())
        .def("enable", &aerotech::PSO::enable, "Enables the PSO Mode on the Axis",
                       py::arg("taskId") = 1)
        .def("disable", &aerotech::PSO::disable, "Disables the PSO Mode on the Axis",
                        py::arg("taskId") = 1)
        .def("enableWindow", &aerotech::PSO::enableWindow, "Enables the PSO Windw Mode on the Axis",
                       py::arg("taskId") = 1)
        .def("disableWindow", &aerotech::PSO::disableWindow, "Disables the PSO Window Mode on the Axis",
                        py::arg("taskId") = 1)

        .def("reset",  &aerotech::PSO::reset, "Resets the counter registers on the Axis",
                       py::arg("hard"), py::arg("taskId") = 1)
        .def("arm", &aerotech::PSO::arm, "Arms the PSO firing on the Axis",
                    py::arg("taskId") = 1)
        .def("disarm", &aerotech::PSO::disarm, "Disarms the PSO firing on the Axis",
                       py::arg("taskId") = 1)
        .def("setFireContiniously", &aerotech::PSO::setFireContiniously, "Arms the PSO to fire continiously",
                                    py::arg("taskId") = 1)
        .def("setFireDistance", &aerotech::PSO::clearFireDistance, "Clear the fire-distance for PSO to fire",
                                py::arg("taskId") = 1)
        .def("clearFireDistance", &aerotech::PSO::setFireDistance, "Sets the PSO to fire at a fixed distance intveral",
                                py::arg("distance"), py::arg("taskId") = 1)
        .def("setOutput", &aerotech::PSO::setOutput, "Sets the output options for the PSO",
                         py::arg("pin") = 1, py::arg("mode") = 1, py::arg("taskId") = 1)
        .def("setPulseDelayOnly", &aerotech::PSO::setPulseDelayOnly, "Sets the PSO to use the pulse delay only mode",
                                  py::arg("totalTime"), py::arg("onTime"), py::arg("delayTime"), py::arg("taskId") = 1)
        .def("setPulseCyclesOnly", &aerotech::PSO::setPulseCyclesOnly, "Sets the PSO to use the pulse cycles only mode",
                                   py::arg("totalTime"), py::arg("onTime"), py::arg("cycles"), py::arg("taskId") = 1)
        .def("setEncoderAxis", &aerotech::PSO::setEncoderAxis, "Set the primary and secondary encoder axes to use as input for the axis",
                               py::arg("encoderId"), py::arg("encoderId2") = -1, py::arg("encoderId3") = -1, py::arg("invert") = false, py::arg("taskId") = 1)
        .def("setWindowMask", &aerotech::PSO::setWindowMask, "Sets the PSO to use a Window Mask for firing by passing an array",
                              py::arg("mask"), py::arg("edgeMode"), py::arg("arrayIdx"), py::arg("hard"), py::arg("taskId") = 1
                              )
        .def("setWindowRange", &aerotech::PSO::setWindowRange, "Sets the PSO fixed window range to fire in",
                              py::arg("low"), py::arg("high"), py::arg("taskId") = 1)
        .def_property_readonly("axis", &aerotech::PSO::axis, "The axis the PSO is operating on")
        .def_property_readonly("isArmed", &aerotech::PSO::isArmed, "Armed status of the PSO")
        .def_property_readonly("windowArrayIndexLocation", &aerotech::PSO::windowArrayIndexLocation, "Index of location axis")
        .def_property_readonly("axis", &aerotech::PSO::isArmed, "Armed status of the PSO")
        .def("psoCounter", &aerotech::PSO::psoCounter, "The current value of the PSO Counter for the Window",
                           py::arg("windowNumber"));


    py::class_<aerotech::Axis, std::shared_ptr<aerotech::Axis>>(m, "Axis")

        .def(py::init<aerotech::A3200Controller::Ptr, uint16_t>())
        .def("enable", &aerotech::Axis::enable,
                       py::arg("taskId") = 1)
        .def("disable", &aerotech::Axis::disable, "Disables the axis",
                        py::arg("taskId") = 1)
        .def("abort",  &aerotech::Axis::abort, "Aborts the current motion on the axis")
        .def("acknowledgeFault", &aerotech::Axis::acknowlegeFault,
                                 py::arg("taskId") = 1 )
            .def("setLabel", &aerotech::Axis::setLabel)
            .def("getLabel", &aerotech::Axis::label)
        .def("home", &aerotech::Axis::home, "Homes the axis",
                     py::arg("taskId") = 1)
        .def("move", &aerotech::Axis::move, "Moves the axis an aboslute or incremental distance",
                     py::arg("position"), py::arg("speed") = -1.0, py::arg("taskId") = 1)

        .def("moveRel", &aerotech::Axis::moveRel,
                        py::arg("distance"), py::arg("speed") = -1.0, py::arg("taskId") = 1)

        .def("wait", &aerotech::Axis::wait,
                      py::arg("timeout") = -1.0, py::arg("inPosition") = false)

        .def("setRampRate", &aerotech::Axis::setRampRate, "Sets the ramp rate, when a LINEAR ramp mode is used for the axis",
                            py::arg("rampRate"), py::arg("taskId") = 1)
        .def("setRampMode", &aerotech::Axis::setRampMode, "Sets the ramp mode for the axis",
                            py::arg("mpde"), py::arg("taskId") = 1)

        .def("setDigitalPin", &aerotech::Axis::setDigitalPin, "Sets the digital output of the pin",
                              py::arg("pin"), py::arg("state"), py::arg("taskId") = 1)
        .def("digitalPin", &aerotech::Axis::digitalPin, "Gets the digital input signal from a selected pin",
                           py::arg("pin"), py::arg("taskId") = 1)
        .def("setAnalogPin", &aerotech::Axis::setAnalogPin, "Sets the analog output of the pin",
                             py::arg("pin"), py::arg("value"), py::arg("taskId") = 1)
        .def("analogPin", &aerotech::Axis::analogPin, "Gets the analog input signal from a selected pin",
                          py::arg("pin"), py::arg("taskId") = 1)

        .def_property_readonly("controller", &aerotech::Axis::controller, "The controller assigned to the axis")
        .def_property_readonly("isEnabled", &aerotech::Axis::isEnabled)
        .def_property_readonly("isHomed", &aerotech::Axis::isHomed)
        .def_property_readonly("isBlocked", &aerotech::Axis::isBlocked)
        .def_property_readonly("id", &aerotech::Axis::id, "Corresponding axis id")
        .def_property_readonly("maxCurrent",  &aerotech::Axis::maxCurrent, "The max current consumption set for the axis")
        .def_property_readonly("position",  &aerotech::Axis::position, "The current position of the axis")
        .def("velocity", &aerotech::Axis::velocity, "Obtains the instantaneous (averaged) velocity of the axis",
                         py::arg("average") = false)
        .def("current",  &aerotech::Axis::current, "Obtains the instantaneous (averaged) current from the motor ",
                         py::arg("average") = false)
        .def_property("defaultSpeed",
                      &Axis::defaultSpeed, &Axis::setDefaultSpeed, py::return_value_policy::copy)
        .def_property("label",
                      &Axis::label, &Axis::setLabel, py::return_value_policy::copy);


    // use for c+11 compilers -  see https://pybind11.readthedocs.io/en/stable/classes.html
    //template <typename... Args>
    //using overload_cast_ = pybind11::detail::overload_cast_impl<Args...>;

    py::class_<aerotech::A3200Controller, std::shared_ptr<aerotech::A3200Controller>>(m, "Controller")
        .def(py::init())
        .def("connect", &aerotech::A3200Controller::connect)
        .def("reset",  &aerotech::A3200Controller::reset)
        .def("disconnect", &aerotech::A3200Controller::disconnect)
        .def("isConnected", &aerotech::A3200Controller::isConnected)
        .def("dwell", &aerotech::A3200Controller::dwell, "Dwells or waits for a period of time (seconds)")
        .def("getLastError", &aerotech::A3200Controller::getErrorString, "Get the last error string reported by the controller")
        .def("getLastErrorCode", &aerotech::A3200Controller::getErrorCode, "Get the last error code reported by the controller")
        .def("runScript", &aerotech::A3200Controller::runScript, "Runs a script from a specified filename",
                          py::arg("filename"), py::arg("taskId"))
        .def("runCommand", &aerotech::A3200Controller::runCommand, "Executes a command",
                           py::arg("command"), py::arg("taskId"))
        .def("startTaskQueue", &aerotech::A3200Controller::startTaskQueue, "Starts a task queue to buffer commands for asynchronous behavior",
                           py::arg("taskId") = 1)
        .def("blockTaskQueue", &aerotech::A3200Controller::blockUntilQueueComplete, "Blocks the Python interpreter until the task queue is complete",
                              py::arg("pollingTime") = 100 , py::arg("taskId") = 1)
        .def("endTaskQueue", &aerotech::A3200Controller::endTaskQueue, "Ends a task queue to buffer commands for asynchronous behavior",
                             py::arg("hard"), py::arg("timeout") = 1000 , py::arg("taskId") = 1)
        .def_property_readonly("queueCount", &aerotech::A3200Controller::queueCount, "The current queue count")
        .def_property_readonly("maximumQueueSize", &aerotech::A3200Controller::maxQueueSize, "Maximum queue size")
        .def("stopProgram", &aerotech::A3200Controller::stopProgram, "Instantly stops the program and optionally waits for axes to complete movement",
                             py::arg("timeout"), py::arg("taskId"))
        .def("criticalStart", &aerotech::A3200Controller::criticalStart, "Start a critical operation period",
                              py::arg("time") = -1.0)
        .def("criticalEnd", &aerotech::A3200Controller::criticalEnd, "Ends a critical operation period")
        .def("acknowledgeAll", &aerotech::A3200Controller::acknowledgeAll, "Acknowledge all fault alerts across all axes",
                               py::arg("taskId") = 1)
        .def("set", &aerotech::A3200Controller::setMotionMode, "Sets the motion mode (absolute/incremental) for all linear motions across all axes",
                              py::arg("mode"), py::arg("taskId") = 1)
        .def("setDefaultSpeed", &aerotech::A3200Controller::setDefaultSpeed, "Sets the default speed of the axis")
        .def("enable", py::overload_cast<const std::vector<Axis::Ptr> &, uint8_t>(&aerotech::A3200Controller::enable), "Enables the axes for motion",
                       py::arg("axes"), py::arg("taskId") = 1)
        .def("enable", py::overload_cast<Axis::Ptr, uint8_t>(&aerotech::A3200Controller::enable), "Enables the axis for motion",
                       py::arg("axis"), py::arg("taskId") = 1)
        .def("disable", py::overload_cast<const std::vector<Axis::Ptr> &, uint8_t>(&aerotech::A3200Controller::disable), "Disables the axes for motion",
                        py::arg("axes"), py::arg("taskId") = 1)
        .def("disable", py::overload_cast<Axis::Ptr, uint8_t>(&aerotech::A3200Controller::disable), "Disable the axis for motion",
                        py::arg("axis"), py::arg("taskId") = 1)
        .def("abort", py::overload_cast<const std::vector<Axis::Ptr> &>(&aerotech::A3200Controller::abort), "Aborts all motion across all axes",
                      py::arg("axes"))
        .def("abort", py::overload_cast<Axis::Ptr>(&aerotech::A3200Controller::abort), "Aborts motion along an axis",
                      py::arg("axis"))
        .def("home", py::overload_cast<const std::vector<Axis::Ptr> &, uint8_t>(&aerotech::A3200Controller::home), "Home all the axes instantaenously",
                     py::arg("axes"), py::arg("taskId") = 1)
        .def("home", py::overload_cast<Axis::Ptr, uint8_t>(&aerotech::A3200Controller::home), "Home the selected axis",
                     py::arg("axis"), py::arg("taskId") = 1)
        .def("move", py::overload_cast<const std::vector<Axis::Ptr> &, std::vector<double>, const double , uint8_t>(&aerotech::A3200Controller::move), "Move all the axes",
                     py::arg("axes"), py::arg("position"), py::arg("velocity") = -1.0, py::arg("taskId") = 1)
        .def("move", py::overload_cast<Axis::Ptr, double, const double, uint8_t>(&aerotech::A3200Controller::move), "Move an axis",
                     py::arg("axis"), py::arg("position"), py::arg("velocity") = -1.0, py::arg("taskId") = 1)
        .def("wait", py::overload_cast<const std::vector<Axis::Ptr> &, uint64_t, bool>(&aerotech::A3200Controller::wait), "Suspend operations and wait for motion to complete across all selected axes",
                     py::arg("axes"), py::arg("timeout"), py::arg("inPosition") = false)
        .def("wait", py::overload_cast<Axis::Ptr, uint64_t, bool>(&aerotech::A3200Controller::wait), "Suspend operations and wait until completion on the selected axis",
                     py::arg("axis"), py::arg("timeout"), py::arg("inPosition") = false)
        .def("setRampRate", py::overload_cast<const std::vector<Axis::Ptr> &, const double, const uint8_t>(&aerotech::A3200Controller::setRampRate), "Sets the ramp rate, when a LINEAR ramp mode is used for the axes",
                            py::arg("axes"), py::arg("rampRate"), py::arg("taskId") = 1)
        .def("setRampRate", py::overload_cast<Axis::Ptr, const double, const uint8_t>(&aerotech::A3200Controller::setRampRate), "Sets the ramp rate, when a LINEAR ramp mode is used for the axis",
                            py::arg("axis"), py::arg("rampRate"), py::arg("taskId") = 1)
        .def("setRampMode", py::overload_cast<const std::vector<Axis::Ptr> &, const RampMode, const uint8_t>(&aerotech::A3200Controller::setRampMode), "Sets the ramp mode for the axes",
                            py::arg("axess"), py::arg("rampMode"), py::arg("taskId") = 1)
        .def("setRampMode", py::overload_cast<Axis::Ptr, const RampMode, const uint8_t>(&aerotech::A3200Controller::setRampMode), "Sets the ramp mode for the axis",
                            py::arg("axis"), py::arg("rampMode"), py::arg("taskId") = 1)
        .def("position", (std::vector<double>(aerotech::A3200Controller::*)(const std::vector<Axis::Ptr> &) const)(&aerotech::A3200Controller::position), "Obtains the instantaneous position of the axes",
                          py::arg("axes"))
        .def("position", (double(A3200Controller::*)(Axis::Ptr) const)(&aerotech::A3200Controller::position), "Obtains the instantaneous position of the axis",
                          py::arg("axis"))
        .def("velocity", (std::vector<double>(aerotech::A3200Controller::*)(const std::vector<Axis::Ptr> &, bool) const)(&aerotech::A3200Controller::velocity), "Obtains the instantaneous (averaged) velocity of the axes",
                          py::arg("axes"), py::arg("average") = false)
        .def("velocity", (double(A3200Controller::*)(Axis::Ptr, bool) const)(&aerotech::A3200Controller::velocity), "Obtains the instantaneous (averaged) velocity of the axis",
                          py::arg("axis"), py::arg("average") = false)
        .def("setGlobalVariable", (void(aerotech::A3200Controller::*)(uint32_t, std::vector<double> &))(&aerotech::A3200Controller::setGlobalVariable), "Set the global variable at index (idx) with value",
                                   py::arg("idx"), py::arg("value"))
        .def("setGlobalVariable", (void(aerotech::A3200Controller::*)(uint32_t, double))(&aerotech::A3200Controller::setGlobalVariable), "Set the global variable at index (idx) with value",
                                   py::arg("idx"), py::arg("value"))
        .def("getSingleDataSignal", &aerotech::A3200Controller::getSingleDataSignal, "Captures immediately a single data signal for an axis",
                                    py::arg("axis"), py::arg("signal"))
        .def("getDataSignal",(Eigen::MatrixXd(aerotech::A3200Controller::*)(Axis::Ptr, std::vector<uint32_t>, uint32_t, uint32_t))( &aerotech::A3200Controller::getDataSignal) , py::return_value_policy::copy,
                              "Captures multiple signals and their occurances from a signle axis, for a sample period",
                              py::arg("axis"), py::arg("signals"), py::arg("numPoints"), py::arg("samplePeriod") = 1);

#endif

#ifdef PROJECT_VERSION
    m.attr("__version__") = "PROJECT_VERSION";
#else
    m.attr("__version__") = "dev";
#endif

}

