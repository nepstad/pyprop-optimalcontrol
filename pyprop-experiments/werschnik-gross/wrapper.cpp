
// Boost Includes ==============================================================
#include <boost/python.hpp>
#include <boost/cstdint.hpp>

// Includes ====================================================================
#include <potential.cpp>

// Using =======================================================================
using namespace boost::python;

// Module ======================================================================
BOOST_PYTHON_MODULE(libpotential)
{
    class_< DynamicPotentialEvaluator<QuantumDotPotential<1>,1> >("QuantumDotPotential_1", init<  >())
        .def(init< const DynamicPotentialEvaluator<QuantumDotPotential<1>,1>& >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<QuantumDotPotential<1>,1>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<QuantumDotPotential<1>,1>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<QuantumDotPotential<1>,1>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<QuantumDotPotential<1>,1>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<QuantumDotPotential<1>,1>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<QuantumDotPotential<1>,1>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<QuantumDotPotential<1>,1>::CalculateExpectationValue)
    ;

    class_< DynamicPotentialEvaluator<DipoleOperator<1>,1> >("DipoleOperator_1", init<  >())
        .def(init< const DynamicPotentialEvaluator<DipoleOperator<1>,1>& >())
        .def("ApplyConfigSection", &DynamicPotentialEvaluator<DipoleOperator<1>,1>::ApplyConfigSection)
        .def("ApplyPotential", &DynamicPotentialEvaluator<DipoleOperator<1>,1>::ApplyPotential)
        .def("MultiplyPotential", &DynamicPotentialEvaluator<DipoleOperator<1>,1>::MultiplyPotential)
        .def("UpdateStaticPotential", &DynamicPotentialEvaluator<DipoleOperator<1>,1>::UpdateStaticPotential)
        .def("GetPotential", &DynamicPotentialEvaluator<DipoleOperator<1>,1>::GetPotential)
        .def("UpdatePotentialData", &DynamicPotentialEvaluator<DipoleOperator<1>,1>::UpdatePotentialData)
        .def("CalculateExpectationValue", &DynamicPotentialEvaluator<DipoleOperator<1>,1>::CalculateExpectationValue)
    ;

}

