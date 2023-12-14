//
// Declare functions to convert "clifford" python package CGA3 objects to/from cga3<double> objects.
//

#ifndef DISPARITY_TO_3D_CLIFFORDINTERFACE_H
#define DISPARITY_TO_3D_CLIFFORDINTERFACE_H

#pragma once

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "cga3d.h"

namespace py = pybind11;
using namespace pybind11::literals;

/// @brief A singleton class to encapsulate the interface between the C++ cga3 class and the Python clifford classes
class CliffordInterface
{
public:

	/// @brief Convert a vector of cga3 multivectors to a clifford ConformalMVArray
	template<class CGA3_ELEMENT>
	py::object to_ConformalMVArray(const std::vector<CGA3_ELEMENT> &elements);

	/// @brief Convert a clifford ConformalMVArray to a vector of cga3 multivectors
	template <class CGA3_ELEMENT>
	void from_ConformalMVArray(const py::sequence &conformal_mv_list, std::vector<CGA3_ELEMENT> &transformed);

	template <class CGA3_ELEMENT>
	class from_MultiVector
	{
	public:
		CGA3_ELEMENT operator()(const py::object& mv_) = delete;
	};

	/// @brief Functions to convert from cga3 objects to a clifford.MultiVector - not a complete set
	py::object to_MultiVector(const cga3<double>::Vector &v);
	py::object to_MultiVector(const cga3<double>::PointPair &pp);
	py::object to_MultiVector(const cga3<double>::Line &line);
	py::object to_MultiVector(const cga3<double>::Rotor &rotor);
	py::object to_MultiVector(const cga3<double>::Translator &translator);
	py::object to_MultiVector(const cga3<double>::Motor &motor);
	py::object to_MultiVector(const cga3<double>::Sphere &sphere);
	py::object to_MultiVector(const cga3<double>::Circle &circle);

	// Construct-on-first-use pattern
	static CliffordInterface& get()
	{
		static auto* ci = new CliffordInterface();
		return *ci;
	}

protected:
	CliffordInterface()
	{
		// Initialise references to the clifford classes
		py::module_ clifford = py::module::import("clifford");
		m_MultiVector = clifford.attr("MultiVector");

		py::module_ clifford_g3c = py::module::import("clifford.g3c");
		m_clifford_g3c_layout = clifford_g3c.attr("layout");

		py::module_ clifford_tools_g3c = py::module::import("clifford.tools.g3c");
		m_ConformalMVArray = clifford_tools_g3c.attr("ConformalMVArray");

	};

	/// Helper functions
	/// @brief Convert from a clifford.MultiVector to a cga3 Vector
	static cga3<double>::Vector from_MultiVector_to_Vector(const py::object &mv_);
	/// @brief Convert from a clifford.MultiVector to a cga3 Line
	static cga3<double>::PointPair from_MultiVector_to_PointPair(const py::object &mv_);
	/// @brief Convert from a clifford.MultiVector to a cga3 Line
	static cga3<double>::Line from_MultiVector_to_Line(const py::object &mv_);
	/// @brief Convert from a clifford.MultiVector to a cga3 Rotor
	static cga3<double>::Rotor from_MultiVector_to_Rotor(const py::object &mv_);
	/// @brief Convert from a clifford.MultiVector to a cga3 Translator
	static cga3<double>::Translator from_MultiVector_to_Translator(const py::object &mv_);
	/// @brief Convert from a clifford.MultiVector to a cga3 Motor
	static cga3<double>::Motor from_MultiVector_to_Motor(const py::object &mv_);
	/// @brief Convert from a clifford.MultiVector to a cga3 Sphere
	static cga3<double>::Sphere from_MultiVector_to_Sphere(const py::object &mv_);
	/// @brief Convert from a clifford.MultiVector to a cga3 Circle
	static cga3<double>::Circle from_MultiVector_to_Circle(const py::object &mv_);

	/// @brief A reference to the clifford.MultiVector() class
	py::object m_MultiVector;

	/// @brief A reference to the clifford.g3c.layout object
	py::object m_clifford_g3c_layout;

	/// @brief A reference to the clifford.tools.g3c.ConformalMVArray() class
	py::object m_ConformalMVArray;

};

/// Inline definitions
template <class CGA3_ELEMENT>
py::object CliffordInterface::to_ConformalMVArray(const std::vector<CGA3_ELEMENT> &elements)
{
	py::list clifford_mvs;
	for (const auto& v : elements)
	{
		clifford_mvs.append(to_MultiVector(v));
	}
	return m_ConformalMVArray(clifford_mvs);
}

template <class CGA3_ELEMENT>
void CliffordInterface::from_ConformalMVArray(const py::sequence &conformal_mv_list,
                                              std::vector<CGA3_ELEMENT> &transformed)
{
	transformed.reserve(conformal_mv_list.size());
	for (const auto& pyobj : conformal_mv_list)
	{
		transformed.emplace_back(CliffordInterface::from_MultiVector<CGA3_ELEMENT>()(pyobj));
	}
}


// Explicit instantiations
template<>
inline cga3<double>::Vector CliffordInterface::from_MultiVector<cga3<double>::Vector>::operator()(const py::object &mv_)
{
	return from_MultiVector_to_Vector(mv_);
}

template<>
inline cga3<double>::PointPair CliffordInterface::from_MultiVector<cga3<double>::PointPair>::operator()(const py::object &mv_)
{
	return from_MultiVector_to_PointPair(mv_);
}

template<>
inline cga3<double>::Line CliffordInterface::from_MultiVector<cga3<double>::Line>::operator()(const py::object &mv_)
{
	return from_MultiVector_to_Line(mv_);
}

template<>
inline cga3<double>::Rotor CliffordInterface::from_MultiVector<cga3<double>::Rotor>::operator()(const py::object &mv_)
{
	return from_MultiVector_to_Rotor(mv_);
}

template<>
inline cga3<double>::Translator CliffordInterface::from_MultiVector<cga3<double>::Translator>::operator()(const py::object &mv_)
{
	return from_MultiVector_to_Translator(mv_);
}

template<>
inline cga3<double>::Motor CliffordInterface::from_MultiVector<cga3<double>::Motor>::operator()(const py::object &mv_)
{
	return from_MultiVector_to_Motor(mv_);
}

template<>
inline cga3<double>::Sphere CliffordInterface::from_MultiVector<cga3<double>::Sphere>::operator()(const py::object &mv_)
{
	return from_MultiVector_to_Sphere(mv_);
}

template<>
inline cga3<double>::Circle CliffordInterface::from_MultiVector<cga3<double>::Circle>::operator()(const py::object &mv_)
{
	return from_MultiVector_to_Circle(mv_);
}


template py::object CliffordInterface::to_ConformalMVArray<cga3<double>::Vector>(const std::vector<cga3<double>::Vector> &elements);
template py::object CliffordInterface::to_ConformalMVArray<cga3<double>::Line>(const std::vector<cga3<double>::Line> &elements);
template py::object CliffordInterface::to_ConformalMVArray<cga3<double>::PointPair>(const std::vector<cga3<double>::PointPair> &elements);


#endif //DISPARITY_TO_3D_CLIFFORDINTERFACE_H
