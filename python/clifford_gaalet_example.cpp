//
// A simple example module showing how to interface the clifford python package with gaalet.
//

#include "pybind11/pybind11.h"
#include "CliffordInterface.h"

namespace py = pybind11;

// Functions (and classes) can be templated, so the element type of the algebra can be decided at compile time.
// However, when interfacing to the clifford python package, the double type must be used.
template <class ELEMENT_TYPE>
void lines_from_points_from_origin(
				std::vector<typename cga3<ELEMENT_TYPE>::Vector> &scene_points,
				std::vector<typename cga3<ELEMENT_TYPE>::Line> &lines_joining_pts_to_origin)
{
	const auto& layout = cga3<ELEMENT_TYPE>::get();
	lines_joining_pts_to_origin.clear();
	for (const auto& pt: scene_points)
	{
		lines_joining_pts_to_origin.emplace_back(eval(cga3<ELEMENT_TYPE>::normalise(layout.eo ^ pt ^ layout.einf)));
	}

}

/// Create lines from the origin to each of the given conformal points
py::object lines_from_origin_to_points(py::sequence &input_points_)
{
	std::vector<typename cga3<double>::Vector> input_points;
	CliffordInterface::get().from_ConformalMVArray(input_points_, input_points);
	std::vector<typename cga3<double>::Line> lines;
	lines_from_points_from_origin<double>(input_points, lines);
	return CliffordInterface::get().to_ConformalMVArray(lines);
}

const char* _doc_lines_from_origin_to_points = R"pbdoc(
		Create a line from the origin to each point
		@param input_points_:  The clifford cga3 Multivector points
		@return: A clifford.g3c.ConformalMVArray() with the lines
	)pbdoc";


PYBIND11_MODULE(clifford_gaalet_example, m)
{
	using namespace pybind11::literals;

	m.doc() = "This module provides an example how to interface the clifford python package with gaalet";

	m.attr("__copyright__") = "Copyright (c) 2023, Felix & Paul Studios Inc.";
	m.attr("__version__") = MODULE_VERSION;

	// Add function and class definitions
	m.def("lines_from_origin_to_points", lines_from_origin_to_points, _doc_lines_from_origin_to_points);
}
