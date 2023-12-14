//
// Define functions to convert "clifford" python package CGA3 objects to cga3<double> objects.
//

#include "pybind11/pybind11.h"
#include "CliffordInterface.h"

/// Convert from a clifford.MultiVector to a cga3 Vector
cga3<double>::Vector CliffordInterface::from_MultiVector_to_Vector(const py::object &mv_)
{
	const auto& mv = py::array_t<double>(mv_.attr("as_array")());
	assert(mv.ndim() == 1 && mv.shape(0) == 32);
	const auto& layout = cga3<double>::get();
	cga3<double>::Vector result = eval(layout.e1 * mv.at(1) + layout.e2 * mv.at(2) + layout.e3 * mv.at(3) +
	                                   layout.ep * mv.at(4) + layout.em * mv.at(5));
	return result;
}

/// Convert from a clifford.MultiVector to a cga3 PointPair
cga3<double>::PointPair CliffordInterface::from_MultiVector_to_PointPair(const py::object &mv_)
{
	const auto& mv = py::array_t<double>(mv_.attr("as_array")());
	assert(mv.ndim() == 1 && mv.shape(0) == 32);

	const auto& layout = cga3<double>::get();
	cga3<double>::PointPair cumulative_sum = eval(mv.at(6) * (layout.e1 ^ layout.e2));
	cumulative_sum = cumulative_sum + eval(mv.at(7) * (layout.e1 ^ layout.e3));
	cumulative_sum = cumulative_sum + eval(mv.at(8) * (layout.e1 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(9) * (layout.e1 ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(10) * (layout.e2 ^ layout.e3));
	cumulative_sum = cumulative_sum + eval(mv.at(11) * (layout.e2 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(12) * (layout.e2 ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(13) * (layout.e3 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(14) * (layout.e3 ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(15) * (layout.ep ^ layout.em));
	return cumulative_sum;
}

/// Convert from a clifford.MultiVector to a cga3 Line
cga3<double>::Line CliffordInterface::from_MultiVector_to_Line(const py::object &mv_)
{
	const auto& mv = py::array_t<double>(mv_.attr("as_array")());
	assert(mv.ndim() == 1 && mv.shape(0) == 32);

	const auto& layout = cga3<double>::get();
	cga3<double>::Line cumulative_sum = eval(mv.at(17) * (layout.e1 ^ layout.e2 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(18) * (layout.e1 ^ layout.e2 ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(19) * (layout.e1 ^ layout.e3 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(20) * (layout.e1 ^ layout.e3 ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(21) * (layout.e1 ^ layout.ep ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(22) * (layout.e2 ^ layout.e3 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(23) * (layout.e2 ^ layout.e3 ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(24) * (layout.e2 ^ layout.ep ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(25) * (layout.e3 ^ layout.ep ^ layout.em));
	return cumulative_sum;
}

cga3<double>::Rotor CliffordInterface::from_MultiVector_to_Rotor(const py::object &mv_)
{
	const auto& mv = py::array_t<double>(mv_.attr("as_array")());
	assert(mv.ndim() == 1 && mv.shape(0) == 32);

	const auto& layout = cga3<double>::get();
	cga3<double>::Rotor cumulative_sum = eval(mv.at(0) * layout.one);
	cumulative_sum = cumulative_sum + eval(mv.at(6) * (layout.e1 ^ layout.e2));
	cumulative_sum = cumulative_sum + eval(mv.at(7) * (layout.e1 ^ layout.e3));
	cumulative_sum = cumulative_sum + eval(mv.at(10) * (layout.e2 ^ layout.e3));
	return cumulative_sum;
}

cga3<double>::Translator CliffordInterface::from_MultiVector_to_Translator(const py::object &mv_)
{
	const auto& mv = py::array_t<double>(mv_.attr("as_array")());
	assert(mv.ndim() == 1 && mv.shape(0) == 32);

	const auto& layout = cga3<double>::get();
	cga3<double>::Translator cumulative_sum = eval(mv.at(0) * layout.one);
	cumulative_sum = cumulative_sum + eval(mv.at(8) * (layout.e1 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(9) * (layout.e1 ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(11) * (layout.e2 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(12) * (layout.e2 ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(13) * (layout.e3 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(14) * (layout.e3 ^ layout.em));
	return cumulative_sum;
}

cga3<double>::Motor CliffordInterface::from_MultiVector_to_Motor(const py::object &mv_)
{
	const auto& mv = py::array_t<double>(mv_.attr("as_array")());
	assert(mv.ndim() == 1 && mv.shape(0) == 32);

	const auto& layout = cga3<double>::get();
	cga3<double>::Motor cumulative_sum = eval(mv.at(0) * layout.one);
	cumulative_sum = cumulative_sum + eval(mv.at(6) * (layout.e1 ^ layout.e2));
	cumulative_sum = cumulative_sum + eval(mv.at(7) * (layout.e1 ^ layout.e3));
	cumulative_sum = cumulative_sum + eval(mv.at(8) * (layout.e1 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(9) * (layout.e1 ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(10) * (layout.e2 ^ layout.e3));
	cumulative_sum = cumulative_sum + eval(mv.at(11) * (layout.e2 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(12) * (layout.e2 ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(13) * (layout.e3 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(14) * (layout.e3 ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(26) * (layout.e1 ^ layout.e2 ^ layout.e3 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(27) * (layout.e1 ^ layout.e2 ^ layout.e3 ^ layout.em));
	return cumulative_sum;
}

/// Convert from a clifford.MultiVector to a cga3 Sphere
cga3<double>::Sphere CliffordInterface::from_MultiVector_to_Sphere(const py::object &mv_)
{
	const auto& mv = py::array_t<double>(mv_.attr("as_array")());
	assert(mv.ndim() == 1 && mv.shape(0) == 32);

	const auto& layout = cga3<double>::get();
	cga3<double>::Sphere cumulative_sum = eval(mv.at(26) * (layout.e1 ^ layout.e2 ^ layout.e3 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(27) * (layout.e1 ^ layout.e2 ^ layout.e3 ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(28) * (layout.e1 ^ layout.e2 ^ layout.ep ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(29) * (layout.e1 ^ layout.e3 ^ layout.ep ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(30) * (layout.e2 ^ layout.e3 ^ layout.ep ^ layout.em));
	return cumulative_sum;
}


/// Convert from a clifford.MultiVector to a cga3 Circle
cga3<double>::Circle CliffordInterface::from_MultiVector_to_Circle(const py::object &mv_)
{
	const auto& mv = py::array_t<double>(mv_.attr("as_array")());
	assert(mv.ndim() == 1 && mv.shape(0) == 32);

	const auto& layout = cga3<double>::get();
	cga3<double>::Circle cumulative_sum = eval(mv.at(16) * (layout.e1 ^ layout.e2 ^ layout.e3));
	cumulative_sum = cumulative_sum + eval(mv.at(17) * (layout.e1 ^ layout.e2 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(18) * (layout.e1 ^ layout.e2 ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(19) * (layout.e1 ^ layout.e3 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(20) * (layout.e1 ^ layout.e3 ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(21) * (layout.e1 ^ layout.ep ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(22) * (layout.e2 ^ layout.e3 ^ layout.ep));
	cumulative_sum = cumulative_sum + eval(mv.at(23) * (layout.e2 ^ layout.e3 ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(24) * (layout.e2 ^ layout.ep ^ layout.em));
	cumulative_sum = cumulative_sum + eval(mv.at(25) * (layout.e3 ^ layout.ep ^ layout.em));
	return cumulative_sum;
}


py::object CliffordInterface::to_MultiVector(const cga3<double>::Vector &v)
{
#define EXTRACT_BLADE_VALUE(NUM) \
   v.template element<vector_conf(NUM)>()

	py::array_t<double> mv = py::array_t<double>(32);
	memset((void*) mv.data(), 0, sizeof(double) * 32);
	mv.mutable_at(1) = EXTRACT_BLADE_VALUE(1);
	mv.mutable_at(2) = EXTRACT_BLADE_VALUE(2);
	mv.mutable_at(3) = EXTRACT_BLADE_VALUE(3);
	mv.mutable_at(4) = EXTRACT_BLADE_VALUE(4);
	mv.mutable_at(5) = EXTRACT_BLADE_VALUE(5);
	return m_MultiVector(m_clifford_g3c_layout, mv);
#undef EXTRACT_BLADE_VALUE
}


py::object CliffordInterface::to_MultiVector(const cga3<double>::PointPair &pp)
{
#define EXTRACT_BLADE_VALUE(CL) \
   pp.template element<CL>()

	py::array_t<double> mv = py::array_t<double>(32);
	memset((void*) mv.data(), 0, sizeof(double) * 32);
	mv.mutable_at(6) = EXTRACT_BLADE_VALUE(blade_conf(1, 2));
	mv.mutable_at(7) = EXTRACT_BLADE_VALUE(blade_conf(1, 3));
	mv.mutable_at(8) = EXTRACT_BLADE_VALUE(blade_conf(1, 4));
	mv.mutable_at(9) = EXTRACT_BLADE_VALUE(blade_conf(1, 5));
	mv.mutable_at(10) = EXTRACT_BLADE_VALUE(blade_conf(2, 3));
	mv.mutable_at(11) = EXTRACT_BLADE_VALUE(blade_conf(2, 4));
	mv.mutable_at(12) = EXTRACT_BLADE_VALUE(blade_conf(2, 5));
	mv.mutable_at(13) = EXTRACT_BLADE_VALUE(blade_conf(3, 4));
	mv.mutable_at(14) = EXTRACT_BLADE_VALUE(blade_conf(3, 5));
	mv.mutable_at(15) = EXTRACT_BLADE_VALUE(blade_conf(4, 5));

	return m_MultiVector(m_clifford_g3c_layout, mv);
#undef EXTRACT_BLADE_VALUE
}

py::object CliffordInterface::to_MultiVector(const cga3<double>::Line &line)
{
#define EXTRACT_BLADE_VALUE(CL) \
   line.template element<CL>()

	py::array_t<double> mv = py::array_t<double>(32);
	memset((void*) mv.data(), 0, sizeof(double) * 32);
	mv.mutable_at(17) = EXTRACT_BLADE_VALUE(blade_conf(1, 2, 4));
	mv.mutable_at(18) = EXTRACT_BLADE_VALUE(blade_conf(1, 2, 5));
	mv.mutable_at(19) = EXTRACT_BLADE_VALUE(blade_conf(1, 3, 4));
	mv.mutable_at(20) = EXTRACT_BLADE_VALUE(blade_conf(1, 3, 5));
	mv.mutable_at(21) = EXTRACT_BLADE_VALUE(blade_conf(1, 4, 5));
	mv.mutable_at(22) = EXTRACT_BLADE_VALUE(blade_conf(2, 3, 4));
	mv.mutable_at(23) = EXTRACT_BLADE_VALUE(blade_conf(2, 3, 5));
	mv.mutable_at(24) = EXTRACT_BLADE_VALUE(blade_conf(2, 4, 5));
	mv.mutable_at(25) = EXTRACT_BLADE_VALUE(blade_conf(3, 4, 5));

	return m_MultiVector(m_clifford_g3c_layout, mv);
#undef EXTRACT_BLADE_VALUE
}

py::object CliffordInterface::to_MultiVector(const cga3<double>::Rotor &rotor)
{
#define EXTRACT_BLADE_VALUE(CL) \
   rotor.template element<CL>()

	py::array_t<double> mv = py::array_t<double>(32);
	memset((void*) mv.data(), 0, sizeof(double) * 32);
	mv.mutable_at(0) = EXTRACT_BLADE_VALUE(0);
	mv.mutable_at(6) = EXTRACT_BLADE_VALUE(blade_conf(1, 2));
	mv.mutable_at(7) = EXTRACT_BLADE_VALUE(blade_conf(1, 3));
	mv.mutable_at(10) = EXTRACT_BLADE_VALUE(blade_conf(2, 3));
	return m_MultiVector(m_clifford_g3c_layout, mv);

#undef EXTRACT_BLADE_VALUE
}


py::object CliffordInterface::to_MultiVector(const cga3<double>::Translator &translator)
{
#define EXTRACT_BLADE_VALUE(CL) \
   translator.template element<CL>()

	py::array_t<double> mv = py::array_t<double>(32);
	memset((void*) mv.data(), 0, sizeof(double) * 32);
	mv.mutable_at(0) = EXTRACT_BLADE_VALUE(0);
	mv.mutable_at(8) = EXTRACT_BLADE_VALUE(blade_conf(1, 4));
	mv.mutable_at(9) = EXTRACT_BLADE_VALUE(blade_conf(1, 5));
	mv.mutable_at(11) = EXTRACT_BLADE_VALUE(blade_conf(2, 4));
	mv.mutable_at(12) = EXTRACT_BLADE_VALUE(blade_conf(2, 5));
	mv.mutable_at(13) = EXTRACT_BLADE_VALUE(blade_conf(3, 4));
	mv.mutable_at(14) = EXTRACT_BLADE_VALUE(blade_conf(3, 5));
	return m_MultiVector(m_clifford_g3c_layout, mv);
#undef EXTRACT_BLADE_VALUE
}

py::object CliffordInterface::to_MultiVector(const cga3<double>::Motor &motor)
{
#define EXTRACT_BLADE_VALUE(CL) \
   motor.template element<CL>()

	py::array_t<double> mv = py::array_t<double>(32);
	memset((void*) mv.data(), 0, sizeof(double) * 32);
	mv.mutable_at(0) = EXTRACT_BLADE_VALUE(0);
	mv.mutable_at(6) = EXTRACT_BLADE_VALUE(blade_conf(1, 2));
	mv.mutable_at(7) = EXTRACT_BLADE_VALUE(blade_conf(1, 3));
	mv.mutable_at(8) = EXTRACT_BLADE_VALUE(blade_conf(1, 4));
	mv.mutable_at(9) = EXTRACT_BLADE_VALUE(blade_conf(1, 5));
	mv.mutable_at(10) = EXTRACT_BLADE_VALUE(blade_conf(2, 3));
	mv.mutable_at(11) = EXTRACT_BLADE_VALUE(blade_conf(2, 4));
	mv.mutable_at(12) = EXTRACT_BLADE_VALUE(blade_conf(2, 5));
	mv.mutable_at(13) = EXTRACT_BLADE_VALUE(blade_conf(3, 4));
	mv.mutable_at(14) = EXTRACT_BLADE_VALUE(blade_conf(3, 5));
	mv.mutable_at(26) = EXTRACT_BLADE_VALUE(blade4_conf(1, 2, 3, 4));
	mv.mutable_at(27) = EXTRACT_BLADE_VALUE(blade4_conf(1, 2, 3, 5));
	return m_MultiVector(m_clifford_g3c_layout, mv);
#undef EXTRACT_BLADE_VALUE
}

py::object CliffordInterface::to_MultiVector(const cga3<double>::Sphere &s)
{
#define EXTRACT_BLADE_VALUE(CL) \
   s.template element<CL>()

	py::array_t<double> mv = py::array_t<double>(32);
	memset((void*) mv.data(), 0, sizeof(double) * 32);
	mv.mutable_at(26) = EXTRACT_BLADE_VALUE(blade4_conf(1, 2, 3, 4));
	mv.mutable_at(27) = EXTRACT_BLADE_VALUE(blade4_conf(1, 2, 3, 5));
	mv.mutable_at(28) = EXTRACT_BLADE_VALUE(blade4_conf(1, 2, 4, 5));
	mv.mutable_at(29) = EXTRACT_BLADE_VALUE(blade4_conf(1, 3, 4, 5));
	mv.mutable_at(30) = EXTRACT_BLADE_VALUE(blade4_conf(2, 3, 4, 5));
	return m_MultiVector(m_clifford_g3c_layout, mv);
#undef EXTRACT_BLADE_VALUE
}

py::object CliffordInterface::to_MultiVector(const cga3<double>::Circle &circle)
{
#define EXTRACT_BLADE_VALUE(CL) \
   circle.template element<CL>()

	py::array_t<double> mv = py::array_t<double>(32);
	memset((void*) mv.data(), 0, sizeof(double) * 32);
	mv.mutable_at(16) = EXTRACT_BLADE_VALUE(blade_conf(1, 2, 3));
	mv.mutable_at(17) = EXTRACT_BLADE_VALUE(blade_conf(1, 2, 4));
	mv.mutable_at(18) = EXTRACT_BLADE_VALUE(blade_conf(1, 2, 5));
	mv.mutable_at(19) = EXTRACT_BLADE_VALUE(blade_conf(1, 3, 4));
	mv.mutable_at(20) = EXTRACT_BLADE_VALUE(blade_conf(1, 3, 5));
	mv.mutable_at(21) = EXTRACT_BLADE_VALUE(blade_conf(1, 4, 5));
	mv.mutable_at(22) = EXTRACT_BLADE_VALUE(blade_conf(2, 3, 4));
	mv.mutable_at(23) = EXTRACT_BLADE_VALUE(blade_conf(2, 3, 5));
	mv.mutable_at(24) = EXTRACT_BLADE_VALUE(blade_conf(2, 4, 5));
	mv.mutable_at(25) = EXTRACT_BLADE_VALUE(blade_conf(3, 4, 5));

	return m_MultiVector(m_clifford_g3c_layout, mv);
#undef EXTRACT_BLADE_VALUE
}
