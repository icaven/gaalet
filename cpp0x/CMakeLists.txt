cmake_minimum_required(VERSION 3.15)
project(cpp0x)

set(CMAKE_CXX_STANDARD 14)

include_directories(../include/cpp0x)
include_directories(../examples/CarDynamics)
include_directories(../examples/CarDynamicsCGA)
include_directories(../examples/ClosedLoopKinematics)
include_directories(../examples/PGA3)

add_subdirectory (../examples)

#
#add_executable(cpp0x
#        benchmarks/benchmarkInverse.cpp
#        benchmarks/Horizon.cpp
#        benchmarks/VectorAdd.cpp
#        benchmarks/VectorAddTbb.cpp
#        benchmarks/VectorMultiply.cpp
#        examples/CarDynamics/CarDynamics.h
#        examples/CarDynamics/MagicFormula2004.h
#        examples/CarDynamics/mainCarDynamics.cpp
#        examples/CarDynamics/RungeKuttaClassic.h
#        examples/CarDynamics/TupleExpressions.h
#        examples/CarDynamicsCGA/CarDynamicsCGA.h
#        examples/CarDynamicsCGA/MagicFormula2004.h
#        examples/CarDynamicsCGA/mainCarDynamicsCGA.cpp
#        examples/CarDynamicsCGA/RungeKuttaClassic.h
#        examples/CarDynamicsCGA/TupleExpressions.h
#        examples/ClosedLoopKinematics/cga_osg.h
#        examples/ClosedLoopKinematics/ClosedLoopKinematics.cpp
#        examples/ElasticCoupling/ElasticCoupling.cpp
#        examples/MonteCarloPi/MonteCarloPi.cpp
#        examples/PGA3/Icosahedron/Icosahedron.cpp
#        examples/PGA3/PointsAndLines/PointsAndLines.cpp
#        examples/PGA3/OSG_Utilities.cpp
#        examples/PGA3/OSG_Utilities.h
#        tests/Basic.cpp
#        tests/CameraRotation.cpp
#        tests/ConfModel.cpp
#        tests/debug_segfault.cpp
#        tests/Degenerate.cpp
#        tests/Dual.cpp
#        tests/ElementTraits.cpp
#        tests/Exponential.cpp
#        tests/ExpressionForward.cpp
#        tests/GradientSearch.cpp
#        tests/Inverse.cpp
#        tests/MetricCombination.cpp
#        tests/MotorAlgebra.cpp
#        tests/Multivector.cpp
#        tests/NumericVersor.cpp
#        tests/Pga3.cpp
#        tests/PointVelocity.cpp
#        tests/Reverse.cpp
#        tests/Rotation.cpp
#        tests/SandwichProduct.cpp
#        tests/Scalar.cpp
#        tests/ScalarMultiplication.cpp
#        )
