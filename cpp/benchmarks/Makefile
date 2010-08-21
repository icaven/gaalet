INCLUDE_DIR = ../../include/cpp
INCLUDE = -I$(INCLUDE_DIR)

DEPENDS = $(INCLUDE_DIR)/*.h Makefile

benchmarkVectorAdd: benchmarkVectorAdd.cpp $(DEPENDS)
	g++45 -o benchmarkVectorAdd benchmarkVectorAdd.cpp -O3 $(INCLUDE)

benchmarkVectorMultiply: benchmarkVectorMultiply.cpp $(DEPENDS)
	g++45 -o benchmarkVectorMultiply benchmarkVectorMultiply.cpp -O3 $(INCLUDE)

