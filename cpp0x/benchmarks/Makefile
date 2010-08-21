INCLUDE_DIR = ../../include/cpp0x
INCLUDE = -I$(INCLUDE_DIR)

DEPENDS = $(INCLUDE_DIR)/*.h Makefile

benchmarkVectorAdd: benchmarkVectorAdd.cpp $(DEPENDS)
	g++45 -o benchmarkVectorAdd benchmarkVectorAdd.cpp -std=c++0x -O3 $(INCLUDE)

benchmarkVectorMultiply: benchmarkVectorMultiply.cpp $(DEPENDS)
	g++45 -o benchmarkVectorMultiply benchmarkVectorMultiply.cpp -std=c++0x -O3 $(INCLUDE)

benchmarkVectorAdd_fast: benchmarkVectorAdd_fast.cpp $(DEPENDS)
	g++45 -o benchmarkVectorAdd_fast benchmarkVectorAdd_fast.cpp -std=c++0x -O3 $(INCLUDE)

benchmarkVectorMultiply_fast: benchmarkVectorMultiply_fast.cpp $(DEPENDS)
	g++45 -o benchmarkVectorMultiply_fast benchmarkVectorMultiply_fast.cpp -std=c++0x -O3 $(INCLUDE)

benchmarkVectorAdd_omp: benchmarkVectorAdd.cpp $(DEPENDS)
	g++45 -o benchmarkVectorAdd_omp benchmarkVectorAdd.cpp -std=c++0x -O3 -fopenmp -lgomp $(INCLUDE)

benchmarkInverse: benchmarkInverse.cpp $(DEPENDS)
	g++45 -o benchmarkInverse benchmarkInverse.cpp -std=c++0x -O3 $(INCLUDE)

