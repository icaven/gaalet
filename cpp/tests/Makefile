INCLUDE_DIR = ../../include/cpp
INCLUDE = -I$(INCLUDE_DIR)

DEPENDS = $(INCLUDE_DIR) Makefile

main: main.cpp $(DEPENDS) 
	g++45 -o main main.cpp $(INCLUDE)

testMultiplicationList: testMultiplicationList.cpp $(DEPENDS) 
	g++45 -o testMultiplicationList testMultiplicationList.cpp $(INCLUDE)

benchmarkVectorAdd: benchmarkVectorAdd.cpp $(DEPENDS) 
	g++45 -o benchmarkVectorAdd benchmarkVectorAdd.cpp -O3 $(INCLUDE)

benchmarkVectorAdd_fast: benchmarkVectorAdd_fast.cpp $(DEPENDS) 
	g++45 -o benchmarkVectorAdd_fast benchmarkVectorAdd_fast.cpp -O3 $(INCLUDE)

testInsertElementMEList: testInsertElementMEList.cpp $(DEPENDS) 
	g++45 -o testInsertElementMEList testInsertElementMEList.cpp $(INCLUDE)

testReverse: testReverse.cpp $(DEPENDS)
	g++45 -o testReverse testReverse.cpp $(INCLUDE)

testInverse: testInverse.cpp $(DEPENDS)
	g++45 -o testInverse testInverse.cpp $(INCLUDE)

testConfModel: testConfModel.cpp $(DEPENDS)
	g++45 -o testConfModel testConfModel.cpp $(INCLUDE)

testMonteCarloPi: testMonteCarloPi.cpp $(DEPENDS)
	g++45 -o testMonteCarloPi testMonteCarloPi.cpp -O3 -fno-signed-zeros $(INCLUDE)

testSpringCGA: testSpringCGA.cpp $(DEPENDS)
	g++45 -O3 -o testSpringCGA testSpringCGA.cpp -losg -losgViewer -Wfatal-errors $(INCLUDE)

testConfModelCuda: testConfModelCuda.cu $(DEPENDS)
	nvcc -o testConfModelCuda testConfModelCuda.cu $(INCLUDE)

testMonteCarloPiCuda: testMonteCarloPiCuda.cu $(DEPENDS)
	nvcc -o testMonteCarloPiCuda testMonteCarloPiCuda.cu -O3 -arch=sm_11 $(INCLUDE)

testRotationCuda: testRotationCuda.cu $(DEPENDS)
	nvcc -o testRotationCuda testRotationCuda.cu -O3 $(INCLUDE)

