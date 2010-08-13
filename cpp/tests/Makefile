INCLUDE_DIR = ../../include
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

testSpringCGA: testSpringCGA.cpp $(DEPENDS)
	g++45 -O3 -o testSpringCGA testSpringCGA.cpp -losg -losgViewer -Wfatal-errors $(INCLUDE)


