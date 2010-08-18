INCLUDE_DIR = ../../include/cpp0x
INCLUDE = -I$(INCLUDE_DIR)

DEPENDS = $(INCLUDE_DIR)/*.h Makefile

main: main.cpp $(DEPENDS)
	g++45 -o main main.cpp -std=c++0x $(INCLUDE)

testMultiplicationList: testMultiplicationList.cpp $(DEPENDS)
	g++45 -o testMultiplicationList testMultiplicationList.cpp -std=c++0x $(INCLUDE)

testInsertElementMEList: testInsertElementMEList.cpp $(DEPENDS)
	g++45 -o testInsertElementMEList testInsertElementMEList.cpp -std=c++0x $(INCLUDE)

testReverse: testReverse.cpp $(DEPENDS)
	g++45 -o testReverse testReverse.cpp -std=c++0x $(INCLUDE)

testInverse: testInverse.cpp $(DEPENDS)
	g++45 -o testInverse testInverse.cpp -std=c++0x $(INCLUDE)

testConfModel: testConfModel.cpp $(DEPENDS)
	g++45 -o testConfModel testConfModel.cpp -std=c++0x $(INCLUDE)

testSpringCGA: testSpringCGA.cpp $(DEPENDS)
	g++45 -O3 -o testSpringCGA testSpringCGA.cpp -losg -losgViewer -std=c++0x -Wfatal-errors $(INCLUDE) #-Winline

testNumericVersor: testNumericVersor.cpp $(DEPENDS)
	g++45 -o testNumericVersor testNumericVersor.cpp -std=c++0x $(INCLUDE)

testMonteCarloPi: testMonteCarloPi.cpp $(DEPENDS)
	g++45 -o testMonteCarloPi testMonteCarloPi.cpp -O3 -std=c++0x -fno-signed-zeros $(INCLUDE)

testHorizon: testHorizon.cpp $(DEPENDS)
	g++45 -o testHorizon testHorizon.cpp -O3 -std=c++0x -fno-signed-zeros $(INCLUDE)

testHorizon_fast: testHorizon_fast.cpp $(DEPENDS)
	g++45 -o testHorizon_fast testHorizon_fast.cpp -O3 -std=c++0x -fno-signed-zeros $(INCLUDE)

