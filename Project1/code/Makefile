Point:testPoint.o Point.o
	g++ -o testPoint testPoint.o Point.o
Point.o:Point.cpp Point.h
	g++ -c Point.cpp
testPoint.o:testPoint.cpp
	g++ -c testPoint.cpp
AdamBashforce:Point.o testAB.o
	g++ -o testAB testAB.o Point.o
testAB.o:testAB.cpp
	g++ -c testAB.cpp
AdamMoulton:Point.o testAM.o
	g++ -o testAM testAM.o Point.o
testAM.o:testAM.cpp
	g++ -c testAM.cpp
BDF:Point.o testBDF.o
	g++ -o testBDF testBDF.o Point.o
testBDF.o:testBDF.cpp
	g++ -c testBDF.cpp
RungeKutta:Point.o testRK.o
	g++ -o testRK testRK.o Point.o
testRK.o:testRK.cpp
	g++ -c testRK.cpp
testall:testall.o Point.o
	g++ -o testall testall.o Point.o
testall.o:testall.cpp
	g++ -c testall.cpp
run:ODESolver.o Point.o testAM.o testRK.o testBDF.o testAB.o
	g++ -o ODESolver ODESolver.o Point.o
	g++ -o testAM testAM.o Point.o
	g++ -o testAB testAB.o Point.o
	g++ -o testBDF testBDF.o Point.o
	g++ -o testRK testRK.o Point.o
ODESolver.o:ODESolver.cpp
	g++ -c ODESolver.cpp
clean:
	rm *.o
	rm testAB testAM testBDF testRK ODESolver Trajectory.m
