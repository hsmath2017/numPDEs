Mesh:testMesh.o
	g++ -o testMesh testMesh.o
testMesh.o:testMesh.cpp
	g++ -c testMesh.cpp
Jacob:testWeightedJacob.o
	g++ -o testJacob testWeightedJacob.o
testWeightedJacob.o:testWeightedJacob.cpp
	g++ -c testWeightedJacob.cpp
VCycle:testVCycle.o
	g++ -o testVCycle testVCycle.o
testVCycle.o:testVCycle.cpp
	g++ -c testVCycle.cpp
FMG:testFMG.o
	g++ -o testFMG testFMG.o
testFMG.o:testFMG.cpp
	g++ -c testFMG.cpp
run:testMesh.o testWeightedJacob.o testVCycle.o testFMG.o main.cpp
	g++ -o testWeightedJacob testWeightedJacob.o
	g++ -o testVCycle testVCycle.o
	g++ -o testFMG testFMG.o
	g++ -o testMesh testMesh.o
	g++ -o main main.cpp
clean:
	rm *.o
	rm testMesh Result.m testVCycle main testFMG testWeightedJacob
