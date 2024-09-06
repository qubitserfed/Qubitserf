test:
	g++ -g -c -std=c++17 src/linear_algebra.cpp 			-o build/linear_algebra_dev.o
	g++ -g -c -std=c++17 src/combinatorics.cpp 				-o build/combinatorics_dev.o
	g++ -g -c -std=c++17 src/quantum_utilities.cpp 			-o build/quantum_utilities_dev.o
	g++ -g -c -std=c++17 src/brouwer_zimmerman.cpp 			-o build/brouwer_zimmerman_dev.o
	g++ -g -c -std=c++17 src/interface.cpp 					-o build/interface_dev.o
	g++ -g -c -std=c++17 src/testing.cpp 					-o build/testing_dev.o
	g++ -g -std=c++17 build/linear_algebra_dev.o build/combinatorics_dev.o build/quantum_utilities_dev.o build/brouwer_zimmerman_dev.o build/testing_dev.o -o build/test
	rm build/*.o

interface:
	g++ -c -Ofast -std=c++17 src/linear_algebra.cpp    		-o build/linear_algebra.o
	g++ -c -Ofast -std=c++17 src/combinatorics.cpp 	   		-o build/combinatorics.o
	g++ -c -Ofast -std=c++17 src/quantum_utilities.cpp 		-o build/quantum_utilities.o
	g++ -c -Ofast -std=c++17 src/brouwer_zimmerman.cpp 		-o build/brouwer_zimmerman.o
	g++ -c -Ofast -std=c++17 src/interface.cpp 		   		-o build/interface.o

	g++ -g -std=c++17 build/linear_algebra.o build/combinatorics.o build/quantum_utilities.o build/brouwer_zimmerman.o build/interface.o -o build/interface
	rm build/*.o
