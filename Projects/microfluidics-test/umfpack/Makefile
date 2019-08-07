clean:
	rm *.o
	rm ldivide

ldivide: ldivide.o
	g++ -O2 -fopenmp -oldivide ldivide.o -L${SUITESPARSE_ROOT}/UMFPACK/Lib -lumfpack -lamd\
		-L${SUITESPARSE_ROOT}/SuiteSparse_config -lsuitesparseconfig\
		-L${SUITESPARSE_ROOT}/CHOLMOD/Lib -lcholmod\
		-L${SUITESPARSE_ROOT}/CCOLAMD/Lib -lccolamd\
		-L${SUITESPARSE_ROOT}/COLAMD/Lib -lcolamd\
		-L${SUITESPARSE_ROOT}/CAMD/Lib -lcamd\
		-L${SUITESPARSE_ROOT}/lib -lmetis\
		-L$(MKLROOT)/lib/intel64\
		-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

ldivide.o: ldivide.cpp
	g++ -O2 -std=c++14 -c ldivide.cpp -fopenmp -Wno-unused-result -I${SUITESPARSE_ROOT}/UMFPACK/Include
