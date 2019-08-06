mumps-seq-omp: mumps-seq-omp.o
	g++ -omumps-seq-omp -fopenmp mumps-seq-omp.o\
		-L$(MUMPS_ROOT)/lib -ldmumps -lmumps_common\
		-lmetis\
		-L$(MUMPS_ROOT)/PORD/lib/ -lpord\
		-lesmumps -lscotch -lscotcherr\
		-L$(MUMPS_ROOT)/libseq -lmpiseq\
		-L$(MKLROOT)/lib/intel64\
		-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread\
		-lgfortran

mumps-seq-omp.o: mumps-seq-omp.cpp
	g++ -std=c++14 -c -omumps-seq-omp.o -I$(MUMPS_ROOT)/include -I$(MUMPS_ROOT)/libseq mumps-seq-omp.cpp

clean:
	rm *.o
	rm mumps-seq-omp
