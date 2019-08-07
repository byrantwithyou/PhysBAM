mumps-seq-omp: mumps-seq-omp.o
	g++ -O2 -omumps-seq-omp -fopenmp mumps-seq-omp.o\
		-L$(MUMPS_ROOT)/lib -ldmumps -lmumps_common\
		-lmetis\
		-L$(MUMPS_ROOT)/PORD/lib/ -lpord\
		-lesmumps -lscotch -lscotcherr\
		-L$(MUMPS_ROOT)/libseq -lmpiseq\
		-L$(MKLROOT)/lib/intel64\
		-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread\
		-lgfortran

mumps-seq-omp.o: mumps-seq-omp.cpp
	g++ -std=c++14 -O2 -c -fopenmp -DBLR_MT -omumps-seq-omp.o -I$(MUMPS_ROOT)/include -I$(MUMPS_ROOT)/libseq mumps-seq-omp.cpp\
		-Wno-unused-result

clean:
	rm *.o
	rm mumps-seq-omp
	rm mumps-par

mumps-par: mumps-par.o
	mpif90 -O2 -omumps-par -fopenmp mumps-par.o\
		-L$(MUMPS_ROOT)/lib -ldmumps -lmumps_common\
		-lmetis\
		-L$(MUMPS_ROOT)/PORD/lib/ -lpord\
		-lesmumps -lscotch -lscotcherr\
		-L$(MKLROOT)/lib/intel64\
		-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread\
		-lstdc++ -lmpi_cxx

mumps-par.o: mumps-par.cpp
	mpic++ -std=c++14 -O2 -c -fopenmp -DBLR_MT -omumps-par.o -I$(MUMPS_ROOT)/include mumps-par.cpp\
		-Wno-unused-result
