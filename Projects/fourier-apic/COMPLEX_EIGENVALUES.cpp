#include <Core/Matrices/MATRIX_MXN.h>
using namespace PhysBAM;
typedef double T;
#ifdef USE_LAPACK
#define lapack_complex_double std::complex<double>
#include <lapacke.h>
void Compute_Eigenvalues(MATRIX_MXN<std::complex<T> > M,ARRAY<std::complex<T> >& eig)
{
    PHYSBAM_ASSERT(M.m==M.n);
    char jobvl = 'N';
    char jobvr = 'N';
    int n = M.m, size_work=10*n;
    ARRAY<std::complex<T> > work(size_work);
    ARRAY<T> rwork(2*n);
    int info = 0;

    LAPACK_zgeev(&jobvl, &jobvr, &n, &M(0,0), &n, &eig(0), 0, &n, 0, &n,
        &work(0), &size_work, &rwork(0), &info);

    if(info) printf("ZGEEV FAILED: %i\n",info);
    PHYSBAM_ASSERT(info==0);
}
#else
void Compute_Eigenvalues(MATRIX_MXN<std::complex<T> > M,ARRAY<std::complex<T> >& eig)
{
    PHYSBAM_FATAL_ERROR("Requires USE_LAPACK");
}
#endif
