//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Tools/Krylov_Solvers/GMRES.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <Tools/Log/LOG.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Vectors/VECTOR.h>

using namespace PhysBAM;

typedef float RW;
typedef double T;

RANDOM_NUMBERS<T> ran;

int ok=0;

template<int n>
void Test()
{
    typedef VECTOR<T,n> TV;
    typedef KRYLOV_VECTOR_WRAPPER<T,TV> T_VEC;
    typedef MATRIX_SYSTEM<MATRIX<T,n>,T,T_VEC> T_MAT;

    MATRIX<T,n> mat;
    ran.Fill_Uniform(mat,-1,1);
    T_MAT M(mat);

    T_VEC x,b,z;
    ran.Fill_Uniform(x.v,-1,1);
    b.v=mat*x.v;

    GMRES<T> gmres;

    ARRAY<KRYLOV_VECTOR_BASE<T>*> av;
    gmres.Solve(M,z,b,av,1e-10,0,10);

    if((x.v-z.v).Magnitude()>1e-12)
        LOG::printf("%P * %P = %P   (vs %P)  (%P)   %P\n",mat,x.v,b.v,z.v,x.v-z.v,(x.v-z.v).Magnitude());
    else ok++;
}

int main(int argc, char* argv[])
{
    for(int i=0;i<20;i++){
        Test<1>();
        Test<2>();
        Test<3>();
        Test<4>();
        Test<5>();
        Test<6>();}

    LOG::printf("OK: %i\n",ok);
    return 0;
}
