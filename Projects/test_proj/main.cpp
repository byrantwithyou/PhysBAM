//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <lapacke.h>

using namespace PhysBAM;

typedef float RW;
typedef double T;
typedef VECTOR<T,3> TV;
typedef VECTOR<int,TV::m> TV_INT;

int main(int argc, char* argv[])
{
    PARSE_ARGS parse_args(argc,argv);
    parse_args.Parse();
    
    MATRIX_MXN<T> M(2,2);

    M(0,0)=1;
    M(0,1)=2;
    M(1,0)=2;
    M(1,1)=0;

    ARRAY<int> idiv(M.m);
    
    lapack_int ret_f = LAPACKE_dsytrf( LAPACK_ROW_MAJOR, 'U', M.m, M.x.Get_Array_Pointer(), M.m, idiv.Get_Array_Pointer() );
    lapack_int ret_i = LAPACKE_dsytri( LAPACK_ROW_MAJOR, 'U', M.m, M.x.Get_Array_Pointer(), M.m, idiv.Get_Array_Pointer() );
    for(int r=0;r<M.m;r++)
        for(int c=0;c<r;c++)
            M(r,c)=M(c,r);
    
    printf("ret: %i %i\n", ret_f, ret_i);

    LOG::printf("%P\n",M);
    
    // double M[3][3] = { {1 , 2 , 3},
    //                    {4 , 5 , 6},
    //                    {7 , 8 , 9}}
    // double pivotArray[3]; //since our matrix has three rows
    // int errorHandler;
    // double lapackWorkspace[9];

    // // dgetrf(M,N,A,LDA,IPIV,INFO) means invert LDA columns of an M by N matrix 
    // // called A, sending the pivot indices to IPIV, and spitting error 
    // // information to INFO.
    // // also don't forget (like I did) that when you pass a two-dimensional array
    // // to a function you need to specify the number of "rows"
    // dgetrf_(3,3,M[3][],3,pivotArray[3],&errorHandler);
    // //some sort of error check

    // dgetri_(3,M[3][],3,pivotArray[3],9,lapackWorkspace,&errorHandler);
    // //another error check
    
    return 0;
}

