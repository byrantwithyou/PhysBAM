//#####################################################################
// Copyright 2007-2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//###################################################################
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/PRIMITIVE_MATRICES.h>
#include <Core/Utilities/PROCESS_UTILITIES.h>
#include <Tools/Tensors/PRIMITIVE_TENSORS.h>
#include <Tools/Tensors/TENSOR.h>
#include "TENSOR_TESTS.h"
using namespace PhysBAM;
//#####################################################################
// Function Run_Test
//#####################################################################
template<class T> typename TENSOR_TESTS<T>::TEST_RESULT TENSOR_TESTS<T>::
Run_Test(int n)
{
    status=true;
    Generate_Tests();
    return status?success:failure;
}
//#####################################################################
// Function Number_Of_Tests
//#####################################################################
template<class T> int TENSOR_TESTS<T>::
Number_Of_Tests() const
{
    return 1;
}
template<class T,int mx,int m,int n,int p,int q>
struct GENERATE_TESTS4
{
    static void F(TENSOR_TESTS<T>& tt)
    {
        GENERATE_TESTS4<T,mx,m,n,p,q-1>::F(tt);
        tt.template Test_4<m,n,p,q>();
    }
};
template<class T,int mx,int m,int n,int p>
struct GENERATE_TESTS4<T,mx,m,n,p,0>
{
    static void F(TENSOR_TESTS<T>& tt){}
};
template<class T,int mx,int m,int n,int p>
struct GENERATE_TESTS3
{
    static void F(TENSOR_TESTS<T>& tt)
    {
        GENERATE_TESTS3<T,mx,m,n,p-1>::F(tt);
        GENERATE_TESTS4<T,mx,m,n,p,mx>::F(tt);
        tt.template Test_3<m,n,p>();
    }
};
template<class T,int mx,int m,int n>
struct GENERATE_TESTS3<T,mx,m,n,0>
{
    static void F(TENSOR_TESTS<T>& tt){}
};
template<class T,int mx,int m,int n>
struct GENERATE_TESTS2
{
    static void F(TENSOR_TESTS<T>& tt)
    {
        GENERATE_TESTS2<T,mx,m,n-1>::F(tt);
        GENERATE_TESTS3<T,mx,m,n,mx>::F(tt);
        tt.template Test_2<m,n>();
    }
};
template<class T,int mx,int m>
struct GENERATE_TESTS2<T,mx,m,0>
{
    static void F(TENSOR_TESTS<T>& tt){}
};
template<class T,int mx,int m>
struct GENERATE_TESTS1
{
    static void F(TENSOR_TESTS<T>& tt)
    {
        GENERATE_TESTS1<T,mx,m-1>::F(tt);
        GENERATE_TESTS2<T,mx,m,mx>::F(tt);
        tt.template Test_1<m>();
    }
};
template<class T,int mx>
struct GENERATE_TESTS1<T,mx,0>
{
    static void F(TENSOR_TESTS<T>& tt){}
};
//#####################################################################
// Function Test_Contract
//#####################################################################
template<class T> void TENSOR_TESTS<T>::
Generate_Tests()
{
    GENERATE_TESTS1<T,3,3>::F(*this);
    Test_0();
}
TENSOR_TESTS<double> tensor_tests;
