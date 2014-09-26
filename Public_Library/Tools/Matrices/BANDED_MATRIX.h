//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BANDED_MATRIX
//#####################################################################
#ifndef __BANDED_MATRIX__
#define __BANDED_MATRIX__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T,int w>
class BANDED_MATRIX
{
    typedef VECTOR<T,2> TV;
public:
    typedef VECTOR<T,w> ROW;
    ARRAY<ROW> A;
    int diagonal_column;

    BANDED_MATRIX(const int m=0)
        :A(m)
    {}

    int Size() const
    {return A.Size();}

    void Resize(const int m_new)
    {A.Resize(m_new);}

//#####################################################################
    VECTOR<T,2> Givens(ROW& x,ROW& y); // make y(0)=0
    template<class T2> void Givens_Shift_Diagonal(ARRAY<T2>& u);
    template<class T2> void Backsolve(ARRAY<T2>& u) const;
    template<class T2> void QR_Solve(ARRAY<T2>& u);
    template<class T2> void QR_Solve(ARRAY<ARRAY<T2>>& u);
//#####################################################################
};
}
#endif
