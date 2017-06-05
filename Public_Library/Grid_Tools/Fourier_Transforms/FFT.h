//#####################################################################
// Copyright 2002-2007, Ron Fedkiw, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FFT
//#####################################################################
#ifndef __FFT__
#define __FFT__


#include <Core/Utilities/TYPE_UTILITIES.h>
#include <Grid_Tools/Arrays/GRIDS_UNIFORM_ARRAYS_FORWARD.h>
#include <Grid_Tools/Grids/GRID.h>
namespace PhysBAM{

template<class T> struct FFT_INDEX_HELPER{typedef int TYPE;};
template<class T,int d> struct FFT_INDEX_HELPER<VECTOR<T,d> >{typedef VECTOR<int,d> TYPE;};

// If TV is scalar, then arrays are normal arrays
// If TV is vector, then arrays are nd arrays.
template<class TV>
struct FFT
{
    typedef typename SCALAR_POLICY<TV>::TYPE T;
    typedef typename FFT_INDEX_HELPER<TV>::TYPE INDEX;
    typedef std::complex<T> C;
    void* plan[2][2];
    INDEX counts[2][2];

    FFT():plan{}{}
    FFT(const FFT&)=delete;
    FFT(FFT&&);
    FFT& operator=(const FFT&)=delete;
    FFT& operator=(FFT&&);
    ~FFT();

    // destroys input
    void Transform(ARRAY<T,INDEX>& in,ARRAY<C,INDEX>& out);
    void Transform(ARRAY<C,INDEX>& in,ARRAY<C,INDEX>& out);
    void Inverse_Transform(ARRAY<C,INDEX>& in,ARRAY<T,INDEX>& out);
    void Inverse_Transform(ARRAY<C,INDEX>& in,ARRAY<C,INDEX>& out);
};
}
#endif
