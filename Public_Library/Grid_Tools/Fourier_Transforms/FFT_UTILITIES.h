//#####################################################################
// Copyright 2002-2006, Ron Fedkiw, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FFT_UTILITIES
//#####################################################################
#ifndef __FFT_UTILITIES__
#define __FFT_UTILITIES__

#include <Core/Arrays/ARRAY.h>
#include <complex>
namespace PhysBAM{
template<class T,int d> void
Enforce_Real_Valued_Symmetry(ARRAY<std::complex<T>,VECTOR<int,d> >& u_hat);
template<class T,int d> void
First_Derivatives(const ARRAY<std::complex<T>,VECTOR<int,d> >& u_hat,
    VECTOR<ARRAY<std::complex<T>,VECTOR<int,d> >,d>& du_hat,
    const VECTOR<T,d>& domain_size);
template<class T,int d> void
Make_Divergence_Free(VECTOR<ARRAY<std::complex<T>,VECTOR<int,d> >,d>& u_hat,
    const VECTOR<T,d>& domain_size);
template<class T,int d> void
Filter_High_Frequencies(ARRAY<std::complex<T>,VECTOR<int,d> >& u_hat,T scale=(T)1);
//#####################################################################
}
#endif
