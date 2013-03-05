//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL
//#####################################################################
#ifndef __SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL__
#define __SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL__
namespace PhysBAM{
template<class TV>
class SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL();
    ~SURFACE_RECONSTRUCTION_ANISOTROPIC_KERNAL();

    void Compute_Kernal_Centers_And_Transformation_And_Density(const ARRAY_VIEW<TV>& X,const ARRAY_VIEW<T>& m,const T h,const T r,const T lambda,const int N_eps,const T kr,const T ks,const T kn,ARRAY<TV>& Xbar,ARRAY<MATRIX<T,TV::m> >& G,ARRAY<T>& density) const;
    void Build_Scalar_Field(const ARRAY<TV>& Xbar,const ARRAY_VIEW<T>& m,const ARRAY_VIEW<T>& density,const ARRAY<MATRIX<T,TV::m> >& G,const GRID<TV>& grid,ARRAY<T,TV_INT>& phi) const;
//#####################################################################
};
}
#endif
