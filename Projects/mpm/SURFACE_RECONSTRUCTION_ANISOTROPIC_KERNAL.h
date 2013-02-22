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

    static void Compute_Kernal_Centers_And_Transformation(const GEOMETRY_PARTICLES<TV>& particles,ARRAY<TV>& Xbar,ARRAY<MATRIX<T,TV::m> >& G);
//#####################################################################
};
}
#endif
