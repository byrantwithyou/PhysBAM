//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON
//#####################################################################
#ifndef __SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON__
#define __SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON__
namespace PhysBAM{
template<class TV>
class SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    T particle_radii;
    T radius_of_neighborhood;

    SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON(){}
    ~SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON(){}

    void Initialize(const T average_particle_spacing);
    void Build_Scalar_Field(const ARRAY_VIEW<TV>& X,const GRID<TV>& grid,ARRAY<T,TV_INT>& phi,const int smooth=0) const;
//#####################################################################
};
}
#endif
