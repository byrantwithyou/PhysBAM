//#####################################################################
// Copyright 2014.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION_UPWIND 
//#####################################################################
#ifndef __ADVECTION_UPWIND__
#define __ADVECTION_UPWIND__

#include <Tools/Advection/ADVECTION.h>
#include <Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
#include <Tools/Vectors/VECTOR.h>
#include <boost/function.hpp>
namespace PhysBAM{

template<class TV> class LEVELSET;

template<class TV,class T2,class Z_INTERP,class T_AVERAGING=AVERAGING_UNIFORM<TV>,
    class U_INTERP=LINEAR_INTERPOLATION_UNIFORM<TV,typename TV::SCALAR>,
    class T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV> >
class ADVECTION_UPWIND:public ADVECTION<TV,T2,T_FACE_LOOKUP>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef INTERPOLATION_UNIFORM<TV,T2,T_FACE_LOOKUP> BASE;

    const LEVELSET<TV>& levelset;
    const T_FACE_LOOKUP& face_velocities;
    T max_in,max_out,time,dt;
    boost::function<T2(const TV& X,T time)> bc_Z;

    ADVECTION_UPWIND(const LEVELSET<TV>& levelset,const T_FACE_LOOKUP& face_velocities,T max_in,T max_out,
        boost::function<T2(const TV& X,T time)> bc_Z,T time);
    virtual ~ADVECTION_UPWIND();

    void Update_Advection_Equation_Cell_Lookup(const GRID<TV>& grid,ARRAY<T2,TV_INT>& Z,const ARRAY<T2,TV_INT>& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,BOUNDARY<TV,T2>& boundary,const T dt,const T time,
        const ARRAY<T2,TV_INT>* Z_min_ghost=0,const ARRAY<T2,TV_INT>* Z_max_ghost=0,ARRAY<T2,TV_INT>* Z_min=0,ARRAY<T2,TV_INT>* Z_max=0) PHYSBAM_OVERRIDE;

    T2 Clamped_To_Array(const GRID<TV>& grid,const ARRAYS_ND_BASE<T2,TV_INT>& Z,const TV& X) const;
//#####################################################################
};
}
#endif
