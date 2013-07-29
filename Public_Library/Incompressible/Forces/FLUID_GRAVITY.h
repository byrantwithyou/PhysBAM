//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_GRAVITY
//#####################################################################
#ifndef __FLUID_GRAVITY__
#define __FLUID_GRAVITY__

#include <Tools/Grids_Uniform/FACE_INDEX.h>
#include <Incompressible/Forces/INCOMPRESSIBLE_FLUIDS_FORCES.h>
namespace PhysBAM{

template<class TV>
class FLUID_GRAVITY:public INCOMPRESSIBLE_FLUIDS_FORCES<TV>
{
    
    typedef typename TV::SCALAR T;
public:
    T gravity;
    TV gravity_direction;

    FLUID_GRAVITY(T gravity=(T)9.8,TV gravity_direction=-TV::Axis_Vector(TV::m==1?0:1))
        :gravity(gravity),gravity_direction(gravity_direction)
    {}

    virtual ~FLUID_GRAVITY()
    {}

//#####################################################################
    void Add_Explicit_Forces(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE
    {for(int axis=0;axis<TV::m;axis++) if(gravity_direction[axis]) face_velocities.Component(axis)+=dt*gravity*gravity_direction[axis];}
    void Add_Implicit_Forces_Projection(const GRID<TV>& grid,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities_ghost,ARRAY<T,FACE_INDEX<TV::m> >& face_velocities,const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Initialize_Grids(const GRID<TV>& grid) PHYSBAM_OVERRIDE {}
    T CFL(const GRID<TV>& grid,const ARRAY<T,FACE_INDEX<TV::m> >& face_velocities) PHYSBAM_OVERRIDE
    {return abs(gravity)*(gravity_direction/grid.DX()).Sum_Abs();}
//#####################################################################
};
}
#endif
