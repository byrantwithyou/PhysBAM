//#####################################################################
// Copyright 2009, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Jerry Talton, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_ADVECTION_MULTIPLE
//##################################################################### 
#ifndef __LEVELSET_ADVECTION_MULTIPLE__
#define __LEVELSET_ADVECTION_MULTIPLE__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <PhysBAM_Dynamics/Level_Sets/LEVELSET_ADVECTION.h>

namespace PhysBAM {
    
template<class T_GRID> class LEVELSET_MULTIPLE;

template<class T_GRID>
class LEVELSET_ADVECTION_MULTIPLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,TV_INT> T_ARRAYS_SCALAR;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
public:
    LEVELSET_MULTIPLE<T_GRID>& levelsets;
    ARRAY<LEVELSET_ADVECTION<TV> > levelset_advections;

    LEVELSET_ADVECTION_MULTIPLE(LEVELSET_MULTIPLE<T_GRID>& _levelsets)
        :levelsets(_levelsets)
    {}

    void Initialize()
    {for(int i=0;i<levelsets.levelsets.m;i++) levelset_advections.Append(LEVELSET_ADVECTION<TV>(levelsets.levelsets(i)));}

    void Set_Custom_Advection(ADVECTION<T,T,T_GRID>& advection_input)
    {for(int i=0;i<levelset_advections.m;i++)levelset_advections(i).Set_Custom_Advection(advection_input);}
    
    void Euler_Step(const T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,const int number_of_ghost_cells)
    {for(int i=0;i<levelset_advections.m;i++) levelset_advections(i).Euler_Step(face_velocities,dt,time,number_of_ghost_cells);}
    
    void Euler_Step(const ARRAY<TV,TV_INT>& velocity,const T dt,const T time,const int number_of_ghost_cells)
    {for(int i=0;i<levelset_advections.m;i++) levelset_advections(i).Euler_Step(velocity,dt,time,number_of_ghost_cells);}
    
    void Reinitialize()
    {for(int i=0;i<levelset_advections.m;i++) levelset_advections(i).Reinitialize();}
};

}
#endif
