//#####################################################################
// Copyright 2009, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Jerry Talton, Elliot English.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_ADVECTION_MULTIPLE
//##################################################################### 
#ifndef __LEVELSET_ADVECTION_MULTIPLE__
#define __LEVELSET_ADVECTION_MULTIPLE__

#include <Geometry/Level_Sets/LEVELSET_POLICY.h>
#include <Dynamics/Level_Sets/LEVELSET_ADVECTION.h>

namespace PhysBAM {
    
template<class TV> class LEVELSET_MULTIPLE;

template<class TV>
class LEVELSET_ADVECTION_MULTIPLE
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef ARRAY<T,FACE_INDEX<TV::m> > T_FACE_ARRAYS_SCALAR;
public:
    LEVELSET_MULTIPLE<TV>& levelsets;
    ARRAY<LEVELSET_ADVECTION<TV> > levelset_advections;

    LEVELSET_ADVECTION_MULTIPLE(LEVELSET_MULTIPLE<TV>& _levelsets)
        :levelsets(_levelsets)
    {}

    void Initialize()
    {levelset_advections.Remove_All();for(int i=0;i<levelsets.levelsets.m;i++) levelset_advections.Append(LEVELSET_ADVECTION<TV>(levelsets.levelsets(i)));}

    void Set_Custom_Advection(ADVECTION<TV,T>& advection_input)
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
