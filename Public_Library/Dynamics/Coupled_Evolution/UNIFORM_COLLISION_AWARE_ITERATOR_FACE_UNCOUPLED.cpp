//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED
//##################################################################### 
#include <Core/Data_Structures/PAIR.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_ID.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<TV>::
UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info,const int number_of_ghost_cells_input,bool use_outside_input,
    const T_REGION& region_type_input,const int side_input,int axis_input)
    :BASE(info.grid,number_of_ghost_cells_input,region_type_input,side_input,axis_input),outside_fluid(*info.outside_fluid),
    collision_index(0),collision_face_info(info.collision_face_info),scan_end(INT_MIN),use_outside(use_outside_input)
{
    index(TV::m-1)-=2;
    Next_Fluid();
}
template<class TV> UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<TV>::
~UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED()
{
}
//#####################################################################
// Function Next_Helper
//#####################################################################
template<class TV> void UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<TV>::
Next_Helper()
{
    BASE::Next();
    for(int c;(c=Compare_Collision_Index())<=0;collision_index++) if(!c) BASE::Next();

    scan_end=grid.counts(TV::m-1)+(TV::m==axis)+number_of_ghost_cells+(region_type==GRID<TV>::WHOLE_REGION?1:0);
    if(collision_index<collision_face_info.Size()){
        const COLLISION_FACE_INFO<TV>& cfi=collision_face_info(collision_index);
        if(axis==cfi.axis && index.Remove_Index(TV::m-1)==cfi.index.Remove_Index(TV::m-1))
            scan_end=cfi.index(TV::m-1);}
}
//#####################################################################
// Function Compare_Collision_Index
//#####################################################################
template<class TV> int UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<TV>::
Compare_Collision_Index() const
{
    if(collision_index>=collision_face_info.Size()) return 1;
    const COLLISION_FACE_INFO<TV>& cfi=collision_face_info(collision_index);
    if(cfi.axis<axis) return -1;
    if(cfi.axis>axis) return 1;
    for(int i=0;i<TV::m;i++){
        if(cfi.index(i)<index(i)) return -1;
        if(cfi.index(i)>index(i)) return 1;}
    return 0;
}
//#####################################################################
namespace PhysBAM{
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<VECTOR<float,1> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<VECTOR<float,2> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<VECTOR<float,3> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<VECTOR<double,1> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<VECTOR<double,2> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<VECTOR<double,3> >;
}
