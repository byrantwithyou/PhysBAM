//#####################################################################
// Copyright 2009, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED
//##################################################################### 
#include <Core/Data_Structures/PAIR.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_ID.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<TV>::
UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED(const UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info,const int number_of_ghost_cells_input,const T_REGION& region_type_input,
    const int side_input,int axis_input)
    :grid(info.grid),collision_index(-1),collision_face_info(info.collision_face_info)
{
    Next();
}
//#####################################################################
namespace PhysBAM{
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<VECTOR<float,1> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<VECTOR<float,2> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<VECTOR<float,3> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<VECTOR<double,1> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<VECTOR<double,2> >;
template class UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<VECTOR<double,3> >;
}
