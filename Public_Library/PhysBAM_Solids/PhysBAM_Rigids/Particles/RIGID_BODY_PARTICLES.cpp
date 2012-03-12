//#####################################################################
// Copyright 2006-2009, Craig Schroeder, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_0X0.h>
#include <PhysBAM_Tools/Matrices/MATRIX_1X1.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGIDS_PARTICLES_FORWARD.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_BODY_PARTICLES<TV>::
RIGID_BODY_PARTICLES()
    :angular_momentum(0,0),mass(0,0),inertia_tensor(0,0),kinematic(0,0)
{
    array_collection->Add_Array(ATTRIBUTE_ID_ANGULAR_MOMENTUM,&angular_momentum);
    array_collection->Add_Array(ATTRIBUTE_ID_RIGID_MASS,&mass);
    array_collection->Add_Array(ATTRIBUTE_ID_RIGID_INERTIA_TENSOR,&inertia_tensor);
    array_collection->Add_Array(ATTRIBUTE_ID_KINEMATIC,&kinematic);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_BODY_PARTICLES<TV>::
~RIGID_BODY_PARTICLES()
{
    Delete_All_Particles();
}
//#####################################################################
template class RIGID_BODY_PARTICLES<VECTOR<float,1> >;
template class RIGID_BODY_PARTICLES<VECTOR<float,2> >;
template class RIGID_BODY_PARTICLES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_BODY_PARTICLES<VECTOR<double,1> >;
template class RIGID_BODY_PARTICLES<VECTOR<double,2> >;
template class RIGID_BODY_PARTICLES<VECTOR<double,3> >;
#endif
}
