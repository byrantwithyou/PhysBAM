//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Function Add_Velocity_Independent_Forces
//#####################################################################
template<class TV> void DEFORMABLE_GRAVITY<TV>::
Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const
{
    if(!gravity) return;
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();F(p)+=gravity*particles.mass(p)*downward_direction;}
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_GRAVITY<TV>::
Potential_Energy(int p,const T time) const
{
    TV acceleration=gravity*downward_direction;
    return -particles.mass(p)*TV::Dot_Product(particles.X(p),acceleration);
}
//#####################################################################
// Function Potential_Energy
//#####################################################################
template<class TV> typename TV::SCALAR DEFORMABLE_GRAVITY<TV>::
Potential_Energy(const T time) const
{
    T potential_energy=0;
    for(ELEMENT_ITERATOR iterator(force_particles);iterator.Valid();iterator.Next()){int p=iterator.Data();
        potential_energy+=Potential_Energy(p,time);}
    return potential_energy;
}
//#####################################################################
// Function Update_Position_Based_State
//#####################################################################
template<class TV> void DEFORMABLE_GRAVITY<TV>::
Update_Position_Based_State(const T time,const bool is_position_update)
{
}
//#####################################################################
template class DEFORMABLE_GRAVITY<VECTOR<float,1> >;
template class DEFORMABLE_GRAVITY<VECTOR<float,2> >;
template class DEFORMABLE_GRAVITY<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DEFORMABLE_GRAVITY<VECTOR<double,1> >;
template class DEFORMABLE_GRAVITY<VECTOR<double,2> >;
template class DEFORMABLE_GRAVITY<VECTOR<double,3> >;
#endif
