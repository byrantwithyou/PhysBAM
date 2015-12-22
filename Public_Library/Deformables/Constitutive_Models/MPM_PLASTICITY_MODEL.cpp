//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Deformables/Constitutive_Models/MPM_PLASTICITY_MODEL.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_PLASTICITY_MODEL<TV>::
MPM_PLASTICITY_MODEL(MPM_PARTICLES<TV>& particles)
    :particles(particles),use_implicit(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_PLASTICITY_MODEL<TV>::
~MPM_PLASTICITY_MODEL()
{
}
//#####################################################################
// Function Updated_Particle
//#####################################################################
template<class TV> bool MPM_PLASTICITY_MODEL<TV>::
Update_Particle(int p) const
{
    MATRIX<T,TV::m> Fe=particles.F(p),U,V;
    DIAGONAL_MATRIX<T,TV::m> singular_values;
    Fe.Fast_Singular_Value_Decomposition(U,singular_values,V);
    TV strain;
    if(Compute(strain,0,0,0,0,singular_values.x,true,p)){
        singular_values.x=strain;
        particles.F(p)=(U*singular_values).Times_Transpose(V);
        particles.Fp(p)=V*singular_values.Inverse().Times_Transpose(U)*Fe*particles.Fp(p);
        return true;}
    return false;
}
template class MPM_PLASTICITY_MODEL<VECTOR<float,2> >;
template class MPM_PLASTICITY_MODEL<VECTOR<float,3> >;
template class MPM_PLASTICITY_MODEL<VECTOR<double,2> >;
template class MPM_PLASTICITY_MODEL<VECTOR<double,3> >;
}
