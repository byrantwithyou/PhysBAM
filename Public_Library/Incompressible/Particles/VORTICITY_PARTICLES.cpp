//#####################################################################
// Copyright 2004-2009, Michael Lentine, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VORTICITY_PARTICLES
//#####################################################################
#include <Core/Vectors/VECTOR.h>
#include <Incompressible/Particles/VORTICITY_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> VORTICITY_PARTICLES<TV>::
VORTICITY_PARTICLES()
    :vorticity(0,0),radius(0,0)
{
    this->Add_Array("vorticity",&vorticity);
    this->Add_Array("radius",&radius);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> VORTICITY_PARTICLES<TV>::
~VORTICITY_PARTICLES()
{}
static int Initialize_Vorticity_Particles()
{
    return 0;
}
int initialize_vorticity_particles=Initialize_Vorticity_Particles();
//#####################################################################
template class VORTICITY_PARTICLES<VECTOR<float,1> >;
template class VORTICITY_PARTICLES<VECTOR<float,2> >;
template class VORTICITY_PARTICLES<VECTOR<float,3> >;
template class VORTICITY_PARTICLES<VECTOR<double,1> >;
template class VORTICITY_PARTICLES<VECTOR<double,2> >;
template class VORTICITY_PARTICLES<VECTOR<double,3> >;
}
