//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function Annotate
//#####################################################################
#ifndef __Annotate__
#define __Annotate__

#include <PhysBAM_Tools/Math_Tools/ONE.h>
#include <PhysBAM_Tools/Math_Tools/ZERO.h>
#include <PhysBAM_Tools/Particles/PARTICLES_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_OBJECT_FORWARD.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_FORWARD.h>
#include "QUANTITY.h"
#include <boost/type_traits/is_integral.hpp>
namespace PhysBAM{

template<class TV> class SOLIDS_PARAMETERS;

namespace UNITS{

static const SIMPLE_UNIT one(0,0,0);
static const SIMPLE_UNIT length(1,0,0);
static const SIMPLE_UNIT mass(0,1,0);
static const SIMPLE_UNIT time(0,0,1);
static const SIMPLE_UNIT velocity(1,0,-1);
static const SIMPLE_UNIT inertia_tensor(2,1,0);

void Annotate(const ZERO,const SIMPLE_UNIT& unit){}
void Annotate(const ONE,const SIMPLE_UNIT& unit){}

template<class T> void Annotate(const QUANTITY<T>& x,const SIMPLE_UNIT& unit)
{Unify(x.unit,unit);}

template<class T_ARRAY> void Annotate(const T_ARRAY& array)
{for(int i=0;i<array.Size();i++) Annotate(array(i));}

template<class T_ARRAY> void Annotate(const T_ARRAY& array,const SIMPLE_UNIT& unit)
{for(int i=0;i<array.Size();i++) Annotate(array(i),unit);}

template<class T> void Annotate(const DIAGONAL_MATRIX<T,2>& matrix,const SIMPLE_UNIT& unit)
{Annotate(matrix.x11,unit);Annotate(matrix.x22,unit);}

template<class T> void Annotate(const DIAGONAL_MATRIX<T,3>& matrix,const SIMPLE_UNIT& unit)
{Annotate(matrix.x11,unit);Annotate(matrix.x22,unit);Annotate(matrix.x33,unit);}

template<class TV,bool world_space> void Annotate(const RIGID_BODY_MASS<TV,world_space>& x)
{Annotate(x.mass,mass);Annotate(x.inertia_tensor,inertia_tensor);}

template<class T> void Annotate(const PARTICLE_ATTRIBUTE<T>&)
{}

template<class TV> void Annotate(const PARTICLE_ATTRIBUTE<RIGID_BODY_MASS<TV> >& m)
{Annotate(m.array);}

template<class TV> void Annotate(const PARTICLE_POSITION_ATTRIBUTE<TV>& X)
{Annotate(X.array,length);}

template<class TV> void Annotate(const PARTICLE_VELOCITY_ATTRIBUTE<TV>& V)
{Annotate(V.array,velocity);}

template<class TV> void Annotate(const PARTICLE_MASS_ATTRIBUTE<TV>& m)
{Annotate(m.array,mass);Annotate(m.one_over_mass,mass.Inverse());Annotate(m.effective_mass,mass);Annotate(m.one_over_effective_mass,mass.Inverse());}

template<class TV> void Annotate(const SOLIDS_PARTICLES<TV>& particles)
{Annotate(particles.X);Annotate(particles.V);Annotate(particles.mass);}

template<class TV> void Annotate(const DEFORMABLE_OBJECT<TV>& deformable_object)
{Annotate(deformable_object.particles);}

template<class TV> void Annotate(const SOLIDS_PARAMETERS<TV>& solids_parameters)
{Annotate(solids_parameters.deformable_object);}

}
}
#endif
