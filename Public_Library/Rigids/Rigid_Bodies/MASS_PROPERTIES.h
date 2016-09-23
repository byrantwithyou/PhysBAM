//#####################################################################
// Copyright 2003-2007, Ron Fedkiw, Eran Guendelman, Don Hatch, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MASS_PROPERTIES
//#####################################################################
#ifndef __MASS_PROPERTIES__
#define __MASS_PROPERTIES__

#include <Core/Matrices/SYMMETRIC_MATRIX.h>
#include <Core/Utilities/NONCOPYABLE.h>
#include <Core/Utilities/TYPE_UTILITIES.h>
#include <Tools/Particles/PARTICLES.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{

template <class TV> class FRAME;

template<class TV,int d=TV::m-1>
class MASS_PROPERTIES:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d>::OBJECT T_SIMPLICIAL_OBJECT;
    struct DUMMY_IMPLEMENTATION{};struct NORMAL_IMPLEMENTATION{};
    typedef typename conditional<d==0,DUMMY_IMPLEMENTATION,NORMAL_IMPLEMENTATION>::type ACCESS_IMPLEMENTATION;
private:
    const T_SIMPLICIAL_OBJECT& object;
    T mass,density;
    bool use_mass;
    T volume;
    TV center;
    SYMMETRIC_MATRIX<T,TV::SPIN::m> inertia_tensor_over_density;
public:

    MASS_PROPERTIES(const T_SIMPLICIAL_OBJECT& object,const bool thin_shell);

    void Set_Mass(const T mass_input)
    {mass=mass_input;use_mass=true;}

    void Set_Density(const T density_input)
    {density=density_input;use_mass=false;}

    T Volume() const
    {return volume;}

    const TV& Center() const
    {return center;}

    T Mass() const
    {return use_mass?mass:density*volume;}

    T Density() const
    {return use_mass?mass/volume:density;}

    SYMMETRIC_MATRIX<T,TV::SPIN::m> Inertia_Tensor() const
    {return Density()*inertia_tensor_over_density;}

//#####################################################################
    static T Thin_Shell_Volume(const T_SIMPLICIAL_OBJECT& object);
    void Transform_To_Object_Frame(FRAME<TV>& frame,DIAGONAL_MATRIX<T,TV::SPIN::m>& object_space_inertia_tensor) const;
    void Transform_To_Object_Frame(FRAME<TV>& frame,DIAGONAL_MATRIX<T,TV::SPIN::m>& object_space_inertia_tensor,GEOMETRY_PARTICLES<TV>& point_cloud) const;
private:
    template<bool thin_shell> void Compute_Properties(MASS_PROPERTIES<TV,d>&,NORMAL_IMPLEMENTATION);
    template<bool thin_shell> void Compute_Properties(MASS_PROPERTIES<TV,d>&,DUMMY_IMPLEMENTATION);
//#####################################################################
};
}
#endif
