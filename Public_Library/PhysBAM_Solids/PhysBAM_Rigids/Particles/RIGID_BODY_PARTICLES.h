//#####################################################################
// Copyright 2006-2009, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_PARTICLES
//#####################################################################
#ifndef __RIGID_BODY_PARTICLES__
#define __RIGID_BODY_PARTICLES__

#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Particles/PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/RIGID_GEOMETRY_PARTICLES.h>
namespace PhysBAM{

template<class TV>
class RIGID_BODY_PARTICLES:public CLONEABLE<RIGID_BODY_PARTICLES<TV>,RIGID_GEOMETRY_PARTICLES<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<RIGID_BODY_PARTICLES<TV>,RIGID_GEOMETRY_PARTICLES<TV> > BASE;
public:
    using BASE::rigid_geometry;using BASE::Delete_All_Particles;using BASE::Add_Array;

    ARRAY_VIEW<typename TV::SPIN> angular_momentum;
    ARRAY_VIEW<T> mass;
    ARRAY_VIEW<DIAGONAL_MATRIX<T,TV::SPIN::m> > inertia_tensor;
    ARRAY_VIEW<bool> kinematic;

    RIGID_BODY_PARTICLES();
    virtual ~RIGID_BODY_PARTICLES();

private:
    void Delete_Particle(const int index)
    {PHYSBAM_FATAL_ERROR();}

    void Delete_Particles_On_Deletion_List(const bool preserve_order=false)
    {PHYSBAM_FATAL_ERROR();}

    void Add_To_Deletion_List(const int index)
    {PHYSBAM_FATAL_ERROR();}
public:
    void Remove_Body(const int p)
    {BASE::Remove_Geometry(p);}

//#####################################################################
};
}
#endif
