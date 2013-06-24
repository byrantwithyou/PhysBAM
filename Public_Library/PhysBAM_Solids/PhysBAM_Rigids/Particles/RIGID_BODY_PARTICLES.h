//#####################################################################
// Copyright 2006-2009, Michael Lentine, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_PARTICLES
//#####################################################################
#ifndef __RIGID_BODY_PARTICLES__
#define __RIGID_BODY_PARTICLES__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Clone/CLONEABLE.h>
#include <PhysBAM_Tools/Data_Structures/ELEMENT_ID.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <PhysBAM_Tools/Particles/PARTICLES.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
namespace PhysBAM{

template<class TV> class RIGID_BODY;
template<class TV> class STRUCTURE;
template<class TV> class COLLISION_BODY_COLLECTION;
template<class TV,class ID> class STRUCTURE_LIST;

template<class TV>
class RIGID_BODY_PARTICLES:public CLONEABLE<RIGID_BODY_PARTICLES<TV>,PARTICLES<TV> >
{
    typedef typename TV::SCALAR T;
    typedef CLONEABLE<RIGID_BODY_PARTICLES<TV>,PARTICLES<TV> > BASE;
public:
    using BASE::Add_Array;

    ARRAY<RIGID_BODY<TV>*> rigid_body;
    ARRAY_VIEW<FRAME<TV> > frame;
    ARRAY_VIEW<TWIST<TV> > twist;
    ARRAY_VIEW<VECTOR<int,3> > structure_ids;
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

//#####################################################################
    void Resize(const int new_size);
    void Remove_Body(const int p);
    void Clean_Memory();
    void Delete_All_Particles();
    void Clone_Helper(const RIGID_BODY_PARTICLES& particles);
//#####################################################################
};
}
#endif
