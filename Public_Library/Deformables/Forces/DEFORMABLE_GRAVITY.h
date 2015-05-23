//#####################################################################
// Copyright 2002-2009, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_GRAVITY
//#####################################################################
#ifndef __DEFORMABLE_GRAVITY__
#define __DEFORMABLE_GRAVITY__

#include <Deformables/Forces/POINTWISE_DEFORMABLE_FORCE.h>
namespace PhysBAM{

template<class TV>
class DEFORMABLE_GRAVITY:public POINTWISE_DEFORMABLE_FORCE<TV>
{
    typedef typename TV::SCALAR T;
public:
    typedef POINTWISE_DEFORMABLE_FORCE<TV> BASE;
    using BASE::particles;using BASE::influenced_particles;using BASE::mpi_solids;using BASE::force_particles;using BASE::is_simulated;
    typedef typename FORCE_ELEMENTS::ITERATOR ELEMENT_ITERATOR;

    TV gravity;
public:

    DEFORMABLE_GRAVITY(DEFORMABLE_PARTICLES<TV>& particles_input,ARRAY<int>* influenced_particles_input,const TV& gravity_input=-(T)9.8*TV::Axis_Vector(1-(TV::dimension==1)))
        :POINTWISE_DEFORMABLE_FORCE<TV>(particles_input,influenced_particles_input),gravity(gravity_input)
    {}

    DEFORMABLE_GRAVITY(DEFORMABLE_PARTICLES<TV>& particles_input,const bool influence_all_particles_input,const TV& gravity_input=-(T)9.8*TV::Axis_Vector(1-(TV::dimension==1)))
        :POINTWISE_DEFORMABLE_FORCE<TV>(particles_input,influence_all_particles_input),gravity(gravity_input)
    {}

    template<class T_MESH>
    DEFORMABLE_GRAVITY(DEFORMABLE_PARTICLES<TV>& particles_input,const T_MESH& mesh,const TV& gravity_input=-(T)9.8*TV::Axis_Vector(1-(TV::dimension==1)))
        :POINTWISE_DEFORMABLE_FORCE<TV>(particles_input,mesh),gravity(gravity_input)
    {
        mesh.elements.Flattened().Get_Unique(*influenced_particles);
    }

    virtual ~DEFORMABLE_GRAVITY()
    {}

    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE
    {}

    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE
    {}

    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE
    {}

    int Velocity_Dependent_Forces_Size() const PHYSBAM_OVERRIDE
    {return 0;}

    void Add_Velocity_Dependent_Forces_First_Half(ARRAY_VIEW<const TV> V,ARRAY_VIEW<T> aggregate,const T time) const PHYSBAM_OVERRIDE
    {}

    void Add_Velocity_Dependent_Forces_Second_Half(ARRAY_VIEW<const T> aggregate,ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE
    {}

//#####################################################################
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,const T time) const PHYSBAM_OVERRIDE;
    T Potential_Energy(int p,const T time) const;
    T Potential_Energy(const T time) const PHYSBAM_OVERRIDE;
    void Update_Position_Based_State(const T time,const bool is_position_update,const bool update_hessian) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
