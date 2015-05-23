//#####################################################################
// Copyright 2007, Nipun Kwatra, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_FLUID_FORCES  
//#####################################################################
#ifndef __EULER_FLUID_FORCES__
#define __EULER_FLUID_FORCES__

#include <Solids/Forces_And_Torques/SOLIDS_FORCES.h>
#include <cfloat>
namespace PhysBAM{

template<class TV> class FLUID_COLLISION_BODY_INACCURATE_UNION;
template<class TV> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;

template<class TV>
class EULER_FLUID_FORCES:public SOLIDS_FORCES<TV>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef SOLIDS_FORCES<TV> BASE;
    typedef typename BASE::RIGID_FREQUENCY_DATA RIGID_FREQUENCY_DATA;
    typedef typename BASE::DEFORMABLE_FREQUENCY_DATA DEFORMABLE_FREQUENCY_DATA;

    using BASE::particles;

    GRID<TV> grid;
    const ARRAY<T,FACE_INDEX<TV::m> >& pressure_at_faces;
    const ARRAY<bool,FACE_INDEX<TV::m> >& solid_fluid_face;
    const ARRAY<bool,TV_INT>& cells_inside_fluid;
    const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>* collision_bodies_affecting_fluid;

public:
    EULER_FLUID_FORCES(const GRID<TV>& grid_input,const ARRAY<T,FACE_INDEX<TV::m> >& pressure_at_faces_input,
        const ARRAY<bool,FACE_INDEX<TV::m> >& solid_fluid_face_input,const ARRAY<bool,TV_INT>& cells_inside_fluid_input,
        const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>* collision_bodies_affecting_fluid_input,DEFORMABLE_PARTICLES<TV>& particles_input,
        RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input);

    virtual ~EULER_FLUID_FORCES();

    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE
    {}

    T CFL_Strain_Rate() const PHYSBAM_OVERRIDE
    {return FLT_MAX;} // TODO

    void Initialize_CFL(ARRAY_VIEW<DEFORMABLE_FREQUENCY_DATA> frequency,ARRAY_VIEW<RIGID_FREQUENCY_DATA> rigid_frequency) PHYSBAM_OVERRIDE
    {}

    void Add_Velocity_Dependent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE
    {}

    void Add_Implicit_Velocity_Independent_Forces(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > rigid_V,ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE
    {}

    void Enforce_Definiteness(const bool enforce_definiteness_input) PHYSBAM_OVERRIDE
    {}

    void Update_Mpi(const ARRAY<bool>& particle_is_simulated,const ARRAY<bool>& rigid_particle_is_simulated,MPI_SOLIDS<TV>* mpi_solids) PHYSBAM_OVERRIDE
    {}

//#####################################################################
    void Add_Velocity_Independent_Forces(ARRAY_VIEW<TV> F,ARRAY_VIEW<TWIST<TV> > rigid_F,const T time) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
