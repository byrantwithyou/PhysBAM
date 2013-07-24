//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Sergey Levine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_BODY_BINDING
//##################################################################### 
#ifndef __RIGID_BODY_BINDING__
#define __RIGID_BODY_BINDING__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/SPARSE_UNION_FIND.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Bindings/BINDING.h>
namespace PhysBAM{

template<class TV>
class RIGID_BODY_BINDING:public BINDING<TV>
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
public:
    typedef BINDING<TV> BASE;
    using BASE::particles;using BASE::particle_index;

    RIGID_BODY_COLLECTION<TV>* rigid_body_collection;
    int rigid_body_particles_index;
    TV object_space_position;

    RIGID_BODY_BINDING(DEFORMABLE_PARTICLES<TV>& particles_input);
    RIGID_BODY_BINDING(DEFORMABLE_PARTICLES<TV>& particles_input,const int particle_index_input,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,
        const int rigid_body_particles_index_input,const TV& object_space_position_input);
    virtual ~RIGID_BODY_BINDING();

    static RIGID_BODY_BINDING* Create(GEOMETRY_PARTICLES<TV>& particles);

    virtual int Name() const PHYSBAM_OVERRIDE {return Static_Name();}
    static int Static_Name()
    {return 4;}

    // TODO: This may not handle kinematic/static bodies
    void Add_Dependencies(SEGMENT_MESH& dependency_mesh) const PHYSBAM_OVERRIDE;

    RIGID_BODY<TV>& Rigid_Body() const
    {return rigid_body_collection->Rigid_Body(rigid_body_particles_index);}

    int Id_Number() const
    {return rigid_body_particles_index;}

    TV Embedded_Position() const PHYSBAM_OVERRIDE;
    TV Embedded_Position(ARRAY_VIEW<const TV> X) const PHYSBAM_OVERRIDE;
    TV Embedded_Velocity() const PHYSBAM_OVERRIDE;
    TV Embedded_Velocity(ARRAY_VIEW<const TV> V) const PHYSBAM_OVERRIDE;
    TV Embedded_Velocity(ARRAY_VIEW<const TV> V,ARRAY_VIEW<const TWIST<TV> > twist) const PHYSBAM_OVERRIDE;
    TV Embedded_Acceleration(ARRAY_VIEW<const TV> F,ARRAY_VIEW<const TWIST<TV> > wrench_input) const PHYSBAM_OVERRIDE; // TODO: consider including centrifugal acceleration
    T One_Over_Effective_Mass(const TV& direction) const  PHYSBAM_OVERRIDE;// assumes direction is normalized
    T One_Over_Effective_Mass() const  PHYSBAM_OVERRIDE; // return a lower bound for effective mass over all directions
    void Apply_Impulse(const TV& impulse) PHYSBAM_OVERRIDE;
    void Apply_Impulse(const TV& impulse,ARRAY_VIEW<TV> V) const PHYSBAM_OVERRIDE;
    void Apply_Impulse(const TV& impulse,ARRAY_VIEW<TV> V,ARRAY_VIEW<TWIST<TV> > rigid_V) const PHYSBAM_OVERRIDE;
    void Apply_Push(const TV& impulse) PHYSBAM_OVERRIDE;
    void Apply_Displacement_To_Parents_Based_On_Embedding(const TV& dX,const ARRAY<bool>* skip_particle) PHYSBAM_OVERRIDE;
    void Apply_Velocity_Change_To_Parents_Based_On_Embedding(const TV& dV,const ARRAY<bool>* skip_particle) PHYSBAM_OVERRIDE;
    void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,const TV& force) const PHYSBAM_OVERRIDE;
    void Distribute_Force_To_Parents(ARRAY_VIEW<TV> F_full,ARRAY_VIEW<TWIST<TV> > wrench_full,const TV& force) const PHYSBAM_OVERRIDE;
    void Distribute_Mass_To_Parents(ARRAY_VIEW<T> mass_full) const PHYSBAM_OVERRIDE;
    SYMMETRIC_MATRIX<T,TV::m> Impulse_Factor() const PHYSBAM_OVERRIDE;
    void Parents(ARRAY<int>& parents) const PHYSBAM_OVERRIDE;
    void Weights(ARRAY<T>& weights) const PHYSBAM_OVERRIDE;
private:
    void Read_Helper(TYPED_ISTREAM& input) PHYSBAM_OVERRIDE;
    void Write_Helper(TYPED_OSTREAM& output) const PHYSBAM_OVERRIDE;
//#####################################################################
}; 
}
#endif
