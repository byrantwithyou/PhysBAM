//#####################################################################
// Copyright 2006, Tamar Shinar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_PARTICLE_EXAMPLE
//#####################################################################
#ifndef __RIGID_PARTICLE_EXAMPLE__
#define __RIGID_PARTICLE_EXAMPLE__

#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
namespace PhysBAM{

template<class T,class RW>
class RIGID_PARTICLE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;typedef VECTOR<T,2> TV;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;

    RIGID_PARTICLE_EXAMPLE()
        :BASE(0,fluids_parameters.NONE)
    {
        last_frame=2400;
        restart=false;restart_frame=0;  
        solids_parameters.cfl=(T)1;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
        output_directory="Rigid_Particle/output";
        solids_parameters.use_constant_mass=true;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    }

    ~RIGID_PARTICLE_EXAMPLE()
    {}

    // Unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
virtual void Get_Initial_Data()
{
    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
    PARTICLES<TV>& particles=deformable_object.particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=deformable_object.rigid_body_particles;
    ARRAY<RIGID_BODY<TV>*>& rigid_bodies=solids_parameters.rigid_body_parameters.list.rigid_bodies;
    SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);

    particles.Store_Velocity();

    // add particles and set up segmented curve
    for(int i=1;i<=20;i++){
        particles.array_collection->Add_Element();
        particles.mass(i)=1;
        particles.X(i)=TV((T)i,0);
        if(i>1) segmented_curve.mesh.elements.Append(i,i-1);}
    deformable_object.Add_Structure(&segmented_curve);

    // rigid bodies
    for(int i=1;i<=3;i++){
        solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<T>(data_directory+"/Rigid_Bodies_2D/square",(T)1);
        rigid_bodies(i)->frame.t=TV((T)5*i+(T).5,0);
        rigid_bodies(i)->frame.r=COMPLEX<T>::Unit_Polar(0);
        rigid_bodies(i)->Set_Coefficient_Of_Friction(0);
        T mass_scale_factor=10/rigid_bodies(i)->mass.mass;
        rigid_bodies(i)->mass*=mass_scale_factor;
        rigid_body_particles.array_collection->Add_Element();
        rigid_body_particles.mass(i)=rigid_bodies(i)->mass;
        rigid_body_particles.id(i)=i;
        solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,5*i,rigid_body_particles,i,TV((T)-.5,0)));
        solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,5*i+1,rigid_body_particles,i,TV((T).5,0)));
        rigid_body_particles.Set_State(i,rigid_bodies(i)->frame);
        particles.array_collection->Add_Element();
        particles.mass(i+20)=rigid_body_particles.mass(i).mass;
        particles.X(i+20)=TV((T)5*i+(T).5,0);
        solid_body_collection.deformable_body_collection.binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,20+i,rigid_body_particles,i,TV()));}
    segmented_curve.Update_Number_Nodes();

    // correct mass
    deformable_object.particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.binding_list);

    // collisions
    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    Get_Initial_Data();

    DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
    PARTICLES<TV>& particles=deformable_object.particles;
    SEGMENTED_CURVE<TV>& segmented_curve=deformable_object.template Find_Structure<SEGMENTED_CURVE<TV>&>();

    PARTICLE_SUBSET<PARTICLES<TV> > gravity_particles(particles);
    for(int i=21;i<=23;i++) gravity_particles.active_indices.Append(i);
    gravity_particles.Update_Subset_Index_From_Element_Index();

    solid_body_collection.Add_Force(new GRAVITY<TV>(gravity_particles));
    solid_body_collection.Add_Force(Create_Edge_Springs(segmented_curve,(T)1e2,(T)1));
    //solid_body_collection.Add_Force(Create_Bending_Elements(segmented_curve,(T)1e1,(T)1));

    solid_body_collection.Update_Fragments();

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Initialize_Bodies();
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T time) PHYSBAM_OVERRIDE
{
    V(1)=V(20)=TV();
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T time) PHYSBAM_OVERRIDE
{
    V(1)=V(20)=TV();
}
//#####################################################################
};
}
#endif
