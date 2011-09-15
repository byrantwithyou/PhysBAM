//#####################################################################
// Copyright 2006-2007, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CLOTH_SPHERES_EXAMPLE
//#####################################################################
#ifndef __CLOTH_SPHERES_EXAMPLE__
#define __CLOTH_SPHERES_EXAMPLE__

#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BINDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDING.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
namespace PhysBAM{

template<class T_input>
class CLOTH_SPHERES_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
public:
    typedef T_input T;typedef VECTOR<T,3> TV;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::frame_rate;using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::solids_parameters;using BASE::fluids_parameters;
    using BASE::data_directory;using BASE::solid_body_collection;using BASE::parse_args;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual

    SOLIDS_STANDARD_TESTS<TV> tests;
    T sphere_scale,sphere_x_position,aspect_ratio,side_length;
    int number_side_panels;
    int number_of_master_particles;
    ARRAY<T> binding_stiffness;
    ARRAY<COLLISION_GEOMETRY_ID> collision_body_ids;

    CLOTH_SPHERES_EXAMPLE(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),sphere_scale((T).5),sphere_x_position((T).7),aspect_ratio((T)1.0),side_length(2),
        number_side_panels(150)
    {
    }

    ~CLOTH_SPHERES_EXAMPLE()
    {}

    // Unused callbacks
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options()
{
    BASE::Register_Options();
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options()
{
    BASE::Parse_Options();
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
virtual void Get_Initial_Data()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,
        RIGID_BODY_STATE<TV>(FRAME<TV>(TV(),ROTATION<TV>::From_Euler_Angles(0,0,0))));
    TRIANGULATED_SURFACE<T>& master_triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
    number_of_master_particles=particles.array_collection->Size();

    for(int i=2;i<=3;i++){
        EMBEDDING<TV>& embedding=*new EMBEDDING<TV>(particles);
        int particle_offset=particles.array_collection->Size();
        particles.array_collection->Add_Elements(number_of_master_particles);
        for(int p=1;p<=number_of_master_particles;p++) particles.array_collection->Copy_Element(*particles.array_collection,p,p+particle_offset);
        embedding.material_surface_mesh.Initialize_Mesh_With_Particle_Offset(master_triangulated_surface.mesh,particle_offset);
        if(i==3) for(int t=1;t<=embedding.material_surface_mesh.elements.m;t++){ // invert 2nd side of drifted surface
            VECTOR<int,3>& element=embedding.material_surface_mesh.elements(t);element=VECTOR<int,3>(element.y,element.x,element.z);}
        deformable_body_collection.deformable_geometry.Add_Structure(&embedding);
    }

    for(int i=1;i<=2;i++){
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("sphere",sphere_scale,(T).1);
        rigid_body.X()=TV::Axis_Vector(2)*(T)(2*i-3)*sphere_x_position;
        rigid_body.Is_Kinematic()=true;}

    LOG::cout<<"deformable_body_collection.deformable_geometry.structures.m="<<deformable_body_collection.deformable_geometry.structures.m<<std::endl;

    // correct number nodes
    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    //solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(&master_triangulated_surface);

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    last_frame=240;
    restart=false;restart_frame=0;  
    solids_parameters.cfl=(T)5;
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
    output_directory="Cloth_Spheres/output_hires_new";
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
    solids_parameters.triangle_collision_parameters.perform_self_collision=true;

    Get_Initial_Data();

    TRIANGULATED_SURFACE<T>& master_triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>(1);
    solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,master_triangulated_surface.mesh,0));
    solid_body_collection.Add_Force(Create_Edge_Springs(master_triangulated_surface,(T)1e2,(T)2));
    solid_body_collection.Add_Force(Create_Altitude_Springs(master_triangulated_surface,(T)1e2,(T)2));
    solid_body_collection.Add_Force(Create_Bending_Elements(master_triangulated_surface,(T)3e-3));

    for(int i=1;i<=number_of_master_particles;i++){
        soft_bindings.Add_Binding(VECTOR<int,2>(number_of_master_particles+i,i),false);
        soft_bindings.Add_Binding(VECTOR<int,2>(2*number_of_master_particles+i,i),false);}
    soft_bindings.Initialize_Binding_Mesh();
    BINDING_SPRINGS<TV>& binding_springs=*Create_Edge_Binding_Springs(particles,*soft_bindings.binding_mesh,(T)1e3,(T)1);
    solid_body_collection.Add_Force(&binding_springs);
    binding_stiffness.Resize(soft_bindings.bindings.m);ARRAYS_COMPUTATIONS::Fill(binding_stiffness,(T)5e2);
    binding_springs.Set_Stiffness(binding_stiffness);
    binding_springs.Set_Overdamping_Fraction((T)1);

    solid_body_collection.Update_Simulated_Particles();
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Initialize_Bodies();
    deformable_body_collection.collisions.Use_Structure_Collide_Collision_Body();
    deformable_body_collection.collisions.structure_collide_collision_body(2).Insert(COLLISION_GEOMETRY_ID(2));
    deformable_body_collection.collisions.structure_collide_collision_body(3).Insert(COLLISION_GEOMETRY_ID(1));
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    T effective_time=(time<(T)3.5)?time:(T)7-time;
    T velocity=(time<(T)3.5)?(T).15:(T)-.15;
    T distance=abs(velocity)*effective_time;
    frame.t=TV::Axis_Vector(2)*(T)(sphere_x_position-distance)*(T)(2*Value(id)-3);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    T velocity=(time<(T)3.5)?(T).15:(T)-.15;
    twist.linear=TV::Axis_Vector(2)*velocity*(T)(3-2*Value(id));
    return true;
}
//#####################################################################
// Function Postprocess_Solids_Substep
//#####################################################################
void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    T critical_distance=(T).01;

    BINDING_SPRINGS<TV>& binding_springs=solid_body_collection.template Find_Force<BINDING_SPRINGS<TV>&>();
    binding_stiffness.Resize(soft_bindings.bindings.m);
    for(int b=1;b<=soft_bindings.bindings.m;b++){
        TV& X1=particles.X(soft_bindings.bindings(b)(1));
        TV& X2=particles.X(soft_bindings.bindings(b)(2));
        T distance=(X1-X2).Magnitude();
        T ratio=distance/critical_distance;
        if(ratio>1)
            binding_stiffness(b)=(T)5e2/ratio;
        else
            binding_stiffness(b)=(T)5e2;}
    binding_springs.Set_Stiffness(binding_stiffness);
    binding_springs.Set_Overdamping_Fraction((T)1);
}

//#####################################################################
};
}
#endif
