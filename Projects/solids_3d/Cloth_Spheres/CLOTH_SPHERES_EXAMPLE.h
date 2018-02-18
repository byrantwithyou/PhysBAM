//#####################################################################
// Copyright 2006-2007, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CLOTH_SPHERES_EXAMPLE
//#####################################################################
#ifndef __CLOTH_SPHERES_EXAMPLE__
#define __CLOTH_SPHERES_EXAMPLE__

#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/BINDING_SPRINGS.h>
#include <Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <Deformables/Fracture/EMBEDDING.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
namespace PhysBAM{

template<class T_input>
class CLOTH_SPHERES_EXAMPLE:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
public:
    typedef T_input T;typedef VECTOR<T,3> TV;
    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::frame_rate;using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::solids_parameters;
    using BASE::data_directory;using BASE::solid_body_collection;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual
    using BASE::user_last_frame;
    
    SOLIDS_STANDARD_TESTS<TV> tests;
    T sphere_scale,sphere_x_position,aspect_ratio,side_length;
    int number_side_panels;
    int number_of_master_particles;
    ARRAY<T> binding_stiffness;
    ARRAY<COLLISION_GEOMETRY_ID> collision_body_ids;

    CLOTH_SPHERES_EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args),tests(stream_type_input,data_directory,solid_body_collection),sphere_scale((T).5),sphere_x_position((T).7),aspect_ratio((T)1.0),side_length(2),
        number_side_panels(150)
    {
        parse_args.Parse();
        tests.data_directory=data_directory;
    }

    virtual ~CLOTH_SPHERES_EXAMPLE()
    {}

    // Unused callbacks
    void Update_Solids_Parameters(const T time) override {}
    void Preprocess_Frame(const int frame) override {}
    void Postprocess_Frame(const int frame) override {}
    void Apply_Constraints(const T dt,const T time) override {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) override {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) override {}
    void Update_Time_Varying_Material_Properties(const T time) override {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) override {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) override {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() override {}
    void Preprocess_Solids_Substep(const T time,const int substep) override {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
virtual void Get_Initial_Data()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,
        RIGID_BODY_STATE<TV>(FRAME<TV>(TV(),ROTATION<TV>::From_Euler_Angles(0,0,0))));
    TRIANGULATED_SURFACE<T>& master_triangulated_surface=deformable_body_collection.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
    number_of_master_particles=particles.Size();

    for(int i=2;i<=3;i++){
        EMBEDDING<TV>& embedding=*new EMBEDDING<TV>(particles);
        int particle_offset=particles.Size();
        particles.Add_Elements(number_of_master_particles);
        for(int p=0;p<number_of_master_particles;p++) particles.Copy_Element(particles,p,p+particle_offset);
        embedding.material_surface_mesh.Initialize_Mesh_With_Particle_Offset(master_triangulated_surface.mesh,particle_offset);
        if(i==3) for(int t=0;t<embedding.material_surface_mesh.elements.m;t++){ // invert 2nd side of drifted surface
            VECTOR<int,3>& element=embedding.material_surface_mesh.elements(t);element=VECTOR<int,3>(element.y,element.x,element.z);}
        deformable_body_collection.Add_Structure(&embedding);
    }

    for(int i=0;i<2;i++){
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("sphere",sphere_scale,(T).1);
        rigid_body.Frame().t=TV::Axis_Vector(2)*(T)(2*i-3)*sphere_x_position;
        rigid_body.Is_Kinematic()=true;}

    LOG::cout<<"deformable_body_collection.structures.m="<<deformable_body_collection.structures.m<<std::endl;

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();

    // add structures and rigid bodies to collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.structures);
    //solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(&master_triangulated_surface);

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() override
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    if(!user_last_frame) last_frame=240;
    restart=false;restart_frame=0;  
    solids_parameters.cfl=(T)5;
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
    if(!this->user_output_directory)
        output_directory="Cloth_Spheres/output_hires_new";
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
    solids_parameters.triangle_collision_parameters.perform_self_collision=true;

    Get_Initial_Data();

    TRIANGULATED_SURFACE<T>& master_triangulated_surface=deformable_body_collection.template Find_Structure<TRIANGULATED_SURFACE<T>&>(1);
    solid_body_collection.Add_Force(new GRAVITY<TV>(particles,rigid_body_collection,master_triangulated_surface.mesh,0));
    solid_body_collection.Add_Force(Create_Edge_Springs(master_triangulated_surface,(T)1e2,(T)2));
    solid_body_collection.Add_Force(Create_Altitude_Springs(master_triangulated_surface,(T)1e2,(T)2));
    solid_body_collection.Add_Force(Create_Bending_Elements(master_triangulated_surface,(T)3e-3));

    for(int i=0;i<number_of_master_particles;i++){
        soft_bindings.Add_Binding(VECTOR<int,2>(number_of_master_particles+i,i),false);
        soft_bindings.Add_Binding(VECTOR<int,2>(2*number_of_master_particles+i,i),false);}
    soft_bindings.Initialize_Binding_Mesh();
    BINDING_SPRINGS<TV>& binding_springs=*Create_Edge_Binding_Springs(particles,*soft_bindings.binding_mesh,(T)1e3,(T)1);
    solid_body_collection.Add_Force(&binding_springs);
    binding_stiffness.Resize(soft_bindings.bindings.m);binding_stiffness.Fill((T)5e2);
    binding_springs.Set_Stiffness(binding_stiffness);
    binding_springs.Set_Overdamping_Fraction((T)1);

    solid_body_collection.Update_Simulated_Particles();
    SOLIDS_EXAMPLE<TV>::Initialize_Bodies();
    deformable_body_collection.collisions.Use_Structure_Collide_Collision_Body();
    deformable_body_collection.collisions.structure_collide_collision_body(2).Insert(COLLISION_GEOMETRY_ID(2));
    deformable_body_collection.collisions.structure_collide_collision_body(3).Insert(COLLISION_GEOMETRY_ID(1));
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) override
{
    T effective_time=(time<(T)3.5)?time:(T)7-time;
    T velocity=(time<(T)3.5)?(T).15:(T)-.15;
    T distance=abs(velocity)*effective_time;
    frame.t=TV::Axis_Vector(2)*(T)(sphere_x_position-distance)*(T)(2*Value(id)-3);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) override
{
    T velocity=(time<(T)3.5)?(T).15:(T)-.15;
    twist.linear=TV::Axis_Vector(2)*velocity*(T)(3-2*Value(id));
    return true;
}
//#####################################################################
// Function Postprocess_Solids_Substep
//#####################################################################
void Postprocess_Solids_Substep(const T time,const int substep) override
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    T critical_distance=(T).01;

    BINDING_SPRINGS<TV>& binding_springs=solid_body_collection.template Find_Force<BINDING_SPRINGS<TV>&>();
    binding_stiffness.Resize(soft_bindings.bindings.m);
    for(int b=0;b<soft_bindings.bindings.m;b++){
        TV& X0=particles.X(soft_bindings.bindings(b)(1));
        TV& X1=particles.X(soft_bindings.bindings(b)(2));
        T distance=(X0-X1).Magnitude();
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
