//#####################################################################
// Copyright 2009, Michael Lentine
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPHERE_POSTPROCESS
//#####################################################################
//   1. Add passively advected spheres to an existing fluid simulation
//   2. Add passively advected spheres to an existing fluid simulation
//#####################################################################
#ifndef __SPHERE_POSTPROCESS__
#define __SPHERE_POSTPROCESS__

#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/SMOKE_STANDARD_TESTS_3D.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/THIN_SHELLS_FLUID_COUPLING_UTILITIES.h>
#include <PhysBAM_Dynamics/Coupled_Evolution/SOLID_FLUID_COUPLED_EVOLUTION.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class SPHERE_POSTPROCESS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;typedef VECTOR<T,3> TV;typedef GRID<TV> T_GRID;typedef VECTOR<int,3> TV_INT;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID> BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;
    using BASE::data_directory;using BASE::stream_type;using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;using BASE::resolution;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual
    using BASE::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization;using BASE::Add_Volumetric_Body_To_Fluid_Simulation;

    SOLIDS_STANDARD_TESTS<TV> tests;
    int num_particles_x,num_particles_y,num_particles_z;
    int frame_1,frame_n;
    std::string input_directory;
    T_GRID grid;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ARRAY<TV,VECTOR<int,3> > v_array;
    T frame_time;
    ARRAY<FRAME<TV> > frame_start;
    ARRAY<PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>*,VECTOR<int,3> > particle_array;
    ARRAY<TV> particle_flat;

    SPHERE_POSTPROCESS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),particle_flat(0)
    {
    }

    // Unused callbacks
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_String_Argument("-input_directory","Specify the input directory");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
    frame_1=1;frame_n=1000;
    last_frame=frame_n-frame_1+1;
    num_particles_x=10*resolution;num_particles_y=30*resolution;num_particles_z=10*resolution;

    if(parse_args->Is_Value_Set("-input_directory")) input_directory=parse_args->Get_String_Value("-input_directory");

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;

    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction=(T).1;
    solids_parameters.rigid_body_evolution_parameters.rigid_geometry_evolution_parameters.use_kinematic_keyframes=false;

    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/grid",grid);

    output_directory=STRING_UTILITIES::string_sprintf("Sphere_Postprocess/Test_%d_%d_%d_%d",test_number,num_particles_x,num_particles_y,num_particles_z);
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==1) frame.t=frame_start(id).t+v_array(grid.Clamp_To_Cell(frame.t))*(time-frame_time);
    else if(test_number==2) frame.t=particle_flat(id);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==1){
        twist.linear=v_array(grid.Clamp_To_Cell(solid_body_collection.rigid_body_collection.rigid_body_particle.X(id)));
        return true;}
    return false;
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE 
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particle;
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame-frame_1+1);
    if(test_number==1){
        frame_time=(float)frame/frame_rate;
        FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/mac_velocities."+f,face_velocities);
        v_array.Resize(face_velocities.Component(1).Domain_Indices(),false,false);
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            v_array(cell_index)=TV((face_velocities(1,cell_index)+face_velocities(1,TV_INT(cell_index.x+1,cell_index.y,cell_index.z)))/(T)2,
                (face_velocities(1,cell_index)+face_velocities(1,TV_INT(cell_index.x,cell_index.y+1,cell_index.z)))/(T)2,
                (face_velocities(1,cell_index)+face_velocities(1,TV_INT(cell_index.x,cell_index.y,cell_index.z+1)))/(T)2);
        }
        for(int id=0;id<rigid_body_particles.array_collection->Size();id++) frame_start(id)=solid_body_collection.rigid_body_collection.Rigid_Body(id).Frame();}
    else if (test_number==2){
        FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/removed_positive_particles."+f,particle_array);
        int index=1;
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
            if(!particle_array(cell_index)) continue;
            for(int i=1;i<=particle_array(cell_index)->array_collection->number;i++) particle_flat(index++)=particle_array(cell_index)->X(i);}
        if(frame>0) rigid_body_particles.Resize(index-1);
        particle_flat.Resize(rigid_body_particles.array_collection->Size());
        solid_body_collection.Update_Simulated_Particles();
    }
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE 
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particle;
    if(test_number==1) for(int id=0;id<rigid_body_particles.array_collection->Size();id++) rigid_body_particles.X(id)=frame_start(id).t+v_array(grid.Clamp_To_Cell(frame_start(id).t))*(1./frame_rate);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    T xwidth=grid.domain.max_corner.x-grid.domain.min_corner.x,ywidth=grid.domain.max_corner.y-grid.domain.min_corner.y,zwidth=grid.domain.max_corner.z-grid.domain.min_corner.z;
    for(int i=0;i<num_particles_x;i++) for(int j=0;j<num_particles_z;j++) for(int k=0;k<num_particles_y;k++){
        TV position=TV((float)i/(float)(num_particles_x+1)*xwidth+grid.domain.min_corner.x,(float)k/(float)(num_particles_y+1)*ywidth+grid.domain.min_corner.y,(float)j/(float)(num_particles_z+1)*zwidth+grid.domain.min_corner.z);
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("sphere",.1,(T).1);
        rigid_body.X()=position;rigid_body.Is_Kinematic()=true;
    }

    frame_start.Resize(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size());
    particle_flat.Resize(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size());
    Preprocess_Frame(0);
    
    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
};
}

#endif
