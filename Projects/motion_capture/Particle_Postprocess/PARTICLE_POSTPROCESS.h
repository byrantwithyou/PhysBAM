//#####################################################################
// Copyright 2009, Michael Lentine
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARTICLE_POSTPROCESS
//#####################################################################
//   1. Add passively advected particles to an existing fluid simulation
//#####################################################################
#ifndef __PARTICLE_POSTPROCESS__
#define __PARTICLE_POSTPROCESS__

#include <PhysBAM_Tools/Arrays/IDENTITY_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
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
class PARTICLE_POSTPROCESS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;typedef VECTOR<T,3> TV;typedef GRID<TV> T_GRID;typedef VECTOR<int,3> TV_INT;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<T_GRID> BASE;
    using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;
    using BASE::data_directory;using BASE::stream_type;using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;using BASE::resolution;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes;using BASE::Set_External_Positions; // silence -Woverloaded-virtual
    using BASE::Initialize_Solid_Fluid_Coupling_Before_Grid_Initialization;using BASE::Add_Volumetric_Body_To_Fluid_Simulation;

    COLLISION_GEOMETRY_COLLECTION<TV> local_collision_body_list;
    DEFORMABLE_BODY_COLLECTION<TV> local_deformable_object;
    TRIANGULATED_SURFACE<T> *surface;
    SMOKE_STANDARD_TESTS_3D<T_GRID> smoke_tests;
    SOLIDS_STANDARD_TESTS<TV> solids_tests;
    int num_particles_x,num_particles_y,num_particles_z;
    int frame_1,frame_n;
    std::string input_directory,embedded_directory;

    PARTICLE_POSTPROCESS(const STREAM_TYPE stream_type)
        :BASE(stream_type,1,fluids_parameters.WATER),local_deformable_object(solid_body_collection.example_forces_and_velocities,local_collision_body_list),
        surface(0),smoke_tests(*this,fluids_parameters,fluid_collection.incompressible_fluid_collection,solid_body_collection.rigid_body_collection),
        solids_tests(*this,solid_body_collection)
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

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_String_Argument("-input_directory","Specify the input directory");
    parse_args->Add_String_Argument("-embedded_directory","Specify the input directory");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    smoke_tests.Initialize(Smoke_Test_Number(test_number),resolution);
    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
    frame_1=1;frame_n=1000;
    last_frame=frame_n-frame_1+1;
    num_particles_x=10*resolution;num_particles_y=30*resolution;num_particles_z=10*resolution;

    if(parse_args->Is_Value_Set("-input_directory")) input_directory=parse_args->Get_String_Value("-input_directory");
    if(parse_args->Is_Value_Set("-embedded_directory")) embedded_directory=parse_args->Get_String_Value("-embedded_directory");

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;

    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction=(T).1;
    solids_parameters.rigid_body_evolution_parameters.rigid_geometry_evolution_parameters.use_kinematic_keyframes=false;

    *fluids_parameters.grid=smoke_tests.grid;
    fluids_parameters.solid_affects_fluid=false;
    fluids_parameters.solve_neumann_regions=false;
    fluids_parameters.fluid_affects_solid=false;
    fluids_parameters.density=(T)100000;
    fluids_parameters.incompressible_iterations=200;
    //fluids_parameters.number_particles_per_cell=0;
    fluids_parameters.write_particles=true;
    fluids_parameters.write_removed_positive_particles=true;
    fluids_parameters.use_removed_positive_particles=true;
    fluids_parameters.reincorporate_removed_particle_velocity=false;
    fluids_parameters.removed_positive_particle_buoyancy_constant=0;
    fluids_parameters.domain_walls[3][2]=fluids_parameters.domain_walls[3][1]=false;
    fluids_parameters.domain_walls[2][2]=fluids_parameters.domain_walls[2][1]=fluids_parameters.domain_walls[1][1]=fluids_parameters.domain_walls[1][2]=true;

    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/grid",*fluids_parameters.grid);

    output_directory=STRING_UTILITIES::string_sprintf("Particle_Postprocess/Test_Prune_%d_%d_%d_%d",test_number,num_particles_x,num_particles_y,num_particles_z);
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE 
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame-frame_1+1);
    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/mac_velocities."+f,fluid_collection.incompressible_fluid_collection.face_velocities);
    local_deformable_object.Read(stream_type,embedded_directory+"/",frame,-1,true,solids_parameters.write_from_every_process);
    surface=local_deformable_object.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>*>();
    surface->Update_Triangle_List();surface->mesh.Initialize_Adjacent_Elements();surface->Initialize_Hierarchy();surface->Update_Bounding_Box();
    if(frame>1) return;

    T xwidth=fluids_parameters.grid->domain.max_corner.x-fluids_parameters.grid->domain.min_corner.x,ywidth=fluids_parameters.grid->domain.max_corner.y-fluids_parameters.grid->domain.min_corner.y,zwidth=fluids_parameters.grid->domain.max_corner.z-fluids_parameters.grid->domain.min_corner.z;
    for(int i=0;i<num_particles_x;i++) for(int j=0;j<num_particles_z;j++) for(int k=1;k<=(frame==1?num_particles_y:1);k++){
        TV position=TV((float)i/(float)(num_particles_x+1)*xwidth+fluids_parameters.grid->domain.min_corner.x,(float)k/(float)(num_particles_y+1)*ywidth+fluids_parameters.grid->domain.min_corner.y,(float)j/(float)(num_particles_z+1)*zwidth+fluids_parameters.grid->domain.min_corner.z);
        TV_INT index=fluids_parameters.grid->Clamp_To_Cell(position);
        if(!fluids_parameters.particle_levelset_evolution->particle_levelset.removed_positive_particles(index))
            fluids_parameters.particle_levelset_evolution->particle_levelset.removed_positive_particles(index)=new PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>();
        int p=fluids_parameters.particle_levelset_evolution->particle_levelset.removed_positive_particles(index)->array_collection->Add_Element();
        fluids_parameters.particle_levelset_evolution->particle_levelset.removed_positive_particles(index)->X(p)=position;}
} 
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE 
{
    for(CELL_ITERATOR iterator(*fluids_parameters.grid);iterator.Valid();iterator.Next()){TV_INT index=iterator.Cell_Index();
        if(!fluids_parameters.particle_levelset_evolution->particle_levelset.removed_positive_particles(index)) continue;
        for(int i=1;i<=fluids_parameters.particle_levelset_evolution->particle_levelset.removed_positive_particles(index)->array_collection->number;i++)
            if(surface->Inside(fluids_parameters.particle_levelset_evolution->particle_levelset.removed_positive_particles(index)->X(i)))
                fluids_parameters.particle_levelset_evolution->particle_levelset.removed_positive_particles(index)->array_collection->Delete_Element(i);}
}
//#####################################################################
// Function Smoke_Test_Number
//#####################################################################
static int Smoke_Test_Number(const int test_number)
{
    switch(test_number){
        case 1:
            return 1;
        default:
            return 1;}
}
//#####################################################################
// Function Initialize_Advection
//#####################################################################
void Initialize_Advection() PHYSBAM_OVERRIDE
{
    fluids_parameters.Use_No_Fluid_Coupling_Defaults();
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    T_GRID& grid=*fluids_parameters.grid;
    ARRAY<T,VECTOR<int,3> >& phi=fluids_parameters.particle_levelset_evolution->phi;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()) phi(iterator.Cell_Index())=-1;
}
//#####################################################################
// Function Initialize_Velocities
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame_1);
    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/mac_velocities."+f,fluid_collection.incompressible_fluid_collection.face_velocities);
    BASE::Initialize_Velocities();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
    solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
};
}
#endif
