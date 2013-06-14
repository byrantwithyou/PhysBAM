//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// 1 cloth bending example - torus invertible test
//#####################################################################
#ifndef __HAIR_TESTS__
#define __HAIR_TESTS__
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_TET_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <fstream>
namespace PhysBAM{

template<class T_input>
class HAIR_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;
    using BASE::stream_type;using BASE::solid_body_collection;using BASE::parse_args;using BASE::test_number;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes; // silence -Woverloaded-virtual
    
    std::ofstream volume_output_file;
    SOLIDS_STANDARD_TESTS<TV> tests;
    TETRAHEDRALIZED_VOLUME<T>* volume;

    int variant;

    HAIR_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(stream_type,data_directory,solid_body_collection),variant(1)
    {
    }

    virtual ~HAIR_TESTS() 
    {
        volume_output_file.close();
    }

    // Unused callbacks
//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options()
{
    BASE::Register_Options();
    parse_args->Add("-variant",&variant,"value","variant");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options()
{
    BASE::Parse_Options();
    tests.data_directory=data_directory;
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    //SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    frame_rate=24;
    last_frame=240;
    solids_parameters.verbose_dt=true;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    
    switch(test_number){
        case 1:
            solids_parameters.cfl=(T)10;
            last_frame=(int)(10*frame_rate);
            break;
        case 2:
            solids_parameters.cfl=(T)10;
            solids_parameters.deformable_object_collision_parameters.collide_with_interior=true;
            last_frame=(int)(5*frame_rate);
            break;
    }
    output_directory=STRING_UTILITIES::string_sprintf("Hair_Tests/Test_%d",test_number);
    LOG::cout<<"output_directory="<<output_directory<<std::endl;
    volume_output_file.open("/home/mlentine/volumes.txt");
    
    // make geometry
    switch(test_number){
        case 1:{
            TRIANGULATED_SURFACE<T>& surface=tests.Create_Triangulated_Object(data_directory+"/Triangulated_Surfaces/noodle_half_torus.tri",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)1.6,0))),false,true); //,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(),ROTATION<TV>((T).8,TV(1,0,0))));
            volume=TETRAHEDRALIZED_VOLUME<T>::Create(particles);
            volume->mesh.Initialize_Bending_Tetrahedrons(surface.mesh);
            for(int t=0;t<volume->mesh.elements.m;t++){
                VECTOR<int,4>& nodes=volume->mesh.elements(t);
                if(TETRAHEDRON<T>(particles.X.Subset(nodes)).Signed_Volume()<0) exchange(nodes[2],nodes[3]);}
            //deformable_body_collection.Add_Structure(volume);
            tests.Add_Ground();
        }break;
        case 2: {
            volume=&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)1.6,0))),false,true,1000);
            for(int t=0;t<volume->mesh.elements.m;t++){
                VECTOR<int,4>& nodes=volume->mesh.elements(t);
                if(TETRAHEDRON<T>(particles.X.Subset(nodes)).Signed_Volume()<0) exchange(nodes[2],nodes[3]);}
            deformable_body_collection.Add_Structure(volume);
            tests.Add_Ground();
        }break;
    }

    // setup collisions
    deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    bool strain_limit=true;
    // make forces
    switch(test_number){
        case 1:{ // initialize very weak forces
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            for(int i=0;i<deformable_body_collection.structures.m;i++){
                if(TRIANGULATED_SURFACE<T>* surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(deformable_body_collection.structures(i))){
                    solid_body_collection.Add_Force(Create_Edge_Springs(*surface,1/(1+sqrt((T)2)),(T)3,strain_limit));
                    solid_body_collection.Add_Force(Create_Bending_Springs(*surface,1/(1+sqrt((T)2)),(T)3,strain_limit));}}
            if(variant==2) solid_body_collection.Add_Force(Create_Altitude_Springs(*volume,(T)10/(1+sqrt((T)2)),(T)3,true,(T).1,strain_limit));
            if(variant==3) solid_body_collection.Add_Force(Create_Tet_Springs(*volume,10/(1+sqrt((T)2)),(T)3,false,(T).1,strain_limit));
        }break;
        case 2:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            /*for(int i=0;i<deformable_body_collection.structures.m;i++){
                if(TETRAHEDRALIZED_VOLUME<T>* volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.structures(i))){
                    solid_body_collection.Add_Force(Create_Edge_Springs(*volume,1/(1+sqrt((T)2)),(T)3));
                    solid_body_collection.Add_Force(Create_Bending_Springs(volume->Get_Boundary_Object(),1/(1+sqrt((T)2)),(T)3));}}
            */
            solid_body_collection.Add_Force(Create_Edge_Springs(*volume,0/(1+sqrt((T)2)),(T)3));
            //solid_body_collection.Add_Force(Create_Bending_Springs(volume->Get_Boundary_Object(),0/(1+sqrt((T)2)),(T)3));
            if(variant==2) solid_body_collection.Add_Force(Create_Altitude_Springs(*volume,0/(1+sqrt((T)2)),(T)3,false,(T).1,false));
            if(variant==3) solid_body_collection.Add_Force(Create_Tet_Springs(*volume,0/(1+sqrt((T)2)),(T)3,false,(T).1,false));
        }break;
    }
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    switch(test_number) {
        case 1:{
            if(time>(T)2){ // up stiffness!
                LINEAR_SPRINGS<TV>* linear_springs=solid_body_collection.template Find_Force<LINEAR_SPRINGS<TV>*>();
                linear_springs->Set_Stiffness((T)10000);
                linear_springs->Set_Overdamping_Fraction((T)3);
                if(LINEAR_ALTITUDE_SPRINGS_3D<T>* alt_springs=solid_body_collection.template Find_Force<LINEAR_ALTITUDE_SPRINGS_3D<T>*>()){
                    alt_springs->Set_Stiffness((T)10000);alt_springs->Set_Overdamping_Fraction((T)3);}
                if(LINEAR_TET_SPRINGS<T>* tet_springs=solid_body_collection.template Find_Force<LINEAR_TET_SPRINGS<T>*>()){
                    tet_springs->Set_Stiffness((T)10000);tet_springs->Set_Overdamping_Fraction((T)3);}
                TRIANGLE_BENDING_SPRINGS<T>* bending_springs=solid_body_collection.template Find_Force<TRIANGLE_BENDING_SPRINGS<T>*>();
                bending_springs->Set_Stiffness((T)10000);
                bending_springs->Set_Overdamping_Fraction((T)3);
            }
        }break;
        case 2:{
            if(time>(T)2){ // up stiffness!
                LINEAR_SPRINGS<TV>* linear_springs=solid_body_collection.template Find_Force<LINEAR_SPRINGS<TV>*>();
                linear_springs->Set_Stiffness((T)100);
                linear_springs->Set_Overdamping_Fraction((T)3);
                if(LINEAR_ALTITUDE_SPRINGS_3D<T>* alt_springs=solid_body_collection.template Find_Force<LINEAR_ALTITUDE_SPRINGS_3D<T>*>()){
                    alt_springs->Set_Stiffness((T)100);alt_springs->Set_Overdamping_Fraction((T)3);}
                if(LINEAR_TET_SPRINGS<T>* tet_springs=solid_body_collection.template Find_Force<LINEAR_TET_SPRINGS<T>*>()){
                    tet_springs->Set_Stiffness((T)100);tet_springs->Set_Overdamping_Fraction((T)3);}
                //TRIANGLE_BENDING_SPRINGS<T>* bending_springs=solid_body_collection.template Find_Force<TRIANGLE_BENDING_SPRINGS<T>*>();
                //bending_springs->Set_Stiffness((T)100);
                //bending_springs->Set_Overdamping_Fraction((T)3);
            }
            for(int i=0;i<deformable_body_collection.structures.m;i++){
                if(TETRAHEDRALIZED_VOLUME<T>* volume=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.structures(i))){
                    volume_output_file<<time<<"\t"<<volume->Total_Volume()<<std::endl;}}
        }break;
    }
}
//#####################################################################
};
}
#endif
