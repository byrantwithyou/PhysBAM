//#####################################################################
// Copyright 2007, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __BODY_TEST__
#define __BODY_TEST__
#include <Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <Rigids/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Deformables/Forces/BINDING_SPRINGS.h>
#include <Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_3D.h>
#include <Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <Deformables/Forces/LINEAR_SPRINGS.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
namespace PhysBAM{

template<class T_input>
class BODY_TEST:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::solids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;
    using BASE::stream_type;using BASE::solid_body_collection;

    SOLIDS_STANDARD_TESTS<TV> tests;
    RIGID_BODY<TV>* body0;
    RIGID_BODY<TV>* ground;

    BODY_TEST(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args),tests(stream_type_input,data_directory,solid_body_collection)
    {
        parse_args.Parse();
        tests.data_directory=data_directory;
    }
    void Update_Solids_Parameters(const T time) override
    {
        COLLISION_GEOMETRY_ID body1_collision_geometry_id=
            solid_body_collection.rigid_body_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(body0->particle_index);
        COLLISION_GEOMETRY_ID ground_collision_geometry_id=
            solid_body_collection.rigid_body_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(ground->particle_index);
        if(time*24.0 > 10 && body1_collision_geometry_id && time*24.0 < 25){
            LOG::cout<<"removing body "<<std::endl;
            solid_body_collection.collision_body_list.Remove_Body(body1_collision_geometry_id);}
        else if (time*24.0 > 15 && time*24.0 < 20 &&  ground_collision_geometry_id)
            solid_body_collection.collision_body_list.Remove_Body(ground_collision_geometry_id);
        else if (time*24.0 > 20 && !ground_collision_geometry_id)
            solid_body_collection.collision_body_list.Add_Body(new RIGID_COLLISION_GEOMETRY<TV>(*ground),ground->particle_index,true);
        else if (time*24.0 > 25 && !body1_collision_geometry_id)
            solid_body_collection.collision_body_list.Add_Body(new RIGID_COLLISION_GEOMETRY<TV>(*body0),body0->particle_index,true);
        else
            LOG::cout<<"updat esolids param"<<std::endl;
    }


    void Get_Initial_Data()
    {
        // deformable bodies
        DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
        tests.Create_Cloth_Panel(20,1,1,0);
        //deformable_body_collection.Add_Structure(&cloth);
        tests.Add_Ground();
        body0=&tests.Add_Rigid_Body("sphere",(T).25,(T)0);body0->Frame().t.z=.5;
        ground=&solid_body_collection.rigid_body_collection.Rigid_Body(1);
        deformable_body_collection.particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
        deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.structures);
    }
    void After_Initialization() override {BASE::After_Initialization();}

    void Initialize_Bodies() override
    {
        DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;

        output_directory="Body_Test/output";
        last_frame=30;
        solids_parameters.cfl=(T)4;

        Get_Initial_Data();
        
        // make forces
        ARRAY<int> referenced_nodes;
        for(int i=0;i<deformable_body_collection.structures.m;i++){
            if(TRIANGULATED_SURFACE<T>* surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(deformable_body_collection.structures(i))){
            solid_body_collection.Add_Force(Create_Edge_Springs(*surface,(T)100,(T)3));
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));}}
    }

//#####################################################################
};
}
#endif
