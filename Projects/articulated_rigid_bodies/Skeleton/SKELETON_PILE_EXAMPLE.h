//#####################################################################
// Copyright 2007, Geoffrey Irving, Tamar Shinar
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SKELETON_PILE_EXAMPLE
//#####################################################################
#ifndef SKELETON_PILE_EXAMPLE__
#define SKELETON_PILE_EXAMPLE__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Parsing/STRING_UTILITIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include "../VISIBLE_HUMAN.h"
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NR3.h>
namespace PhysBAM{

template<class T>
class SKELETON_PILE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >
{
    typedef UNIFORM_GRID_ITERATOR_CELL<TV> CELL_ITERATOR;typedef UNIFORM_GRID_ITERATOR_FACE_3D<T> FACE_ITERATOR;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;typedef VECTOR<T,3> TV;
    using BASE::fluids_parameters;using BASE::output_directory;using BASE::solids_parameters;using BASE::write_last_frame;using BASE::data_directory;using BASE::stream_type;
    using BASE::restart;using BASE::initial_time;using BASE::first_frame;using BASE::last_frame;using BASE::restart_frame;using BASE::frame_rate;

    ARTICULATED_RIGID_BODY<TV>& arb;
    SOLIDS_STANDARD_TESTS<TV> tests;
    VISIBLE_HUMAN<T> skeleton;
    int num_skeletons;
    BOX_3D<T> skeleton_bounding_box;

    SKELETON_PILE_EXAMPLE(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),arb(*new ARTICULATED_RIGID_BODY<TV>(solid_body_collection.deformable_object.particles,solids_parameters.rigid_body_parameters.list)),
        tests(*this,solids_parameters),skeleton(stream_type,&arb,solid_body_collection.deformable_object.rigid_body_particles,data_directory,FRAME<TV>(),true),num_skeletons(1),
        skeleton_bounding_box(-(T)0.00641811,(T)0.600359,-(T)0.00898567,(T)0.368807,-(T)0.00963866,(T)1.88996)
    {
        solids_parameters.rigid_body_parameters.Set_Articulated_Rigid_Body(&arb);
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.perform_self_collision=false;

        output_directory=STRING_UTILITIES::string_sprintf("Skeleton_Pile/%s",output_directory.c_str());
        frame_rate=48;last_frame=720;

        arb.Set_Use_Shock_Propagation(false);
    }

    ~SKELETON_PILE_EXAMPLE()
    {}

    // Unused callbacks
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(TWIST<TV>& wrench,const T time,const int id) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Initialize
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    tests.Add_Ground((T).5,0);

    Make_Skeletons();

    ARRAY<RIGID_BODY<TV>*>& rigid_bodies=solids_parameters.rigid_body_parameters.list.rigid_bodies;
    LOG::cout<<"rigid_bodies.m="<<rigid_bodies.m<<std::endl;
    for(int i=1;i<=rigid_bodies.m;i++) if(!rigid_bodies(i)->is_static && !rigid_bodies(i)->is_kinematic)
        rigid_bodies(i)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Initialize_Bodies();
}
//#####################################################################
// Function Make_Skeletons
//#####################################################################
void Make_Skeletons()
{
    RANDOM_NR3 random;random.Set_Seed(65539);
    VECTOR<int,3> counts(2,5,2);
    GRID<TV> grid(counts,BOX_3D<T>());
    ARRAY<FRAME<TV> ,VECTOR<int,3> > frames(grid);
    TV max_corner((T).9,(T).7,(T).9);T base_height=3,max_angle=(T)pi;
    TV box_center=skeleton_bounding_box.Center();
    for(;;max_corner+=(T).1){
        LOG::cout<<"trying cell_size = "<<max_corner<<std::endl;
        grid=GRID_3D<T>(counts,BOX_3D<T>(TV(),TV(counts)*max_corner)+TV(0,base_height,0),true);
        for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){VECTOR<int,3> I=iterator.Cell_Index();
            frames(I).r=QUATERNION<T>::From_Rotation_Vector(max_angle*random.Get_Vector_In_Unit_Sphere<TV>());
            frames(I).t=iterator.Location()-box_center;}
        int intersections;
        for(int attempts=1;attempts<=200;attempts++){intersections=0;
            for(FACE_ITERATOR iterator(grid,0,GRID<TV>::INTERIOR_REGION);iterator.Valid();iterator.Next()){
                VECTOR<int,3> I=iterator.First_Cell_Index(),J=iterator.Second_Cell_Index();
                if(ORIENTED_BOX<TV>(skeleton_bounding_box,frames(I).Inverse_Times(frames(J))).Intersection(skeleton_bounding_box)){
                    frames(I).r=QUATERNION<T>::From_Rotation_Vector(max_angle*random.Get_Vector_In_Unit_Sphere<TV>());
                    frames(J).r=QUATERNION<T>::From_Rotation_Vector(max_angle*random.Get_Vector_In_Unit_Sphere<TV>());
                    intersections++;}}
            if(!intersections) break;}
        LOG::cout<<"intersections = "<<intersections<<std::endl;
        if(!intersections){LOG::cout<<"final cell size = "<<max_corner<<std::endl;break;}}

    const ARRAY<bool> bone_filter;
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){VECTOR<int,3> I=iterator.Cell_Index();
        skeleton.transform=frames(I);
        skeleton.Make_Bones(bone_filter);skeleton.Make_Joints();skeleton.Add_Joints_To_ARB();}

// compute the skeleton bounding box
//    BOX_3D<T> bounding_box=BOX_3D<T>::Empty_Box();
//    for(int i=1;i<=skeleton.bones.m;i++){skeleton.bones(i)->Update_Bounding_Box();bounding_box.Enlarge_To_Include_Box(skeleton.bones(i)->axis_aligned_bounding_box);}
//    LOG::cout<<"Axis aligned bounding box="<<bounding_box<<std::endl;
}
//#####################################################################
};
}

#endif
