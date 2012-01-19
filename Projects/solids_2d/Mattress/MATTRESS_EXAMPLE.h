//#####################################################################
// Copyright 2002-2004, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATTRESS_EXAMPLE
//#####################################################################
#ifndef __MATTRESS_EXAMPLE__
#define __MATTRESS_EXAMPLE__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_SPLINE_MODEL_2D.h>
#include <Forces_And_Torques/DIAGONALIZED_FINITE_VOLUME_2D.h>
namespace PhysBAM{

template <class T,class RW>
class MATTRESS_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
    typedef VECTOR<T,2> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::data_directory;using BASE::fluids_parameters;using BASE::solids_parameters;
    using BASE::verbose_dt;

    T initial_height;
    T initial_orientation;
    TV initial_velocity;
    T initial_angular_velocity;
    GRID<TV> mattress_grid;
    bool fix_boundary;
    
    MATTRESS_EXAMPLE()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(0,fluids_parameters.NONE),initial_height(4),initial_orientation(0),
        initial_velocity(0,0),initial_angular_velocity(3),
        mattress_grid(3,5,-.25,.25,-.5,.5),
        fix_boundary(false)
    {
        if(fix_boundary){initial_angular_velocity=0;initial_velocity=TV(0,0);}
        last_frame=10*24;
        solids_parameters.cfl=(T).5;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        output_directory="Mattress/output";
        solids_parameters.collide_with_interior=true;
        verbose_dt=true;
    }

    ~MATTRESS_EXAMPLE()
    {}

    // unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data() PHYSBAM_OVERRIDE
{
    DEFORMABLE_OBJECT_2D<T>& deformable_object=solid_body_collection.deformable_object;
    PARTICLES<TV>& particles=deformable_object.particles;
    TRIANGULATED_AREA<T>& triangulated_area=*TRIANGULATED_AREA<T>::Create(particles);
    deformable_object.Add_Structure(&triangulated_area);

    triangulated_area.Initialize_Square_Mesh_And_Particles(mattress_grid);

    triangulated_area.Set_Density(1000);
    triangulated_area.Set_Mass_Of_Particles(false);

    triangulated_area.Update_Bounding_Box();
    TV center(triangulated_area.bounding_box->Center());T bottom=triangulated_area.bounding_box->ymin;
    for(int i=0;i<triangulated_area.particles.array_size;i++){
        particles.X(i)=center+MATRIX<T,2>::Rotation_Matrix(initial_orientation)*(particles.X(i)-center);
        particles.V(i)=initial_velocity+TV::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i).y+=initial_height-bottom;}

    triangulated_area.Check_Signed_Areas_And_Make_Consistent(true);

    int index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    Get_Initial_Data();
    DEFORMABLE_OBJECT_2D<T>& deformable_object=solid_body_collection.deformable_object;
    TRIANGULATED_AREA<T>& triangulated_area=deformable_object.template Find_Structure<TRIANGULATED_AREA<T>&>();

    solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_object.particles));
    solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(triangulated_area,new DIAGONALIZED_SPLINE_MODEL_2D<T>((T)5e4,(T).3,(T).5,7,(T).01)));
    //deformable_object.Add_Edge_Springs(deformable_object.triangulated_area->triangle_mesh,3000,(T).9);
    //deformable_object.Add_Altitude_Springs(deformable_object.triangulated_area->triangle_mesh,300);
    //deformable_object.Add_Neo_Hookean_Elasticity(*deformable_object.triangulated_area,(T)5e4);
    //deformable_object.Add_Linear_Finite_Volume(*deformable_object.triangulated_area,3e4);

    solid_body_collection.Update_Fragments();

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Initialize_Bodies();
}
//#####################################################################
// Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T time) PHYSBAM_OVERRIDE
{
    switch(id_number){
    case 1:
        if(fix_boundary)for(int j=1,index=0;j<=mattress_grid.n;j++)for(int k=0;k<mattress_grid.m;k++){index++;V(index)=TV(0,0);V(V.m+1-index)=TV(0,0);}
        break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
// Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T time) PHYSBAM_OVERRIDE
{
    switch(id_number){
    case 1:
        if(fix_boundary)for(int j=1,index=0;j<=mattress_grid.n;j++)for(int k=0;k<mattress_grid.m;k++){index++;V(index)=TV();V(V.m+1-index)=TV();}
        break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
};
}
#endif
