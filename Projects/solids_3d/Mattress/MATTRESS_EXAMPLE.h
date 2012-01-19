//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATTRESS_EXAMPLE
//#####################################################################
#ifndef __MATTRESS_EXAMPLE__
#define __MATTRESS_EXAMPLE__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <Constitutive_Models/DIAGONALIZED_NEO_HOOKEAN_3D.h>
#include <Forces_And_Torques/DIAGONALIZED_FINITE_VOLUME_3D.h>
namespace PhysBAM{

template <class T,class RW>
class MATTRESS_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;
    using BASE::output_directory;using BASE::data_directory;using BASE::fluids_parameters;
    using BASE::solids_parameters;

    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    VECTOR_3D<T> attachment_velocity;
    VECTOR_3D<T> external_force_rate;
    int m_input,n_input,mn_input;
    T xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input; 
    GRID<TV> mattress_grid;
    bool cube_mesh;
    bool fix_boundary;
    T fiber_tension;
    ARRAY<VECTOR_3D<T> > fiber_direction;
    
    MATTRESS_EXAMPLE()
        :BASE(0,fluids_parameters.NONE),initial_height(4),initial_orientation(0,VECTOR_3D<T>(1,0,0)),
        initial_velocity(0,0,0),initial_angular_velocity(5,2,0),
        m_input(3),n_input(5),mn_input(10),
        xmin_input(-.25),xmax_input(.25),ymin_input(-.5),ymax_input(.5),zmin_input(-1),zmax_input(1),
        mattress_grid(m_input,n_input,mn_input,xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input),
        cube_mesh(true),fix_boundary(true)
    {
        if(fix_boundary){initial_angular_velocity=VECTOR_3D<T>(0,0,0);initial_velocity=VECTOR_3D<T>(0,0,0);}
        last_frame=10*24;
        restart=false;restart_frame=0;   
        solids_parameters.cfl=(T)10;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-5;
        output_directory="Mattress/output";
        if(fix_boundary) solids_parameters.perform_collision_body_collisions=false;
        attachment_velocity=VECTOR_3D<T>(0,0,0);
        solids_parameters.check_initial_mesh_for_self_intersection=false;
        solids_parameters.perform_self_collision=false;
        solids_parameters.collide_with_interior=false;
        solids_parameters.deformable_body_parameters.print_residuals=true;
        solids_parameters.deformable_body_parameters.print_diagnostics=true;
        solids_parameters.backward_euler=true;
    }

    ~MATTRESS_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    Get_Initial_Data();
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume;
    DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);

    solid_body_collection.Add_Force(Create_Body_Forces<T>(tetrahedralized_volume));
    if(!solids_parameters.backward_euler) solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tetrahedralized_volume,new DIAGONALIZED_NEO_HOOKEAN_3D<T>((T)1e5,(T).45)));
    else{
        solid_body_collection.Add_Force(Create_Backward_Euler_Diagonalized_Finite_Volume(tetrahedralized_volume,
            new DIAGONALIZED_NEO_HOOKEAN_3D<T>((T)1e5,(T).45,(T).25,(T).005),false,(T).1,true,false,true));
        tetrahedralized_volume.particles.mass.Compute_One_Over_Mass();}
}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data() PHYSBAM_OVERRIDE
{
    solids_parameters.deformable_body_parameters.list.Add_Deformable_Object();
    solids_parameters.deformable_body_parameters.list(1).Allocate_Tetrahedralized_Volume();
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume;
    TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.tetrahedron_mesh;

    if(!cube_mesh) tetrahedralized_volume.Initialize_Octahedron_Mesh_And_Particles(mattress_grid);
    else tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(mattress_grid);

    tetrahedralized_volume.Set_Density(1000);
    tetrahedralized_volume.Set_Mass_Of_Particles(true);

    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<tetrahedralized_volume.particles.array_size;i++){
        tetrahedralized_volume.particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}

    fiber_direction.Resize(tetrahedron_mesh.tetrahedrons.m);
    Update_Time_Varying_Material_Properties(0,1);
}
//#####################################################################
// Function Initialize_Collision_Bodies
//#####################################################################
void Initialize_Collision_Bodies()
{  
    if(fix_boundary) return; 
    solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("../../Public_Data/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(solids_parameters.rigid_body_parameters.list.rigid_bodies.m)->coefficient_of_friction=(T).3;
    solids_parameters.collision_body_list.Append_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR_3D<T> >& V,const T time){
    switch(id_number){
    case 1:
        if(fix_boundary)for(int j=1,index=0;j<=n_input;j++)for(int k=0;k<m_input;k++){index++;V(index)=attachment_velocity;V(V.m+1-index)=-attachment_velocity;}
        break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
// Update_Time_Varying_Material_Properties
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {
    for(int t=0;t<fiber_direction.m;t++)fiber_direction(t)=VECTOR_3D<T>(0,(T)1,0);
    fiber_tension=(T)5e4*time;
}
//#####################################################################
// Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> >& V,const T time){
    switch(id_number){
    case 1:
        if(fix_boundary)for(int j=1,index=0;j<=n_input;j++)for(int k=0;k<m_input;k++){index++;V(index)=VECTOR_3D<T>(0,0,0);V(V.m+1-index)=VECTOR_3D<T>(0,0,0);}
        break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
};
}
#endif
