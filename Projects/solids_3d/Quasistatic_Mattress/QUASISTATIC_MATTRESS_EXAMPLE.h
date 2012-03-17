//#####################################################################
// Copyright 2004-2006, Mike Rodgers, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUASISTATIC_MATTRESS_EXAMPLE
//#####################################################################
#ifndef __QUASISTATIC_MATTRESS_EXAMPLE__
#define __QUASISTATIC_MATTRESS_EXAMPLE__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COLLISION_PENALTY_FORCES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_NEO_HOOKEAN_3D.h>
#include <Forces_And_Torques/DIAGONALIZED_FINITE_VOLUME_3D.h>
namespace PhysBAM{

template<class T,class RW>
class QUASISTATIC_MATTRESS_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::solids_parameters;using BASE::fluids_parameters;
    using BASE::output_directory;

    T initial_height;
    QUATERNION<T> initial_orientation;
    ARRAY<VECTOR<T,3> > X_save;
    VECTOR<T,3> attachment_velocity;
    VECTOR<T,3> attachment_angular_velocity;
    int m_input,n_input,mn_input;
    T xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input; 
    GRID<TV> mattress_grid;
    bool cube_mesh;
    bool static_example;
    VECTOR<T,3> end1_center_of_mass,end2_center_of_mass;
    VECTOR<T,3> ball_velocity,initial_ball_position;
    COLLISION_PENALTY_FORCES<T>* penalty_forces;
    
    QUASISTATIC_MATTRESS_EXAMPLE()
        :BASE(0,fluids_parameters.NONE),initial_height(4),initial_orientation(0,VECTOR<T,3>(1,0,0)),
        m_input(6),n_input(10),mn_input(20),
        xmin_input(-.25),xmax_input(.25),ymin_input(-.5),ymax_input(.5),zmin_input(-1),zmax_input(1),
        mattress_grid(m_input,n_input,mn_input,xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input),
        cube_mesh(true),static_example(false)
    {
        solids_parameters.quasistatic=true;

        last_frame=10*24;
        restart=false;restart_frame=0;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
        solids_parameters.implicit_solve_parameters.cg_iterations=900;
        solids_parameters.newton_iterations=10;
        solids_parameters.newton_tolerance=(T)1e-2;
        output_directory="Quasistatic_Mattress/output";
        attachment_velocity=VECTOR<T,3>(0,0,(T).1);
        attachment_angular_velocity=VECTOR<T,3>(0,0,(T)0);
        solids_parameters.use_partially_converged_result=true;
        solids_parameters.collision_tolerance=(T)1e-4;
        solids_parameters.perform_self_collision=false;
        solids_parameters.perform_collision_body_collisions=false;
        ball_velocity=VECTOR<T,3>((T)0,(T).5,(T)0);
        initial_ball_position=VECTOR<T,3>((T)0,(T)3.2,(T)0);
        if(static_example) solids_parameters.newton_iterations=1;
    }

    ~QUASISTATIC_MATTRESS_EXAMPLE()
    {delete penalty_forces;}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    Get_Initial_Data();

    //deformable bodies
    //deformable_body_parameters.list.deformable_objects(1)->Add_Body_Forces(*deformable_body_parameters.list(1).tetrahedralized_volume,gravity,gravity_direction);
    //body_forces=(BODY_FORCES_3D<T>*)deformable_body_parameters.list.deformable_objects(1)->solids_forces(1);
    solid_body_collection.Add_Force(Create_Quasistatic_Diagonalized_Finite_Volume(
        solid_body_collection.deformable_object.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(),new DIAGONALIZED_NEO_HOOKEAN_3D<T>((T)1e5,(T).45)));
    //deformable_body_parameters.list.deformable_objects(1)->Add_Quasistatic_Diagonalized_Linear_Elasticity(*deformable_body_parameters.list(1).tetrahedralized_volume,(T)5e4);
    //deformable_body_parameters.list.deformable_objects(1)->Add_Quasistatic_Diagonalized_Linear_Finite_Volume(*deformable_body_parameters.list(1).tetrahedralized_volume,(T)5e4);
    //deformable_body_parameters.list.deformable_objects(1)->Add_Edge_Springs(deformable_body_parameters.list(1).tetrahedralized_volume->tetrahedron_mesh,3000,(T).9);
    //deformable_body_parameters.list.deformable_objects(1)->Add_Quasistatic_Neo_Hookean_Elasticity(*deformable_body_parameters.list(1).tetrahedralized_volume,(T)5e5);
    //deformable_body_parameters.list.deformable_objects(1)->Add_Linear_Finite_Volume(*deformable_body_parameters.list(1).tetrahedralized_volume,1e5);
    //deformable_body_parameters.list.deformable_objects(1)->Add_Diagonalized_Spline_Model(*deformable_body_parameters.list(1).tetrahedralized_volume,5e4,.475,.05,0,1);

    //if(solids_parameters.perform_collision_body_collisions){
    //    //add prototype penalty forces
    //    penalty_forces=new COLLISION_PENALTY_FORCES<T>(solids_parameters.deformable_body_parameters.list.deformable_objects(1)->tetrahedralized_volume->particles);
    //    penalty_forces->Set_Collision_Body_List(solids_parameters.collision_body_list);
    //    penalty_forces->Set_Boundary_Only_Collisions(solids_parameters.deformable_body_parameters.list.deformable_objects(1)->tetrahedralized_volume->tetrahedron_mesh);
    //    solids_parameters.deformable_body_parameters.list.deformable_objects(1)->solids_forces.Append(penalty_forces);}

    solid_body_collection.Update_Fragments();
}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    DEFORMABLE_OBJECT_3D<T>& deformable_object=solid_body_collection.deformable_object;
    deformable_object.Add_Structure(TETRAHEDRALIZED_VOLUME<T>::Create(deformable_object.particles));
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_object.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
    TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.mesh;
    DEFORMABLE_PARTICLES<T,VECTOR<T,3> >& particles=deformable_object.particles;

    if(!cube_mesh) tetrahedralized_volume.Initialize_Octahedron_Mesh_And_Particles(mattress_grid);
    else tetrahedralized_volume.Initialize_Cube_Mesh_And_Particles(mattress_grid);
    std::cout << "total vertices = " << particles.Size() << std::endl;std::cout << "total tetrahedra = " << tetrahedron_mesh.elements.m << std::endl;

    particles.Store_Velocity(false);particles.Store_Mass();
    tetrahedralized_volume.Set_Density(1000);tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);

    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR<T,3> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<tetrahedralized_volume.particles.array_size;i++){
        tetrahedralized_volume.particles.X(i)=center+initial_orientation.Rotate(tetrahedralized_volume.particles.X(i)-center);
        tetrahedralized_volume.particles.X(i).y+=initial_height-bottom;}

    X_save.Resize(tetrahedralized_volume.particles.X.array.m);ARRAY<VECTOR<T,3> >::copy(tetrahedralized_volume.particles.X.array,X_save);

    int index=0;
    for(int j=0;j<n_input;j++)for(int k=0;k<m_input;k++){
        index++;end1_center_of_mass+=X_save(index);end2_center_of_mass+=X_save(X_save.m+1-index);}
    end1_center_of_mass=((T)1/(T)index)*end1_center_of_mass;end2_center_of_mass=((T)1/(T)index)*end2_center_of_mass;

    //rigid bodies
    if(solids_parameters.perform_collision_body_collisions){
        int index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<T>("../../Public_Data/Rigid_Bodies/sphere",(T).75);
        solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T)0);
        solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);}
    else{solids_parameters.rigid_body_evolution_parameters.write_rigid_bodies=false;}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<VECTOR<T,3> > X,const T time,const int id_number) PHYSBAM_OVERRIDE
{
    switch(id_number){
    case 1:
        for(int j=1,index=0;j<=n_input;j++)for(int k=0;k<m_input;k++){
            index++;if(static_example){X(index)=X_save(index)+attachment_velocity;
            X(X.m+1-index)=X_save(X.m+1-index)-attachment_velocity;}
            else{
                QUATERNION<T> orientation,opposite_orientation;
                T magnitude=attachment_angular_velocity.Magnitude();
                if(magnitude>1e-6){
                    orientation=QUATERNION<T>(time*magnitude,attachment_angular_velocity);
                    opposite_orientation=QUATERNION<T>(-time*magnitude,attachment_angular_velocity);}
                X(index)=end1_center_of_mass+orientation.Rotate(X_save(index)-end1_center_of_mass)+time*attachment_velocity;
                X(X.m+1-index)=end2_center_of_mass+opposite_orientation.Rotate(X_save(X.m+1-index)-end2_center_of_mass)-time*attachment_velocity;}}
        break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<VECTOR<T,3> > X,const T time,const int id_number) PHYSBAM_OVERRIDE
{
    switch(id_number){
    case 1:
        for(int j=1,index=0;j<=n_input;j++)for(int k=0;k<m_input;k++){index++;X(index)=VECTOR<T,3>(0,0,0);X(X.m+1-index)=VECTOR<T,3>(0,0,0);}
        break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
// Function Update_Collision_Body_Positions_And_Velocities(const T time)
//#####################################################################
void Update_Collision_Body_Positions_And_Velocities(const T time) PHYSBAM_OVERRIDE
{   
    //if(time>(T)4) return;
    if(solids_parameters.perform_collision_body_collisions) solids_parameters.rigid_body_parameters.list.rigid_bodies(1)->frame.t=time*ball_velocity+initial_ball_position;
}
//#####################################################################
};
}
#endif
