//#####################################################################
// Copyright 2004, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FACE_SMASH
//#####################################################################
#ifndef __FACE_SMASH__
#define __FACE_SMASH__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <Constitutive_Models/DIAGONALIZED_NEO_HOOKEAN_3D.h>
#include <Constitutive_Models/DIAGONALIZED_SPLINE_MODEL_3D.h>
#include <Forces_And_Torques/DIAGONALIZED_FINITE_VOLUME_3D.h>
namespace PhysBAM{

template<class T,class RW>
class FACE_SMASH:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::frame_rate;using BASE::output_directory;using BASE::data_directory;using BASE::fluids_parameters;
    using BASE::solids_parameters;

    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity;
    std::string input_file;
    int force,tweak;
    int strength;
    VECTOR_3D<T> gravity;

    FACE_SMASH()
        :BASE(0,fluids_parameters.NONE),initial_height(1),initial_orientation(pi/2,VECTOR_3D<T>(-1,0,0)),strength(1)
    {   
        frame_rate*=5;
        last_frame=100;
        restart=false;restart_frame=38;

        force=3;
        tweak=3;

        output_directory=STRING_UTILITIES::string_sprintf("Face/output%d_%d",force,tweak);
        if(strength!=1) output_directory+=STRING_UTILITIES::string_sprintf("_%d",strength);
        input_file=data_directory+"/Face_Data/Eftychis_700k/Front_300k/face_simulation.tet";
        solids_parameters.perform_self_collision=true;
        solids_parameters.collide_with_interior=true;

/*
        solids_parameters.collisions_repulsion_thickness=(T)1e-3;
        solids_parameters.collisions_repulsion_clamp_fraction=(T).9;
        solids_parameters.collision_repulsion_spring_multiplier=100;
*/

        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        solids_parameters.asynchronous_average_steps_between_self_collisions=.5;
        gravity=solids_parameters.gravity*solids_parameters.gravity_direction;
    }

    ~FACE_SMASH()
    {}

    void Postprocess_Substep(const T time,const int substep)
    {std::cout << "minimum volume = " << solids_parameters.deformable_body_parameters.list(1).tetrahedralized_volume->Minimum_Signed_Volume() << std::endl;}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    Get_Initial_Data();

    DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*deformable_object.tetrahedralized_volume;
    solid_body_collection.Add_Force(Create_Body_Forces(tetrahedralized_volume,solids_parameters.gravity));
    T max_strain=0; 
    
    if(force==1){ // diagonalized finite volume
        if(tweak==1){solids_parameters.cfl=10;max_strain=(T).1;}
        if(tweak==2){solids_parameters.cfl=20;max_strain=(T).1;}
        if(tweak==3){solids_parameters.cfl=10;max_strain=(T).15;}
        solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tetrahedralized_volume,new DIAGONALIZED_NEO_HOOKEAN_3D<T>(strength*(T)20000,(T).45,(T).01,(T).25),true,max_strain));}
    if(force==2){ // semi-implicit diagonalized finite volume
        solids_parameters.semi_implicit=true;
        if(tweak==1){solids_parameters.cfl=10;max_strain=(T).1;}
        if(tweak==3){solids_parameters.cfl=10;max_strain=(T).15;}
        if(tweak==4){solids_parameters.cfl=20;max_strain=(T).15;}
        if(tweak==5){solids_parameters.cfl=20;max_strain=(T).2;}
        if(tweak==6){solids_parameters.cfl=30;max_strain=(T).3;}
        if(tweak==7){solids_parameters.cfl=200;max_strain=(T)2;}
        if(tweak==8){solids_parameters.cfl=2000;max_strain=(T)20;}
        solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tetrahedralized_volume,new DIAGONALIZED_NEO_HOOKEAN_3D<T>(strength*(T)20000,(T).45,(T).01,(T).25),true,max_strain));}
    if(force==3){ // asynchronous semi-implicit diagonalized finite volume
        solids_parameters.asynchronous=true;
        if(tweak==1){solids_parameters.cfl=1;max_strain=(T).05;}
        if(tweak==2){solids_parameters.cfl=.7;max_strain=(T).03;}
        if(tweak==3){solids_parameters.cfl=1;max_strain=(T).05;strength=10;solids_parameters.asynchronous_average_steps_between_self_collisions=2;}
        //deformable_object.Add_Diagonalized_Neo_Hookean_Elasticity(tetrahedralized_volume,strength*(T)20000,(T).45,(T).02,(T).25,true,max_strain);}
        solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(tetrahedralized_volume,new DIAGONALIZED_SPLINE_MODEL_3D<T>(strength*(T)200000,(T).45,(T).5,(T)5,(T).05),true,max_strain));}
}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Tetrahedralized_Volume();
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(index).tetrahedralized_volume;
    TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.tetrahedron_mesh;
    PARTICLES<T,VECTOR_3D<T> >& particles=tetrahedralized_volume.particles;

    FILE_UTILITIES::Read_From_File<RW>(input_file,tetrahedralized_volume);
    std::cout << "total vertices = " << particles.array_collection->Size() << std::endl;std::cout << "total tetrahedra = " << tetrahedron_mesh.tetrahedrons.m << std::endl;
    particles.Update_Velocity();particles.Store_Mass(); // in case they're not stored in the file!

    tetrahedralized_volume.Rescale(10);
    tetrahedralized_volume.Set_Density(1000);tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<particles.array_collection->Size();i++){
        particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.X(i).y+=initial_height-bottom;}
    tetrahedralized_volume.Update_Bounding_Box();
    std::cout<<"bounding box = "<<*tetrahedralized_volume.bounding_box<<std::endl;

    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Add_External_Impulses
//#####################################################################
void Add_External_Impulses(ARRAY<VECTOR_3D<T> >& V,const T time,const T dt)
{
    V+=dt*gravity; // since semi-implicit doesn't notice body forces yet 
}
//#####################################################################
// Function Add_External_Impulse
//#####################################################################
void Add_External_Impulse(ARRAY<VECTOR_3D<T> >& V,const int node,const T time,const T dt)
{
    V(node)+=dt*gravity; // since semi-implicit doesn't notice body forces yet 
}
//#####################################################################
};
}
#endif
