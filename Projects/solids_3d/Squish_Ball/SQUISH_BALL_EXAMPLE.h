//#####################################################################
// Copyright 2003-2005, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SQUISH_BALL_EXAMPLE
//#####################################################################
#ifndef __SQUISH_BALL_EXAMPLE__
#define __SQUISH_BALL_EXAMPLE__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template<class T,class RW>
class SQUISH_BALL_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;

    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity,initial_angular_velocity,fiber_direction;
    std::string input_file;
    ARRAY<T> muscle_activations;
    bool use_force_length;

    SQUISH_BALL_EXAMPLE()
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.NONE),initial_height(1),initial_orientation(0,VECTOR_3D<T>(1,0,0)),initial_velocity(0,0,0),initial_angular_velocity(0,0,0),fiber_direction((T).5*sqrt((T)2),(T).5*sqrt((T)2),0)
    {   
        last_frame=100;
        restart=true;restart_frame=20;
        solids_parameters.cfl=(T)50;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;

        output_directory=STRING_UTILITIES::string_sprintf("Squish_Ball/output");
        input_file=data_directory+"/Tetrahedralized_Volumes/sphere.tet";
        solids_parameters.perform_self_collision=false;
        solids_parameters.collide_with_interior=true;
        use_force_length=false;

    }

    ~SQUISH_BALL_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    Get_Initial_Data();

    DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*deformable_object.tetrahedralized_volume;
    solids_parameters.deformable_body_parameters.list(1).Add_Body_Forces(tetrahedralized_volume);
    solids_parameters.deformable_body_parameters.list(1).Add_Edge_Springs(tetrahedralized_volume.tetrahedron_mesh,3000,(T)1);
    solids_parameters.deformable_body_parameters.list(1).Add_Altitude_Springs(tetrahedralized_volume.tetrahedron_mesh,300,(T)1);

    //set up muscle related arrays
    ARRAY<ARRAY<VECTOR_3D<T> > > fiber_directions(1);fiber_directions(1).Resize(tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m);
    //initialize fiber direction
    for(int t=0;t<tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m;t++) fiber_directions(1)(t)=fiber_direction;
    ARRAY<ARRAY<int> > muscle_tets(1);ARRAY<ARRAY<T> > muscle_densities(1);
    muscle_tets(1).Resize(tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m);
    muscle_densities(1).Resize(tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m);
    muscle_activations.Resize(1);muscle_activations(1)=(T)0;
    //set up densities for one activation region
    for(int t=0;t<tetrahedralized_volume.tetrahedron_mesh.tetrahedrons.m;t++){muscle_tets(1)(t)=t;muscle_densities(1)(t)=(T)1;}
    solids_parameters.deformable_body_parameters.list.deformable_objects(1)->Add_Diagonalized_Fiber_Tension(tetrahedralized_volume,muscle_tets,fiber_directions,muscle_densities,muscle_activations,use_force_length);

    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Initialize_Bodies();
}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Tetrahedralized_Volume();
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(index).tetrahedralized_volume;
    TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.tetrahedron_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=tetrahedralized_volume.particles;

    FILE_UTILITIES::Read_From_File<RW>(input_file,tetrahedralized_volume);
    std::cout << "total vertices = " << particles.array_collection->Size() << std::endl;std::cout << "total tetrahedra = " << tetrahedron_mesh.tetrahedrons.m << std::endl;
    particles.Update_Velocity();particles.Store_Mass(); // in case they're not stored in the file!

    tetrahedralized_volume.Set_Density(1000);tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);
    tetrahedralized_volume.Update_Bounding_Box();
    VECTOR_3D<T> center(tetrahedralized_volume.bounding_box->Center());T bottom=tetrahedralized_volume.bounding_box->ymin;
    for(int i=0;i<particles.array_collection->Size();i++){
        particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.X(i).y+=initial_height-bottom;}

    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Update_Time_Varrying_Material_Properties
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time){
    if(time>1&&time<2) if(!use_force_length) muscle_activations(1)=(T).25;else muscle_activations(1)=(T)1;
    else muscle_activations(1)=(T)0;
}
//#####################################################################
};
}
#endif
