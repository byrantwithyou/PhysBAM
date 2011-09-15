//#####################################################################
// Copyright 2004, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MULTIPLE_TORI
//#####################################################################
#ifndef __MULTIPLE_TORI__
#define __MULTIPLE_TORI__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Solids_And_Fluids/SOLIDS_3D_EXAMPLE.h>
namespace PhysBAM{

template<class T,class RW>
class MULTIPLE_TORI:public SOLIDS_3D_EXAMPLE<T,RW>
{
public:
    int m,n,mn;
    T base_height;
    T max_angular_velocity;

    QUATERNION<T> initial_orientation1,initial_orientation2;
    VECTOR_3D<T> initial_velocity1,initial_angular_velocity1,initial_velocity2,initial_angular_velocity2;
    std::string input_file;

    MULTIPLE_TORI()
        :m(2),n(11),mn(2),base_height(1.5),max_angular_velocity(2)
    {   
        collisions_repulsion_thickness=(T)1e-2;
        collisions_repulsion_clamp_fraction=(T).9;
        collision_repulsion_spring_multiplier=100;

        last_frame=10*24;
        //last_frame=17;//10*24;
        restart=false;restart_frame=16;   
        cfl=(T)10;
        cg_tolerance=(T)1e-2;
        output_directory="Multiple_Tori/output";
        perform_self_collision=true;
        solids_parameters.collide_with_interior=true;
    }

    ~MULTIPLE_TORI()
    {}

    void Postprocess_Substep(const T time,const int substep)
    {for(int o=1;o<=solids_parameters.deformable_body_parameters.list.deformable_objects.m;o++)
        std::cout<<"minimum volume for torus "<<o<<" = "<<solids_parameters.deformable_body_parameters.list(o).tetrahedralized_volume->Minimum_Signed_Volume()<<std::endl;}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    TETRAHEDRON_MESH torus_mesh;
    PARTICLES<T,VECTOR_3D<T> > torus_particles;
    TETRAHEDRALIZED_VOLUME<T> torus_volume(torus_mesh,torus_particles);
    FILE_UTILITIES::Read_From_File<RW>(data_directory+"/Tetrahedralized_Volumes/adaptive_torus_float.tet",torus_volume);
    std::cout<<"torus vertices = "<<torus_particles.array_collection->Size()<<"\ntorus tetrahedra = "<<torus_mesh.tetrahedrons.m<<"\n";

    torus_volume.Update_Bounding_Box();
    VECTOR_3D<T> center=torus_volume.bounding_box->Center();
    T size=torus_volume.bounding_box->Size().Max();
    std::cout<<"torus size = "<<size<<"\n";
    GRID<TV> grid(m,n,mn,0,(m-1)*size,base_height,base_height+(n-1)*size,0,(mn-1)*size);
    RANDOM_NUMBERS random;random.Set_Seed(12321);

    for(int i=1;i<=m;i++)for(int j=1;j<=n;j++)for(int ij=1;ij<=mn;ij++){
        int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Tetrahedralized_Volume();
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(index).tetrahedralized_volume;
        TETRAHEDRON_MESH& tetrahedron_mesh=tetrahedralized_volume.tetrahedron_mesh;
        PARTICLES<T,VECTOR_3D<T> >& particles=tetrahedralized_volume.particles;
        tetrahedron_mesh.Initialize_Tetrahedron_Mesh(torus_mesh);particles.Initialize_Particles(torus_particles);
        particles.Update_Velocity();particles.Store_Mass();
        tetrahedralized_volume.Set_Density(1000);tetrahedralized_volume.Set_Mass_Of_Particles(solids_parameters.use_constant_mass);

        VECTOR_3D<T> new_center=grid.Node(i,j,ij);
        std::cout<<i<<" "<<j<<" "<<ij<<" "<<grid.dy<<"\n";
        VECTOR_3D<T> angular_velocity=random.Get_Uniform_Number((T)0,max_angular_velocity)*random.Get_Direction();
        QUATERNION<T> orientation;
        do{orientation.s=random.Get_Uniform_Number(0,1);orientation.v=random.Get_Uniform_Vector(VECTOR_3D<T>(-1,-1,-1),VECTOR_3D<T>(1,1,1));}while(orientation.Magnitude()>1);
        orientation.Normalize();
        std::cout<<"Adding torus "<<index<<" at center "<<new_center<<", orientation "<<orientation<<", angular velocity "<<angular_velocity<<"\n";
        for(int p=1;p<=particles.array_collection->Size();p++){
            VECTOR_3D<T> dX=particles.X(p)-center;
            particles.V(p)=VECTOR_3D<T>::Cross_Product(angular_velocity,dX);
            particles.X(p)=new_center+orientation.Rotate(dX);}}
         
    int index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground",1,true,true,false);
    solids_parameters.rigid_body_parameters.list(index).Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list(index).is_static=true;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    SOLIDS_3D_EXAMPLE<T,RW>::Initialize_Bodies();
    for(int o=1;o<=solids_parameters.deformable_body_parameters.list.deformable_objects.m;o++){
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*solids_parameters.deformable_body_parameters.list(o).tetrahedralized_volume;
        solids_parameters.deformable_body_parameters.list(o).Add_Body_Forces(tetrahedralized_volume);
        solids_parameters.deformable_body_parameters.list(o).Add_Diagonalized_Neo_Hookean_Elasticity(tetrahedralized_volume,(T)2e5,(T).45,(T).01);}
}
//#####################################################################
};
}
#endif
