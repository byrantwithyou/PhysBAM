//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KINEMATIC_DEFORMABLE_AND_RIGID
//#####################################################################
#ifndef __KINEMATIC_DEFORMABLE_AND_RIGID__
#define __KINEMATIC_DEFORMABLE_AND_RIGID__

#include <Geometry/LEVELSET_IMPLICIT_SURFACE.h>
#include <Rigid_Bodies/DEFORMABLE_OBJECT_RIGID_BODY_WRAPPER_3D.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template <class T,class RW>
class KINEMATIC_DEFORMABLE_AND_RIGID:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::verbose_dt;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::initial_time;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;

    int example_number;
    GRID<TV> cloth_grid;
    DEFORMABLE_OBJECT_RIGID_BODY_WRAPPER_3D<T>* stub_rigid_body;
    std::string presimulation_output_directory;

    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}

    KINEMATIC_DEFORMABLE_AND_RIGID(const int example_number_input=0)
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.NONE),example_number(example_number_input),stub_rigid_body(0)
    {
        last_frame=800;
        restart=false;restart_frame=0;   
        solids_parameters.cfl=(T).5;
        output_directory="Kinematic_Deformable_And_Rigid/output";
        presimulation_output_directory="Kinematic_Deformable_And_Rigid/deformable_output/";
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.perform_self_collision=false;
        verbose_dt=true;
        if(example_number==0){
//            solids_parameters.gravity=0;
        }
        else if(example_number==1){
            output_directory=presimulation_output_directory;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=false;
            solids_parameters.gravity_direction.y=1;
        }
    }

    ~KINEMATIC_DEFORMABLE_AND_RIGID()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    int id;
    if(example_number==0){
        T epsilon=0.5;

        solids_parameters.deformable_body_parameters.list.template Read_Static_Variables<RW>(presimulation_output_directory);
        solids_parameters.deformable_body_parameters.list.template Read_Dynamic_Variables<RW>(presimulation_output_directory,0);

        DEFORMABLE_OBJECT_3D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
        deformable_object.simulate=false;
        
        {
            stub_rigid_body=new DEFORMABLE_OBJECT_RIGID_BODY_WRAPPER_3D<T>();
            stub_rigid_body->Initialize_Triangulated_Surface(*deformable_object.triangulated_surface);
            id=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body(stub_rigid_body,0,0,0);
            stub_rigid_body->Set_Name("deforming body");
            stub_rigid_body->Set_Coefficient_Of_Restitution(epsilon);
            stub_rigid_body->is_static=true;

            // Artificially enlarge bounding box to avoid small timesteps
            if(stub_rigid_body->triangulated_surface->bounding_box){
                BOX_3D<T>& box=*stub_rigid_body->triangulated_surface->bounding_box;
                if(box.Size().y<1){box.ymin-=1;}
                stub_rigid_body->Update_Bounding_Box();}
        }

        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere",2);
        solids_parameters.rigid_body_parameters.list(id)->Set_Name("first object");
        solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(epsilon);
        solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_3D<T>(1,5,0);
        solids_parameters.rigid_body_parameters.list(id)->velocity=VECTOR_3D<T>(0,0,0);
        solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

        solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
    }
    else if(example_number==1){
        cloth_grid=GRID_2D<T>(20,20,-10,10,-10,10);
        id=Add_Deformable_Object(cloth_grid,MATRIX<T,4>::Rotation_Matrix_X_Axis(pi/2));
        DEFORMABLE_OBJECT_3D<T>* deformable_object=solids_parameters.deformable_body_parameters.list.deformable_objects(id);
        T edge_stiffness_scaling=200,altitude_stiffness_scaling=100;
        deformable_object->triangulated_surface->Set_Density(1);
        deformable_object->triangulated_surface->Set_Mass_Of_Particles();
        deformable_object->Add_Body_Forces(*deformable_object->triangulated_surface,solids_parameters.gravity,solids_parameters.gravity_direction);
        deformable_object->Add_Edge_Springs(deformable_object->triangulated_surface->triangle_mesh,edge_stiffness_scaling/(1+sqrt((T)2)),4);
        deformable_object->Add_Altitude_Springs(deformable_object->triangulated_surface->triangle_mesh,altitude_stiffness_scaling*4/(1+sqrt((T)2)),8);
        deformable_object->Add_Bending_Elements(deformable_object->triangulated_surface->triangle_mesh,10);
    }
    else if(example_number==2){
        T epsilon=0.7;

        LEVELSET_IMPLICIT_SURFACE<T>* implicit_surface=LEVELSET_IMPLICIT_SURFACE<T>::Create();
        implicit_surface->levelset.grid=GRID_3D<T>(61,31,31,-8,8,-4,4,-4,4);
        implicit_surface->levelset.phi.Resize(implicit_surface->levelset.grid);
        ARRAY<VECTOR_3D<T> ,VECTOR<int,3> >* velocity_field=new ARRAY<VECTOR_3D<T> ,VECTOR<int,3> >(implicit_surface->levelset.grid);
        stub_rigid_body=new DEFORMABLE_OBJECT_RIGID_BODY_WRAPPER_3D<T>();
        stub_rigid_body->Initialize_Implicit_Surface(*implicit_surface);
        stub_rigid_body->Set_Velocity_Field(&implicit_surface->levelset.grid,velocity_field);
        Update_Deforming_Volume(0);
        int implicit_surface_id=solids_parameters.rigid_body_parameters.list.implicit_surface_list.Add_Element(implicit_surface);
        id=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body(stub_rigid_body,0,implicit_surface_id,0);
        stub_rigid_body->Set_Name("deforming body");
        stub_rigid_body->Set_Coefficient_Of_Restitution(epsilon);
        stub_rigid_body->is_static=true;

        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere");
        solids_parameters.rigid_body_parameters.list(id)->Set_Name("first object");
        solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(epsilon);
        solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_3D<T>(-2,5,0);
        solids_parameters.rigid_body_parameters.list(id)->velocity=VECTOR_3D<T>(0,0,0);
        solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("../../Public_Data/Rigid_Bodies/ground");
        solids_parameters.rigid_body_parameters.list(id)->orientation=QUATERNION<T>(pi/2,VECTOR_3D<T>(0,0,1));
        solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
        solids_parameters.rigid_body_parameters.list(id)->is_static=true;

        solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
    }
}
//#####################################################################
// Function Add_Circle_Deformable_Object
//#####################################################################
int Add_Deformable_Object(const GRID<TV>& cloth_grid,const MATRIX<T,4>& transform)
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Triangulated_Surface();
    TRIANGULATED_SURFACE<T>& triangulated_surface=*solids_parameters.deformable_body_parameters.list(index).triangulated_surface;
    TRIANGLE_MESH& triangle_mesh=triangulated_surface.triangle_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_3D<T> >& particles=triangulated_surface.particles;

    triangle_mesh.Initialize_Herring_Bone_Mesh(cloth_grid.m,cloth_grid.n);
    particles.Add_Elements(triangle_mesh.number_nodes);
    for(int i=0;i<cloth_grid.m;i++) for(int j=0;j<cloth_grid.n;j++){
        int node=i+cloth_grid.m*(j-1);particles.X(node)=transform*VECTOR_3D<T>(cloth_grid.X(i,j));particles.V(node)=VECTOR_3D<T>();}

    return index;
}
//#####################################################################
// Function Update_Deforming_Volume
//#####################################################################
void Update_Deforming_Volume(const T time)
{
    LEVELSET_3D<GRID<TV> >& levelset=((LEVELSET_IMPLICIT_SURFACE<T>*)stub_rigid_body->implicit_surface)->levelset;
//    VECTOR_3D<T> center(0,0);T radius=0.4;
//    for(int i=0;i<levelset.grid.m;i++) for(int j=0;j<levelset.grid.n;j++)
//        levelset.phi(i,j)=(levelset.grid.X(i,j)-center).Magnitude()-radius;
    VECTOR_3D<T> velocity(0,0,0),corner=VECTOR_3D<T>(-4,-1,-1)+time*velocity,size(8,2,2);
    VECTOR_3D<T> angular_velocity(0,0,-0.3);
    BOX_3D<T> box(corner,corner+size);
    MATRIX<T,3> rotation_matrix=MATRIX<T,3>::Rotation_Matrix(-time*angular_velocity);
    for(int i=0;i<levelset.grid.m;i++) for(int j=0;j<levelset.grid.n;j++) for(int k=0;k<levelset.grid.mn;k++){
        levelset.phi(i,j,k)=box.Signed_Distance(rotation_matrix*levelset.grid.X(i,j,k));
        if(stub_rigid_body->velocity_field) (*stub_rigid_body->velocity_field)(i,j,k)=velocity+VECTOR_3D<T>::Cross_Product(angular_velocity,levelset.grid.X(i,j,k));}
    stub_rigid_body->implicit_surface->Update_Box();
    stub_rigid_body->implicit_surface->Compute_Cell_Minimum_And_Maximum();
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(example_number==2 && id==int(1)){
//        VECTOR_2D<T> velocity(1,0);
//        state.position=time*velocity;
//        state.orientation=0.1*time;
    }
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    Get_Initial_Data();
    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Initialize_Bodies();
}
//#####################################################################
// Function Solids_Example_Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR_3D<T> >& V,const T time) PHYSBAM_OVERRIDE
{
    Zero_Out_Enslaved_Velocity_Nodes(V,time,id);
}
//#####################################################################
// Function Solids_Example_Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_3D<T> >& V,const T time) PHYSBAM_OVERRIDE
{
    if(example_number==1){
        for(int i=0;i<cloth_grid.m;i++){V(i)=V(i+cloth_grid.m*(cloth_grid.n-1))=VECTOR_3D<T>();}
        for(int j=0;j<cloth_grid.n;j++){V(1+cloth_grid.m*(j-1))=V(cloth_grid.m*j)=VECTOR_3D<T>();}}
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    std::cout << "In Update_Solids_Parameters (time=" << time << ")" << std::endl;
    int frame=(int)((time-initial_time)*frame_rate);
    if(example_number==0){
        solids_parameters.deformable_body_parameters.list.template Read_Dynamic_Variables<RW>(presimulation_output_directory,frame);
        stub_rigid_body->triangulated_surface->Refresh_Auxiliary_Structures();
        stub_rigid_body->Update_Bounding_Box();
    }
    else if(example_number==2){
        Update_Deforming_Volume(time);

        std::string filename=STRING_UTILITIES::string_sprintf("%s/levelset.%d",output_directory.c_str(),frame);
        LEVELSET_3D<GRID<TV> >& levelset=((LEVELSET_IMPLICIT_SURFACE<T>*)stub_rigid_body->implicit_surface)->levelset;
        FILE_UTILITIES::Write_To_File<RW>(filename,levelset);

        if(stub_rigid_body->velocity_field){
            filename=STRING_UTILITIES::string_sprintf("%s/velocities.%d",output_directory.c_str(),frame);
            FILE_UTILITIES::Write_To_File<RW>(filename,*stub_rigid_body->velocity_field);}
    }
}
//#####################################################################
};
}
#endif
