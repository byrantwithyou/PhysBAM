//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KINEMATIC_DEFORMABLE_AND_RIGID
//#####################################################################
#ifndef __KINEMATIC_DEFORMABLE_AND_RIGID__
#define __KINEMATIC_DEFORMABLE_AND_RIGID__

#include <Geometry/LEVELSET_IMPLICIT_CURVE.h>
#include <Rigid_Bodies/DEFORMABLE_OBJECT_RIGID_BODY_WRAPPER_2D.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_2D.h>
namespace PhysBAM{

template <class T,class RW>
class KINEMATIC_DEFORMABLE_AND_RIGID:public SOLIDS_FLUIDS_EXAMPLE_2D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart_frame;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::data_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::fluids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::verbose_dt;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::initial_time;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::frame_rate;

    int example_number;
    VECTOR_2D<T> center;
    T radius,speed;
    DEFORMABLE_OBJECT_RIGID_BODY_WRAPPER_2D<T>* stub_rigid_body;
    std::string presimulation_output_directory;

    KINEMATIC_DEFORMABLE_AND_RIGID(const int example_number_input=0)
        :SOLIDS_FLUIDS_EXAMPLE_2D<RW>(fluids_parameters.NONE),example_number(example_number_input),stub_rigid_body(0)
    {
        last_frame=800;
        restart=false;restart_frame=0;   
        solids_parameters.cfl=(T).5;
        output_directory="Kinematic_Deformable_And_Rigid/output";
        presimulation_output_directory="Kinematic_Deformable_And_Rigid/deformable_output/";
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        verbose_dt=true;
        if(example_number==1){
            output_directory=presimulation_output_directory;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=false;
        }
    }

    ~KINEMATIC_DEFORMABLE_AND_RIGID()
    {}

    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    int id;
    if(example_number==0){
        T epsilon=0.2;

        solids_parameters.deformable_body_parameters.list.template Read_Static_Variables<RW>(presimulation_output_directory);
        solids_parameters.deformable_body_parameters.list.template Read_Dynamic_Variables<RW>(presimulation_output_directory,0);

//        center=VECTOR_2D<T>(0,2);radius=1;speed=1;
//        id=Add_Circle_Deformable_Object(30);
//        DEFORMABLE_OBJECT_2D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(id);
        DEFORMABLE_OBJECT_2D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
        deformable_object.simulate=false;
        
        {
            stub_rigid_body=new DEFORMABLE_OBJECT_RIGID_BODY_WRAPPER_2D<T>();
            stub_rigid_body->Initialize_Segmented_Curve(*deformable_object.segmented_curve);
            id=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body(stub_rigid_body,0,0,0);
            //int segmented_curve_id=solids_parameters.rigid_body_parameters.list.segmented_curve_list.Add_Element(deformable_object.segmented_curve);
            //id=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body(rigid_body,segmented_curve_id,0,0);
            stub_rigid_body->Set_Name("deforming body");
            stub_rigid_body->Set_Coefficient_Of_Restitution(epsilon);
            stub_rigid_body->is_static=true;
        }

        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/circle");
        solids_parameters.rigid_body_parameters.list(id)->Set_Name("first object");
        solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(epsilon);
        solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_2D<T>(1,5);
        solids_parameters.rigid_body_parameters.list(id)->velocity=VECTOR_2D<T>(0,0);
        solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D(),solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

#if 0
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("../../Public_Data/Rigid_Bodies_2D/ground");
        solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
        solids_parameters.rigid_body_parameters.list(id)->is_static=true;
#endif

        solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
    }
    else if(example_number==1){
        id=Add_Deformable_Object(80,VECTOR_2D<T>(-10,0),VECTOR_2D<T>(10,0));
        DEFORMABLE_OBJECT_2D<T>* deformable_object=SOLIDS_FLUIDS_EXAMPLE_2D<RW>::solids_parameters.deformable_body_parameters.list.deformable_objects(id);
        deformable_object->segmented_curve->Set_Density(1);
        deformable_object->segmented_curve->Set_Mass_Of_Particles();
        deformable_object->Add_Body_Forces(*deformable_object->segmented_curve,solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D());
        deformable_object->Add_Edge_Springs(deformable_object->segmented_curve->segment_mesh,3000,(T).9);
    }
    else if(example_number==2){
        T epsilon=0.7;

        LEVELSET_IMPLICIT_CURVE<T>* implicit_curve=LEVELSET_IMPLICIT_CURVE<T>::Create();
        implicit_curve->levelset.grid=GRID_2D<T>(101,51,-8,8,-4,4);
        implicit_curve->levelset.phi.Resize(implicit_curve->levelset.grid);
        ARRAY<VECTOR_2D<T> ,VECTOR<int,2> >* velocity_field=new ARRAY<VECTOR_2D<T> ,VECTOR<int,2> >(implicit_curve->levelset.grid);
        stub_rigid_body=new DEFORMABLE_OBJECT_RIGID_BODY_WRAPPER_2D<T>();
        stub_rigid_body->Initialize_Implicit_Curve(*implicit_curve);
        stub_rigid_body->Set_Velocity_Field(&implicit_curve->levelset.grid,velocity_field);
        Update_Deforming_Volume(0);
        //int segmented_curve_id=solids_parameters.rigid_body_parameters.list.segmented_curve_list.Add_Element(deformable_object.segmented_curve);
        int implicit_curve_id=solids_parameters.rigid_body_parameters.list.implicit_curve_list.Add_Element(implicit_curve);
        id=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body(stub_rigid_body,0,implicit_curve_id,0);
        //id=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body(stub_rigid_body,0,0,0);
        stub_rigid_body->Set_Name("deforming body");
        stub_rigid_body->Set_Coefficient_Of_Restitution(epsilon);
        stub_rigid_body->rigid_body_collection.rigid_body_particle.kinematic(rigid_body->particle_index)=true;

        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/circle");
        solids_parameters.rigid_body_parameters.list(id)->Set_Name("first object");
        solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Restitution(epsilon);
        solids_parameters.rigid_body_parameters.list(id)->position=VECTOR_2D<T>(-2,5);
        solids_parameters.rigid_body_parameters.list(id)->velocity=VECTOR_2D<T>(0,0);
        solids_parameters.rigid_body_parameters.list(id)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D(),solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

#if 1
        id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("../../Public_Data/Rigid_Bodies_2D/ground");
        solids_parameters.rigid_body_parameters.list(id)->orientation=pi/2;
        solids_parameters.rigid_body_parameters.list(id)->Set_Coefficient_Of_Friction((T).3);
        solids_parameters.rigid_body_parameters.list(id)->is_static=true;
#endif

        solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
    }
}
//#####################################################################
// Function Add_Circle_Deformable_Object
//#####################################################################
int Add_Deformable_Object(const int number_of_vertices,const VECTOR_2D<T>& start_position,const VECTOR_2D<T>& end_position)
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Segmented_Curve();
    SEGMENTED_CURVE_2D<T>& segmented_curve=*solids_parameters.deformable_body_parameters.list(index).segmented_curve;
    SEGMENT_MESH& segment_mesh=segmented_curve.segment_mesh;
    PARTICLES<T,VECTOR_2D<T> >& particles=segmented_curve.particles;
    segment_mesh.Initialize_Straight_Mesh(number_of_vertices);
    for(int i=0;i<number_of_vertices;i++){
        int index=particles.array_collection->Add_Element();assert(index==i);
        particles.X(i)=start_position+((T)(i-1)/(number_of_vertices-1))*(end_position-start_position);
        particles.V(i)=VECTOR_2D<T>(0,0);}
    return index;
}
//#####################################################################
// Function Add_Circle_Deformable_Object
//#####################################################################
int Add_Circle_Deformable_Object(const int number_of_vertices)
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Segmented_Curve();
    SEGMENTED_CURVE_2D<T>& segmented_curve=*solids_parameters.deformable_body_parameters.list(index).segmented_curve;
    SEGMENT_MESH& segment_mesh=segmented_curve.segment_mesh;
    PARTICLES<T,VECTOR_2D<T> >& particles=segmented_curve.particles;
    segment_mesh.Initialize_Straight_Mesh(number_of_vertices,true);
    for(int i=0;i<number_of_vertices;i++){
        int index=particles.array_collection->Add_Element();assert(index==i);
        T angle=(i-1)*(T)2*pi/number_of_vertices;
        VECTOR_2D<T> radial_direction(cos(angle+.5*pi),sin(angle+.5*pi));
        particles.X(i)=center+radius*radial_direction;
        particles.V(i)=speed*radial_direction;}
    return index;
}
//#####################################################################
// Function Update_Deforming_Volume
//#####################################################################
void Update_Deforming_Volume(const T time)
{
    LEVELSET_2D<T>& levelset=((LEVELSET_IMPLICIT_CURVE<T>*)stub_rigid_body->implicit_curve)->levelset;
//    VECTOR_2D<T> center(0,0);T radius=0.4;
//    for(int i=0;i<levelset.grid.m;i++) for(int j=0;j<levelset.grid.n;j++)
//        levelset.phi(i,j)=(levelset.grid.X(i,j)-center).Magnitude()-radius;
    VECTOR_2D<T> velocity(0,0),corner=VECTOR_2D<T>(-4,-1)+time*velocity,size(8,2);
    T angular_velocity=-0.3;
    BOX_2D<T> box(corner.x,corner.x+size.x,corner.y,corner.y+size.y);
    T orientation=angular_velocity*time;MATRIX<T,2> rotation_matrix=MATRIX<T,2>::Rotation_Matrix(-orientation);
    for(int i=0;i<levelset.grid.m;i++) for(int j=0;j<levelset.grid.n;j++){
        levelset.phi(i,j)=box.Signed_Distance(rotation_matrix*levelset.grid.X(i,j));
#if 1
        if(stub_rigid_body->velocity_field){
            VECTOR_2D<T> arm=levelset.grid.X(i,j);
            (*stub_rigid_body->velocity_field)(i,j)=velocity+angular_velocity*VECTOR_2D<T>(-arm.y,arm.x);}
#endif
    }
    stub_rigid_body->implicit_curve->Update_Box();
    stub_rigid_body->implicit_curve->Compute_Cell_Minimum_And_Maximum();
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(example_number==2 && id==int(1)){
//        VECTOR_2D<T> velocity(1,0);twist.t=time*velocity;twist.r=0.1*time;
    }
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    Get_Initial_Data();
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Initialize_Bodies();
}
//#####################################################################
// Function Solids_Example_Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR_2D<T> >& V,const T time) PHYSBAM_OVERRIDE
{
    Zero_Out_Enslaved_Velocity_Nodes(V,time,fragment_id);
}
//#####################################################################
// Function Solids_Example_Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_2D<T> >& V,const T time) PHYSBAM_OVERRIDE
{
    if(example_number==1) V(1)=V(V.m)=VECTOR_2D<T>(0,0);
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    int frame=(int)((time-initial_time)*frame_rate);
    if(example_number==0){
        std::cout << "In Update_Solids_Parameters (time=" << time << ")" << std::endl;

#if 1
        solids_parameters.deformable_body_parameters.list.template Read_Dynamic_Variables<RW>(presimulation_output_directory,frame);
#else
        DEFORMABLE_OBJECT_2D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);
        PARTICLES<T,VECTOR_2D<T> >& particles=deformable_object.particles;
        for(int i=0;i<particles.array_collection->Size();i++){
            particles.X(i)=center+radius*particles.V(i).Normalized()+particles.V(i)*time;}

        SEGMENTED_CURVE_2D<T>* segmented_curve=solids_parameters.deformable_body_parameters.list(1).segmented_curve;
#endif

        stub_rigid_body->segmented_curve->Refresh_Auxiliary_Structures();
        stub_rigid_body->Update_Bounding_Box();
    }
    else if(example_number==2){
        std::cout << "In Update_Solids_Parameters (time=" << time << ")" << std::endl;
        Update_Deforming_Volume(time);

        std::string filename=STRING_UTILITIES::string_sprintf("%s/%d/levelset",output_directory.c_str(),frame);
        LEVELSET_2D<T>& levelset=((LEVELSET_IMPLICIT_CURVE<T>*)stub_rigid_body->implicit_curve)->levelset;
        FILE_UTILITIES::Write_To_File<RW>(filename,levelset);

        if(stub_rigid_body->velocity_field){
            filename=STRING_UTILITIES::string_sprintf("%s/%d/velocities",output_directory.c_str(),frame);
            FILE_UTILITIES::Write_To_File<RW>(filename,*stub_rigid_body->velocity_field);}
    }
}
//#####################################################################
};
}
#endif
