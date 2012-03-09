//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEFORMABLE_AND_RIGID
//#####################################################################
#ifndef __DEFORMABLE_AND_RIGID__
#define __DEFORMABLE_AND_RIGID__

#include <Rigid_Bodies/RIGID_BODY_IMPULSE_ACCUMULATOR_2D.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_2D.h>
namespace PhysBAM{

template <class T,class RW>
class DEFORMABLE_AND_RIGID:public SOLIDS_FLUIDS_EXAMPLE_2D<RW>
{
public:
    T initial_height;
    T initial_orientation;
    VECTOR_2D<T> initial_velocity;
    T initial_angular_velocity;
    T side_length;
    int m,n;
    GRID<TV> mattress_grid;
    bool fix_boundary;
    
    DEFORMABLE_AND_RIGID()
        :initial_height(4),initial_orientation(pi/2),
        initial_velocity(0,0),initial_angular_velocity(0),
        side_length(1),m(5),n(40),
        mattress_grid(m+1,n+1,0,side_length,0,side_length*n/m),
        fix_boundary(true)
    {
        if(fix_boundary){initial_angular_velocity=0;initial_velocity=VECTOR_2D<T>(0,0);}
        last_frame=10*24;
        restart=false;restart_frame=0;   
        solids_parameters.cfl=(T).5;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        output_directory="Deformable_And_Rigid/output";
        solids_parameters.collide_with_interior=false;
        verbose_dt=true;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    }

    ~DEFORMABLE_AND_RIGID()
    {}

//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    solids_parameters.deformable_body_parameters.list.Add_Deformable_Triangulated_Area();
    TRIANGULATED_AREA<T>& triangulated_area=*solids_parameters.deformable_body_parameters.list(1).triangulated_area;
    TRIANGLE_MESH& triangle_mesh=triangulated_area.triangle_mesh;
    DEFORMABLE_PARTICLES<T,VECTOR_2D<T> >& particles=solids_parameters.deformable_body_parameters.list(1).particles;

    triangulated_area.Initialize_Square_Mesh_And_Particles(mattress_grid);

    triangulated_area.Set_Density(1000);
    triangulated_area.Set_Mass_Of_Particles(true);

    std::cout << "DEFORMABLE_PARTICLES MASS " << particles.mass(1) << std::endl;

    triangulated_area.Update_Bounding_Box();
    VECTOR_2D<T> center(triangulated_area.bounding_box->Center());T bottom=triangulated_area.bounding_box->ymin;
    for(int i=0;i<triangulated_area.particles.array_size;i++){
        triangulated_area.particles.X(i)=center+MATRIX<T,2>::Rotation_Matrix(initial_orientation)*(triangulated_area.particles.X(i)-center);
        VECTOR_2D<T> radial = triangulated_area.particles.X(i)-center;
        T temp=radial.y;radial.y = -radial.x;radial.x=temp;
        triangulated_area.particles.V(i)=initial_velocity+initial_angular_velocity*(radial);
        triangulated_area.particles.X(i).y+=initial_height-bottom;}

    triangulated_area.Check_Signed_Areas_And_Make_Consistent(true);

    int index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies_2D/circle");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Mass(1000);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->position=VECTOR_2D<T>(0,10);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->velocity=VECTOR_2D<T>(0,0);
//    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->position=VECTOR_2D<T>(-3,5);
//    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->velocity=VECTOR_2D<T>(10,0);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction.Vector_2D(),solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

#if 0
    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("../../Public_Data/Rigid_Bodies_2D/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;
#endif

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Initialize_Bodies();
    DEFORMABLE_OBJECT_2D<T>& deformable_object=solids_parameters.deformable_body_parameters.list(1);

    solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_object.particles));
    solid_body_collection.Add_Force(Create_Diagonalized_Finite_Volume(triangulated_area,new DIAGONALIZED_SPLINE_MODEL_2D<T>(*deformable_object.triangulated_area,(T)5e4,(T).3,(T).5,7,(T).01)));
    //deformable_object.Add_Edge_Springs(deformable_object.triangulated_area->triangle_mesh,3000,(T).9);
    //deformable_object.Add_Altitude_Springs(deformable_object.triangulated_area->triangle_mesh,300);
    //deformable_object.Add_Neo_Hookean_Elasticity(*deformable_object.triangulated_area,(T)5e4);
    //deformable_object.Add_Linear_Finite_Volume(*deformable_object.triangulated_area,3e4);
}
//#####################################################################
// Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY<VECTOR_2D<T> >& V,const T time){
    switch(id_number){
    case 1:
        if(fix_boundary) for(int i=0;i<mattress_grid.m;i++) V(i)=V(i+mattress_grid.m*(mattress_grid.n-1))=VECTOR_2D<T>();
        break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
// Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY<VECTOR_2D<T> >& V,const T time){
    switch(id_number){
    case 1:
        if(fix_boundary) for(int i=0;i<mattress_grid.m;i++) V(i)=V(i+mattress_grid.m*(mattress_grid.n-1))=VECTOR_2D<T>();
        break;
    default:std::cout<<"Unrecognized deformable object id number"<<std::endl;exit(1);}
}
//#####################################################################
};
}
#endif
