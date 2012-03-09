//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MELTING_CIRCLE
//#####################################################################
#ifndef __MELTING_CIRCLE__
#define __MELTING_CIRCLE__

#include <Geometry/EMBEDDED_TRIANGULATED_AREA.h>
#include <Level_Sets/LEVELSET_TRIANGULATED_OBJECT.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_3D.h>
#include <Solids_And_Fluids/MELTING_EXAMPLE_2D.h>
namespace PhysBAM{

template <class T,class RW>
class MELTING_CIRCLE:public MELTING_EXAMPLE_2D<T,RW>
{
public:
    typedef MELTING_EXAMPLE_2D<T,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::frame_rate;using BASE::restart_frame;using BASE::output_directory;using BASE::data_directory;using BASE::fluids_parameters;
    using BASE::solids_parameters;using BASE::verbose_dt;using MELTING_EXAMPLE_2D<T,RW>::melting_parameters;

    CIRCLE<T> circle;
    T initial_height;
    T initial_orientation;
    VECTOR_2D<T> initial_velocity;
    T initial_angular_velocity;
    T side_length;
    int m,n;
    
    MELTING_CIRCLE()
        :MELTING_EXAMPLE_2D<T,RW>(0,fluids_parameters.NONE),initial_height(0),initial_orientation(0),initial_velocity(0,0),initial_angular_velocity(1),
        side_length(3),m(3),n(3)
    {
        frame_rate=24*10;
        last_frame=10*(int)frame_rate;
        restart=false;restart_frame=205;
        solids_parameters.cfl=(T).5;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        melting_parameters.maximum_depth=3;

        melting_parameters.refine_near_interface=true;
        //melting_parameters.successive_refinement_multiplier=(T)1.5;
        //melting_parameters.expansive_refinement_threshold=(T).4;
        //melting_parameters.compressive_refinement_threshold=(T).4;

        melting_parameters.write_overlay_levelset=false;        
        
        output_directory="Melting_Circle/output";
        solids_parameters.collide_with_interior=false;
        verbose_dt=true;
    }

    ~MELTING_CIRCLE()
    {}

//#####################################################################
// Function Initialize_Deformable_And_Rigid_Bodies
//#####################################################################
void Initialize_Deformable_And_Rigid_Bodies()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Embedded_Triangulated_Area();
    Add_Melting_Object(melting_parameters.DEFORMABLE,index);

    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("../../Public_Data/Rigid_Bodies_2D/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Grid(const int object,RED_GREEN_GRID_2D<T>& grid)
{
    assert(object==1);T y=2;
    grid.Initialize(GRID_2D<T>(4*m+1,4*n+1,0,side_length,y,y+side_length*n/m),melting_parameters.maximum_depth);
    circle.center=grid.uniform_grid.Domain().Center();    
    circle.radius=1.2;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi(const int object,ARRAY<T>& phi)
{
    assert(object==1);
    RED_GREEN_GRID_2D<T>& grid=melting_parameters.levelsets(1)->grid;
    ARRAY<VECTOR_2D<T> >& node_locations=grid.Node_Locations();
    
    for(int p=0;p<phi.m;p++) phi(p)=circle.Signed_Distance(node_locations(p));    
}
//#####################################################################
// Function Initialize_Levelset_Velocity
//#####################################################################
void Initialize_Levelset_Velocity(const int object,ARRAY<VECTOR_2D<T> >& V)
{
    assert(object==1);
    RED_GREEN_GRID_2D<T>& grid=melting_parameters.levelsets(1)->grid;
    ARRAY<VECTOR_2D<T> >& node_locations=grid.Node_Locations();
    
    //for(int p=0;p<V.m;p++) V(p)=VECTOR_2D<T>(1,0);    
    for(int p=0;p<V.m;p++) V(p)=VECTOR_2D<T>(node_locations(p).y-circle.center.y,0);
}
//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Particle_Positions_And_Velocities(const int object)
{
    TRIANGULATED_AREA<T>& triangulated_area=solids_parameters.deformable_body_parameters.list(object).embedded_triangulated_area->triangulated_area;
    DEFORMABLE_PARTICLES<T,VECTOR_2D<T> >& particles=solids_parameters.deformable_body_parameters.list(object).particles;

    particles.Update_Velocity();
    triangulated_area.Update_Bounding_Box();
    VECTOR_2D<T> center(triangulated_area.bounding_box->Center());
    for(int i=0;i<particles.array_collection->Size();i++){
        particles.X(i)=center+MATRIX<T,2>::Rotation_Matrix(initial_orientation)*(particles.X(i)-center);
        particles.V(i)=initial_velocity+initial_angular_velocity*(particles.X(i)-center).Rotate_Counterclockwise_90();
        particles.X(i).y+=initial_height;}
}
//#####################################################################
// Function Initialize_Forces
//#####################################################################
void Initialize_Forces()
{
    TRIANGULATED_AREA<T>& triangulated_area=solids_parameters.deformable_body_parameters.list(1).embedded_triangulated_area->triangulated_area;

    solids_parameters.deformable_body_parameters.list(1).particles.Store_Mass();
    triangulated_area.Set_Density(1000);
    triangulated_area.Set_Mass_Of_Particles(false);

    solids_parameters.deformable_body_parameters.list(1).Delete_Forces();
    solids_parameters.deformable_body_parameters.list(1).Add_Force(Create_Body_Forces<T>(triangulated_area,solids_parameters.gravity,solids_parameters.gravity_direction));
    solids_parameters.deformable_body_parameters.list(1).Add_Force(Create_Diagonalized_Finite_Volume(triangulated_area,new DIAGONALIZED_SPLINE_MODEL_2D<T>((T)5e4,(T).3,(T).5,7,(T).01)));
}
//#####################################################################
// Function Postprocess_Solids_Substep
//#####################################################################
void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE
{
    Update_Solids_Topology_For_Melting(1/frame_rate,time); // ridiculous hack to test if things work
}
//#####################################################################
};
}
#endif

