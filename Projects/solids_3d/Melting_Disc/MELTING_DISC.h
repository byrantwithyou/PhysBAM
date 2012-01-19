//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MELTING_DISC
//#####################################################################
#ifndef __MELTING_DISC__
#define __MELTING_DISC__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <Level_Sets/LEVELSET_TRIANGULATED_OBJECT.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_2D.h>
#include <Solids_And_Fluids/MELTING_EXAMPLE_S3D.h>
namespace PhysBAM{

template <class T,class RW>
class MELTING_DISC:public MELTING_EXAMPLE_S3D<T,RW>
{
public:
    typedef MELTING_EXAMPLE_S3D<T,RW> BASE;
    using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::data_directory;using BASE::verbose_dt;
    using BASE::fluids_parameters;using BASE::solids_parameters;using MELTING_EXAMPLE_S3D<T,RW>::melting_parameters;

    CIRCLE<T> circle;
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity;
    VECTOR_3D<T> initial_angular_velocity;
    T side_length;
    int m,n;
    
    MELTING_DISC()
        :MELTING_EXAMPLE_S3D<T,RW>(fluids_parameters.NONE),initial_height(0),initial_orientation(),initial_velocity(),initial_angular_velocity(),//.5,1,0),
        side_length(3),m(5),n(5)
    {
        frame_rate=24;
        last_frame=10*(int)frame_rate;
        restart=false;restart_frame=205;
        solids_parameters.cfl=(T).5;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        melting_parameters.maximum_depth=4;

        melting_parameters.refine_near_interface=true;
        //melting_parameters.successive_refinement_multiplier=(T)1.5;
        //melting_parameters.expansive_refinement_threshold=(T).4;
        //melting_parameters.compressive_refinement_threshold=(T).4;

        melting_parameters.write_overlay_levelset=false;        
        solids_parameters.perform_self_collision=false;
        
        output_directory="Melting_Disc/output";
        solids_parameters.collide_with_interior=false;
        verbose_dt=true;
    }

    ~MELTING_DISC()
    {}

//#####################################################################
// Function Initialize_Deformable_And_Rigid_Bodies
//#####################################################################
void Initialize_Deformable_And_Rigid_Bodies()
{
    int index=solids_parameters.deformable_body_parameters.list.Add_Deformable_Embedded_Triangulated_Surface();
    Add_Melting_Object(melting_parameters.DEFORMABLE,index);

    index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("../../Public_Data/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T)3);
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
    circle.center=grid.uniform_grid.Domain().Center();circle.radius=.9;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi(const int object,ARRAY<T>& phi)
{
    assert(object==1);
    RED_GREEN_GRID_2D<T>& grid=melting_parameters.levelsets(1)->grid;
    ARRAY<VECTOR_2D<T> >& node_locations=grid.Node_Locations();
    
    for(int p=0;p<phi.m;p++) phi(p)=circle.Phi(node_locations(p));    
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
    TRIANGULATED_SURFACE<T>& triangulated_surface=solids_parameters.deformable_body_parameters.list(object).embedded_triangulated_surface->triangulated_surface;
    PARTICLES<T,VECTOR_3D<T> >& particles=solids_parameters.deformable_body_parameters.list(object).particles;

    particles.Update_Velocity();
    triangulated_surface.Update_Bounding_Box();
    VECTOR_3D<T> center(triangulated_surface.bounding_box->Center());
    for(int i=1;i<=particles.array_collection->Size();i++){
        particles.X(i)=center+initial_orientation.Rotate(particles.X(i)-center);
        particles.V(i)=initial_velocity+VECTOR_3D<T>::Cross_Product(initial_angular_velocity,particles.X(i)-center);
        particles.X(i).y+=initial_height;}
}
//#####################################################################
// Function Initialize_Forces
//#####################################################################
void Initialize_Forces()
{
    TRIANGULATED_SURFACE<T>& triangulated_surface=solids_parameters.deformable_body_parameters.list(1).embedded_triangulated_surface->triangulated_surface;

    solids_parameters.deformable_body_parameters.list(1).particles.Store_Mass();
    triangulated_surface.Set_Density(1000);
    triangulated_surface.Set_Mass_Of_Particles(false);

    solids_parameters.deformable_body_parameters.list(1).Delete_Forces();
    solids_parameters.deformable_body_parameters.list(1).Add_Body_Forces(triangulated_surface,solids_parameters.gravity,solids_parameters.gravity_direction);
    solids_parameters.deformable_body_parameters.list(1).Add_Diagonalized_Linear_Finite_Volume(triangulated_surface,(T)5e5,(T).3,(T).01);
    solids_parameters.deformable_body_parameters.list(1).Add_Bending_Elements(triangulated_surface.triangle_mesh,(T)1e3,(T)1e2);
}
//#####################################################################
// Function Postprocess_Solids_Substep
//#####################################################################
void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE
{
    //Update_Solids_Topology_For_Melting(1/frame_rate,time); // ridiculous hack to test if things work
}
//#####################################################################
};
}
#endif

