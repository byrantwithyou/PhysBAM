//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MELTING_RIGID_SPHERE
//#####################################################################
#ifndef __MELTING_RIGID_SPHERE__
#define __MELTING_RIGID_SPHERE__

#include <PhysBAM_Tools/Matrices/QUATERNION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Level_Sets/LEVELSET_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_3D.h>
#include <Solids_And_Fluids/MELTING_EXAMPLE_3D.h>
namespace PhysBAM{

template <class T,class RW>
class MELTING_RIGID_SPHERE:public MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::frame_rate;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::restart_frame;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::data_directory;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::verbose_dt;using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::fluids_parameters;
    using SOLIDS_FLUIDS_EXAMPLE_3D<RW>::solids_parameters;using MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >::melting_parameters;

    SPHERE<T> sphere;
    T initial_height;
    QUATERNION<T> initial_orientation;
    VECTOR_3D<T> initial_velocity;
    VECTOR_3D<T> initial_angular_velocity;
    T side_length;
    int m,n,mn;
    
    MELTING_RIGID_SPHERE()
        :MELTING_EXAMPLE_3D<T,RW,SOLIDS_FLUIDS_EXAMPLE_3D<RW> >(fluids_parameters.NONE),initial_height(0),initial_orientation(),initial_velocity(0,0,0),initial_angular_velocity(0,0,0),
        side_length(2),m(2),n(2),mn(2)
    {
        frame_rate=24;
        last_frame=10*(int)frame_rate;
        last_frame=200;
        restart=false;restart_frame=205;
        solids_parameters.cfl=(T).5;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        melting_parameters.maximum_depth=4;

        melting_parameters.refine_near_interface=true;
        //successive_refinement_multiplier=(T)1.5;
        //expansive_refinement_threshold=(T).4;
        //compressive_refinement_threshold=(T).4;

        melting_parameters.write_overlay_levelset=false;
        solids_parameters.perform_self_collision=false;
        
        output_directory="Melting_Rigid_Sphere/output";
        solids_parameters.collide_with_interior=false;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        verbose_dt=true;
    }

    ~MELTING_RIGID_SPHERE()
    {}

//#####################################################################
// Function Initialize_Deformable_And_Rigid_Bodies
//#####################################################################
void Initialize_Deformable_And_Rigid_Bodies()
{
    Add_Melting_Object(melting_parameters.RIGID,0);
    Add_Melting_Object(melting_parameters.RIGID,0);

    int index=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>("/n/bcc/data/losasso/PhysBAM-Full/Public_Data/Rigid_Bodies/ground");
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->Set_Coefficient_Of_Friction((T).3);
    solids_parameters.rigid_body_parameters.list.rigid_bodies(index)->is_static=true;

//    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
}
//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Grid(const int object,RED_GREEN_GRID_3D<T>& grid)
{
    T y=.5;
    grid.Initialize(GRID_3D<T>(4*m+1,4*n+1,4*mn+1,0,side_length,y,y+side_length*n/m,0,side_length*mn/m),melting_parameters.maximum_depth);
    sphere.center=grid.uniform_grid.Domain().Center();
    sphere.radius=.9;
}
//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi(const int object,ARRAY<T>& phi)
{
    RED_GREEN_GRID_3D<T>& grid=melting_parameters.levelsets(object)->grid;
    ARRAY<VECTOR_3D<T> >& node_locations=grid.Node_Locations();
    
    for(int p=1;p<=phi.m;p++) phi(p)=sphere.Phi(node_locations(p));    
    //ARRAY<T>::copy(-1,phi);
}
//#####################################################################
// Function Initialize_Levelset_Velocity
//#####################################################################
void Initialize_Levelset_Velocity(const int object,ARRAY<VECTOR_3D<T> >& V)
{
    RED_GREEN_GRID_3D<T>& grid=melting_parameters.levelsets(object)->grid;
    ARRAY<VECTOR_3D<T> >& node_locations=grid.Node_Locations();
    
//    melting_parameters.use_constant_melting_speed=true;
//    melting_parameters.constant_melting_speed=.1;
    
    //for(int p=1;p<=V.m;p++) V(p)=VECTOR_3D<T>(node_locations(p).y-sphere.center.y,0,0);
    for(int p=1;p<=V.m;p++) V(p)=VECTOR_3D<T>();//VECTOR_3D<T>(1*max(1-2*fabs(node_locations(p).y-sphere.center.y),0.),0,0);
}
//#####################################################################
// Function Initialize_Particle_Positions_And_Velocities
//#####################################################################
void Initialize_Particle_Positions_And_Velocities(const int object)
{
    RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(melting_parameters.body_index(object));
    rigid_body.position.y+=initial_height+(object-1)*2.1;
    rigid_body.velocity=initial_velocity;
    rigid_body.angular_velocity=initial_angular_velocity;
}
//#####################################################################
// Function Initialize_Forces
//#####################################################################
void Initialize_Forces()
{
    for(int object=1;object<=melting_parameters.body_index.m;object++){
        int index=melting_parameters.body_index(object);if(!index)return;
        RIGID_BODY<TV>& rigid_body=*solids_parameters.rigid_body_parameters.list(index);
        rigid_body.Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);}
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    Update_Solids_Topology_For_Melting(.5/frame_rate,time); // ridiculous hack to test if things work
}
//#####################################################################
};
}
#endif

