//#####################################################################
// Copyright 2003-2004, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNCE_EXAMPLE
//##################################################################### 
#ifndef __BOUNCE_EXAMPLE__
#define __BOUNCE_EXAMPLE__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RGB_COLORS.h>
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_3D.h>
namespace PhysBAM{

template<class T,class RW>
class BOUNCE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_3D<RW>
{
public:
    BOUNCE_EXAMPLE()
        :SOLIDS_FLUIDS_EXAMPLE_3D<RW>(fluids_parameters.NONE)
    {
        last_frame=240;
        output_directory="Bounce_Example/output";
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    }

    ~BOUNCE_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies()
{
    //VECTOR_3D<double> ground_color=RGB_COLORS<T>().Ground_Tan(),color1=RGB_COLORS<T>().Light_Steel_Blue();

    RIGID_BODY_LIST_3D<T>& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    int id=0;

    // SPHERE
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere");
    rigid_body_particles.Rigid_Body(id).position=VECTOR_3D<T>(-3,5,0);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Restitution(1.0);
    rigid_body_particles.Rigid_Body(id).Set_Name("sphere (coeff 1)");

    // SPHERE
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere");
    rigid_body_particles.Rigid_Body(id).position=VECTOR_3D<T>(0,5,0);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Restitution(0.5);
    rigid_body_particles.Rigid_Body(id).Set_Name("sphere (coeff 0.5)");

    // SPHERE
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/sphere");
    rigid_body_particles.Rigid_Body(id).position=VECTOR_3D<T>(3,5,0);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Restitution(0.0);
    rigid_body_particles.Rigid_Body(id).Set_Name("sphere (coeff 0)");

    // PLANE -- no gravity
    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground");
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Restitution(1.0);
    rigid_body_particles.Rigid_Body(id).Set_Name("ground");
    rigid_body_particles.Rigid_Body(id).is_static=true;
    rigid_body_particles.Rigid_Body(id).add_to_spatial_partition=false;

    for(int i=1;i<=rigid_body_list.Number_Of_Elements();i++) if(!rigid_body_particles.Rigid_Body(i).is_static)
        rigid_body_particles.Rigid_Body(i).Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    SOLIDS_FLUIDS_EXAMPLE_3D<RW>::Initialize_Bodies();
}
//#####################################################################
};
}
#endif
