//#####################################################################
// Copyright 2003-2004, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PLANK_TEST
//##################################################################### 
#ifndef __PLANK_TEST__
#define __PLANK_TEST__

#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RGB_COLORS.h>
namespace PhysBAM{

template<class T,class RW>
class PLANK_TEST:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW> BASE;
    using BASE::last_frame;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::data_directory;using BASE::fluids_parameters;
    using BASE::solids_parameters;using BASE::verbose_dt;

    int parameter;

    PLANK_TEST()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(0,fluids_parameters.NONE),parameter(0)
    {
        last_frame=480;
        output_directory="Plank_Test/output";

        solids_parameters.rigid_body_parameters.artificial_maximum_speed=20;
        solids_parameters.perform_self_collision=false;

        if (parameter == 0) {

        }
        else if (parameter == 1) {
            solids_parameters.rigid_body_parameters.use_triangle_hierarchy = true;
            solids_parameters.rigid_body_parameters.use_edge_intersection = true;
            solids_parameters.rigid_body_parameters.use_triangle_hierarchy_center_phi_test = true;
        }

        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    }

    // unused callbacks
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    RIGID_BODY_LIST<T,TV>& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    int id=0;

    T baseboxsize=2;
    T dropboxsize=1.2;
    T plankscale=1;
    T dropheight=70;
    T smallboxsize=0.5;
    T smallboxmass=1;
    T offset = 0.05;
    std::string boxfile = (parameter == 0) ? "subdivided_box" : "box";
    std::string plankfile = (parameter == 0) ? "plank" : "unsubdivided_plank";

    T stack_epsilon = 0.3;
    T stack_mu = 0.5;

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+boxfile,smallboxsize);
    rigid_body_particles.Rigid_Body(id).frame.t=TV(-0.65, baseboxsize+4, 0);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Restitution(stack_epsilon);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Friction(stack_mu);
    rigid_body_particles.Rigid_Body(id).Set_Mass(smallboxmass);
    rigid_body_particles.Rigid_Body(id).Set_Name("stack box 1a");

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+boxfile,smallboxsize);
    rigid_body_particles.Rigid_Body(id).frame.t=TV(0.65, baseboxsize+5, 0);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Restitution(stack_epsilon);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Friction(stack_mu);
    rigid_body_particles.Rigid_Body(id).Set_Mass(smallboxmass);
    rigid_body_particles.Rigid_Body(id).Set_Name("stack box 1b");

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+boxfile,smallboxsize);
    rigid_body_particles.Rigid_Body(id).frame.t=TV(offset, baseboxsize+7, 0);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Restitution(stack_epsilon);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Friction(stack_mu);
    rigid_body_particles.Rigid_Body(id).Set_Mass(smallboxmass);
    rigid_body_particles.Rigid_Body(id).Set_Name("stack box 2");

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+boxfile,smallboxsize);
    rigid_body_particles.Rigid_Body(id).frame.t=TV(0, baseboxsize+9, offset);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Restitution(stack_epsilon);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Friction(stack_mu);
    rigid_body_particles.Rigid_Body(id).Set_Mass(smallboxmass);
    rigid_body_particles.Rigid_Body(id).Set_Name("stack box 3");

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+boxfile,smallboxsize);
    rigid_body_particles.Rigid_Body(id).frame.t=TV(-offset, baseboxsize+14, -offset);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Restitution(stack_epsilon);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Friction(stack_mu);
    rigid_body_particles.Rigid_Body(id).Set_Mass(smallboxmass);
    rigid_body_particles.Rigid_Body(id).Set_Name("stack box 4");

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+boxfile,smallboxsize);
    rigid_body_particles.Rigid_Body(id).frame.t=TV(0, baseboxsize+17, 0);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Restitution(stack_epsilon);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Friction(stack_mu);
    rigid_body_particles.Rigid_Body(id).Set_Mass(smallboxmass);
    rigid_body_particles.Rigid_Body(id).Set_Name("stack box 5");

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+boxfile,smallboxsize);
    rigid_body_particles.Rigid_Body(id).frame.t=TV(0, baseboxsize+19, 0);
    rigid_body_particles.Rigid_Body(id).frame.r=QUATERNION<T>(pi/4,TV(0,1,0));
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Restitution(stack_epsilon);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Friction(stack_mu);
    rigid_body_particles.Rigid_Body(id).Set_Mass(smallboxmass);
    rigid_body_particles.Rigid_Body(id).Set_Name("stack box 6");

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+boxfile);
    rigid_body_particles.Rigid_Body(id).frame.t=TV(baseboxsize+plankscale*5-dropboxsize,2*baseboxsize+0.5*plankscale+dropboxsize+dropheight,0);
    rigid_body_particles.Rigid_Body(id).frame.r=QUATERNION<T>(pi/2,TV(0,1,0));
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Restitution(0.05);
    rigid_body_particles.Rigid_Body(id).Set_Mass(100);
    rigid_body_particles.Rigid_Body(id).Set_Name("drop box");

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+plankfile,plankscale);
    rigid_body_particles.Rigid_Body(id).frame.t=TV(baseboxsize,2*baseboxsize+0.25,0);
    rigid_body_particles.Rigid_Body(id).frame.r=QUATERNION<T>(pi/2,TV(0,1,0));
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Restitution(0.1);
    rigid_body_particles.Rigid_Body(id).Set_Mass(10);
    rigid_body_particles.Rigid_Body(id).Set_Name("plank");

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+boxfile,baseboxsize);
    rigid_body_particles.Rigid_Body(id).frame.t=TV(0,baseboxsize,0);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Restitution(0.1);
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Friction(0.3);
    rigid_body_particles.Rigid_Body(id).Set_Name("base box");
    rigid_body_particles.Rigid_Body(id).is_static = true;

    id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground");
    rigid_body_particles.Rigid_Body(id).Set_Name("ground");
    rigid_body_particles.Rigid_Body(id).is_static=true;
    rigid_body_particles.Rigid_Body(id).add_to_spatial_partition=false;
    rigid_body_particles.Rigid_Body(id).Set_Coefficient_Of_Restitution(1);

    for(int i=0;i<rigid_body_list.Number_Of_Elements();i++) if(!rigid_body_particles.Rigid_Body(i).is_static)
        rigid_body_particles.Rigid_Body(i).Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Initialize_Bodies();
}
//#####################################################################
};
}
#endif
