//#####################################################################
// Copyright 2005, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARB_PARAMETERS
//#####################################################################
#ifndef __ARB_PARAMETERS__
#define __ARB_PARAMETERS__

#include <Tools/Parsing/PARAMETER_LIST.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{
namespace ARB_PARAMETERS{

template<class T>
inline void Read_Common_Parameters(const std::string& filename,SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T,2> > >& example,PARAMETER_LIST& parameter_list)
{
    typedef VECTOR<T,2> TV;
    SOLIDS_PARAMETERS<TV>& solids_parameters=example.solids_parameters;
    ARTICULATED_RIGID_BODY<TV>& arb=example.solid_body_collection.rigid_body_collection.articulated_rigid_body;

    parameter_list.Set_Verbose(true);
    if(FILE_UTILITIES::File_Exists(filename)){
        LOG::cout << "Reading parameter file '" << filename << "'" << std::endl;
        parameter_list.Read(filename);}

    parameter_list.Get_Parameter_In_Place("frame_rate",example.frame_rate);
    parameter_list.Get_Parameter_In_Place("last_frame",example.last_frame);

    parameter_list.Get_Parameter_In_Place("solids_parameters.cfl",example.solids_parameters.cfl);
    parameter_list.Get_Parameter_In_Place("example.frame_rate",example.frame_rate);
    
    parameter_list.Get_Parameter_In_Place("arb.iterative_tolerance",arb.iterative_tolerance);
    parameter_list.Get_Parameter_In_Place("arb.max_iterations",arb.max_iterations);
    parameter_list.Get_Parameter_In_Place("arb.poststabilization_iterations",arb.poststabilization_iterations);
    parameter_list.Get_Parameter_In_Place("arb.use_epsilon_scale",arb.use_epsilon_scale);
    parameter_list.Get_Parameter_In_Place("arb.contact_level_iterations",arb.contact_level_iterations);
    parameter_list.Get_Parameter_In_Place("arb.shock_propagation_level_iterations",arb.shock_propagation_level_iterations);
    parameter_list.Get_Parameter_In_Place("arb.actuation_iterations",arb.actuation_iterations);
    parameter_list.Get_Parameter_In_Place("arb.max_line_search_iterations",arb.max_line_search_iterations);
    parameter_list.Get_Parameter_In_Place("arb.line_search_interval_tolerance",arb.line_search_interval_tolerance);
    parameter_list.Get_Parameter_In_Place("arb.use_shock_propagation",arb.use_shock_propagation);
    parameter_list.Get_Parameter_In_Place("arb.do_final_pass",arb.do_final_pass);
    parameter_list.Get_Parameter_In_Place("arb.enforce_nonnegative_activations",arb.enforce_nonnegative_activations);
    parameter_list.Get_Parameter_In_Place("arb.clamp_negative_activations",arb.clamp_negative_activations);
    parameter_list.Get_Parameter_In_Place("arb.activation_optimization_iterations",arb.activation_optimization_iterations);
    parameter_list.Get_Parameter_In_Place("arb.global_post_stabilization",arb.global_post_stabilization);
    parameter_list.Get_Parameter_In_Place("arb.verbose",arb.verbose);
    parameter_list.Get_Parameter_In_Place("arb.use_angular_damping",arb.use_angular_damping);

    if(parameter_list.Get_Parameter("arb.use_pd_actuators",false)){LOG::cout << "Using PD actuators" << std::endl;arb.Use_PD_Actuators();}
    if(parameter_list.Get_Parameter("arb.use_muscle_actuators",false)){LOG::cout << "Using muscle actuators" << std::endl;arb.Use_Muscle_Actuators();}
    if(parameter_list.Get_Parameter("arb.use_no_actuators",false)){LOG::cout << "Using no actuators" << std::endl;arb.Use_No_Actuators();}

    parameter_list.Get_Parameter_In_Place("solids_parameters.rigid_body_collision_parameters.use_epsilon_scaling",solids_parameters.rigid_body_collision_parameters.use_epsilon_scaling);
    parameter_list.Get_Parameter_In_Place("solids_parameters.rigid_body_collision_parameters.use_epsilon_scaling_for_level",solids_parameters.rigid_body_collision_parameters.use_epsilon_scaling_for_level);
    parameter_list.Get_Parameter_In_Place("solids_parameters.rigid_body_collision_parameters.contact_iterations",solids_parameters.rigid_body_collision_parameters.contact_iterations);
    parameter_list.Get_Parameter_In_Place("solids_parameters.rigid_body_collision_parameters.use_shock_propagation",solids_parameters.rigid_body_collision_parameters.use_shock_propagation);
    parameter_list.Get_Parameter_In_Place("solids_parameters.rigid_body_collision_parameters.use_push_out",example.solids_parameters.rigid_body_collision_parameters.use_push_out);
    parameter_list.Get_Parameter_In_Place("solids_parameters.rigid_body_collision_parameters.collision_iterations",solids_parameters.rigid_body_collision_parameters.collision_iterations);
    parameter_list.Get_Parameter_In_Place("solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies",solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies);
    parameter_list.Get_Parameter_In_Place("solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction",solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction);
}

}
}
#endif
