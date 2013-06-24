//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOLIDS_PARAMETERS
//#####################################################################
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Vectors/VECTOR_3D.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <climits>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOLIDS_PARAMETERS<TV>::
SOLIDS_PARAMETERS()
    :triangle_collision_parameters(*new TRIANGLE_COLLISION_PARAMETERS<TV>),implicit_solve_parameters(*new IMPLICIT_SOLVE_PARAMETERS<TV>),
    rigid_body_collision_parameters(*new RIGID_BODY_COLLISION_PARAMETERS<TV>),rigid_body_evolution_parameters(*new RIGID_BODY_EVOLUTION_PARAMETERS<TV>),
    deformable_object_collision_parameters(*new DEFORMABLE_OBJECT_COLLISION_PARAMETERS<TV>),write_deformable_body(true),verbose(true),verbose_dt(false),cfl((T)10),min_dt((T)0),
    fracture(false),write_static_variables_every_frame(false),newton_tolerance((T)1e-3),newton_iterations(1),use_partially_converged_result(true),
    write_from_every_process(true),enforce_repulsions_in_cg(true),use_post_cg_constraints(true),use_rigid_deformable_contact(false),rigid_cluster_fracture_frequency(INT_MAX),
    use_trapezoidal_rule_for_velocities(true),enforce_poststabilization_in_cg(true),
    no_contact_friction(false),use_projections_in_position_update(false),allow_altitude_spring_change_between_updates(true)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOLIDS_PARAMETERS<TV>::
~SOLIDS_PARAMETERS()
{
    delete &triangle_collision_parameters;
    delete &implicit_solve_parameters;
    delete &rigid_body_collision_parameters;
    delete &rigid_body_evolution_parameters;
    delete &deformable_object_collision_parameters;
}
//#####################################################################
template class SOLIDS_PARAMETERS<VECTOR<float,1> >;
template class SOLIDS_PARAMETERS<VECTOR<float,2> >;
template class SOLIDS_PARAMETERS<VECTOR<float,3> >;
template class SOLIDS_PARAMETERS<VECTOR<double,1> >;
template class SOLIDS_PARAMETERS<VECTOR<double,2> >;
template class SOLIDS_PARAMETERS<VECTOR<double,3> >;
}
