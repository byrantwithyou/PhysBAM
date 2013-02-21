//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_SIMULATION
//#####################################################################
#ifndef __MPM_SIMULATION__
#define __MPM_SIMULATION__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include "MPM_CONSTITUTIVE_MODEL.h"
#include "MPM_CUBIC_B_SPLINE.h"
#include "MPM_PARTICLES.h"
namespace PhysBAM{
template<class TV>
class MPM_SIMULATION
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    static const bool PROFILING=true;
public:
    enum WORKAROUND{basis_function_order=3,IN=basis_function_order+1};
    //#################################################################
    // need external input
    //#################################################################
    MPM_PARTICLES<TV> particles; // Resize(), X, V, mass, Fe, Fp
    GRID<TV> grid;    
    bool dump_matrix,test_system;
    T min_mass,min_rho;
    T dt;
    T mu0,lambda0,ym0,pr0;
    T xi;
    bool use_plasticity_yield;
    T yield_min;
    T yield_max;
    bool use_plasticity_clamp;
    T clamp_min;
    T clamp_max;
    bool use_gravity;
    T FLIP_alpha;
    T friction_coefficient;
    T ground_level;
    ARRAY<RANGE<TV> > dirichlet_box;
    ARRAY<TV> dirichlet_velocity;
    ARRAY<SPHERE<TV> > rigid_ball;
    ARRAY<TV> rigid_ball_velocity;
    //#################################################################
    TV gravity_constant;
    ARRAY<bool> valid;
    ARRAY<T> mu,lambda;
    ARRAY<T> Je,Jp;
    ARRAY<MATRIX<T,TV::m> > Re,Se;
    ARRAY<TV_INT> influence_corner;
    ARRAY<ARRAY<T,TV_INT> > weight;
    ARRAY<ARRAY<TV,TV_INT> > grad_weight;
    MPM_CUBIC_B_SPLINE<TV,basis_function_order> grid_basis_function; 
    MPM_CONSTITUTIVE_MODEL<TV> constitutive_model;
    ARRAY<T,TV_INT> node_mass;
    ARRAY<TV,TV_INT> node_V;
    ARRAY<TV,TV_INT> node_V_star;
    ARRAY<TV,TV_INT> node_V_old;
    ARRAY<TV,TV_INT> node_force;
    ARRAY<TV,TV_INT> node_external_force;
    int frame;

    MPM_SIMULATION();
    ~MPM_SIMULATION();

    void Initialize();
    void Advance_One_Time_Step_Forward_Euler();
    void Advance_One_Time_Step_Backward_Euler();
protected:
    void Build_Weights_And_Grad_Weights();
    void Build_Helper_Structures_For_Constitutive_Model();
    void Rasterize_Particle_Data_To_The_Grid();
    void Compute_Particle_Volumes_And_Densities();
    void Compute_Grid_Forces();
    void Apply_Gravity_To_Grid_Forces();
    void Update_Velocities_On_Grid();
    void Grid_Based_Body_Collisions();
    void Solve_The_Linear_System_Explicit();
    void Solve_The_Linear_System();
    void Update_Deformation_Gradient();
    void Update_Particle_Velocities();
    void Particle_Based_Body_Collisions();
    void Update_Particle_Positions();
    void Update_Dirichlet_Box_Positions();
    void Update_Colliding_Object_Positions();
private:
    T Get_Maximum_Node_Velocity() const;
    T Get_Maximum_Particle_Velocity() const;
//#####################################################################
};
}
#endif
