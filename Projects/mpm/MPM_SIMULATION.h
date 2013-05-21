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
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include "MPM_CONSTITUTIVE_MODEL.h"
#include "MPM_LINEAR_BASIS.h"
#include "MPM_CUBIC_B_SPLINE.h"
#include "MPM_PARTICLES.h"
namespace PhysBAM{
template<class TV>
class MPM_SIMULATION
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    // enum WORKAROUND{basis_function_order=1,IN=basis_function_order+1};
    enum WORKAROUND{basis_function_order=3,IN=basis_function_order+1};

    //#################################################################
    // need external input
    //#################################################################
    bool PROFILING;
    MPM_PARTICLES<TV> particles; // Resize(), X, V, mass, Fe, Fp
    GRID<TV> grid;    
    bool dump_matrix,test_system;
    T min_mass,min_rho;
    T dt;
    bool use_gravity;
    T FLIP_alpha;
    T friction_coefficient;
    ARRAY<RANGE<TV> > dirichlet_box;
    ARRAY<TV> dirichlet_velocity;
    ARRAY<SPHERE<TV> > rigid_ball;
    ARRAY<TV> rigid_ball_velocity;
    T xi;
    //#################################################################
    TV gravity_constant;
    ARRAY<T> Je;
    ARRAY<MATRIX<T,TV::m> > Re,Se;
    ARRAY<TV_INT> influence_corner;
    ARRAY<ARRAY<T,TV_INT> > weight;
    ARRAY<ARRAY<TV,TV_INT> > grad_weight;
    // MPM_LINEAR_BASIS<TV,basis_function_order> grid_basis_function; 
    MPM_CUBIC_B_SPLINE<TV,basis_function_order> grid_basis_function; 
    MPM_CONSTITUTIVE_MODEL<TV> constitutive_model;
    ARRAY<T,TV_INT> node_mass;
    ARRAY<TV,TV_INT> node_V;
    ARRAY<TV,TV_INT> node_V_star;
    ARRAY<TV,TV_INT> node_V_old;
    ARRAY<TV,TV_INT> node_force;
    ARRAY<TV_INT> nodes_need_projection;
    int frame;

    MPM_SIMULATION();
    ~MPM_SIMULATION();

    void Initialize();
    void Advance_One_Time_Step_Forward_Euler();
    void Advance_One_Time_Step_Backward_Euler();

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

    T Get_Maximum_Node_Velocity() const;
    T Get_Maximum_Particle_Velocity() const;
    TV Get_Total_Momentum_On_Nodes() const;
    TV Get_Total_Momentum_On_Particles() const;
//#####################################################################
};
}
#endif
