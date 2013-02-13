//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_SIMULATION
//#####################################################################
#ifndef __MPM_SIMULATION__
#define __MPM_SIMULATION__
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX.h>
#include "MPM_PARTICLES.h"
#include "MPM_CONSTITUTIVE_MODEL.h"
#include "MPM_CUBIC_B_SPLINE.h"

namespace PhysBAM{

template<class TV>
class MPM_SIMULATION
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    static const int basis_function_order=3;
    static const int IN=basis_function_order+1;
    int IN_POWER;
public:
    MPM_PARTICLES<TV> particles;
    int N_particles;
    T mu0,lambda0;
    T xi;
    bool use_plasticity;
    T theta_c;
    T theta_s;
    T FLIP_alpha;

    ARRAY<T> mu,lambda;
    ARRAY<T> Je,Jp;
    ARRAY<MATRIX<T,TV::m> > Ue,Ve,Re,Se;
    ARRAY<DIAGONAL_MATRIX<T,TV::m> > SIGMAe;
    ARRAY<TV_INT> influence_corner;
    ARRAY<ARRAY<T> > weight;
    ARRAY<ARRAY<TV> > grad_weight;
    MPM_CUBIC_B_SPLINE<TV,basis_function_order> grid_basis_function; 
    MPM_CONSTITUTIVE_MODEL<TV> constitutive_model;

    GRID<TV> grid;
    int N_nodes;
    ARRAY<T> node_mass;
    ARRAY<TV> node_V;
    ARRAY<TV> node_V_star;
    ARRAY<TV> node_V_old; // used for FLIP
    ARRAY<TV> node_force;

    T dt;
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
    void Update_Velocities_On_Grid();
    void Grid_Based_Body_Collisions();
    void Solve_The_Linear_System();
    void Update_Deformation_Gradient();
    void Update_Particle_Velocities();
    void Particle_Based_Body_Collisions();
    void Update_Particle_Positions();
//#####################################################################
};
}
#endif
