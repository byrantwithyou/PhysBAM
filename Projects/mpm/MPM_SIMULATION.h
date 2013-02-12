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
    static const int influenced_node_N=pow(IN,TV::m);
public:
    MPM_PARTICLES<TV> particles;
    int N_particles;
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

    T dt;
    int frame;

    MPM_SIMULATION();
    ~MPM_SIMULATION();

    int Initialize(/*to be determined*/);
    void Advance_One_Time_Step_Forward_Euler();

protected:
    void Build_Weights_And_Grad_Weights();
    void Rasterize_Particle_Data_To_Grid();
//#####################################################################
};
}
#endif
