//#####################################################################
// Copyright 2013
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_PROJECTION
// This scheme always assume the outmost boundary of the grid
// is absolutely dirichlet. i.e., the grid domain is larger enough.
//#####################################################################
#ifndef __MPM_PROJECTION__
#define __MPM_PROJECTION__
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_ITERATOR.h>
#include "MPM_SIMULATION.h"
namespace PhysBAM{
template<class TV>
class MPM_PROJECTION
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    MPM_SIMULATION<TV>& sim;

    bool independent_from_MPM;
    ARRAY<TV> particle_X;
    ARRAY<TV> particle_V;

    GRID<TV> mac_grid; // One cell wider than MPM grid so that MPM grid v lives on cell centers of mac_grid.
    ARRAY<bool,TV_INT> cell_dirichlet;
    ARRAY<bool,TV_INT> cell_neumann;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ARRAY<T,TV_INT> div_u;
    ARRAY<T,TV_INT> pressure;

    MPM_PROJECTION(MPM_SIMULATION<TV>& sim_in,bool independent=true);
    ~MPM_PROJECTION();

    void Reinitialize();
    void Load_New_Particle_Data(ARRAY<TV>& particle_X_in,ARRAY<TV>& particle_V_in);
    void Identify_Dirichlet_Cells();
    void Identify_Neumann_Cells();
    void Generate_Face_Velocities();
    void Build_Velocity_Divergence();
    void Solve_For_Pressure(const T dt,const T rho);
    void Do_Projection(const T dt,const T rho);
    void Send_Velocities_Back_To_MPM_Grid();
};
}
#endif
