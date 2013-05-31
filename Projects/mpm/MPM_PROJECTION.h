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
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include "MPM_SIMULATION.h"
namespace PhysBAM{
template<class TV>
class MPM_PROJECTION
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    MPM_SIMULATION<TV>& sim;

    GRID<TV> mac_grid;
    ARRAY<bool,TV_INT> cell_dirichlet;
    ARRAY<bool,TV_INT> cell_neumann;
    HASHTABLE<TV_INT,bool> nodes_non_dirichlet_cells;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_masses;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_momenta;
    ARRAY<T,TV_INT> div_u;
    ARRAY<T,TV_INT> pressure;

    MPM_PROJECTION(MPM_SIMULATION<TV>& sim_in);
    ~MPM_PROJECTION();

    void Reinitialize();                                   // step 1
    void Identify_Dirichlet_Cells();                       // step 2 
    void Identify_Neumann_Cells();                         // step 3
    void Identify_Nodes_Of_Non_Dirichlet_Cells();          // step 4
    void Velocities_Corners_To_Faces_MPM_Style();          // step 5
    void Build_Velocity_Divergence();                      // step 6 
    void Solve_For_Pressure();                             // step 7
    void Do_Projection();                                  // step 8
    void Velocities_Faces_To_Corners_MPM_Style();          // step 9

    void Fix_RHS_Neumann_Cells(ARRAY<T,TV_INT>& rhs);      // called by Solve_For_Pressure
    
    TV Get_Total_Momentum_On_Faces() const;
};
}
#endif
