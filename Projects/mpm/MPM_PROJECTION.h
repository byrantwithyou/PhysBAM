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
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
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

    GRID<TV> mac_grid;
    ARRAY<bool,TV_INT> cell_dirichlet;
    ARRAY<bool,TV_INT> cell_neumann;
    ARRAY<int,TV_INT> neumann_cell_normal_axis; // +-1 +-2 +-3
    HASHTABLE<TV_INT,bool> nodes_non_dirichlet_cells;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities_old;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_masses;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_volumes;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_densities;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_momenta;
    ARRAY<T,TV_INT> div_u;
    T max_div;
    ARRAY<T,TV_INT> pressure_unknown;
    ARRAY<T,TV_INT> pressure_rasterized;

    MPM_PROJECTION(MPM_SIMULATION<TV>& sim_in);
    ~MPM_PROJECTION();

    void Reinitialize();                                   
    void Identify_Dirichlet_Cells();                       
    void Identify_Neumann_Cells();                         
    void Identify_Nodes_Of_Non_Dirichlet_Cells();         
    void Velocities_Corners_To_Faces_MPM_Style();          
    void Build_Velocity_Divergence();                      
    void Solve_For_Pressure();                            
    void Do_Projection();                                  
    void Velocities_Faces_To_Corners_MPM_Style(T FLIP_alpha=T(0));         

    void Fix_RHS_Neumann_Cells(ARRAY<T,TV_INT>& rhs);      // called by Solve_For_Pressure
    
    TV Get_Total_Momentum_On_Faces() const;
};
}
#endif
