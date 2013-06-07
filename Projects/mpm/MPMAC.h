//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPMAC
//#####################################################################
#ifndef __MPMAC__
#define __MPMAC__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include "MPM_PARTICLES.h"
#include "MPM_LINEAR_BASIS.h"
namespace PhysBAM{
template<class TV>
class MPMAC
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    enum WORKAROUND{basis_function_order=1,IN=basis_function_order+1};

    // set externally
    MPM_PARTICLES<TV> particles;
    GRID<TV> grid;
    T dt;
    bool use_gravity;
    T FLIP_alpha;

    T min_mass;
    int frame;
    GRID<TV> mac_grid;
    GRID<TV> face_grid[TV::m];
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_masses;
    ARRAY<T,FACE_INDEX<TV::dimension> > face_momenta;
    ARRAY<bool,TV_INT> cell_dirichlet;
    ARRAY<bool,TV_INT> cell_neumann;
    ARRAY<T,TV_INT> div_u;
    ARRAY<T,TV_INT> pressure;
    MPM_LINEAR_BASIS<TV,basis_function_order> grid_basis_function;

    ARRAY<TV_INT> influence_corner[TV::m];
    ARRAY<ARRAY<T,TV_INT> > weight[TV::m];

    MPMAC(){};
    ~MPMAC(){};
    
    void Initialize();
    void Reinitialize();
    void Weights();
    void Rasterize();
//#####################################################################
};
}
#endif
