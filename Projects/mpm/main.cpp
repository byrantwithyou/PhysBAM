//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Tools/Read_Write/TYPED_STREAM.h>
#include <PhysBAM_Tools/Read_Write/READ_WRITE_FORWARD.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include "MPM_SIMULATION.h"
using namespace PhysBAM;

int main(int argc,char *argv[])
{
    static const int dimension=2;
    typedef double T;
    typedef float RW;
    typedef VECTOR<T,dimension> TV;
    typedef VECTOR<int,dimension> TV_INT;

    MPM_SIMULATION<TV> sim;

    static const int test=1;
    if(test==1){ // stretched cube
        TV_INT grid_counts(7,7);
        RANGE<TV> grid_box(TV(-3,-3),TV(3,3));
        TV_INT particle_counts(7,4);
        RANGE<TV> particle_box(TV(-0.5,-0.5),TV(0.5,0.5));
        GRID<TV> grid(grid_counts,grid_box);
        sim.particles.Resize(particle_counts.Product());
        sim.particles.Initialize_X_As_A_Grid(particle_counts,particle_box);
        sim.grid=grid;
        sim.dt=1e-3;
        T ym=3000;
        T pr=0.3;
        sim.mu0=ym;
        sim.lambda0=pr;
        
    }

    // sim.Initialize();
    VIEWER_OUTPUT<TV> vo(STREAM_TYPE((RW)0),sim.grid,"output");
    for(int i=0;i<sim.particles.X.m;i++)
        Add_Debug_Particle(sim.particles.X(i),VECTOR<T,3>(0,1,0));
    Flush_Frame<TV>("mpm");
    // Add_Debug_Particle(TV(1,1), VECTOR<T,3>(0,1,0));
    // Flush_Frame<TV>("mpm");    
    // Add_Debug_Particle(TV(2,2), VECTOR<T,3>(0,1,0));
    // Flush_Frame<TV>("mpm");    
    // Add_Debug_Particle(TV(3,3), VECTOR<T,3>(0,1,0));
    // Add_Debug_Particle(TV(3,4), VECTOR<T,3>(0,1,0));
    // Flush_Frame<TV>("mpm");    

    return 0;
}
