//#####################################################################
// Copyright 2013, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef __MPLE_DRIVER__
#define __MPLE_DRIVER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include "MPLE_POINT.h"

namespace PhysBAM{

template<class TV,int w>
class MPLE_DRIVER: public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    
    enum WORKAROUND{ghost=3};

public:
    
    ARRAY<MPLE_POINT<TV,w> > points;   // data
    GRID<TV> grid;                     // grid

    ARRAY<T,TV_INT>* u;                // segmentation function
    ARRAY<T,TV_INT>* u_new;            // new segmentations funtion

    int frames;
    int timesteps;
    T dt,mu,nu;
    
    MPLE_DRIVER():frames(10),timesteps(10),dt(.1),mu(1),nu(.1)
    {
        u=new ARRAY<T,TV_INT>;
        u_new=new ARRAY<T,TV_INT>;
    }

    ~MPLE_DRIVER()
    {
        delete u;
        delete u_new;
    }

    void Initialize()
    {
        u->Resize(grid.Node_Indices(ghost),false);
        u_new->Resize(grid.Node_Indices(ghost),false);

        for(int i=0;i<points.m;i++)
            points(i).Update_Base_And_Weights(grid);
    }

    void Write(const char* title)
    {
        for(int i=0;i<points.m;i++){
            Add_Debug_Particle<TV>(points(i).X,VECTOR<T,3>(1,0,0));
            LOG::cout<<"X="<<points(i).X<<std::endl;}
        Flush_Frame<TV>(title);
    }

    void Diffusion_Timestep()
    {

    }

    void Threshold_Segmentation_Function()
    {

    }

    void Advance_Frame()
    {
        for(int i=0;i<timesteps;i++) Diffusion_Timestep();
        Threshold_Segmentation_Function();
    }

    void Run()
    {
        Initialize();
        for(int k=0;k<frames;k++){
            Advance_Frame();
            char buff[100];
            sprintf(buff,"Frame %d",k);
            Write(buff);}
    }
};
}

#endif
