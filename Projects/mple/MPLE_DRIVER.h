//#####################################################################
// Copyright 2013, Alexey Stomakhin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################

#ifndef __MPLE_DRIVER__
#define __MPLE_DRIVER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/NODE_ITERATOR.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include "MPLE_DOUBLE_WELL.h"
#include "MPLE_ITERATOR.h"
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

    ARRAY<T,TV_INT> u;                // segmentation function
    ARRAY<T,TV_INT> u_new;            // new segmentation funtion
    ARRAY<TV,TV_INT> location;        // flat index to location
    ARRAY<TV_INT,TV_INT> index;       // flat to vector index
    
    int frames;
    int timesteps;
    T mu,nu,epsilon;
    T dt,one_over_dx_squared;
    
    MPLE_DRIVER():frames(1000),timesteps(1),mu(.1),nu(.05){}

    ~MPLE_DRIVER(){}

    void Initialize()
    {
        u.Resize(grid.Node_Indices(ghost),false);
        u_new.Resize(grid.Node_Indices(ghost),false);
        location.Resize(grid.Node_Indices(ghost),false);
        index.Resize(grid.Node_Indices(ghost),false);
        
        int k=0;
        for(NODE_ITERATOR<TV> it(grid,ghost);it.Valid();it.Next(),k++){
            location.array(k)=it.Location();
            index.array(k)=it.Node_Index();}

        for(int i=0;i<u.array.m;i++){
            u_new.array(i)=0;
            u.array(i)=0;}
        
        for(int i=0;i<points.m;i++)
            points(i).Update_Base_And_Weights(grid);

        dt=sqr(grid.dX(0))/(2*(TV::m+1));
        epsilon=3*grid.dX(0);
        one_over_dx_squared=1/sqr(grid.dX(0));
        PHYSBAM_ASSERT(grid.dX.Min()==grid.dX.Max());
    }

    void Write(const char* title)
    {
        for(int i=0;i<points.m;i++)
            Add_Debug_Particle<TV>(points(i).X,VECTOR<T,3>(1,0,0));
        
        for(int i=0;i<location.array.m;i++){
            T value=u.array(i);
            if(value){
                Add_Debug_Particle(location.array(i),value>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(0,.5,1));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,value);}}
        
        Flush_Frame<TV>(title);
    }

    void Diffusion_Step()
    {
        for(NODE_ITERATOR<TV> it(grid);it.Valid();it.Next()){
            const TV_INT& this_node=it.Node_Index();
            T& value=u_new(this_node);
            value=0;
            for(int k=0;k<TV::m;k++){
                value+=u(this_node+TV_INT::Axis_Vector(k));
                value+=u(this_node-TV_INT::Axis_Vector(k));}
            value-=u(this_node)*2*TV::m;
            value*=epsilon*dt*one_over_dx_squared;
            value+=u(this_node);
            value-=(dt/epsilon)*MPLE_DOUBLE_WELL<T>::Gradient(u(this_node));
        }
    }

    void Rasterization_Step()
    {
        for(int i=0;i<points.m;i++){
            for(MPLE_ITERATOR<TV,w> it(points(i));it.Valid();it.Next())
                u_new(it.Node())+=mu*dt*it.Weight()*((1-2*nu)/nu+sqr(2*nu-1)/((1-nu)*nu)*u(it.Node()));}
    }

    void Clamp_Step()
    {
        for(int i=0;i<u.array.m;i++)
            if(u_new.array(i)<0) u_new.array(i)=0;
    }

    void Advance_Frame()
    {
        for(int i=0;i<timesteps;i++){
            Diffusion_Step();
            Rasterization_Step();
            Clamp_Step();
            u.Exchange(u_new);}
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
