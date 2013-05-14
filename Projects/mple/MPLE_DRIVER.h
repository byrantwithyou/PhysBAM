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
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include "MPLE_DOUBLE_WELL.h"
#include "MPLE_ITERATOR.h"
#include "MPLE_POINT.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace PhysBAM{

template<class TV,int w>
class MPLE_DRIVER: public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
    
public:
    
    ARRAY<MPLE_POINT<TV,w> > points;   // data
    GRID<TV> grid;                     // grid

    ARRAY<T,TV_INT> u;                // segmentation function
    ARRAY<T,TV_INT> u_new;            // new segmentation funtion
    ARRAY<TV,TV_INT> location;        // flat index to location
    ARRAY<TV_INT,TV_INT> index;       // flat to vector index
    
    int array_m;
    
    T cfl;
    T spread;
    T frame_dt;
    int frames;
    int timesteps;
    T mu,nu,epsilon;
    T dt,one_over_dx_squared;
    
    MPLE_DRIVER():cfl((T)1),spread((T)1),frame_dt((T)1/24),frames(100),mu(5e-4),nu(.05){}

    ~MPLE_DRIVER(){}

    void Initialize()
    {
        u.Resize(grid.Node_Indices(),false);
        u_new.Resize(grid.Node_Indices(),false);
        location.Resize(grid.Node_Indices(),false);
        index.Resize(grid.Node_Indices(),false);
        array_m=u.array.m;

        int k=0;
        for(NODE_ITERATOR<TV> it(grid);it.Valid();it.Next(),k++){
            location.array(k)=it.Location();
            index.array(k)=it.Node_Index();}

        RANDOM_NUMBERS<T> random;
        random.Set_Seed(0);
        for(int i=0;i<array_m;i++)
            u.array(i)=random.Get_Uniform_Number(0,1);

#pragma omp parallel for        
        for(int i=0;i<points.m;i++)
            points(i).Update_Base_And_Weights(grid);

        dt=cfl*sqr(grid.dX(0))/(2*(TV::m+1));
        timesteps=(int)(frame_dt/dt);
        epsilon=spread*grid.dX(0);
        one_over_dx_squared=1/sqr(grid.dX(0));
        PHYSBAM_ASSERT(grid.dX.Min()==grid.dX.Max());
    }

    void Write(const char* title)
    {
        for(int i=0;i<points.m;i++)
            Add_Debug_Particle<TV>(points(i).X,VECTOR<T,3>(1,0,0));
        
        for(int i=0;i<array_m;i++){
            T value=u.array(i);
            if(value){
                Add_Debug_Particle(location.array(i),value>0?VECTOR<T,3>(0,1,0):VECTOR<T,3>(0,.5,1));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_DISPLAY_SIZE,value);}}
        
        Flush_Frame<TV>(title);
    }

    void Diffusion_Step()
    {
#pragma omp parallel for        
        for(int i=0;i<array_m;i++){
            const TV_INT& this_node=index.array(i);
            T& value=u_new.array(i);
            value=0;
            for(int k=0;k<TV::m;k++){
                if(this_node(k)+1<grid.Domain_Indices().max_corner(k))
                    value+=u(this_node+TV_INT::Axis_Vector(k));
                if(this_node(k)-1>=grid.Domain_Indices().min_corner(k))
                    value+=u(this_node-TV_INT::Axis_Vector(k));}
            value-=u.array(i)*2*TV::m;
            value*=epsilon*dt*one_over_dx_squared;
            value+=u.array(i);
            value-=(dt/epsilon)*MPLE_DOUBLE_WELL<T>::Gradient(u.array(i));}
    }

    void Rasterization_Step()
    {
        for(int i=0;i<points.m;i++){
            for(MPLE_ITERATOR<TV,w> it(points(i));it.Valid();it.Next())
                u_new(it.Node())+=mu*dt*it.Weight()*((1-2*nu)/nu+sqr(2*nu-1)/((1-nu)*nu)*u(it.Node()));}
    }

    void Clamp_Step()
    {
#pragma omp parallel for        
        for(int i=0;i<array_m;i++)
            if(u_new.array(i)<0) u_new.array(i)=0;
    }

    void Advance_Frame(int frame)
    {
        LOG::cout<<"Frame "<<frame<<std::endl;
        for(int i=0;i<timesteps;i++){
            Diffusion_Step();
            Rasterization_Step();
            Clamp_Step();
            u.Exchange(u_new);}
    }

    void Run()
    {
        Initialize();
        Write("Frame 0");
        for(int k=1;k<=frames;k++){
            Advance_Frame(k);
            char buff[100];
            sprintf(buff,"Frame %d",k);
            Write(buff);}
    }
};
}

#endif
