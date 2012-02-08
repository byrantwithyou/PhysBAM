//#####################################################################
// Copyright 2005, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BREAKING_WAVE
//#####################################################################
#ifndef __BREAKING_WAVE__
#define __BREAKING_WAVE__

#include "BOUNDARY_CUSTOM.h"
#include <Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_2D.h>
namespace PhysBAM{

template<class T,class RW=T>
class BREAKING_WAVE:public SOLIDS_FLUIDS_EXAMPLE_2D<RW>
{
public:
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::frame_rate;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart;
    using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::output_directory;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::fluids_parameters;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_frame_title;using SOLIDS_FLUIDS_EXAMPLE_2D<RW>::write_substeps;

    BOUNDARY_CUSTOM<T,T> phi_periodic;
    BOUNDARY_CUSTOM<T,VECTOR_2D<T> > fluid_periodic;

    BOUNDARY<T,T> phi_constant;

    T omega;//2*pi/period;
    T epsilon;
    T depth;

    ARRAY<T,VECTOR<int,1> > initial_depth;

    ARRAY<T,VECTOR<int,2> > object_phi;
    ARRAY<VECTOR_2D<T> ,VECTOR<int,2> > object_V;

    BREAKING_WAVE(std::string desc,T omega_input,T epsilon_input,T depth_input)
        :SOLIDS_FLUIDS_EXAMPLE_2D<RW>(fluids_parameters.WATER),omega(omega_input),epsilon(epsilon_input),depth(depth_input)
    {       
        std::cout<<"Running with Omega="<<omega<<", Epsilon="<<epsilon<<", Depth="<<depth<<std::endl;

        fluids_parameters.phi_boundary=&phi_constant;
        //fluids_parameters.fluid_boundary=&fluid_periodic;
        fluids_parameters.grid.Initialize(300,75,-.5,3.5,0,1);
        first_frame=0;last_frame=2000;
        frame_rate=24;
        restart=false;restart_frame=18;
        fluids_parameters.domain_walls[0][0]=false;fluids_parameters.domain_walls[0][1]=false;fluids_parameters.domain_walls[1][0]=true;fluids_parameters.domain_walls[1][1]=false; 
        fluids_parameters.number_particles_per_cell=32;
        fluids_parameters.particle_half_bandwidth=1; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fluids_parameters.reseeding_frame_rate=1;
        fluids_parameters.bias_towards_negative_particles=true;
        fluids_parameters.viscosity=0;//(T)(.001137*2e7);
        fluids_parameters.incompressible_iterations=200;
        write_frame_title=true;;
        fluids_parameters.write_levelset=true;fluids_parameters.write_velocity=true;fluids_parameters.write_particles=true;fluids_parameters.write_removed_positive_particles=false;fluids_parameters.write_removed_negative_particles=false;
        fluids_parameters.write_debug_data=true;fluids_parameters.write_ghost_values=true;//write_substeps=true;write_frame_title=true;
        output_directory=STRING_UTILITIES::string_sprintf("Breaking_Wave/%s_%.2f_%.2f_%.2f",desc.c_str(),omega,epsilon,depth);
//        output_directory="Breaking_Wave/output";
        fluids_parameters.gravity=9.8;
        //gravity_direction.x=.5;gravity_direction.Normalize();
        fluids_parameters.delete_fluid_inside_objects=true;
        fluids_parameters.enforce_divergence_free_extrapolation=false;
        fluids_parameters.levelset_substeps=2.9;


        initial_depth.Resize(1,fluids_parameters.grid.m);
/*        for(int i=0;i<fluids_parameters.grid.m;i++){
            T x=fluids_parameters.grid.x(i);
            initial_depth(i)=depth+1./(2*pi)*(epsilon*cos(2*pi*x)+.5*sqr(epsilon)*cos(4*pi*x)+3./8.*cube(epsilon)*cos(6*pi*x));}*/
/*        for(int i=0;i<fluids_parameters.grid.m;i++){
            T x=fluids_parameters.grid.x(i);
            if(fabs(x)<.25)
                initial_depth(i)=depth+.1*(1+cos(4*pi*x));
            else
            initial_depth(i)=depth;*/
        for(int i=0;i<fluids_parameters.grid.m;i++){
            initial_depth(i)=depth;
        }

        object_phi.Resize(fluids_parameters.grid,3);object_V.Resize(fluids_parameters.grid,3);
        for(int i=object_phi.m_start;i<=object_phi.m_end;i++)for(int j=object_phi.n_start;j<=object_phi.n_end;j++){
            object_phi(i,j)=1;
/*            T x=fluids_parameters.grid.x(i),y=fluids_parameters.grid.y(j);
            if(1.5<x&&x<2.5)
                object_phi(i,j)=y-.25;
            if(x<1.5)
                object_phi(i,j)=y-.25+.5*(1.5-x);
            if(x>2.5)
            object_phi(i,j)=y-.25+.5*(x-2.5);*/}
        LEVELSET_2D<T> levelset(fluids_parameters.grid,object_phi);
        levelset.Fast_Marching_Method();
    }

//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi() PHYSBAM_OVERRIDE
{
    GRID<TV>& grid=fluids_parameters.grid;
    for(int i=0;i<grid.m;i++)for(int j=0;j<grid.n;j++){
        //fluids_parameters.particle_levelset_evolution.phi(i,j)=grid.y(j)-.5;
        fluids_parameters.particle_levelset_evolution.phi(i,j)=grid.y(j)-initial_depth(i);
    }
}

//#####################################################################
// Function Get_Dirichlet_Boundary_Conditions
//#####################################################################
void Set_Dirichlet_Boundary_Conditions(const T time) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_2D<RW>::Set_Dirichlet_Boundary_Conditions(time);
    PROJECTION_2D<T>& projection=fluids_parameters.incompressible.projection;
    ARRAY<bool,VECTOR<int,2> >& psi_D=fluids_parameters.incompressible.projection.elliptic_solver->psi_D;
    for(int j=0;j<fluids_parameters.p_grid.n;j++){
        psi_D(0,j)=true;
        projection.p(0,j)=max((T)0,(depth-fluids_parameters.p_grid.y(j))*fluids_parameters.gravity);
        //printf("%f\n",projection.p(0,j));
    }
    for(int j=0;j<fluids_parameters.p_grid.n;j++){
        psi_D(fluids_parameters.p_grid.m+1,j)=true;
        projection.p(fluids_parameters.p_grid.m+1,j)=max((T)0,(depth-fluids_parameters.p_grid.y(j))*fluids_parameters.gravity);
        //printf("%f\n",projection.p(fluids_parameters.p_grid.m+1,j));
    }
}

//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Velocities() PHYSBAM_OVERRIDE
{
    /*   for(int i=0;i<fluids_parameters.grid.m;i++) for(int j=0;j<fluids_parameters.grid.n;j++){
        T x=fluids_parameters.grid.x(i),y=fluids_parameters.grid.y(j);
        T first_part=omega*exp(-2*pi*(depth-y));
        if(x>1.5)
            fluids_parameters.incompressible.V(i,j).x=-1./(2*pi)*first_part*(epsilon*cos(2*pi*x)+.5*sqr(epsilon)*cos(4*pi*x)+3./8.*cube(epsilon)*cos(6*pi*x));
        else
            fluids_parameters.incompressible.V(i,j).x=1./(2*pi)*first_part*(epsilon*cos(2*pi*x)+.5*sqr(epsilon)*cos(4*pi*x)+3./8.*cube(epsilon)*cos(6*pi*x));
        fluids_parameters.incompressible.V(i,j).y=1./(2*pi)*first_part*(epsilon*sin(2*pi*x)+.5*sqr(epsilon)*sin(4*pi*x)+3./8.*cube(epsilon)*sin(6*pi*x));
        }*/
       for(int i=0;i<fluids_parameters.grid.m;i++) for(int j=0;j<fluids_parameters.grid.n;j++){
           T x=fluids_parameters.grid.x(i),y=fluids_parameters.grid.y(j);
           T center_x=2;
           if(x>center_x){
               T distance=x-center_x;
               fluids_parameters.incompressible.V(i,j).x=-2*min(distance*4,(T)1);
           }
           else{
               T distance=center_x-x;
               fluids_parameters.incompressible.V(i,j).x=1*min(distance*3,(T)1)*2*(depth-y);
               fluids_parameters.incompressible.V(i,j).y=-.5*min(distance*5,(T)1)*2*(depth-y);
           }
        }
}
/*
//#####################################################################
// Function Extrapolate_Phi_Into_Objects
//#####################################################################
void Extrapolate_Phi_Into_Objects(const T time) PHYSBAM_OVERRIDE
{
    fluids_parameters.Extrapolate_Phi_Into_Object(object_phi);
}
//#####################################################################
// Function Extrapolate_Phi_Into_Objects
//#####################################################################
void Adjust_Phi_With_Objects(const T time)
{
    fluids_parameters.Adjust_Phi_With_Object(object_phi,object_V,time);
}
//#####################################################################
// Function Get_Object_Velocities
//#####################################################################
void Get_Object_Velocities(const T dt,const T time) PHYSBAM_OVERRIDE
{
    fluids_parameters.Extrapolate_Velocity_Into_Object(object_phi,object_V,3,true,time);
}
//#####################################################################*/

};    
}
#endif
