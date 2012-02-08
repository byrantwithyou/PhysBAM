//#####################################################################
// Copyright 2002, Sergey Koltakov.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_VISCOSITY
//##################################################################### 
//
//#####################################################################
// Koltakov - August 12, 2003
//#####################################################################
#ifndef __IMPLICIT_VISCOSITY__
#define __IMPLICIT_VISCOSITY__

#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_CYLINDER.h>
#include "../WATER_FREE_SURFACE_3D_EXAMPLE.h"
namespace PhysBAM{

template<class T,class RW=T>
class IMPLICIT_VISCOSITY : public WATER_FREE_SURFACE_3D_EXAMPLE<T,RW>
{
public:
    IMPLICIT_VISCOSITY() 
    {
        first_frame=0;
        last_frame=5*(3*24);
        restart=false;
        restart_frame=0;
        m=50;n=50;mn=50;
        xmin=(T)0;xmax=(T)1;ymin=(T)0;ymax=(T)1;zmin=(T)0;zmax=(T)1;
        domain_walls[0][0]=true;domain_walls[0][1]=true;domain_walls[1][0]=true;domain_walls[1][1]=false;domain_walls[2][0]=true;domain_walls[2][1]=true;
        fluids_parameters.Initialize_Domain_Boundary_Conditions(); // sets up the proper wall states
        number_particles_per_cell=16;
        use_removed_positive_particles=false;use_removed_negative_particles=false;
        viscosity=(T)1e6*.001137;implicit_viscosity=true;use_central_differencing=false;second_order_pressure=false;variable_viscosity=true;
        write_levelset=true;write_velocity=true;write_removed_positive_particles=false;write_removed_negative_particles=false;
        matlab_directory="implicit_viscosity/matlab";output_directory="implicit_viscosity/output";
        delete_fluid_inside_objects=false;
    }

    ~IMPLICIT_VISCOSITY() 
    {}

//#####################################################################
// Function Initialize_Phi
//#####################################################################
void Initialize_Phi()
{
    for(int i=0;i<m;i++) for(int j=0;j<n;j++) for(int ij=0;ij<mn;ij++) 
        phi(i,j,ij)=minmag(grid.y(j)-(T).5,(T)sqrt(sqr(grid.x(i)-(T).5)+sqr(grid.y(j)-(T).75)+sqr(grid.z(ij)-(T).5))-(T).125);
}
//#####################################################################
// Function Specify_Variable_Viscosity_Values
//#####################################################################
void Specify_Variable_Viscosity_Values(ARRAY<T,VECTOR<int,3> >& viscosity)
{
    for(int i=0;i<m;i++) for(int j=0;j<n;j++) for(int ij=0;ij<n;ij++)
        if(grid.x(i) > .5) viscosity(i,j,ij)=1e6*.001137;
        else viscosity(i,j,ij)=(T)(1e6*.001137*1e-6);
}
};      
}
#endif
