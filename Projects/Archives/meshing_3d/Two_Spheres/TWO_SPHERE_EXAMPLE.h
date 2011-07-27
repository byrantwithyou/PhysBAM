//#####################################################################
// Copyright 2002, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TWO_SPHERE_EXAMPLE
//#####################################################################
//
//#####################################################################
// Molino - October 8, 2002
//#####################################################################
#ifndef __TWO_SPHERE_EXAMPLE__
#define __TWO_SPHERE_EXAMPLE__

#include "../MESHING_EXAMPLE.h"
namespace PhysBAM{

template<class T,class RW=T>
class TWO_SPHERE_EXAMPLE:public MESHING_EXAMPLE<T,RW>
{
public:
    TWO_SPHERE_EXAMPLE()
    {
        //altitude_spring_overdamping_fraction=0; // debugging
    }

    ~TWO_SPHERE_EXAMPLE()
    {}

//#####################################################################
// Function Initialize_Implicit_Surface
//#####################################################################
void Initialize_Implicit_Surface(LEVELSET_IMPLICIT_SURFACE<T>& implicit_surface) PHYSBAM_OVERRIDE
{
    std::fstream input;
    char File_Name[256];sprintf(File_Name,"two_spheres.phi");
    input.open(File_Name,std::ios::in|std::ios::binary);
    if (input.is_open()) {
        implicit_surface.template Read<RW>(input); input.close(); 
    }
    else {
        for(int i=1;i<=implicit_surface.levelset.grid.m;i++) for(int j=1;j<=implicit_surface.levelset.grid.n;j++) for(int ij=1;ij<=implicit_surface.levelset.grid.mn;ij++){
            double x=implicit_surface.levelset.grid.x(i),y=implicit_surface.levelset.grid.y(j),z=implicit_surface.levelset.grid.z(ij);
            VECTOR_3D X(x,y,z);

            implicit_surface.levelset.phi(i,j,ij)=min((X-VECTOR_3D(-.7,0,0)).Magnitude()-1, (X-VECTOR_3D(.5,0,0)).Magnitude()-.64);}
        implicit_surface.levelset.Reinitialize();
    }
}
//#####################################################################
};
}
#endif

