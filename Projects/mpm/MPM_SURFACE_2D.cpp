//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include "MPM_SURFACE_2D.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_SURFACE_2D<TV>::
MPM_SURFACE_2D(const MPM_SIMULATION<TV>& sim):sim(sim)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_SURFACE_2D<TV>::
~MPM_SURFACE_2D()
{}
//#####################################################################
// Function Initialize_With_A_Box
//#####################################################################
template<class TV> void MPM_SURFACE_2D<TV>::
Initialize_With_A_Box(const T h,const RANGE<TV>& box)
{
    int Nx=floor((box.max_corner(0)-box.min_corner(0))/h)+1;
    int Ny=floor((box.max_corner(1)-box.min_corner(1))/h)+1;
    GRID<TV> grid(TV_INT(Nx,Ny),box);
    int N_boundary_nodes=2*(Nx+Ny)-4;
    curve.particles.Add_Elements(N_boundary_nodes);
    curve.mesh.number_nodes=N_boundary_nodes;
    curve.mesh.elements.Exact_Resize(N_boundary_nodes);
    int count=0;
    for(int i=0,j=0;j<Ny-1;j++){
        curve.particles.X(count)=grid.Node(TV_INT(i,j));
        curve.mesh.elements(count).Set(count,count+1);
        count++;}
    for(int j=Ny-1,i=0;i<Nx-1;i++){
        curve.particles.X(count)=grid.Node(TV_INT(i,j));
        curve.mesh.elements(count).Set(count,count+1);
        count++;}
    for(int i=Nx-1,j=Ny-1;j>0;j--){
        curve.particles.X(count)=grid.Node(TV_INT(i,j));
        curve.mesh.elements(count).Set(count,count+1);
        count++;}
    for(int j=0,i=Nx-1;i>1;i--){
        curve.particles.X(count)=grid.Node(TV_INT(i,j));
        curve.mesh.elements(count).Set(count,count+1);
        count++;}
    curve.particles.X(count)=grid.Node(TV_INT(1,0));
    curve.mesh.elements(count).Set(count,0);
    count++;
    PHYSBAM_ASSERT(count==N_boundary_nodes);
    curve.Update_Segment_List();
    Xm=curve.particles.X;
}

//#####################################################################
template class MPM_SURFACE_2D<VECTOR<float,2> >;
template class MPM_SURFACE_2D<VECTOR<double,2> >;
}


