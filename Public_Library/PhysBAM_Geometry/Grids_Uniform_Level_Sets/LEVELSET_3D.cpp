//#####################################################################
// Copyright 2002-2008, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Frank Losasso, Neil Molino, Igor Neverov, Avi Robinson-Mosher,
//     Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/RUNGEKUTTA.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> LEVELSET<VECTOR<T,3> >::
LEVELSET(GRID<TV>& grid_input,ARRAY<T,TV_INT>& phi_input,const int number_of_ghost_cells_input)
    :LEVELSET_BASE<TV>(grid_input,phi_input,number_of_ghost_cells_input)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> LEVELSET<VECTOR<T,3> >::
~LEVELSET()
{}
//#####################################################################
// Function Principal_Curvatures
//#####################################################################
template<class T> VECTOR<T,2> LEVELSET<VECTOR<T,3> >::
Principal_Curvatures(const TV& X) const
{
    TV grad_phi=TV((Phi(TV(X.x+grid.dX.x,X.y,X.z))-Phi(TV(X.x-grid.dX.x,X.y,X.z)))/(2*grid.dX.x),
                                     (Phi(TV(X.x,X.y+grid.dX.y,X.z))-Phi(TV(X.x,X.y-grid.dX.y,X.z)))/(2*grid.dX.y),
                                     (Phi(TV(X.x,X.y,X.z+grid.dX.z))-Phi(TV(X.x,X.y,X.z-grid.dX.z)))/(2*grid.dX.z));
    TV N=grad_phi;T grad_phi_magnitude=N.Normalize();
    SYMMETRIC_MATRIX<T,3> P=(T)1-SYMMETRIC_MATRIX<T,3>::Outer_Product(N),M=SYMMETRIC_MATRIX<T,3>::Conjugate(P,Hessian(X))/grad_phi_magnitude;
    T trace=M.Trace();
    QUADRATIC<T> quadratic(-1,trace,sqr(M(1,1))-M(0,1)*M(1,2)+sqr(M(2,1))-M(0,1)*M(2,3)+sqr(M(2,2))-M(1,2)*M(2,3));
    quadratic.Compute_Roots();
    if(quadratic.roots == 0) (T).5*VECTOR<T,2>(trace,trace);
    else if(quadratic.roots == 1) return VECTOR<T,2>(quadratic.root1,quadratic.root1);
    return VECTOR<T,2>(quadratic.root1,quadratic.root2);
}
//#####################################################################
// Function Approximate_Surface_Area
//#####################################################################
// calculates the approximate perimeter using delta functions
template<class T> T LEVELSET<VECTOR<T,3> >::
Approximate_Surface_Area(const T interface_thickness,const T time) const
{
    int ghost_cells=number_of_ghost_cells;
    ARRAY<T,TV_INT> phi_ghost(grid.Domain_Indices(ghost_cells));boundary->Fill_Ghost_Cells(grid,phi,phi_ghost,0,time,ghost_cells);
    T interface_half_width=interface_thickness*grid.dX.Max()/2,one_over_two_dx=1/(2*grid.dX.x),one_over_two_dy=1/(2*grid.dX.y),one_over_two_dz=1/(2*grid.dX.z),surface_area=0;
    for(int i=0;i<grid.counts.x;i++) for(int j=0;j<grid.counts.y;j++) for(int ij=0;ij<grid.counts.z;ij++)
        surface_area+=(LEVELSET_UTILITIES<T>::Delta(phi_ghost(i,j,ij),interface_half_width)*sqrt(sqr((phi_ghost(i+1,j,ij)-phi_ghost(i-1,j,ij))*one_over_two_dx)+
                                sqr((phi_ghost(i,j+1,ij)-phi_ghost(i,j-1,ij))*one_over_two_dy)+sqr((phi_ghost(i,j,ij+1)-phi_ghost(i,j,ij-1))*one_over_two_dz)));
    return surface_area*grid.dX.x*grid.dX.y*grid.dX.z;
}
//#####################################################################
// Function Calculate_Triangulated_Surface_From_Marching_Tetrahedra
//#####################################################################
// uses levelset grid for tet marching - faster than version below because we don't need to interpolate phi
template<class T> void LEVELSET<VECTOR<T,3> >::
Calculate_Triangulated_Surface_From_Marching_Tetrahedra(TRIANGULATED_SURFACE<T>& triangulated_surface,const bool include_ghost_values) const
{
    triangulated_surface.Clean_Memory();triangulated_surface.mesh.Clean_Memory();triangulated_surface.particles.Clean_Memory();
    int m_start=1,m_end=grid.counts.x,n_start=1,n_end=grid.counts.y,mn_start=1,mn_end=grid.counts.z;
    if(include_ghost_values){m_start=phi.domain.min_corner.x;m_end=phi.domain.max_corner.x;n_start=phi.domain.min_corner.y;n_end=phi.domain.max_corner.y;mn_start=phi.domain.min_corner.z;mn_end=phi.domain.max_corner.z;}
    ARRAY<VECTOR<int,6>,TV_INT> edge(m_start,m_end,n_start,n_end,mn_start,mn_end);edge.Fill(VECTOR<int,6>()-1);
    // create particles
    for(int i=m_start;i<m_end;i++) for(int j=n_start;j<n_end;j++) for(int k=mn_start;k<mn_end;k++){TV_INT index(i,j,k);
        if(i<m_end-1) edge(index)(0)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,index,TV_INT(i+1,j,k));
        if(j<n_end-1) edge(index)(1)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,index,TV_INT(i,j+1,k));
        if(k<mn_end-1)edge(index)(2)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,index,TV_INT(i,j,k+1));
        if((i+j+k)%2 == 0){
            if(j<n_end-1 && k<mn_end-1)edge(index)(3)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,TV_INT(i,j+1,k),TV_INT(i,j,k+1));
            if(i<m_end-1 && k<mn_end-1)edge(index)(4)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,TV_INT(i+1,j,k),TV_INT(i,j,k+1));
            if(i<m_end-1 && j<n_end-1) edge(index)(5)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,TV_INT(i,j+1,k),TV_INT(i+1,j,k));}
        else{
            if(j<n_end-1 && k<mn_end-1)edge(index)(3)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,index,TV_INT(i,j+1,k+1));
            if(i<m_end-1 && k<mn_end-1)edge(index)(4)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,index,TV_INT(i+1,j,k+1));
            if(i<m_end-1 && j<n_end-1) edge(index)(5)=If_Zero_Crossing_Add_Particle_By_Index(triangulated_surface,index,TV_INT(i+1,j+1,k));}}
    // calculate triangles
    for(int i=m_start;i<m_end-1;i++) for(int j=n_start;j<n_end-1;j++) for(int k=mn_start;k<mn_end-1;k++)
        if((i+j+k)%2 == 0){
            Append_Triangles(triangulated_surface,edge(i,j,k)(0),edge(i,j,k)(1),edge(i,j,k)(2),edge(i,j,k)(5),edge(i,j,k)(3),edge(i,j,k)(4),phi(i,j,k)); // bottom left
            Append_Triangles(triangulated_surface,edge(i,j,k+1)(5),edge(i+1,j,k)(3),edge(i+1,j,k+1)(1),edge(i,j,k)(4),edge(i+1,j,k)(2),edge(i,j,k+1)(0),phi(i+1,j+1,k+1)); // bottom right
            Append_Triangles(triangulated_surface,edge(i,j+1,k)(0),edge(i+1,j,k)(1),edge(i+1,j+1,k)(2),edge(i,j,k)(5),edge(i+1,j,k)(3),edge(i,j+1,k)(4),phi(i+1,j+1,k)); // top front
            Append_Triangles(triangulated_surface,edge(i,j+1,k)(2),edge(i,j,k)(3),edge(i,j+1,k)(4),edge(i,j,k+1)(1),edge(i,j,k+1)(5),edge(i,j+1,k+1)(0),phi(i,j+1,k)); // top back
            Append_Triangles(triangulated_surface,edge(i,j,k)(5),edge(i,j+1,k)(4),edge(i,j,k)(3),edge(i+1,j,k)(3),edge(i,j,k+1)(5),edge(i,j,k)(4),phi(i,j+1,k));} // center
        else{
            Append_Triangles(triangulated_surface,edge(i,j,k)(5),edge(i+1,j,k)(1),edge(i+1,j,k)(3),edge(i,j,k)(0),edge(i+1,j,k)(2),edge(i,j,k)(4),phi(i+1,j+1,k)); // bottom front
            Append_Triangles(triangulated_surface,edge(i,j,k)(4),edge(i,j,k)(3),edge(i,j,k)(2),edge(i,j,k+1)(5),edge(i,j,k+1)(1),edge(i,j,k+1)(0),phi(i,j,k)); // bottom back
            Append_Triangles(triangulated_surface,edge(i,j,k)(5),edge(i,j,k)(1),edge(i,j,k)(3),edge(i,j+1,k)(0),edge(i,j+1,k)(2),edge(i,j+1,k)(4),phi(i,j,k)); // top left
            Append_Triangles(triangulated_surface,edge(i,j+1,k)(4),edge(i,j+1,k+1)(0),edge(i,j,k+1)(5),edge(i+1,j+1,k)(2),edge(i+1,j,k+1)(1),edge(i+1,j,k)(3),phi(i,j+1,k+1)); // top right
            Append_Triangles(triangulated_surface,edge(i,j,k)(5),edge(i,j,k)(3),edge(i,j,k)(4),edge(i,j+1,k)(4),edge(i,j,k+1)(5),edge(i+1,j,k)(3),phi(i,j,k));} // center
    triangulated_surface.mesh.number_nodes=triangulated_surface.particles.Size();
    triangulated_surface.Remove_Degenerate_Triangles();
}
//#####################################################################
// Function If_Zero_Crossing_Add_Particle
//#####################################################################
template<class T> int LEVELSET<VECTOR<T,3> >::
If_Zero_Crossing_Add_Particle_By_Index(TRIANGULATED_SURFACE<T>& triangulated_surface,const TV_INT& index1,const TV_INT& index2) const
{
    int index=-1;T phi1=phi(index1),phi2=phi(index2);
    if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
        index=triangulated_surface.particles.Add_Element();
        triangulated_surface.particles.X(index)=LINEAR_INTERPOLATION<T,TV>::Linear(grid.X(index1),grid.X(index2),LEVELSET_UTILITIES<T>::Theta(phi1,phi2));}
    return index;
}
//#####################################################################
// Function Calculate_Triangulated_Surface_From_Marching_Tetrahedra
//#####################################################################
template<class T> void LEVELSET<VECTOR<T,3> >::
Calculate_Triangulated_Surface_From_Marching_Tetrahedra(const GRID<TV>& tet_grid,TRIANGULATED_SURFACE<T>& triangulated_surface) const
{
    assert(tet_grid.domain.min_corner.x >= grid.domain.min_corner.x && tet_grid.domain.max_corner.x < grid.domain.max_corner.x && tet_grid.domain.min_corner.y >= grid.domain.min_corner.y && tet_grid.domain.max_corner.y < grid.domain.max_corner.y && tet_grid.domain.min_corner.z >= grid.domain.min_corner.z &&
               tet_grid.domain.max_corner.z < grid.domain.max_corner.z);
    triangulated_surface.Clean_Memory();triangulated_surface.mesh.Clean_Memory();triangulated_surface.particles.Clean_Memory();
    ARRAY<VECTOR<int,6>,TV_INT> edge(0,tet_grid.counts.x,0,tet_grid.counts.y,0,tet_grid.counts.z);
    // create particles
    int i;for(i=0;i<tet_grid.counts.x;i++) for(int j=0;j<tet_grid.counts.y;j++) for(int k=0;k<tet_grid.counts.z;k++){
        edge(i,j,k)(0)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i,j,k),tet_grid.X(i+1,j,k));
        edge(i,j,k)(1)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i,j,k),tet_grid.X(i,j+1,k));
        edge(i,j,k)(2)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i,j,k),tet_grid.X(i,j,k+1));
        if((i+j+k)%2 == 0){
            edge(i,j,k)(3)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i,j+1,k),tet_grid.X(i,j,k+1));
            edge(i,j,k)(4)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i+1,j,k),tet_grid.X(i,j,k+1));
            edge(i,j,k)(5)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i,j+1,k),tet_grid.X(i+1,j,k));}
        else{
            edge(i,j,k)(3)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i,j,k),tet_grid.X(i,j+1,k+1));
            edge(i,j,k)(4)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i,j,k),tet_grid.X(i+1,j,k+1));
            edge(i,j,k)(5)=If_Zero_Crossing_Add_Particle(triangulated_surface,tet_grid.X(i,j,k),tet_grid.X(i+1,j+1,k));}}
    // calculate triangles
    for(i=0;i<tet_grid.counts.x-1;i++) for(int j=0;j<tet_grid.counts.y-1;j++) for(int k=0;k<tet_grid.counts.z-1;k++)
        if((i+j+k)%2 == 0){
            Append_Triangles(triangulated_surface,edge(i,j,k)(0),edge(i,j,k)(1),edge(i,j,k)(2),edge(i,j,k)(5),edge(i,j,k)(3),edge(i,j,k)(4),Phi(tet_grid.X(i,j,k))); // bottom left
            Append_Triangles(triangulated_surface,edge(i,j,k+1)(5),edge(i+1,j,k)(3),edge(i+1,j,k+1)(1),edge(i,j,k)(4),edge(i+1,j,k)(2),edge(i,j,k+1)(0),Phi(tet_grid.X(i+1,j+1,k+1))); // bottom right
            Append_Triangles(triangulated_surface,edge(i,j+1,k)(0),edge(i+1,j,k)(1),edge(i+1,j+1,k)(2),edge(i,j,k)(5),edge(i+1,j,k)(3),edge(i,j+1,k)(4),Phi(tet_grid.X(i+1,j+1,k))); // top front
            Append_Triangles(triangulated_surface,edge(i,j+1,k)(2),edge(i,j,k)(3),edge(i,j+1,k)(4),edge(i,j,k+1)(1),edge(i,j,k+1)(5),edge(i,j+1,k+1)(0),Phi(tet_grid.X(i,j+1,k))); // top back
            Append_Triangles(triangulated_surface,edge(i,j,k)(5),edge(i,j+1,k)(4),edge(i,j,k)(3),edge(i+1,j,k)(3),edge(i,j,k+1)(5),edge(i,j,k)(4),Phi(tet_grid.X(i,j+1,k)));} // center
        else{
            Append_Triangles(triangulated_surface,edge(i,j,k)(5),edge(i+1,j,k)(1),edge(i+1,j,k)(3),edge(i,j,k)(0),edge(i+1,j,k)(2),edge(i,j,k)(4),Phi(tet_grid.X(i+1,j+1,k))); // bottom front
            Append_Triangles(triangulated_surface,edge(i,j,k)(4),edge(i,j,k)(3),edge(i,j,k)(2),edge(i,j,k+1)(5),edge(i,j,k+1)(1),edge(i,j,k+1)(0),Phi(tet_grid.X(i,j,k))); // bottom back
            Append_Triangles(triangulated_surface,edge(i,j,k)(5),edge(i,j,k)(1),edge(i,j,k)(3),edge(i,j+1,k)(0),edge(i,j+1,k)(2),edge(i,j+1,k)(4),Phi(tet_grid.X(i,j,k))); // top left
            Append_Triangles(triangulated_surface,edge(i,j+1,k)(4),edge(i,j+1,k+1)(0),edge(i,j,k+1)(5),edge(i+1,j+1,k)(2),edge(i+1,j,k+1)(1),edge(i+1,j,k)(3),Phi(tet_grid.X(i,j+1,k+1))); // top right
            Append_Triangles(triangulated_surface,edge(i,j,k)(5),edge(i,j,k)(3),edge(i,j,k)(4),edge(i,j+1,k)(4),edge(i,j,k+1)(5),edge(i+1,j,k)(3),Phi(tet_grid.X(i,j,k)));} // center
    triangulated_surface.mesh.number_nodes=triangulated_surface.particles.Size();
    triangulated_surface.Remove_Degenerate_Triangles();
}
//#####################################################################
// Function If_Zero_Crossing_Add_Particle
//#####################################################################
template<class T> int LEVELSET<VECTOR<T,3> >::
If_Zero_Crossing_Add_Particle(TRIANGULATED_SURFACE<T>& triangulated_surface,const TV& x1,const TV& x2) const
{
    int index=-1;T phi1=Phi(x1),phi2=Phi(x2);
    if(LEVELSET_UTILITIES<T>::Interface(phi1,phi2)){
        index=triangulated_surface.particles.Add_Element();
        triangulated_surface.particles.X(index)=LINEAR_INTERPOLATION<T,TV>::Linear(x1,x2,LEVELSET_UTILITIES<T>::Theta(phi1,phi2));}
    return index;
}
//#####################################################################
// Function Append_Triangles
//#####################################################################
// looking down at the node with phi1, e1-e2-e3 is conunter clockwise
// edges 1,2,3 come out of phi1 - edge4 is opposite edge3 - edge2 is opposite edge6 - edge1 is opposite edge5
template<class T> void LEVELSET<VECTOR<T,3> >::
Append_Triangles(TRIANGULATED_SURFACE<T>& triangulated_surface,const int e1,const int e2,const int e3,const int e4,const int e5,const int e6,const T phi1) const
{
    int number_positive=(e1>=0)+(e2>=0)+(e3>=0)+(e4>=0)+(e5>=0)+(e6>=0);if(number_positive==0) return;assert(number_positive==3 || number_positive==4);
    if(e1>=0 && e2>=0 && e5>=0 && e6>=0){ // 2 triangles
        if(phi1>0){triangulated_surface.mesh.elements.Append(TV_INT(e1,e6,e2));triangulated_surface.mesh.elements.Append(TV_INT(e2,e6,e5));}
        else{triangulated_surface.mesh.elements.Append(TV_INT(e1,e2,e6));triangulated_surface.mesh.elements.Append(TV_INT(e2,e5,e6));}}
    else if(e2>=0 && e3>=0 && e4>=0 && e6>=0){
        if(phi1>0){triangulated_surface.mesh.elements.Append(TV_INT(e2,e4,e3));triangulated_surface.mesh.elements.Append(TV_INT(e3,e4,e6));}
        else{triangulated_surface.mesh.elements.Append(TV_INT(e2,e3,e4));triangulated_surface.mesh.elements.Append(TV_INT(e4,e3,e6));}}
    else if(e1>=0 && e3>=0 && e4>=0 && e5>=0){
        if(phi1>0){triangulated_surface.mesh.elements.Append(TV_INT(e1,e3,e5));triangulated_surface.mesh.elements.Append(TV_INT(e1,e5,e4));}
        else {triangulated_surface.mesh.elements.Append(TV_INT(e1,e5,e3));triangulated_surface.mesh.elements.Append(TV_INT(e1,e4,e5));}}
    else if(e1>=0 && e2>=0 && e3>=0){ // 1 triangle
        if(phi1>0) triangulated_surface.mesh.elements.Append(TV_INT(e1,e3,e2));
        else triangulated_surface.mesh.elements.Append(TV_INT(e1,e2,e3));}
    else if(e1>=0 && e4>=0 && e6>=0){
        if(phi1>0) triangulated_surface.mesh.elements.Append(TV_INT(e1,e6,e4));
        else triangulated_surface.mesh.elements.Append(TV_INT(e1,e4,e6));}
    else if(e3>=0 && e5>=0 && e6>=0){
        if(phi1>0) triangulated_surface.mesh.elements.Append(TV_INT(e3,e5,e6));
        else triangulated_surface.mesh.elements.Append(TV_INT(e3,e6,e5));}
    else if(e4>=0 && e2>=0 && e5>=0){
        if(phi1>0) triangulated_surface.mesh.elements.Append(TV_INT(e4,e5,e2));
        else triangulated_surface.mesh.elements.Append(TV_INT(e4,e2,e5));}
}
//#####################################################################
template class LEVELSET<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LEVELSET<VECTOR<double,3> >;
#endif
