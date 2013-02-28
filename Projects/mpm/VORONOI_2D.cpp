//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include "VORONOI_2D.h"
namespace PhysBAM{
//#####################################################################
// Function Initialize_With_A_Regular_Grid_Of_Particles
//#####################################################################
template<class T> void VORONOI_2D<T>::
Initialize_With_A_Regular_Grid_Of_Particles(const GRID<TV>& particle_grid)
{
    HASHTABLE<TV_INT,int> particle_index;
    int ID=0;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+particle_grid.counts));it.Valid();it.Next()) particle_index.Get_Or_Insert(it.index)=ID++;
    RANGE<TV> segmented_mesh_domain(particle_grid.domain.min_corner-(T)0.5*particle_grid.dX,particle_grid.domain.max_corner+(T)0.5*particle_grid.dX);
    GRID<TV> node_grid(particle_grid.counts+TV_INT(1,1),segmented_mesh_domain);
    HASHTABLE<TV_INT,int> node_index;
    ID=0;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+particle_grid.counts+1));it.Valid();it.Next()) node_index.Get_Or_Insert(it.index)=ID++;
    int N_mesh_nodes=3*node_grid.counts.Product()-node_grid.counts.Sum();
    X.Resize(N_mesh_nodes);
    association.Resize(N_mesh_nodes);
    mesh.Set_Number_Nodes(N_mesh_nodes);
    ID=0;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+node_grid.counts));it.Valid();it.Next()){
        X(ID++)=node_grid.Node(it.index);}
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+particle_grid.counts));it.Valid();it.Next()){
        TV particleX=particle_grid.Node(it.index);
        TV_INT cell=node_grid.Cell(particleX,0);
        association(node_index.Get_Or_Insert(cell)).Append(particle_index.Get_Or_Insert(it.index));
        association(node_index.Get_Or_Insert(cell+TV_INT(1,0))).Append(particle_index.Get_Or_Insert(it.index));
        association(node_index.Get_Or_Insert(cell+TV_INT(0,1))).Append(particle_index.Get_Or_Insert(it.index));
        association(node_index.Get_Or_Insert(cell+TV_INT(1,1))).Append(particle_index.Get_Or_Insert(it.index));}
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+node_grid.counts));it.Valid();it.Next()){
        if(it.index(1)+1<node_grid.counts(1)){
            X(ID++)=node_grid.Node(it.index)+TV((T)0,(T)0.5*node_grid.dX(1));
            mesh.elements.Append(TV_INT(node_index.Get_Or_Insert(it.index),ID-1));
            mesh.elements.Append(TV_INT(ID-1,node_index.Get_Or_Insert(it.index+TV_INT(0,1))));
            int A=node_index.Get_Or_Insert(it.index),B=node_index.Get_Or_Insert(it.index+TV_INT(0,1));
            association(ID-1).Find_Common_Elements(association(A),association(B));}}
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+node_grid.counts));it.Valid();it.Next()){
        if(it.index(0)+1<node_grid.counts(0)){
            X(ID++)=node_grid.Node(it.index)+TV((T)0.5*node_grid.dX(0),(T)0);
            mesh.elements.Append(TV_INT(node_index.Get_Or_Insert(it.index),ID-1));
            mesh.elements.Append(TV_INT(ID-1,node_index.Get_Or_Insert(it.index+TV_INT(1,0))));
            int A=node_index.Get_Or_Insert(it.index),B=node_index.Get_Or_Insert(it.index+TV_INT(1,0));
            association(ID-1).Find_Common_Elements(association(A),association(B));}}
    Xm=X;
}
//#####################################################################
// Function Deform_Mesh_Using_Particle_Deformation
//#####################################################################
template<class T> void VORONOI_2D<T>::
Deform_Mesh_Using_Particle_Deformation(const ARRAY_VIEW<TV>& particle_Xm,const ARRAY_VIEW<TV>& particle_X,const ARRAY_VIEW<MATRIX<T,TV::m> >& particle_Fe,const ARRAY_VIEW<MATRIX<T,TV::m> >& particle_Fp)
{
    ARRAY<TV> particle_b(particle_X.m);
    ARRAY<MATRIX<T,TV::m> > particle_F(particle_X.m);
    for(int i=0;i<particle_b.m;i++){
        particle_F(i)=particle_Fe(i)*particle_Fp(i);
        particle_b(i)=particle_X(i)-particle_F(i)*particle_Xm(i);}
    for(int i=0;i<X.m;i++){
        X(i)=TV();
        for(int p=0;p<association(i).m;p++){
            int particle=association(i)(p);
            X(i)+=particle_F(particle)*Xm(i)+particle_b(particle);}
        X(i)/=association(i).m;}
}
//#####################################################################
template class VORONOI_2D<float>;
template class VORONOI_2D<double>;
}
