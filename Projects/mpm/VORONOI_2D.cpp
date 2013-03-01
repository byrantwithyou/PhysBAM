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
    GRID<TV> node_grid(TV_INT(particle_grid.counts.x*2+1,particle_grid.counts.y*2+1),RANGE<TV>(particle_grid.domain.min_corner-particle_grid.dX/2.0,TV(particle_grid.domain.max_corner+particle_grid.dX/2.0)));
    HASHTABLE<TV_INT,int> particle_index;
    int ID=0;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+particle_grid.counts));it.Valid();it.Next()) particle_index.Get_Or_Insert(it.index)=ID++;
    HASHTABLE<TV_INT,int> node_index;
    ID=0;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+node_grid.counts));it.Valid();it.Next()){
        int n_old_component=(int)(it.index.x%2==1)+(int)(it.index.y%2==1);
        if(n_old_component==0){
            node_index.Get_Or_Insert(it.index)=ID++;
            Xm.Append(node_grid.Node(it.index));
            if(it.index.x==0||it.index.x==node_grid.counts.x-1||it.index.y==0||it.index.y==node_grid.counts.y-1) type.Append(1);
            else type.Append(2);}
        else if(n_old_component==1){
            node_index.Get_Or_Insert(it.index)=ID++;
            Xm.Append(node_grid.Node(it.index));
            type.Append(3);}}
    ID=0;
    elements.Resize(particle_grid.counts.Product());
    for(int i=0;i<node_grid.counts.x-1;i+=2){
        for(int j=0;j<node_grid.counts.y-1;j+=2){
            TV_INT aa(i,j),bb(i+1,j),cc(i+2,j),dd(i+2,j+1),ee(i+2,j+2),ff(i+1,j+2),gg(i,j+2),hh(i,j+1);
            int a=node_index.Get_Or_Insert(aa),b=node_index.Get_Or_Insert(bb),c=node_index.Get_Or_Insert(cc),d=node_index.Get_Or_Insert(dd),e=node_index.Get_Or_Insert(ee),f=node_index.Get_Or_Insert(ff),g=node_index.Get_Or_Insert(gg),h=node_index.Get_Or_Insert(hh);
            elements(ID).Append(a);elements(ID).Append(b);elements(ID).Append(c);elements(ID).Append(d);elements(ID).Append(e);elements(ID).Append(f);elements(ID).Append(g);elements(ID).Append(h);
            ID++;}}
    X=Xm;
    Initialize_Neighbor_Cells();
    Build_Association();
}
//#####################################################################
// Function Initialize_Neighbor_Cells
//#####################################################################
template<class T> void VORONOI_2D<T>::
Initialize_Neighbor_Cells()
{
    HASHTABLE<TV_INT,bool> neighbor_cells_hash;
    for(int c1=0;c1<elements.m;c1++){
        for(int c2=c1+1;c2<elements.m;c2++){
            ARRAY<int> common_nodes;
            common_nodes.Find_Common_Elements(elements(c1),elements(c2));
            for(int v=0;v<common_nodes.m;v++){
                if(type(common_nodes(v))==3){
                    neighbor_cells_hash.Get_Or_Insert(TV_INT(c1,c2))=true;
                    break;}}}}
    for(typename HASHTABLE<TV_INT,bool>::ITERATOR it(neighbor_cells_hash);it.Valid();it.Next()) neighbor_cells.Append(TRIPLE<int,int,bool>(it.Key().x,it.Key().y,true));
}
//#####################################################################
// Function Build_Association
//#####################################################################
template<class T> void VORONOI_2D<T>::
Build_Association()
{
    association.Clean_Memory();
    association.Resize(Xm.m);
    for(int c=0;c<elements.m;c++) for(int v=0;v<elements(c).m;v++) association(elements(c)(v)).Append(c);
}
//#####################################################################
// Function Build_Segments
//#####################################################################
template<class T> void VORONOI_2D<T>::
Build_Segments()
{
    segments.Clean_Memory();
    HASHTABLE<TV_INT,bool> seg;
    for(int e=0;e<elements.m;e++){
        int total_nodes=elements(e).m;
        for(int a=0;a<total_nodes;a++){
            int first=elements(e)(a),second=elements(e)(a+1);
            if(a==total_nodes-1) second=elements(e)(0);
            if(first<second) seg.Get_Or_Insert(TV_INT(first,second))=true;
            else seg.Get_Or_Insert(TV_INT(second,first))=true;}}
    for(typename HASHTABLE<TV_INT,bool>::ITERATOR it(seg);it.Valid();it.Next()) segments.Append(it.Key());
}
//#####################################################################
// Function Build_Boundary_Segments
//#####################################################################
template<class T> void VORONOI_2D<T>::
Build_Boundary_Segments()
{
    boundary_segments.Clean_Memory();
    HASHTABLE<TV_INT,int> seg;
    for(int e=0;e<elements.m;e++){
        int total_nodes=elements(e).m;
        for(int a=0;a<total_nodes;a++){
            int first=elements(e)(a),second=elements(e)(a+1);
            if(a==total_nodes-1) second=elements(e)(0);
            if(first<second) seg.Get_Or_Insert(TV_INT(first,second))++;
            else seg.Get_Or_Insert(TV_INT(second,first))++;}}
    for(typename HASHTABLE<TV_INT,int>::ITERATOR it(seg);it.Valid();it.Next()) if(it.Data()==1) boundary_segments.Append(it.Key());
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
