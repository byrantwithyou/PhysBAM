//#####################################################################
// Copyright 2005-2007, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Read_Write/FILE_UTILITIES.h>
#include <Grid_Tools/Grids/CELL_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_ITERATOR.h>
#include <Grid_Tools/Grids/FACE_RANGE_ITERATOR.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Basic_Geometry/RAY.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Level_Sets/FAST_MARCHING_METHOD_UNIFORM.h>
#include <Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
template<class T,int d>
struct HEAP_ENTRY
{
    VECTOR<int,d> index;
    int* close_ptr;
    T dist;
};
//#####################################################################
// Function Exchange
//#####################################################################
template<class T,int d> void
Exchange(ARRAY<HEAP_ENTRY<T,d> >& heap,int& ia,int& ib)
{
    HEAP_ENTRY<T,d> &ha=heap(ia),&hb=heap(ib);
    exchange(ha,hb);
    *ha.close_ptr=ia;
    *hb.close_ptr=ib;
    exchange(ia,ib);
}
//#####################################################################
// Function Up_Heap
//#####################################################################
template<class T,int d> static void
Up_Heap(ARRAY<HEAP_ENTRY<T,d> >& heap,int index)
{
    while(index>0){
        int parent=(index-1)/2;
        if(heap(index).dist>=heap(parent).dist) break;
        Exchange(heap,index,parent);}
}
//#####################################################################
// Function Down_Heap
//#####################################################################
template<class T,int d> static void
Down_Heap(ARRAY<HEAP_ENTRY<T,d> >& heap,int index)
{
    while(1){
        int next=2*index+1;
        if(next>=heap.m) break;
        if(next+1<heap.m && heap(next).dist>heap(next+1).dist) next++;
        if(heap(index).dist<=heap(next).dist) break;
        Exchange(heap,index,next);}
}
//#####################################################################
// Function Down_Heap
//#####################################################################
template<class T,int d> static void
Make_Heap(ARRAY<HEAP_ENTRY<T,d> >& heap)
{
    for(int i=heap.m/2-1;i>=0;i--) Down_Heap(heap,i);
}
//#####################################################################
// Function Pop_Heap
//#####################################################################
template<class T,int d> static HEAP_ENTRY<T,d>
Pop_Heap(ARRAY<HEAP_ENTRY<T,d> >& heap)
{
    int ia=0,ie=heap.m-1;
    Exchange(heap,ia,ie);
    HEAP_ENTRY<T,d> heap_entry=heap.Pop();
    Down_Heap(heap,0);
    return heap_entry;
}
//#####################################################################
// Function Solve_Quadratic
//#####################################################################
template<class T> static T
Solve_Quadratic(const T phi,const T value_x,const T value_y,const T dx,const T dy)
{
    assert(value_x*value_y>=0);
    if(abs(value_x) >= abs(value_y)+dy) return value_y+LEVELSET_UTILITIES<T>::Sign(phi)*dy;
    if(abs(value_y) >= abs(value_x)+dx) return value_x+LEVELSET_UTILITIES<T>::Sign(phi)*dx;
    T dx2=sqr(dx),dy2=sqr(dy);
    return (dy2*value_x+dx2*value_y+LEVELSET_UTILITIES<T>::Sign(phi)*dx*dy*sqrt(dx2+dy2-sqr(value_x-value_y)))/(dx2+dy2);
}
//#####################################################################
// Function Solve_Close_Point
//#####################################################################
template<class T,int d> static T
Solve_Close_Point(const T phi,const int number_of_axis,const VECTOR<T,d>& value,const VECTOR<T,d>& dx)
{
    assert(number_of_axis);
    if(d==1 || number_of_axis==1) return value[0]+LEVELSET_UTILITIES<T>::Sign(phi)*dx[0];
    if(d==2 || number_of_axis==2) return Solve_Quadratic(phi,value[0],value[1],dx[0],dx[1]);
    assert(d==3); // candidates exist in all three directions (must be in 3d)
    T value_yz=Solve_Quadratic(phi,value[1],value[2],dx[1],dx[2]);
    if(abs(value[0]) >= abs(value_yz)) return value_yz;
    T value_xz=Solve_Quadratic(phi,value[0],value[2],dx[0],dx[2]);
    if(abs(value[1]) >= abs(value_xz)) return value_xz;
    T value_xy=Solve_Quadratic(phi,value[0],value[1],dx[0],dx[1]);
    if(abs(value[2]) >= abs(value_xy)) return value_xy;
    // use the candidates in all three directions
    T dx2=sqr(dx[0]),dy2=sqr(dx[1]),dz2=sqr(dx[2]),dx2dy2=dx2*dy2,dx2dz2=dx2*dz2,dy2dz2=dy2*dz2;
    return (dy2dz2*value[0]+dx2dz2*value[1]+dx2dy2*value[2]+LEVELSET_UTILITIES<T>::Sign(phi)*dx[0]*dx[1]*dx[2]*
        sqrt(dx2dy2+dx2dz2+dy2dz2-dx2*sqr(value[1]-value[2])-dy2*sqr(value[0]-value[2])-dz2*sqr(value[0]-value[1])))/(dx2dy2+dx2dz2+dy2dz2);
}
//#####################################################################
// Function Update_Close_Point
//##################################################################### 
// needs done=0 around the outside of the domain
template<class T,class TV,class TV_INT> static void
Update_Close_Point(const GRID<TV>& grid,ARRAY<T,TV_INT>& phi,const TV_INT& index,
    ARRAY<int,TV_INT>& close_k,
    std::function<bool(FACE_INDEX<TV_INT::m>& face)> Neighbor_Visible)
{
    TV value; // the phi value to use in the given direction
    int number_of_axis=0; // the number of axis that we want to use later
    for(int axis=0;axis<TV::m;axis++){
        T vals[2]={};
        int used=0;
        for(int s=0;s<2;s++){
            TV_INT neighbor(index);
            neighbor(axis)+=2*s-1;
            if(close_k(neighbor)!=-2) continue;
            if(Neighbor_Visible){
                FACE_INDEX<TV::m> face(axis,index);
                face.index(axis)+=s-1;
                if(!Neighbor_Visible(face)) continue;}
            vals[used++]=phi(neighbor);}
        if(used==2) vals[0]=minmag(vals[0],vals[1]);
        if(used) value[number_of_axis++]=vals[0];}
    phi(index)=Solve_Close_Point(phi(index),number_of_axis,value,grid.dX);
}
//#####################################################################
// Function Compute_Initial_Interface_Estimate
//#####################################################################
template<class TV,class TV_INT,class EP> inline auto
Compute_Initial_Interface_Estimate(const LEVELSET<TV>& levelset,
    ARRAY<int,TV_INT>& close_k,const TV_INT& index,const EP& ep)
{
    typedef typename TV::SCALAR T;
    T phi_new=2*levelset.grid.dX.Max();
    T value[3]={0}; // the phi value to use in the given direction
    int number_of_axis=0; // the number of axis that we want to use later
    for(int axis=0;axis<TV::m;axis++){
        int used=0;
        T vals[2]={};
        T dx=levelset.grid.dX(axis),phi_index=levelset.phi(index);
        for(int s=0;s<2;s++){
            TV_INT i(index);
            i(axis)+=2*s-1;
            if(close_k(i)!=-2) continue;
            T phi_i=levelset.phi(i);
            if(LEVELSET_UTILITIES<T>::Interface(phi_index,phi_i))
                vals[used++]=LEVELSET_UTILITIES<T>::Theta(phi_index,phi_i);}
        if(used==2) vals[0]=min(vals[0],vals[1]);
        if(used) value[number_of_axis++]=vals[0]*dx;}
    if(number_of_axis==1) phi_new=value[0];
    else if(number_of_axis==2){
        if(T d2=sqr(value[0])+sqr(value[1])) phi_new=value[0]*value[1]/sqrt(d2);
        else phi_new=0;}
    else{PHYSBAM_ASSERT(TV::m==3); // 2d should never get to this point
        T value_xy=value[0]*value[1],value_xz=value[0]*value[2],value_yz=value[1]*value[2],d2=sqr(value_xy)+sqr(value_xz)+sqr(value_yz);
        if(d2) phi_new=value_xy*value[2]/sqrt(d2);
        else phi_new=min(value[0],value[1],value[2]);}

    if(ep.refine_with_iterative_solver){
        T iterative_tolerance_absolute=ep.iterative_tolerance*levelset.grid.dX.Min();
        TV location=levelset.grid.Center(index),vertex=location;
        T phi_value=levelset.phi(index);
        for(int iterations=0;iterations<ep.iterations;iterations++){
            if(abs(phi_value)>10*levelset.grid.dX.Min())break; // stop if it looks like it's blowing up
            vertex-=phi_value*levelset.Normal(vertex);
            phi_value=levelset.Phi(vertex);
            if(abs(phi_value)<=iterative_tolerance_absolute){
                T new_phi_value=(vertex-location).Magnitude();
                if((new_phi_value-phi_new)/max(iterative_tolerance_absolute,phi_new)<ep.iterative_drift_fraction && new_phi_value>0)
                    phi_new=new_phi_value;
                break;}}}
    return phi_new*LEVELSET_UTILITIES<T>::Sign(levelset.phi(index));
}

//#####################################################################
// Function Add_To_Initial
//##################################################################### 
template<class T,class TV_INT> static void
Add_To_Initial(const ARRAY<T,TV_INT>& phi,ARRAY<int,TV_INT>& close_k,
    ARRAY<HEAP_ENTRY<T,TV_INT::m> >& heap,const TV_INT& index,
    std::function<bool(FACE_INDEX<TV_INT::m>& face)> Neighbor_Visible,
    int process_sign)
{
    // add neighbors to close if not done 
    for(int a=0;a<TV_INT::m;a++){
        for(int s=0;s<2;s++){
            TV_INT neighbor(index);
            neighbor(a)+=2*s-1;
            int& k=close_k(neighbor);
            if(k!=-1) continue;
            if(process_sign && process_sign*phi(neighbor)<0) continue;
            if(Neighbor_Visible){
                FACE_INDEX<TV_INT::m> face(a,index);
                face.index(a)+=s-1;
                if(!Neighbor_Visible(face)) continue;}
            k=heap.Append({neighbor,&k});}}
}
//#####################################################################
// Function Initialize_Interface
//#####################################################################
// pass heap_length by reference
template<class T,class TV,class TV_INT,class EP> inline void
Initialize_Interface(GRID<TV>& grid,ARRAY<T,TV_INT>& phi,ARRAY<int,TV_INT>& close_k,
    ARRAY<HEAP_ENTRY<T,TV_INT::m> >& heap,
    std::function<bool(FACE_INDEX<TV_INT::m>& face)> Neighbor_Visible,
    const ARRAY<TV_INT>* seed_indices,
    const bool add_seed_indices_for_ghost_cells,int ghost_cells,
    const EP& ep,int process_sign,bool correct_interface_phi)
{
    ARRAY<TV_INT> seeds;
    // Use provided seeds
    if(seed_indices){
        seeds=*seed_indices;
        for(int i=0;i<seeds.m;i++)
            close_k(seeds(i))=-2;}

    // Gather additional seeds; done check interior if seeds provided.
    if(!seed_indices || add_seed_indices_for_ghost_cells){
        for(FACE_RANGE_ITERATOR<TV::m> it(grid.Domain_Indices(),ghost_cells,
                0,seed_indices?RF::ghost|RF::skip_outer:RF::skip_outer);it.Valid();it.Next()){
            TV_INT a[2]={it.face.First_Cell_Index(),it.face.Second_Cell_Index()};
            if(Neighbor_Visible && !Neighbor_Visible(it.face)) continue;
            if(!LEVELSET_UTILITIES<T>::Interface(phi(a[0]),phi(a[1]))) continue;
            for(int i=0;i<2;i++)
                if(close_k(a[i])==-1){
                    close_k(a[i])=-2;
                    if(process_sign && process_sign*phi(a[i])<0) continue;
                    seeds.Append(a[i]);}}}

    for(int i=0;i<seeds.m;i++)
        Add_To_Initial(phi,close_k,heap,seeds(i),Neighbor_Visible,process_sign);

    if(correct_interface_phi){
        ARRAY<T> phi_new(seeds.m);
        LEVELSET<TV> levelset(grid,phi);
        for(int i=0;i<seeds.m;i++)
            phi_new(i)=Compute_Initial_Interface_Estimate(levelset,close_k,seeds(i),ep);
        for(int i=0;i<seeds.m;i++)
            phi(seeds(i))=phi_new(i);}

    for(int i=0;i<heap.m;i++){
        PHYSBAM_ASSERT(close_k(heap(i).index)>=0);
        Update_Close_Point(grid,phi,heap(i).index,close_k,Neighbor_Visible);
        heap(i).dist=abs(phi(heap(i).index));}
    Make_Heap(heap);
}
//#####################################################################
// Function Fast_Marching_Method
//#####################################################################
// Updates phi only out ghost cells.
template<class TV> void FAST_MARCHING_METHOD_UNIFORM<TV>::
Fast_Marching_Method(const GRID<TV>& grid_input,int ghost,ARRAY<T,TV_INT>& phi,const T stopping_distance)
{
    GRID<TV> grid(grid_input.Is_MAC_Grid()?grid_input:grid_input.Get_MAC_Grid_At_Regular_Positions());
    PHYSBAM_ASSERT(grid.Is_MAC_Grid());
    RANGE<TV_INT> domain=grid.Domain_Indices(ghost);

    // -1=unvisited, -2=done, -3=sentinal, (>=0)=index in heap
    ARRAY<int,TV_INT> close_k(domain.Thickened(1),true,-1);

    // Pad with a band of invalid cells; this avoids the need to bounds check.
    for(int a=0;a<TV::m;a++)
        for(RANGE_ITERATOR<TV::m-1> it(close_k.domain.Remove_Dimension(a));it.Valid();it.Next()){
            close_k(it.index.Insert(close_k.domain.min_corner(a),a))=-3;
            close_k(it.index.Insert(close_k.domain.max_corner(a)-1,a))=-3;}

    ARRAY<HEAP_ENTRY<T,TV::m> > heap;
    Initialize_Interface(grid,phi,close_k,heap,Neighbor_Visible,seed_indices,
        add_seed_indices_for_ghost_cells,ghost,estimation_parameters,process_sign,
        correct_interface_phi);

    while(heap.m!=0){
        HEAP_ENTRY<T,TV::m> heap_entry=Pop_Heap(heap);
        if(stopping_distance && heap_entry.dist>stopping_distance){ // exit early
            for(CELL_ITERATOR<TV> it(grid,domain);it.Valid();it.Next())
                if(close_k(it.index)>-2)
                    phi(it.index)=LEVELSET_UTILITIES<T>::Sign(phi(it.index))*stopping_distance;
            break;}
        *heap_entry.close_ptr=-2; // flag as done

        for(int axis=0;axis<TV::m;axis++)
            for(int s=0;s<2;s++){
                TV_INT neighbor(heap_entry.index);
                neighbor(axis)+=2*s-1;
                if(close_k(neighbor)<=-2) continue;
                if(Neighbor_Visible){
                    FACE_INDEX<TV::m> face(axis,heap_entry.index);
                    face.index(axis)+=s-1;
                    if(!Neighbor_Visible(face)) continue;}
                Update_Close_Point(grid,phi,neighbor,close_k,Neighbor_Visible);
                int& k=close_k(neighbor);
                if(k<0) k=heap.Append({neighbor,&k,abs(phi(neighbor))}); // add to close
                else heap(k).dist=abs(phi(neighbor));
                Up_Heap(heap,k);}}
}
//#####################################################################
namespace PhysBAM{
template class FAST_MARCHING_METHOD_UNIFORM<VECTOR<float,1> >;
template class FAST_MARCHING_METHOD_UNIFORM<VECTOR<float,2> >;
template class FAST_MARCHING_METHOD_UNIFORM<VECTOR<float,3> >;
template class FAST_MARCHING_METHOD_UNIFORM<VECTOR<double,1> >;
template class FAST_MARCHING_METHOD_UNIFORM<VECTOR<double,2> >;
template class FAST_MARCHING_METHOD_UNIFORM<VECTOR<double,3> >;
}
