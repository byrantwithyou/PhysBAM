//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include "MPM_PARTICLES.h"
#include "MPM_PARTICLES_FORWARD.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_PARTICLES<TV>::
MPM_PARTICLES()
{
    Store_Velocity();
    Store_Mass();
    Add_Array(ATTRIBUTE_ID_MATERIAL_COORDINATE,&Xm);
    Add_Array(ATTRIBUTE_ID_VOLUME,&volume);
    Add_Array(ATTRIBUTE_ID_FE,&Fe);
    Add_Array(ATTRIBUTE_ID_FP,&Fp);
    rand_generator.Set_Seed(0);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_PARTICLES<TV>::
~MPM_PARTICLES()
{}
//#####################################################################
// Function Initialize_X_As_A_Grid
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Initialize_X_As_A_Grid(const VECTOR<int,TV::m>& count,const RANGE<TV>& box)
{
    GRID<TV> grid(count,box);
    ARRAY<TV> sample_X;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+count));it.Valid();it.Next()){
        TV x=grid.X(it.index);
        sample_X.Append(x);}
    Resize(sample_X.m);
    X=sample_X;
    Xm=X;
}
//#####################################################################
// Function Initialize_X_As_A_Randomly_Sampled_Box
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Initialize_X_As_A_Randomly_Sampled_Box(const int N,const RANGE<TV>& box,const T exclude_radius)
{
    
    ARRAY<TV> sample_X;
    for(int n=0;n<N;n++){
        TV x;
        rand_generator.Fill_Uniform(x,box);
        sample_X.Append(x);}
    ARRAY<bool> should_go_away(sample_X.m);
    if(exclude_radius<(T)100){
        T exclude_radius2=exclude_radius*exclude_radius;
        HASHTABLE<TV_INT,ARRAY<int> > buckets;
        for(int i=0;i<sample_X.m;i++)
            buckets.Get_Or_Insert(TV_INT(floor(sample_X(i)/exclude_radius))).Append(i);
        for(int i=0;i<sample_X.m;i++){
            if(!should_go_away(i)){
                TV_INT my_cell=TV_INT(floor(sample_X(i)/exclude_radius));
                TV_INT aaa=TV_INT(floor(sample_X(i)/exclude_radius)-1.0);
                TV_INT bbb=TV_INT(floor(sample_X(i)/exclude_radius)+2.0);
                for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(aaa,bbb));it.Valid();it.Next()){
                     if(const ARRAY<int>* l1=buckets.Get_Pointer(it.index)){
                         for(int j=0;j<l1->m;j++){
                             if((*l1)(j)!=i && !should_go_away((*l1)(j)) && (sample_X(i)-sample_X((*l1)(j))).Magnitude_Squared()<exclude_radius2){
                                 should_go_away((*l1)(j))=true;}}}}}}}
    ARRAY<TV> filted_X;
    for(int i=0;i<sample_X.m;i++)
        if(!should_go_away(i)) filted_X.Append(sample_X(i));
    Resize(filted_X.m);
    X=filted_X;
    Xm=X;
}
//#####################################################################
// Function Add_Randomly_Sampled_Object
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Add_Randomly_Sampled_Implicit_Object(const IMPLICIT_OBJECT<TV>& object,const T exclude_radius)
{
    const_cast<IMPLICIT_OBJECT<TV>&>(object).Update_Box();
    RANGE<TV> bounding_box(object.box);
    TV_INT cells(ceil(bounding_box.Edge_Lengths()/exclude_radius));
    bounding_box.Scale_About_Center(exclude_radius*TV(cells)/bounding_box.Edge_Lengths());
    GRID<TV> grid(cells,bounding_box,true);
    ARRAY<bool,TV_INT> ok(grid.Domain_Indices(1));
    ARRAY<TV_INT> todo;
    for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid);it.Valid();it.Next())
        if(object.Extended_Phi(it.Location())<=0){
            ok(it.index)=true;
            todo.Append(it.index);}
    rand_generator.Random_Shuffle(todo);
    RANGE<TV_INT> neigh_index(TV_INT()-1,TV_INT()+2);
    for(int i=0;i<todo.m;i++){
        TV_INT index=todo(i);
        if(!ok(index)) continue;
        TV new_X;
        for(int j=0;j<1000;j++){
            new_X=rand_generator.Get_Uniform_Vector(grid.Cell_Domain(index));
            if(object.Extended_Phi(new_X)<=0) break;}
        for(RANGE_ITERATOR<TV::m> it(neigh_index);it.Valid();it.Next())
            ok(it.index+index)=false;
        int p=this->Add_Element();
        X(p)=new_X;
        Xm(p)=new_X;}
}
//#####################################################################
// Function Initialize_X_As_A_Ball
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Initialize_X_As_A_Ball(const VECTOR<int,TV::m>& count,const RANGE<TV>& square_box)
{
    GRID<TV> grid(count,square_box);
    TV center=square_box.Center();
    T r=0.5*((-square_box.min_corner+square_box.max_corner)(0));
    ARRAY<TV> sample_X;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+count));it.Valid();it.Next()){
        TV x=grid.X(it.index);
        if((x-center).Magnitude()<=r) sample_X.Append(x);}
    Resize(sample_X.m);
    X=sample_X;
    Xm=X;
}
//#####################################################################
// Function Add_X_As_A_Grid
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Add_X_As_A_Grid(const VECTOR<int,TV::m>& count,const RANGE<TV>& box)
{
    GRID<TV> grid(count,box);
    ARRAY<TV> old_X;
    old_X.Resize(X.m);
    for(int i=0;i<X.m;i++) old_X(i)=X(i);
    ARRAY<TV> sample_X;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+count));it.Valid();it.Next()){
        TV x=grid.X(it.index);
        sample_X.Append(x);}
    old_X.Append_Elements(sample_X);
    Resize(old_X.m);
    X=old_X;
    Xm=X;
}
//#####################################################################
// Function Add_X_As_A_Randomly_Sampled_Box
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Add_X_As_A_Randomly_Sampled_Box(const int N,const RANGE<TV>& box)
{
    ARRAY<TV> sample_X;
    for(int n=0;n<N;n++){
        TV x;
        rand_generator.Fill_Uniform(x,box);
        sample_X.Append(x);}
    ARRAY<TV> old_X;
    old_X.Resize(X.m);
    for(int i=0;i<X.m;i++) old_X(i)=X(i);
    old_X.Append_Elements(sample_X);
    Resize(old_X.m);
    X=old_X;
    Xm=X;
}
//#####################################################################
// Function Add_X_As_A_Ball
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Add_X_As_A_Ball(const VECTOR<int,TV::m>& count,const RANGE<TV>& square_box)
{
    GRID<TV> grid(count,square_box);
    TV center=square_box.Center();
    T r=0.5*((-square_box.min_corner+square_box.max_corner)(0));
    ARRAY<TV> old_X;
    old_X.Resize(X.m);
    for(int i=0;i<X.m;i++) old_X(i)=X(i);
    ARRAY<TV> sample_X;
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),TV_INT()+count));it.Valid();it.Next()){
        TV x=grid.X(it.index);
        if((x-center).Magnitude()<=r) sample_X.Append(x);}
    old_X.Append_Elements(sample_X);
    Resize(old_X.m);
    X=old_X;
    Xm=X;
}
//#####################################################################
// Function Reduce_X_As_A_Box
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Reduce_X_In_A_Box(const RANGE<TV>& box)
{
    ARRAY<TV> old_X;
    old_X.Resize(X.m);
    for(int i=0;i<X.m;i++) old_X(i)=X(i);
    ARRAY<TV> sample_X;
    for(int i=0;i<old_X.m;i++)
        if(!box.Lazy_Inside(old_X(i)))
            sample_X.Append(old_X(i));
    Resize(sample_X.m);
    X=sample_X;
    Xm=X;
}
//#####################################################################
// Function Reduce_X_As_A_Ball
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Reduce_X_As_A_Ball(const RANGE<TV>& square_box)
{
    ARRAY<TV> old_X;
    old_X.Resize(X.m);
    for(int i=0;i<X.m;i++) old_X(i)=X(i);
    TV center=square_box.Center();
    T r=0.5*((-square_box.min_corner+square_box.max_corner)(0));
    ARRAY<TV> sample_X;
    for(int i=0;i<old_X.m;i++)
        if((old_X(i)-center).Magnitude()>r)
            sample_X.Append(old_X(i));
    Resize(sample_X.m);
    X=sample_X;
    Xm=X;
}
//#####################################################################
// Function Reduce_X_Where_Not_In_A_Ball
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Reduce_X_Where_Not_In_A_Ball(const SPHERE<TV>& ball)
{
    ARRAY<TV> old_X;
    old_X.Resize(X.m);
    for(int i=0;i<X.m;i++) old_X(i)=X(i);
    ARRAY<TV> sample_X;
    for(int i=0;i<old_X.m;i++)
        if(ball.Lazy_Inside(old_X(i)))
            sample_X.Append(old_X(i));
    Resize(sample_X.m);
    X=sample_X;
    Xm=X;
}
//#####################################################################
// Function Reduce_X_Where_Not_In_A_Ball_But_In_A_Box
//#####################################################################
template<class TV> void MPM_PARTICLES<TV>::
Reduce_X_Where_Not_In_A_Ball_But_In_A_Box(const SPHERE<TV>& ball,const RANGE<TV>& box)
{
    ARRAY<TV> old_X;
    old_X.Resize(X.m);
    for(int i=0;i<X.m;i++) old_X(i)=X(i);
    ARRAY<TV> sample_X;
    for(int i=0;i<old_X.m;i++){
        if(!ball.Lazy_Inside(old_X(i)) && box.Lazy_Inside(old_X(i))) continue;
        sample_X.Append(old_X(i));}
    Resize(sample_X.m);
    X=sample_X;
    Xm=X;
}
static int Initialize_MPM_Particles()
{
    Register_Attribute_Name(ATTRIBUTE_ID_MATERIAL_COORDINATE,"Xm");
    Register_Attribute_Name(ATTRIBUTE_ID_VOLUME,"volume");
    Register_Attribute_Name(ATTRIBUTE_ID_FE,"Fe");
    Register_Attribute_Name(ATTRIBUTE_ID_FP,"Fp");
    return 0;
}
int initialize_mpm_particles=Initialize_MPM_Particles();
//#####################################################################
template class MPM_PARTICLES<VECTOR<float,2> >;
template class MPM_PARTICLES<VECTOR<float,3> >;
template class MPM_PARTICLES<VECTOR<double,2> >;
template class MPM_PARTICLES<VECTOR<double,3> >;
}
