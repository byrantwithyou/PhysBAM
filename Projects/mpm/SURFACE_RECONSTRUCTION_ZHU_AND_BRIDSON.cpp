//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX.h>
#include "SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON.h"
namespace PhysBAM{
//#####################################################################
// Function - kernel stuff
//#####################################################################
template<class T> T
kernel(const T s_squared)
{
    return max((T)0,cube((T)1-s_squared));
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON<TV>::
Initialize(const T average_particle_spacing)
{
    particle_radii=average_particle_spacing;
    radius_of_neighborhood=(T)2*average_particle_spacing;
}
//#####################################################################
// Function Build_Scalar_Field
//#####################################################################
template<class TV> void SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON<TV>::
Build_Scalar_Field(const ARRAY_VIEW<TV>& X,const GRID<TV>& grid,ARRAY<T,TV_INT>& phi) const
{
    T R2=sqr(radius_of_neighborhood);
    T one_over_R2=Inverse(R2);
    HASHTABLE<TV_INT,ARRAY<int> > buckets;
    for(int i=0;i<X.m;i++)
        buckets.Get_Or_Insert(TV_INT(floor(X(i)/radius_of_neighborhood))).Append(i);
    phi.Resize(RANGE<TV_INT>(TV_INT(),grid.counts));phi.Fill(T(0.01));     
    for(RANGE_ITERATOR<TV::m> it(RANGE<TV_INT>(TV_INT(),grid.counts));it.Valid();it.Next()){
        TV x=grid.Node(it.index);
        ARRAY<int> my_neighbors;
        ARRAY<T> numerators;
        TV_INT mycell=TV_INT(floor(x/radius_of_neighborhood));
        TV_INT aaa=TV_INT(floor(x/radius_of_neighborhood)-1.0);
        TV_INT bbb=TV_INT(floor(x/radius_of_neighborhood)+2.0);
        for(RANGE_ITERATOR<TV::m> it2(RANGE<TV_INT>(aaa,bbb));it2.Valid();it2.Next()){
            if(const ARRAY<int>* l1=buckets.Get_Pointer(it2.index)){
                for(int j=0;j<l1->m;j++){
                    int that_particle=(*l1)(j);
                    T distance2=(x-X(that_particle)).Magnitude_Squared();
                    if(distance2<=R2){
                        my_neighbors.Append(that_particle);
                        numerators.Append(kernel(distance2*one_over_R2));}}}}
        if(my_neighbors.m!=0){
            ARRAY<T> w(my_neighbors.m);
            T one_over_denominator;
            if(numerators.Sum()<(T)1e-8) one_over_denominator=T(0);
            else one_over_denominator=Inverse(numerators.Sum());
            for(int i=0;i<my_neighbors.m;i++) w(i)=numerators(i)*one_over_denominator;
            TV x_average=TV();
            T r_average=T(0);
            for(int i=0;i<my_neighbors.m;i++){
                x_average+=w(i)*X(my_neighbors(i));
                r_average+=w(i)*particle_radii;}
            phi(it.index)=(x-x_average).Magnitude()-r_average;}}
}
//#####################################################################
template class SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON<VECTOR<float,2> >;
template class SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON<VECTOR<double,2> >;
template class SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON<VECTOR<float,3> >;
template class SURFACE_RECONSTRUCTION_ZHU_AND_BRIDSON<VECTOR<double,3> >;
}
