//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays/ARRAY.h>
#include <Core/Log/LOG.h>
#include <Core/Math_Tools/RANGE.h>
#include <OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <OpenGL/OpenGL/OPENGL_GRID_BASED_VECTOR_FIELD_3D.h>
#include <OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
using namespace PhysBAM;
//#####################################################################
// OPENGL_GRID_BASED_VECTOR_FIELD_3D
//#####################################################################
template<class T> OPENGL_GRID_BASED_VECTOR_FIELD_3D<T>::
OPENGL_GRID_BASED_VECTOR_FIELD_3D(STREAM_TYPE stream_type,GRID<TV> &grid, ARRAY<VECTOR<T,3> ,VECTOR<int,3> > &V)
    :OPENGL_VECTOR_FIELD_3D<T>(stream_type,*(new ARRAY<VECTOR<T,3> >),*(new ARRAY<VECTOR<T,3> >),OPENGL_COLOR::Gray(.8f),.25f,true,false,false),grid(grid),V(V)
{
    max_vectors_3d = 100000;
}
//#####################################################################
// ~OPENGL_GRID_BASED_VECTOR_FIELD_3D
//#####################################################################
template<class T> OPENGL_GRID_BASED_VECTOR_FIELD_3D<T>::
~OPENGL_GRID_BASED_VECTOR_FIELD_3D()
{
    delete &vector_field;
    delete &vector_locations;
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > OPENGL_GRID_BASED_VECTOR_FIELD_3D<T>::
Bounding_Box() const
{
    return World_Space_Box(grid.domain);
}
//#####################################################################
// Update
//#####################################################################
template<class T> void OPENGL_GRID_BASED_VECTOR_FIELD_3D<T>::
Update()
{
    vector_field.Resize(0);
    vector_locations.Resize(0);

    if(V.domain.Empty_Half_Open()) return;

    OPENGL_UNIFORM_SLICE<T>* slice=(OPENGL_UNIFORM_SLICE<T>*)this->slice;
    if(!slice || slice->mode==OPENGL_SLICE::NO_SLICE){
        int count=0;int index;
        for(index=0;index<V.array.Size();index++)if(V.array(index)!=VECTOR<T,3>()) count++;
        int num_vectors=min(count,max_vectors_3d);
        vector_field.Resize(num_vectors);
        vector_locations.Resize(num_vectors);
        int m=grid.counts.x,n=grid.counts.y,mn=grid.counts.z;index=0;
        if(count>max_vectors_3d){
            int inc=max(1,(int)ceil(pow((T)count/max_vectors_3d,((T)1/3))));
            m/=inc;n/=inc;mn/=inc;
            for(int i=0;i<m;i++)for(int j=0;j<n;j++)for(int ij=0;ij<mn;ij++){
                VECTOR<int,3> grid_index(i*inc,j*inc,ij*inc);
                if(V(grid_index)!=VECTOR<T,3>()){
                    index++;if(index>num_vectors)break;
                    vector_field(index)=V(grid_index);
                    vector_locations(index)=grid.X(grid_index);}}}
        else for(int i=0;i<m;i++)for(int j=0;j<n;j++)for(int ij=0;ij<mn;ij++)if(V(i,j,ij)!=VECTOR<T,3>()){
            vector_field(index)=V(i,j,ij);vector_locations(index++)=grid.X(TV_INT(i,j,ij));}}
    else{
        VECTOR<int,3> domain_start(V.domain.min_corner.x,V.domain.min_corner.y,V.domain.min_corner.z),domain_end(V.domain.max_corner.x,V.domain.max_corner.y,V.domain.max_corner.z);
        if((slice->mode==OPENGL_SLICE::CELL_SLICE && (!grid.Is_MAC_Grid() || slice->index<domain_start[slice->axis] || slice->index>=domain_end[slice->axis])) ||
            (slice->mode==OPENGL_SLICE::NODE_SLICE && (grid.Is_MAC_Grid() || slice->index<domain_start[slice->axis] || slice->index>=domain_end[slice->axis]))){
            // Currently we don't draw anything if the slice doesn't match where the vector field lives
            return;}

        RANGE<TV_INT> domain=V.domain;
        domain.min_corner(slice->axis)=slice->index;
        domain.max_corner(slice->axis)=slice->index+1;

        int idx=0;
        vector_field.Resize(domain.Size());
        vector_locations.Resize(domain.Size());
        for(RANGE_ITERATOR<TV::m> it(domain);it.Valid();it.Next()){
            vector_field(idx)=V(it.index);
            vector_locations(idx)=grid.X(it.index);
            idx++;}}
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_GRID_BASED_VECTOR_FIELD_3D<T>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION<T>* current_selection) const
{
    // TODO: interpolate to particles
    if(current_selection && current_selection->type==OPENGL_SELECTION<T>::GRID_NODE_3D && !grid.Is_MAC_Grid()){
        VECTOR<int,3> index=((OPENGL_SELECTION_GRID_NODE_3D<T>*)current_selection)->index;
        if(V.Valid_Index(index)) stream<<V(index);}
    if(current_selection && current_selection->type==OPENGL_SELECTION<T>::GRID_CELL_3D && grid.Is_MAC_Grid()){
        VECTOR<int,3> index=((OPENGL_SELECTION_GRID_CELL_3D<T>*)current_selection)->index;
        if(V.Valid_Index(index)) stream<<V(index);}
    stream<<std::endl;
}
//#####################################################################
namespace PhysBAM{
template class OPENGL_GRID_BASED_VECTOR_FIELD_3D<float>;
template class OPENGL_GRID_BASED_VECTOR_FIELD_3D<double>;
}
