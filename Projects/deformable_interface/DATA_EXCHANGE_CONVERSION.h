//#####################################################################
// Copyright 2011, Craig Schroeder.
//#####################################################################
#ifndef __DATA_EXCHANGE_CONVERSION__
#define __DATA_EXCHANGE_CONVERSION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <vector>
#include "deformable_body.h"
using namespace PhysBAM;
template<class T>
void To_Pb(T& o, const T& i){o=i;}

void To_Pb(double& o, float i){o=i;}

template<class T,class U,int d>
void To_Pb(VECTOR<T,d>& o,const data_exchange::fixed_vector<U,d>& v)
{
    for(size_t i=0; i<d; i++) To_Pb(o(i+1), v.data[i]);
}

template<class T,int d>
VECTOR<T,d> To_Pb(const data_exchange::fixed_vector<T,d>& v)
{
    VECTOR<T,d> x;
    To_Pb(x, v);
    return x;
}

template<class T,class U>
void To_Pb(ARRAY_VIEW<T> view,const std::vector<U>& v)
{
    for(size_t i=0; i<v.size(); i++) To_Pb(view(i+1), v[i]);
}

template<class T,class U>
void To_Pb(ARRAY<T>& array,const std::vector<U>& v)
{
    array.Resize(v.size());
    To_Pb(array, v);
}

void Triangle_Mesh_From_Data_Exchange(TRIANGLE_MESH& mesh, const data_exchange::polygon_mesh& poly, ARRAY<int>* parent_index)
{
    for(size_t i=0, k=0; i<poly.polygon_counts.size(); i++){
        int n=poly.polygon_counts[i];
        for(int j=2; j<n; j++){
            mesh.elements.Append(VECTOR<int,3>(poly.polygons[k]+1, poly.polygons[k+j-1]+1, poly.polygons[k+j]+1));
            if(parent_index) parent_index->Append(i);}
        k+=n;}
}

template<class TV>
void Triangulated_Surface_From_Data_Exchange(TRIANGULATED_SURFACE<TV>& surface, const data_exchange::polygon_mesh& poly, const std::vector<data_exchange::vf3>& pos, ARRAY<int>* parent_index)
{
    Triangle_Mesh_From_Data_Exchange(surface.mesh, poly, parent_index);
    surface.particles.array_collection->Add_Elements(pos.size());
    To_Pb(surface.particles.X, pos);
    surface.Update_Number_Nodes();
}

template<class T>
void From_Pb(T& o, const T& i){o=i;}

void From_Pb(float& o, double i){o=i;}

template<class T,class U,int d>
void From_Pb(data_exchange::fixed_vector<U,d>& v, const VECTOR<T,d>& o)
{
    for(size_t i=0; i<d; i++) From_Pb(v.data[i], o(i+1));
}

template<class T,class U,int d>
data_exchange::fixed_vector<U,d>  From_Pb(const VECTOR<T,d>& v)
{
    data_exchange::fixed_vector<T,d> x;
    From_Pb(x, v);
    return x;
}

template<class T,class T_ARRAY,class U>
void From_Pb(std::vector<U>& v, const ARRAY_BASE<T,T_ARRAY>& array)
{
    v.resize(array.Size());
    for(size_t i=0; i<v.size(); i++) From_Pb(v[i], array(i+1));
}

#endif
