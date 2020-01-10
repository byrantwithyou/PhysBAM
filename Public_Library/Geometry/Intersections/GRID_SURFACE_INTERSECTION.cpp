//#####################################################################
// Copyright 2020, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Arrays_Nd/ARRAYS_ND.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Math_Tools/clamp.h>
#include <Grid_Tools/Grids/GRID.h>
#include <Geometry/Geometry_Particles/VIEWER_OUTPUT.h>
#include <Geometry/Intersections/GRID_SURFACE_INTERSECTION.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{

// Node grid, return cut edges.
template<class T,class TV>
void Grid_Surface_Intersection(
    HASHTABLE<EDGE_INDEX<2>,GRID_SURFACE_INTERSECTION_DATA<TV> >& hash,
    const GRID<TV>& grid,const SEGMENTED_CURVE_2D<T>& surface,bool compute_inside)
{
    typedef VECTOR<int,TV::m> TV_INT;
    ARRAY<TV> X(surface.particles.X.m);
    for(int i=0;i<X.m;i++)
        X(i)=(surface.particles.X(i)-grid.domain.min_corner)*grid.one_over_dX;

    for(int e=0;e<surface.mesh.elements.m;e++)
    {
        const TV& X0=X(surface.mesh.elements(e).x),&X1=X(surface.mesh.elements(e).y);
        TV u=X1-X0;
        TV_INT i0(floor(X0)),i1(floor(X1)),off;
        for(int i=0;i<2;i++) off(i)=i1(i)>=i0(i);
        TV_INT dir=off*2-1;
        bool canon_sgn=u.y>=0;
        T last_b = 0;

        for(TV_INT i=i0;i!=i1;)
        {
            TV Z(i+off);
            T cp=(Z-X0).Cross(u).x;
            bool sgn=cp>=0;
            int axis=i.y!=i1.y && (i.x==i1.x || (sgn==canon_sgn)==off.x);
            T b1=clamp((Z(axis)-X0(axis))/u(axis),last_b,(T)1);
            last_b=b1;
            T theta=clamp((b1*u+X0)(1-axis)-i(1-axis),(T)0,(T)1);
            EDGE_INDEX<2> edge(1-axis,i.Add_Axis(axis,off(axis)));
            auto& data=hash.Get_Or_Insert(edge);
            data.cut_elements.Append({e,theta,TV(1-b1,b1)});
            i(axis)+=dir(axis);
        }
    }
    if(compute_inside)
    {
        VECTOR<HASHTABLE<int,ARRAY<PAIR<int,GRID_SURFACE_INTERSECTION_DATA<TV>*> > >,2> scan_hash;
        for(auto& h:hash)
        {
            int a=h.key.axis;
            TV_INT i=h.key.index;
            scan_hash(a).Get_Or_Insert(i(1-a)).Append({i(a),&h.data});
        }
        for(auto& h:scan_hash)
            for(auto& A:h)
            {
                A.data.Sort([](const auto&a,const auto&b){return a.x<b.x;});
                bool in=false;
                for(const auto& p:A.data)
                {
                    p.y->in[0]=in;
                    in^=p.y->cut_elements.m%2;
                    p.y->in[1]=in;
                }
                assert(!in);
            }
    }
}

template<class T,class TV>
void Grid_Surface_Intersection(
    HASHTABLE<EDGE_INDEX<3>,GRID_SURFACE_INTERSECTION_DATA<TV> >& hash,
    const GRID<TV>& grid,const TRIANGULATED_SURFACE<T>& surface,bool compute_inside)
{
    PHYSBAM_NOT_IMPLEMENTED();
}
template void Grid_Surface_Intersection<double,VECTOR<double,2> >(HASHTABLE<EDGE_INDEX<2>,GRID_SURFACE_INTERSECTION_DATA<VECTOR<double,2> > >&,GRID<VECTOR<double,2> > const&,SEGMENTED_CURVE_2D<double> const&,bool);
template void Grid_Surface_Intersection<float,VECTOR<float,2> >(HASHTABLE<EDGE_INDEX<2>,GRID_SURFACE_INTERSECTION_DATA<VECTOR<float,2> > >&,GRID<VECTOR<float,2> > const&,SEGMENTED_CURVE_2D<float> const&,bool);
}
