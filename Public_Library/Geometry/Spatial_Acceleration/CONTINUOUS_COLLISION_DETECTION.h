//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Andrew Selle, Joseph Teran, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __CONTINUOUS_COLLISION_DETECTION__
#define __CONTINUOUS_COLLISION_DETECTION__
#include <Core/Arrays/ARRAY.h>
#include <Core/Arrays/ARRAY_VIEW.h>
#include <Core/Data_Structures/TRIPLE.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
namespace PhysBAM{
template<class TV> class STRUCTURE;
template<int d>
struct CCD_PAIR
{
    int s0,s1;
    VECTOR<int,d+1> f; // pf: 0 1 1 1, ee: 0 0 1 1
    int e0,e1; // element index for face or edge; particle index for particle.
};
template<class TV,int d>
struct CCD_DATA
{
    BOX_HIERARCHY<TV> h;
    ARRAY<TRIPLE<int,VECTOR<int,d>,int> > p; // object index, simplex vertex indices, element index
};
template<class TV>
class CONTINUOUS_COLLISION_DETECTION
{
    typedef typename TV::SCALAR T;
public:
    CCD_DATA<TV,1> c1;
    CCD_DATA<TV,2> c2;
    CCD_DATA<TV,3> c3; // Unused in 2D

    struct UPDATE_DATA
    {
        int a[3],b[3]; // begin,end of data in hierarchies.
        UPDATE_DATA(): a{},b{} {}
    };
    
    ARRAY<UPDATE_DATA> update_data;
    bool stale_leaves=true;
    bool update_hierarchy=true;
    
    void Compute_Pairs_PF(ARRAY<CCD_PAIR<TV::m> >& point_face,T thickness=0)
    {Compute_Pairs_PF(point_face,*this,thickness);}
    void Compute_Pairs_EE(ARRAY<CCD_PAIR<TV::m> >& edge_edge,T thickness=0)
    {Compute_Pairs_EE(edge_edge,*this,thickness);}

    void Compute_Pairs_PF(ARRAY<CCD_PAIR<TV::m> >& point_face,
        CONTINUOUS_COLLISION_DETECTION<TV>& ccd,T thickness=0);
    void Compute_Pairs_EE(ARRAY<CCD_PAIR<TV::m> >& edge_edge,
        CONTINUOUS_COLLISION_DETECTION<TV>& ccd,T thickness=0);

    void Update_Positions(ARRAY_VIEW<const TV> X0,ARRAY_VIEW<const TV> X1);
    void Update_Positions(ARRAY_VIEW<const TV> X0,ARRAY_VIEW<const TV> X1,int b);
    template<int d>
    void Update_Positions(ARRAY_VIEW<const TV> X0,ARRAY_VIEW<const TV> X1,
        CCD_DATA<TV,d>& cd);
    template<int d>
    void Update_Positions(ARRAY_VIEW<const TV> X0,ARRAY_VIEW<const TV> X1,int b,
        CCD_DATA<TV,d>& cd);
    void Update_Positions(ARRAY_VIEW<const TV> X,const FRAME<TV>& f0,
        const FRAME<TV>& f1,int b);
    template<int d>
    void Update_Positions(ARRAY_VIEW<const TV> X,const FRAME<TV>& f0,
        const FRAME<TV>& f1,int b,CCD_DATA<TV,d>& cd);
    void Update_Positions(ARRAY_VIEW<const TV> X,int b);
    template<int d>
    void Update_Positions(ARRAY_VIEW<const TV> X,int b,CCD_DATA<TV,d>& cd);

    // flags: 0=points, 1=edges, 2=faces
    int Add_Structure(STRUCTURE<TV>* s,int flags);
    int Add_Particles(int n);

    void Update_Hierarchy();
    void Update_Boxes();
    void Clean_Memory();
};
}
#endif

