//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Andrew Selle, Joseph Teran, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Matrices/FRAME.h>
#include <Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <Geometry/Spatial_Acceleration/CONTINUOUS_COLLISION_DETECTION.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{
template<class T,class TV,int d,int e> static void
Compute_Pairs_Helper(ARRAY<CCD_PAIR<TV::m> >& pairs,
    const CCD_DATA<TV,d>& d0,const CCD_DATA<TV,e>& d1,T thickness)
{
    struct VISITOR
    {
        ARRAY<CCD_PAIR<TV::m> >& pairs;
        const ARRAY<PAIR<int,VECTOR<int,d> > >& p0;
        const ARRAY<PAIR<int,VECTOR<int,e> > >& p1;
        VISITOR(ARRAY<CCD_PAIR<TV::m> >& pairs,
            const ARRAY<PAIR<int,VECTOR<int,d> > >& p0,
            const ARRAY<PAIR<int,VECTOR<int,e> > >& p1)
            :pairs(pairs),p0(p0),p1(p1)
        {}

        bool Cull_Self(const int self_box_index) const
        {return false;}

        bool Cull(const int self_box_index,const int other_box_index) const
        {return false;}

        void Store(const int i,const int j) const
        {pairs.Append({p0(i).x,p1(j).x,p0(i).y.Append_Elements(p1(j).y)});}
    } visitor(pairs,d0.p,d1.p); 

    d0.h.Intersection_List(d1.h,visitor,thickness);
}
template<class TV,class T> static void PF_Helper(
    CONTINUOUS_COLLISION_DETECTION<TV>& ccd0,
    CONTINUOUS_COLLISION_DETECTION<TV>& ccd1,
    ARRAY<CCD_PAIR<2> >& point_face,T thickness)
{
    Compute_Pairs_Helper(point_face,ccd0.c1,ccd1.c2,thickness);
}
template<class TV,class T> static void PF_Helper(
    CONTINUOUS_COLLISION_DETECTION<TV>& ccd0,
    CONTINUOUS_COLLISION_DETECTION<TV>& ccd1,
    ARRAY<CCD_PAIR<3> >& point_face,T thickness)
{
    Compute_Pairs_Helper(point_face,ccd0.c1,ccd1.c3,thickness);
}
template<class TV,class T> static void EE_Helper(
    CONTINUOUS_COLLISION_DETECTION<TV>& ccd0,
    CONTINUOUS_COLLISION_DETECTION<TV>& ccd1,
    ARRAY<CCD_PAIR<2> >& point_face,T thickness)
{
}
template<class TV,class T> static void EE_Helper(
    CONTINUOUS_COLLISION_DETECTION<TV>& ccd0,
    CONTINUOUS_COLLISION_DETECTION<TV>& ccd1,
    ARRAY<CCD_PAIR<3> >& point_face,T thickness)
{
    Compute_Pairs_Helper(point_face,ccd0.c2,ccd1.c2,thickness);
}
//#####################################################################
// Function Compute_Pairs_PF
//#####################################################################
template<class TV> void CONTINUOUS_COLLISION_DETECTION<TV>::
Compute_Pairs_PF(ARRAY<CCD_PAIR<TV::m> >& point_face,
    CONTINUOUS_COLLISION_DETECTION<TV>& ccd,T thickness)
{
    if(update_hierarchy) Update_Hierarchy();
    if(stale_leaves) Update_Boxes();
    if(ccd.update_hierarchy) ccd.Update_Hierarchy();
    if(ccd.stale_leaves) ccd.Update_Boxes();
    PF_Helper(*this,ccd,point_face,thickness);
}
//#####################################################################
// Function Compute_Pairs_EE
//#####################################################################
template<class TV> void CONTINUOUS_COLLISION_DETECTION<TV>::
Compute_Pairs_EE(ARRAY<CCD_PAIR<TV::m> >& edge_edge,
    CONTINUOUS_COLLISION_DETECTION<TV>& ccd,T thickness)
{
    if(update_hierarchy) Update_Hierarchy();
    if(stale_leaves) Update_Boxes();
    if(ccd.update_hierarchy) ccd.Update_Hierarchy();
    if(ccd.stale_leaves) ccd.Update_Boxes();
    EE_Helper(*this,ccd,edge_edge,thickness);
}
//#####################################################################
// Function Update_Positions
//#####################################################################
template<class TV> template<int d> void CONTINUOUS_COLLISION_DETECTION<TV>::
Update_Positions(ARRAY_VIEW<const TV> X0,ARRAY_VIEW<const TV> X1,
        CCD_DATA<TV,d>& cd)
{
    if(cd.h.box_hierarchy.m<cd.p.m)
        cd.h.box_hierarchy.Resize(cd.p.m);
    for(int i=0;i<cd.p.m;i++){
        auto e=cd.p(i).y;
        cd.h.box_hierarchy(i)=RANGE<TV>::Bounding_Box(X0.Subset(e)).Unite(RANGE<TV>::Bounding_Box(X1.Subset(e)));}
    stale_leaves=true;
}
//#####################################################################
// Function Update_Positions
//#####################################################################
template<class TV> void CONTINUOUS_COLLISION_DETECTION<TV>::
Update_Positions(ARRAY_VIEW<const TV> X0,ARRAY_VIEW<const TV> X1)
{
    Update_Positions(X0,X1,c1);
    Update_Positions(X0,X1,c2);
    Update_Positions(X0,X1,c3);
}
//#####################################################################
// Function Update_Positions
//#####################################################################
template<class TV> template<int d> void CONTINUOUS_COLLISION_DETECTION<TV>::
Update_Positions(ARRAY_VIEW<const TV> X0,ARRAY_VIEW<const TV> X1,int b,
        CCD_DATA<TV,d>& cd)
{
    if(cd.h.box_hierarchy.m<cd.p.m)
        cd.h.box_hierarchy.Resize(cd.p.m);
    for(int i=update_data(b).a[d-1],b=update_data(b).b[d-1];i<b;i++){
        auto e=cd.p(i).y;
        cd.h.box_hierarchy(i)=RANGE<TV>::Bounding_Box(X0.Subset(e)).Unite(RANGE<TV>::Bounding_Box(X1.Subset(e)));}
    stale_leaves=true;
}
//#####################################################################
// Function Update_Positions
//#####################################################################
template<class TV> void CONTINUOUS_COLLISION_DETECTION<TV>::
Update_Positions(ARRAY_VIEW<const TV> X0,ARRAY_VIEW<const TV> X1,int b)
{
    Update_Positions(X0,X1,c1);
    Update_Positions(X0,X1,c2);
    Update_Positions(X0,X1,c3);
}
//#####################################################################
// Function Update_Positions
//#####################################################################
template<class TV> template<int d> void CONTINUOUS_COLLISION_DETECTION<TV>::
Update_Positions(ARRAY_VIEW<const TV> X,const FRAME<TV>& f0,
    const FRAME<TV>& f1,int b,CCD_DATA<TV,d>& cd)
{
    if(cd.h.box_hierarchy.m<cd.p.m)
        cd.h.box_hierarchy.Resize(cd.p.m);
    for(int i=update_data(b).a[d-1],n=update_data(b).b[d-1];i<n;i++){
        VECTOR<TV,d> Y(X.Subset(cd.p(i).y));
        VECTOR<TV,2*d> Z;
        for(int j=0;j<d;j++){
            Z(j)=f0*Y(j);
            Z(j+d)=f1*Y(j);}
        cd.h.box_hierarchy(i)=RANGE<TV>::Bounding_Box(Z);}
    stale_leaves=true;
}
//#####################################################################
// Function Update_Positions
//#####################################################################
template<class TV> void CONTINUOUS_COLLISION_DETECTION<TV>::
Update_Positions(ARRAY_VIEW<const TV> X,const FRAME<TV>& f0,
    const FRAME<TV>& f1,int b)
{
    Update_Positions(X,f0,f1,b,c1);
    Update_Positions(X,f0,f1,b,c2);
    Update_Positions(X,f0,f1,b,c3);
}
//#####################################################################
// Function Update_Positions
//#####################################################################
template<class TV> template<int d> void CONTINUOUS_COLLISION_DETECTION<TV>::
Update_Positions(ARRAY_VIEW<const TV> X,int b,CCD_DATA<TV,d>& cd)
{
    if(cd.h.box_hierarchy.m<cd.p.m)
        cd.h.box_hierarchy.Resize(cd.p.m);
    for(int i=update_data(b).a[d-1],n=update_data(b).b[d-1];i<n;i++)
        cd.h.box_hierarchy(i)=RANGE<TV>::Bounding_Box(X.Subset(cd.p(i).y));
    stale_leaves=true;
}
//#####################################################################
// Function Update_Positions
//#####################################################################
template<class TV> void CONTINUOUS_COLLISION_DETECTION<TV>::
Update_Positions(ARRAY_VIEW<const TV> X,int b)
{
    Update_Positions(X,b,c1);
    Update_Positions(X,b,c2);
    Update_Positions(X,b,c3);
}
template<class T,class TV> static int
Add_Structure_Helper(CONTINUOUS_COLLISION_DETECTION<TV>& ccd,
    STRUCTURE<VECTOR<T,2> >* s,int flags)
{
    ARRAY<VECTOR<int,2> >* a2=0;
    SEGMENTED_CURVE<TV>* sc=dynamic_cast<SEGMENTED_CURVE<TV>*>(s);
    if(!sc)
        if(TRIANGULATED_AREA<T>* ta=dynamic_cast<TRIANGULATED_AREA<T>*>(s))
            sc=&ta->Get_Boundary_Object();
    if(sc) a2=&sc->mesh.elements;
    else return -1;

    int b=ccd.update_data.Add_End();
    auto& u=ccd.update_data(b);
    u.a[0]=ccd.c1.p.m;
    u.a[1]=ccd.c2.p.m;
    u.a[2]=ccd.c3.p.m;
    if(flags&1){
        HASHTABLE<int> hash;
        for(auto i:a2->Flattened())
            if(hash.Set(i))
                ccd.c1.p.Append({b,VECTOR<int,1>(i)});}
    if(flags&4)
        for(auto i:*a2)
            ccd.c2.p.Append({b,i});
    u.b[0]=ccd.c1.p.m;
    u.b[1]=ccd.c2.p.m;
    u.b[2]=ccd.c3.p.m;
    return b;
}
template<class T,class TV> static int
Add_Structure_Helper(CONTINUOUS_COLLISION_DETECTION<TV>& ccd,
    STRUCTURE<VECTOR<T,3> >* s,int flags)
{
    ARRAY<VECTOR<int,2> >* a2=0;
    ARRAY<VECTOR<int,3> >* a3=0;
    if(SEGMENTED_CURVE<TV>* sc=dynamic_cast<SEGMENTED_CURVE<TV>*>(s))
        a2=&sc->mesh.elements;
    else{
        TRIANGULATED_SURFACE<T>* ts=dynamic_cast<TRIANGULATED_SURFACE<T>*>(s);
        if(!ts)
            if(TETRAHEDRALIZED_VOLUME<T>* tv=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(s))
                ts=&tv->Get_Boundary_Object();
        if(ts){
            a2=&ts->mesh.Get_Segment_Mesh().elements;
            a3=&ts->mesh.elements;}
        else return -1;}

    int b=ccd.update_data.Add_End();
    auto& u=ccd.update_data(b);
    u.a[0]=ccd.c1.p.m;
    u.a[1]=ccd.c2.p.m;
    u.a[2]=ccd.c3.p.m;
    if(flags&1){
        HASHTABLE<int> hash;
        for(auto i:(a3?a3->Flattened():a2->Flattened()))
            if(hash.Set(i))
                ccd.c1.p.Append({b,VECTOR<int,1>(i)});}
    if(flags&2)
        for(auto i:*a2)
            ccd.c2.p.Append({b,i});
    if(flags&4)
        for(auto i:*a3)
            ccd.c3.p.Append({b,i});
    u.b[0]=ccd.c1.p.m;
    u.b[1]=ccd.c2.p.m;
    u.b[2]=ccd.c3.p.m;
    return b;
}
//#####################################################################
// Function Add_Structure
//#####################################################################
template<class TV> int CONTINUOUS_COLLISION_DETECTION<TV>::
Add_Structure(STRUCTURE<TV>* s,int flags)
{
    return Add_Structure_Helper(*this,s,flags);
}
//#####################################################################
// Function Add_Particles
//#####################################################################
template<class TV> int CONTINUOUS_COLLISION_DETECTION<TV>::
Add_Particles(int n)
{
    int b=update_data.Add_End();
    auto& u=update_data(b);
    u.a[0]=c1.p.m;
    u.a[1]=c2.p.m;
    u.a[2]=c3.p.m;
    for(int i=0;i<n;i++)
        c1.p.Append({b,VECTOR<int,1>(i)});
    u.b[0]=c1.p.m;
    u.b[1]=c2.p.m;
    u.b[2]=c3.p.m;
    return b;
}
//#####################################################################
// Function Update_Hierarchy
//#####################################################################
template<class TV> void CONTINUOUS_COLLISION_DETECTION<TV>::
Update_Hierarchy()
{
    // TODO: Make sure size is right!!!
    if(c1.p.m){
        c1.h.box_hierarchy.Resize(c1.p.m);
        c1.h.Initialize_Hierarchy_Using_KD_Tree();}
    if(c2.p.m){
        c2.h.box_hierarchy.Resize(c2.p.m);
        c2.h.Initialize_Hierarchy_Using_KD_Tree();}
    if(c3.p.m){
        c3.h.box_hierarchy.Resize(c3.p.m);
        c3.h.Initialize_Hierarchy_Using_KD_Tree();}
    update_hierarchy=false;
    stale_leaves=true;
}
//#####################################################################
// Function Update_Boxes
//#####################################################################
template<class TV> void CONTINUOUS_COLLISION_DETECTION<TV>::
Update_Boxes()
{
    if(c1.p.m) c1.h.Update_Nonleaf_Boxes();
    if(c2.p.m) c2.h.Update_Nonleaf_Boxes();
    if(c3.p.m) c3.h.Update_Nonleaf_Boxes();
    stale_leaves=false;
}
//#####################################################################
// Function Clean_Memor
//#####################################################################
template<class TV> void CONTINUOUS_COLLISION_DETECTION<TV>::
Clean_Memory()
{
    c1.h.Clean_Memory();
    c2.h.Clean_Memory();
    c3.h.Clean_Memory();
    c1.p.Clean_Memory();
    c2.p.Clean_Memory();
    c3.p.Clean_Memory();
}
template class CONTINUOUS_COLLISION_DETECTION<VECTOR<float,2> >;
template class CONTINUOUS_COLLISION_DETECTION<VECTOR<float,3> >;
template class CONTINUOUS_COLLISION_DETECTION<VECTOR<double,2> >;
template class CONTINUOUS_COLLISION_DETECTION<VECTOR<double,3> >;
}
