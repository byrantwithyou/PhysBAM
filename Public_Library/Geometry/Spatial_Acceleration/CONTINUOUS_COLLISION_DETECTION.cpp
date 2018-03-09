//#####################################################################
// Copyright 2004-2006, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Andrew Selle, Joseph Teran, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Data_Structures/HASHTABLE.h>
#include <Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <Geometry/Spatial_Acceleration/BOX_HIERARCHY_DEFINITION.h>
#include <Geometry/Spatial_Acceleration/CONTINUOUS_COLLISION_DETECTION.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{
template<class TV,int d>
void Add_Points(ARRAY_VIEW<const TV> X0,ARRAY_VIEW<const TV> X1,
    const ARRAY<VECTOR<int,d> >& e,int s,
    BOX_HIERARCHY<TV>& h,ARRAY<PAIR<int,int> >& p)
{
    HASHTABLE<int> hash;
    for(auto i:e.Flattened())
        if(hash.Set(i)){
            p.Append({s,i});
            h.box_hierarchy.Append(RANGE<TV>::Bounding_Box(X0(i),X1(i)));}
}
template<class TV,int d>
void Add_Elements(ARRAY_VIEW<const TV> X0,ARRAY_VIEW<const TV> X1,
    const ARRAY<VECTOR<int,d> >& ar,int s,
    BOX_HIERARCHY<TV>& h,ARRAY<PAIR<int,VECTOR<int,d> > >& p)
{
    for(auto i:ar){
        p.Append({s,i});
        h.box_hierarchy.Append(RANGE<TV>::Bounding_Box(X0.Subset(i)).Unite(RANGE<TV>::Bounding_Box(X1.Subset(i))));}
}

template<class T,class TV>
void Continuous_Collision_Detection_Helper(ARRAY_VIEW<const TV> X0,ARRAY_VIEW<const TV> X1,
    ARRAY_VIEW<STRUCTURE<VECTOR<T,2> >*> structures,ARRAY<CCD_PAIR<2> >* point_face,
    ARRAY<CCD_PAIR<2> >* edge_edge,T thickness)
{
    if(!point_face) return;
    BOX_HIERARCHY<TV> h1,h2;
    ARRAY<PAIR<int,int> > p1;
    ARRAY<PAIR<int,VECTOR<int,2> > > p2;
    for(int i=0;i<structures.m;i++){
        SEGMENTED_CURVE<TV>* sc=dynamic_cast<SEGMENTED_CURVE<TV>*>(structures(i));
        if(!sc)
            if(TRIANGULATED_AREA<T>* ta=dynamic_cast<TRIANGULATED_AREA<T>*>(structures(i)))
                sc=&ta->Get_Boundary_Object();
        if(sc){
            Add_Points(X0,X1,sc->mesh.elements,i,h1,p1);
            Add_Elements(X0,X1,sc->mesh.elements,i,h2,p2);}}
    h1.Initialize_Hierarchy_Using_KD_Tree();
    h2.Initialize_Hierarchy_Using_KD_Tree();

    if(point_face){
        h1.Initialize_Hierarchy_Using_KD_Tree();
        h2.Initialize_Hierarchy_Using_KD_Tree();
        ARRAY<ARRAY<int> > intersection_list;
        BOX_VISITOR_TRIVIAL visitor(intersection_list);
        h2.Intersection_List(h1,visitor,thickness);
        for(int i=0;i<intersection_list.m;i++)
            for(auto j:intersection_list(i))
                point_face->Append({p1(j).x,p2(i).x,p2(i).y.Insert(p1(j).y,0)});}
}

template<class T,class TV>
void Continuous_Collision_Detection_Helper(ARRAY_VIEW<const TV> X0,ARRAY_VIEW<const TV> X1,
    ARRAY_VIEW<STRUCTURE<VECTOR<T,3> >*> structures,ARRAY<CCD_PAIR<3> >* point_face,
    ARRAY<CCD_PAIR<3> >* edge_edge,T thickness)
{
    BOX_HIERARCHY<TV> h1,h2,h3;
    ARRAY<PAIR<int,int> > p1;
    ARRAY<PAIR<int,VECTOR<int,2> > > p2;
    ARRAY<PAIR<int,VECTOR<int,3> > > p3;
    for(int i=0;i<structures.m;i++){
        if(SEGMENTED_CURVE<TV>* sc=dynamic_cast<SEGMENTED_CURVE<TV>*>(structures(i))){
            if(point_face) Add_Points(X0,X1,sc->mesh.elements,i,h1,p1);
            if(edge_edge) Add_Elements(X0,X1,sc->mesh.elements,i,h2,p2);}
        else{
            TRIANGULATED_SURFACE<T>* ts=dynamic_cast<TRIANGULATED_SURFACE<T>*>(structures(i));
            if(!ts)
                if(TETRAHEDRALIZED_VOLUME<T>* tv=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(structures(i)))
                    ts=&tv->Get_Boundary_Object();
            if(ts){
                if(point_face) Add_Points(X0,X1,ts->mesh.elements,i,h1,p1);
                if(edge_edge) Add_Elements(X0,X1,ts->mesh.Get_Segment_Mesh().elements,i,h2,p2);
                if(point_face) Add_Elements(X0,X1,ts->mesh.elements,i,h3,p3);}}}

    if(point_face){
        h1.Initialize_Hierarchy_Using_KD_Tree();
        h3.Initialize_Hierarchy_Using_KD_Tree();
        h1.Update_Nonleaf_Boxes();
        h3.Update_Nonleaf_Boxes();
        ARRAY<ARRAY<int> > intersection_list(p3.m);
        BOX_VISITOR_TRIVIAL visitor(intersection_list);
        h3.Intersection_List(h1,visitor,thickness);
        for(int i=0;i<intersection_list.m;i++)
            for(auto j:intersection_list(i))
                point_face->Append({p1(j).x,p3(i).x,p3(i).y.Insert(p1(j).y,0)});}

    if(edge_edge){
        h2.Initialize_Hierarchy_Using_KD_Tree();
        h2.Update_Nonleaf_Boxes();
        ARRAY<ARRAY<int> > intersection_list(p2.m);
        BOX_VISITOR_TRIVIAL visitor(intersection_list);
        h2.Intersection_List(h2,visitor,thickness);
        for(int i=0;i<intersection_list.m;i++)
            for(auto j:intersection_list(i))
                edge_edge->Append({p2(i).x,p2(j).x,p2(i).y.Append_Elements(p2(j).y)});}
}

template<class TV>
void Continuous_Collision_Detection(ARRAY_VIEW<const TV> X0,ARRAY_VIEW<const TV> X1,
    ARRAY_VIEW<STRUCTURE<TV>*> structures,ARRAY<CCD_PAIR<TV::m> >* point_face,
    ARRAY<CCD_PAIR<TV::m> >* edge_edge,typename TV::SCALAR thickness)
{
    Continuous_Collision_Detection_Helper(X0,X1,structures,point_face,edge_edge,thickness);
}
template void Continuous_Collision_Detection<VECTOR<double,2> >(
    ARRAY_VIEW<VECTOR<double,2> const,int>,ARRAY_VIEW<VECTOR<double,2> const,int>,
    ARRAY_VIEW<STRUCTURE<VECTOR<double,2> >*,int>,
    ARRAY<CCD_PAIR<VECTOR<double,2>::m>,int>*,
    ARRAY<CCD_PAIR<VECTOR<double,2>::m>,int>*,VECTOR<double,2>::SCALAR);
template void Continuous_Collision_Detection<VECTOR<double,3> >(
    ARRAY_VIEW<VECTOR<double,3> const,int>,ARRAY_VIEW<VECTOR<double,3> const,int>,
    ARRAY_VIEW<STRUCTURE<VECTOR<double,3> >*,int>,
    ARRAY<CCD_PAIR<VECTOR<double,3>::m>,int>*,
    ARRAY<CCD_PAIR<VECTOR<double,3>::m>,int>*,VECTOR<double,3>::SCALAR);
template void Continuous_Collision_Detection<VECTOR<float,2> >(
    ARRAY_VIEW<VECTOR<float,2> const,int>,ARRAY_VIEW<VECTOR<float,2> const,int>,
    ARRAY_VIEW<STRUCTURE<VECTOR<float,2> >*,int>,
    ARRAY<CCD_PAIR<VECTOR<float,2>::m>,int>*,
    ARRAY<CCD_PAIR<VECTOR<float,2>::m>,int>*,VECTOR<float,2>::SCALAR);
template void Continuous_Collision_Detection<VECTOR<float,3> >(
    ARRAY_VIEW<VECTOR<float,3> const,int>,ARRAY_VIEW<VECTOR<float,3> const,int>,
    ARRAY_VIEW<STRUCTURE<VECTOR<float,3> >*,int>,
    ARRAY<CCD_PAIR<VECTOR<float,3>::m>,int>*,
    ARRAY<CCD_PAIR<VECTOR<float,3>::m>,int>*,VECTOR<float,3>::SCALAR);
}
