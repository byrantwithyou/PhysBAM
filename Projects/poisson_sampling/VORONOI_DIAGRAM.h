#ifndef __VORONOI_DIAGRAM__
#define __VORONOI_DIAGRAM__
//#####################################################################
// Copyright 2016, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VORONOI_DIAGRAM
//##################################################################### 

#include <set>

namespace PhysBAM{

template<class T>
struct VORONOI_DIAGRAM
{
    typedef VECTOR<T,2> TV;
    enum VERTEX_STATE {unknown,in,out};
    enum CELL_STATE {non_incident,incident};
    struct COEDGE;
    struct CELL;

    static char next_cell;
    static char next_vertex;
    static char next_coedge;
    
    struct VERTEX
    {
        TV X;
        VERTEX_STATE state;
        COEDGE* first_coedge; // outgoing
        char name;
        
        VERTEX()
            :state(unknown),first_coedge(0),name(next_vertex++)
        {}

        T Criterion(const TV& A,const TV& B,const TV& C,const TV& D) const;
        T Criterion(const TV& X) const;
        void Compute_Point(const TV& A,const TV& B,const TV& C);
        void Compute_Point();
        void Print() const;
    };

    // Trace regions ccw
    struct COEDGE
    {
        VERTEX* head, *tail;
        COEDGE* pair, *next, *prev;
        CELL* cell;
        char name;

        COEDGE()
            :head(0),tail(0),pair(0),next(0),prev(0),cell(0),name(next_coedge++)
        {}
        void Print() const;
    };

    struct CELL
    {
        TV X;
        CELL_STATE state;
        COEDGE * first_coedge;
        char name;

        CELL()
            :state(non_incident),first_coedge(0),name(next_cell++)
        {}
        void Print() const;
    };

    ARRAY<CELL*> cells;
    std::set<COEDGE*> coedges;
    std::set<VERTEX*> vertices;

    void Insert_Coedge(COEDGE* ce)
    {
        coedges.insert(ce);
    }

    void Remove_Coedge(COEDGE* ce)
    {
        coedges.erase(ce);
    }

    void Update_Coedge(COEDGE* ce)
    {
    }

    void Discover_Inside(ARRAY<COEDGE*>& in,ARRAY<COEDGE*>& adj,
        ARRAY<VERTEX*>& out_v,ARRAY<VERTEX*>& in_v,COEDGE* ce,const TV& new_pt);
    void Insert_Point(COEDGE* start,const TV& new_pt);
    void First_Three_Points(TV A,TV B,TV C);
    void Visualize_State(const char* title) const;
    void Sanity_Checks() const;
    VERTEX* Select_First_In_Vertex(COEDGE* start,const TV& new_pt) const;
    void Print() const;
};



}
#endif
