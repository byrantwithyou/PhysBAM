#ifndef __VORONOI_DIAGRAM__
#define __VORONOI_DIAGRAM__
//#####################################################################
// Copyright 2016, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VORONOI_DIAGRAM
//##################################################################### 
#include <Core/Math_Tools/RANGE.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
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
        int piece;
        
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
        bool outside;
        char name;

        CELL()
            :state(non_incident),first_coedge(0),outside(false),name(next_cell++)
        {}
        void Print() const;
    };

    ARRAY<CELL*> cells;
    std::set<COEDGE*> coedges;
    std::set<VERTEX*> vertices;
    T radius;
    RANGE<TV> bounding_box;
    RANDOM_NUMBERS<T> random;
    enum PIECE_TYPE {unset,empty,no_disc,full_disc,out0,out1,both_out};

    struct PIECE_HELPER
    {
        PIECE_TYPE type;
        TV A,B,C;
        T aux0,aux1,aux2;

        PIECE_HELPER()
            :type(unset),aux0(0),aux1(0),aux2(0)
        {}

        T Compute(COEDGE* ce,T radius,bool clipped);
        TV Choose_Feasible_Point(RANDOM_NUMBERS<T>& random,T radius) const;

        void Print() const;
    };
    
    static const int first_clipped_piece_index=1<<30;
    struct PIECE
    {
        T this_area;
        T subtree_area;
        COEDGE* coedge;
        PIECE_HELPER h;
        
        PIECE()
            :this_area(0),subtree_area(0),coedge(0)
        {}

        void Print() const;
    };
    ARRAY<PIECE> pieces;

    struct CLIPPED_PIECE
    {
        T this_area;
        T subtree_area;
        COEDGE* coedge;
        int num_sub_pieces;
        PIECE_HELPER sub_pieces[4];
        CLIPPED_PIECE()
            :this_area(0),subtree_area(0),coedge(0),num_sub_pieces(0)
        {}

        void Print() const;
    };
    ARRAY<CLIPPED_PIECE> clipped_pieces;

    void Init(const RANGE<TV>& box);
    
    void Update_Piece_Tree(int i,T diff_area);
    void Update_Clipped_Piece_Tree(int i,T diff_area);
    void Insert_Coedge(COEDGE* ce);
    void Remove_Coedge(COEDGE* ce);
    void Update_Coedge(COEDGE* ce);
    void Remove_Piece(int p);
    void Remove_Clipped_Piece(int p);
    void Insert_Clipped_Coedge(COEDGE* ce);
    
    T Compute_Available_Area(TV A,TV B,TV C,bool first_vertex_disk); // ccw order
    int Choose_Piece();
    
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
