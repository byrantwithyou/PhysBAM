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
    enum CELL_TYPE {type_inside,type_outside};
    struct COEDGE;
    struct CELL;
    
    struct VERTEX
    {
        TV X;
        VERTEX_STATE state;
        COEDGE* coedge; // outgoing
        
        VERTEX()
            :state(unknown),coedge(0)
        {}

        T Criterion(const TV& A,const TV& B,const TV& C,const TV& D) const;
        T Criterion(const TV& X) const;
        void Compute_Point(const TV& A,const TV& B,const TV& C);
        void Compute_Point();
    };

    // Trace regions ccw
    struct COEDGE
    {
        VERTEX* head, *tail;
        COEDGE* pair, *next, *prev;
        CELL* cell;
        int piece;
        
        COEDGE()
            :head(0),tail(0),pair(0),next(0),prev(0),cell(0),piece(-1)
        {}
    };

    struct CELL
    {
        TV X;
        CELL_STATE state;
        CELL_TYPE type;
        COEDGE * coedge;
        bool outside;

        CELL()
            :state(non_incident),type(type_inside),coedge(0),outside(false)
        {}
    };

    ARRAY<CELL*> cells;
    T radius;
    RANGE<TV> bounding_box;
    mutable RANDOM_NUMBERS<T> random;
    enum PIECE_TYPE {unset,empty,no_disc,full_disc,out0,out1,both_out};

    struct PIECE_HELPER
    {
        PIECE_TYPE type;
        TV A,B,C;
        T aux0,aux1,aux2;
        T area;

        PIECE_HELPER()
            :type(unset),aux0(0),aux1(0),aux2(0)
        {}

        T Compute(COEDGE* ce,T radius,bool clipped);
        TV Choose_Feasible_Point(RANDOM_NUMBERS<T>& random,T radius) const;
    };
    
    static const int clipped_piece_offset=1<<30;
    struct PIECE
    {
        T area;
        T subtree_area;
        COEDGE* coedge;
        PIECE_HELPER h;
        
        PIECE()
            :area(0),subtree_area(0),coedge(0)
        {}
    };
    ARRAY<PIECE> pieces;

    struct CLIPPED_PIECE
    {
        T area;
        T subtree_area;
        COEDGE* coedge;
        int num_sub_pieces;
        PIECE_HELPER sub_pieces[4];
        CLIPPED_PIECE()
            :area(0),subtree_area(0),coedge(0),num_sub_pieces(0)
        {}
    };
    ARRAY<CLIPPED_PIECE> clipped_pieces;

    VORONOI_DIAGRAM();

    void Init(const RANGE<TV>& box,T radius_input);
    
    void Update_Piece_Tree(int i);
    void Update_Clipped_Piece_Tree(int i);
    void Insert_Coedge(COEDGE* ce);
    void Remove_Coedge(COEDGE* ce);
    void Update_Coedge(COEDGE* ce);
    void Remove_Piece(int p);
    void Remove_Clipped_Piece(int p);
    void Insert_Clipped_Coedge(COEDGE* ce);

    int Choose_Piece() const;
    TV Choose_Feasible_Point(int p) const;

    void Discover_Inside(ARRAY<COEDGE*>& in,ARRAY<COEDGE*>& adj,
        ARRAY<VERTEX*>& in_v,COEDGE* ce,const TV& new_pt);
    void Insert_Point(COEDGE* start,const TV& new_pt);
    void Insert_Point(int p,const TV& new_pt);
    void Visualize_State(const char* title) const;
    void Sanity_Checks() const;
    void Sample_Fully(const RANGE<TV>& box,T radius_input);
    void Get_Samples(ARRAY<TV>& X) const;
};
}
#endif
