//#####################################################################
// Copyright 2005-2010, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BLOCK_UNIFORM
//#####################################################################
#ifndef __BLOCK_UNIFORM__
#define __BLOCK_UNIFORM__

#include <Grid_Tools/Grids/GRID.h>
namespace PhysBAM{

template<class TV>
class BLOCK_UNIFORM
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:
    const GRID<TV>& grid;
    const TV_INT block_index;

    BLOCK_UNIFORM(const GRID<TV>& grid_input,const TV_INT& block_index_input)
        :grid(grid_input),block_index(block_index_input)
    {}

    BLOCK_UNIFORM(const GRID<TV>& grid_input,const TV& X,const int number_of_ghost_cells=3)
        :grid(grid_input),block_index(grid.Block_Index(X,number_of_ghost_cells))
    {}

    TV_INT Block() const
    {return block_index;}

private:
    TV_INT Face_X(const int face_index) const
    {assert(0<=face_index&&face_index<GRID<TV>::number_of_faces_per_block/TV::m);
    static const int lookup[][3]={{-1,-1,-1},{0,-1,-1},{1,-1,-1},{-1,0,-1},{0,0,-1},{1,0,-1},{-1,-1,0},{0,-1,0},{1,-1,0},{-1,0,0},{0,0,0},{1,0,0}};
    TV_INT face;for(int i=0;i<TV::m;i++)face[i]=block_index[i]+lookup[face_index][i];return face;}

    TV_INT Face_Y(const int face_index) const
    {assert(0<=face_index&&face_index<GRID<TV>::number_of_faces_per_block/TV::m&&TV::m>=2);
    static const int lookup[][3]={{-1,-1,-1},{0,-1,-1},{-1,0,-1},{0,0,-1},{-1,1,-1},{0,1,-1},{-1,-1,0},{0,-1,0},{-1,0,0},{0,0,0},{-1,1,0},{0,1,0}};
    TV_INT face;for(int i=0;i<TV::m;i++)face[i]=block_index[i]+lookup[face_index][i];return face;}

    TV_INT Face_Z(const int face_index) const
    {assert(0<=face_index&&face_index<GRID<TV>::number_of_faces_per_block/TV::m&&TV::m==3);
    static const int lookup[][3]={{-1,-1,-1},{0,-1,-1},{-1,0,-1},{0,0,-1},{-1,-1,0},{0,-1,0},{-1,0,0},{0,0,0},{-1,-1,1},{0,-1,1},{-1,0,1},{0,0,1}};
    TV_INT face;for(int i=0;i<TV::m;i++)face[i]=block_index[i]+lookup[face_index][i];return face;}

public:
    template<class T_FACE_LOOKUP>
    typename T_FACE_LOOKUP::ELEMENT Face_X_Value(const T_FACE_LOOKUP& face_value,const int block_face_index) const
    {return face_value(0,Face_X(block_face_index));}

    T& Face_X_Reference(ARRAY<T,FACE_INDEX<TV::m> >& face_value,const int block_face_index) const
    {return face_value(0,Face_X(block_face_index));}

    template<class T_FACE_LOOKUP>
    typename T_FACE_LOOKUP::ELEMENT Face_Y_Value(const T_FACE_LOOKUP& face_value,const int block_face_index) const
    {assert(TV::m>=2);return face_value(1,Face_Y(block_face_index));}

    T& Face_Y_Reference(ARRAY<T,FACE_INDEX<TV::m> >& face_value,const int block_face_index) const
    {assert(TV::m>=2);return face_value(1,Face_Y(block_face_index));}

    template<class T_FACE_LOOKUP>
    typename T_FACE_LOOKUP::ELEMENT Face_Z_Value(const T_FACE_LOOKUP& face_value,const int block_face_index) const
    {assert(TV::m==3);return face_value(2,Face_Z(block_face_index));}

    T& Face_Z_Reference(ARRAY<T,FACE_INDEX<TV::m> >& face_value,const int block_face_index) const
    {assert(TV::m==3);return face_value(2,Face_Z(block_face_index));}

    FACE_INDEX<TV::m> Face_X_Index(const int block_face_index) const
    {return FACE_INDEX<TV::m>(0,Face_X(block_face_index));}
    
    FACE_INDEX<TV::m> Face_Y_Index(const int block_face_index) const
    {assert(TV::m>=2);return FACE_INDEX<TV::m>(1,Face_Y(block_face_index));}
    
    FACE_INDEX<TV::m> Face_Z_Index(const int block_face_index) const
    {assert(TV::m==3);return FACE_INDEX<TV::m>(2,Face_Z(block_face_index));}

    TV_INT Cell(const int cell_index) const
    {assert(0<=cell_index&&cell_index<GRID<TV>::number_of_cells_per_block);
    static const int lookup[][3]={{-1,-1,-1},{0,-1,-1},{-1,0,-1},{0,0,-1},{-1,-1,0},{0,-1,0},{-1,0,0},{0,0,0}};
    TV_INT cell;for(int i=0;i<TV::m;i++)cell[i]=block_index[i]+lookup[cell_index][i];return cell;}

    void All_Cell_Indices(TV_INT cell_indices[GRID<TV>::number_of_cells_per_block]) const
    {for(int i=0;i<GRID<TV>::number_of_cells_per_block;i++)cell_indices[i]=Cell(i);}

    TV Minimum_Corner() const
    {return grid.Center(block_index-TV_INT::All_Ones_Vector());}

    TV Center() const
    {return grid.Node(block_index);}

    TV One_Over_DX() const
    {return grid.one_over_dX;}

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>(Minimum_Corner(),grid.Center(block_index));}

    // TODO: these are not the typical parameters to Inside functions
    bool Inside(const TV& X,const T thickness_multiplier=-(T)1e-3) const
    {return Bounding_Box().Inside(X,thickness_multiplier*grid.dX.Min());}

    bool Lazy_Inside(const TV& X) const
    {return Bounding_Box().Lazy_Inside(X);}

    TV_INT Incident_Face(const int axis,const int face_index)
    {static const int lookup[3][4][3]={{{0,-1,-1},{0,0,-1},{0,-1,0},{0,0,0}},{{-1,0,-1},{0,0,-1},{-1,0,0},{0,0,0}},{{-1,-1,0},{0,-1,0},{-1,0,0},{0,0,0}}};
    TV_INT incident_face;for(int i=0;i<TV::m;i++) incident_face[i]=block_index[i]+lookup[axis][face_index][i];
    return incident_face;}

//#####################################################################
};

// workaround for a memory problem that arises in the optimized version of the 2d float code
template<>
inline VECTOR<int,2> BLOCK_UNIFORM<VECTOR<float,2> >::Face_X(const int face_index) const
{static const VECTOR<int,2> lookup[]={VECTOR<int,2>(-1,-1),VECTOR<int,2>(0,-1),VECTOR<int,2>(1,-1),VECTOR<int,2>(-1,0),VECTOR<int,2>(0,0),VECTOR<int,2>(1,0)};
return block_index+lookup[face_index];}

template<>
inline VECTOR<int,2> BLOCK_UNIFORM<VECTOR<float,2> >::Face_Y(const int face_index) const
{static const VECTOR<int,2> lookup[]={VECTOR<int,2>(-1,-1),VECTOR<int,2>(0,-1),VECTOR<int,2>(-1,0),VECTOR<int,2>(0,0),VECTOR<int,2>(-1,1),VECTOR<int,2>(0,1)};
return block_index+lookup[face_index];}

}
#endif
