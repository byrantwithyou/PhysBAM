//#####################################################################
// Copyright 2005-2006, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __AVERAGING_UNIFORM__
#define __AVERAGING_UNIFORM__

#include <Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
namespace PhysBAM{

template<class TV,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<TV>
class AVERAGING_UNIFORM
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
public:

    typedef T_FACE_LOOKUP FACE_LOOKUP;

    AVERAGING_UNIFORM()
    {}

    ~AVERAGING_UNIFORM()
    {}

    TV Face_To_Cell_Vector(const GRID<TV>& grid,const TV_INT& cell_index,const T_FACE_LOOKUP& u_face) const
    {const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Cell(cell_index);
    TV value;for(int axis=0;axis<TV::m;axis++)
        value[axis]=(T).5*(lookup(axis,grid.First_Face_Index_In_Cell(axis,cell_index))+lookup(axis,grid.Second_Face_Index_In_Cell(axis,cell_index)));
    return value;}

    TV Face_To_Node_Vector(const GRID<TV>& grid,const TV_INT& node_index,const T_FACE_LOOKUP& u_face) const
    {const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Cell(node_index);
    TV value;for(int axis=0;axis<TV::m;axis++)for(int face=0;face<GRID<TV>::number_of_nodes_per_face;face++)
        value[axis]+=lookup(axis,GRID<TV>::Node_Face_Index(axis,node_index,face));
    return value/GRID<TV>::number_of_nodes_per_face;}

    // Note this never has to be replaced because it is phi averaged to the faces, not the velocities
    T Cell_To_Face(const GRID<TV>& grid,const int axis,const TV_INT& face_index,const ARRAY<T,TV_INT>& u_cell) const
    {TV_INT cell1,cell2;grid.Cells_Touching_Face(axis,face_index,cell1,cell2);
    return (T).5*(u_cell(cell1)+u_cell(cell2));}

    // Note this never has to be replaced because it is phi averaged to the faces, not the velocities
    T Cell_To_Face(const GRID<TV>& grid,const FACE_INDEX<TV::m>& face_index,const ARRAY<T,TV_INT>& u_cell) const
    {return Cell_To_Face(grid,face_index.axis,face_index.index,u_cell);}

    TV Face_To_Face_Vector(const GRID<TV>& grid,const int axis,const TV_INT& face_index,const T_FACE_LOOKUP& u_face) const
    {FACE_ITERATOR<TV> iterator(grid,axis,face_index);
    return Average_Face_To_Face_Vector(grid,iterator,u_face);}

    TV Face_To_Face_Vector(const GRID<TV>& grid,const FACE_INDEX<TV::m>& face_index,const T_FACE_LOOKUP& u_face) const
    {return Face_To_Face_Vector(grid,face_index.axis,face_index.index,u_face);}

    template<class T_FACE_LOOKUP_LOOKUP,class T_FACE_ITERATOR>
    static TV Average_Face_To_Face_Vector(const GRID<TV>& grid,const T_FACE_ITERATOR& iterator,const T_FACE_LOOKUP_LOOKUP& u_face)
    {const typename T_FACE_LOOKUP_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Face(iterator.Axis(),iterator.Face_Index());
    return Average_Face_To_Face_Vector_Helper(grid,iterator,lookup);}

    template<class T_FACE_LOOKUP_LOOKUP>
    static VECTOR<T,1> Average_Face_To_Face_Vector_Helper(const GRID<VECTOR<T,1> >& grid,const FACE_ITERATOR<VECTOR<T,1> >& iterator,const T_FACE_LOOKUP_LOOKUP& u_face)
    {return VECTOR<T,1>(u_face(iterator.Axis(),iterator.Face_Index()));}

    template<class T_FACE_LOOKUP_LOOKUP>
    static VECTOR<T,2> Average_Face_To_Face_Vector_Helper(const GRID<VECTOR<T,2> >& grid,const FACE_ITERATOR<VECTOR<T,2> >& iterator,const T_FACE_LOOKUP_LOOKUP& u_face)
    {int axis=iterator.Axis();typename GRID<VECTOR<T,2> >::INDEX face=iterator.Face_Index(),cell1,cell2;iterator.Unordered_Cell_Indices_Touching_Face(cell1,cell2);
    VECTOR<T,2> value;value[axis]=u_face(axis,face);
    int other_axis=1-axis;
    value[other_axis]=(T).25*(u_face(other_axis,grid.First_Face_Index_In_Cell(other_axis,cell1))+u_face(other_axis,grid.Second_Face_Index_In_Cell(other_axis,cell1))+
                              u_face(other_axis,grid.First_Face_Index_In_Cell(other_axis,cell2))+u_face(other_axis,grid.Second_Face_Index_In_Cell(other_axis,cell2)));
    return value;}

    template<class T_FACE_LOOKUP_LOOKUP>
    static VECTOR<T,3> Average_Face_To_Face_Vector_Helper(const GRID<VECTOR<T,3> >& grid,const FACE_ITERATOR<VECTOR<T,3> >& iterator,const T_FACE_LOOKUP_LOOKUP& u_face)
    {static const int axis_to_other_axis[3][2]={{1,2},{0,2},{0,1}};
    int axis=iterator.Axis();typename GRID<VECTOR<T,3> >::INDEX face=iterator.Face_Index(),cell1,cell2;iterator.Unordered_Cell_Indices_Touching_Face(cell1,cell2);
    VECTOR<T,3> value;value[axis]=u_face(axis,face);
    for(int i=0;i<TV::m-1;i++){
        int other_axis=axis_to_other_axis[axis][i];
        value[other_axis]=(T).25*(u_face(other_axis,grid.First_Face_Index_In_Cell(other_axis,cell1))+u_face(other_axis,grid.Second_Face_Index_In_Cell(other_axis,cell1))+
                            u_face(other_axis,grid.First_Face_Index_In_Cell(other_axis,cell2))+u_face(other_axis,grid.Second_Face_Index_In_Cell(other_axis,cell2)));}
    return value;}
//#####################################################################
};
}
#endif
