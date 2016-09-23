//#####################################################################
// Copyright 2005-2006, Geoffrey Irving, Frank Losasso, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __AVERAGING_COLLIDABLE_UNIFORM__
#define __AVERAGING_COLLIDABLE_UNIFORM__

#include <Grid_PDE/Interpolation/AVERAGING_UNIFORM.h>
namespace PhysBAM{

template<class TV> class GRID_BASED_COLLISION_GEOMETRY_UNIFORM;

template<class TV,class T_FACE_LOOKUP>
class AVERAGING_COLLIDABLE_UNIFORM
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    typedef T_FACE_LOOKUP FACE_LOOKUP;

    const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list;
    T default_replacement_value;

    AVERAGING_COLLIDABLE_UNIFORM(const GRID_BASED_COLLISION_GEOMETRY_UNIFORM<TV>& body_list_input,const T default_replacement_value_input)
        :body_list(body_list_input),default_replacement_value(default_replacement_value_input)
    {}

    ~AVERAGING_COLLIDABLE_UNIFORM()
    {}

    static TV Face_To_Node_Vector(const GRID<TV>& grid,const TV_INT& node_index,const T_FACE_LOOKUP& u_face)
    {const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Cell(node_index);
    TV value;for(int axis=0;axis<TV::m;axis++)for(int face=0;face<GRID<TV>::number_of_nodes_per_face;face++)
        value[axis]+=lookup(axis,GRID<TV>::Node_Face_Index(axis,node_index,face));
    return value/GRID<TV>::number_of_nodes_per_face;}

    static TV Face_To_Cell_Vector(const GRID<TV>& grid,const TV_INT& cell_index,const T_FACE_LOOKUP& u_face)
    {const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Cell(cell_index);lookup.Set_Reference_Point(grid.X(cell_index));
    TV value;for(int axis=0;axis<TV::m;axis++)
        value[axis]=(T).5*(lookup(axis,grid.First_Face_Index_In_Cell(axis,cell_index))+lookup(axis,grid.Second_Face_Index_In_Cell(axis,cell_index)));
    return value;}

    T Cell_To_Face(const GRID<TV>& grid,const int axis,const TV_INT& face_index,const ARRAY<T,TV_INT>& u_cell) const // this never needs to set starting points doesn't use velocities
    {FACE_ITERATOR<TV> iterator(grid,axis,face_index);
    TV_INT cell1,cell2;grid.Cells_Touching_Face(axis,face_index,cell1,cell2);
    if(body_list.Occupied_Face_Center(iterator)){
        T cell1_value=default_replacement_value,cell2_value=default_replacement_value;
        if(body_list.Cell_Center_Visible_From_Face(cell1,axis,face_index)) cell1_value=u_cell(cell1);
        if(body_list.Cell_Center_Visible_From_Face(cell2,axis,face_index)) cell2_value=u_cell(cell2);
        return (T).5*(cell1_value+cell2_value);}
    else return (T).5*(u_cell(cell1)+u_cell(cell2));}

    T Cell_To_Face(const GRID<TV>& grid,const int axis,const TV_INT& face_index,const ARRAY<TV,TV_INT>& u_cell) const // this never needs to set starting points doesn't use velocities
    {FACE_ITERATOR<TV> iterator(grid,axis,face_index);
    TV_INT cell1,cell2;grid.Cells_Touching_Face(axis,face_index,cell1,cell2);
    if(body_list.Occupied_Face_Center(iterator)){
        T cell1_value=default_replacement_value,cell2_value=default_replacement_value;
        if(body_list.Cell_Center_Visible_From_Face(cell1,axis,face_index)) cell1_value=u_cell(cell1)[axis];
        if(body_list.Cell_Center_Visible_From_Face(cell2,axis,face_index)) cell2_value=u_cell(cell2)[axis];
        return (T).5*(cell1_value+cell2_value);}
    else return (T).5*(u_cell(cell1)[axis]+u_cell(cell2)[axis]);}

    TV Face_To_Face_Vector(const GRID<TV>& grid,const int axis,const TV_INT& face_index,const T_FACE_LOOKUP& u_face) const
    {FACE_ITERATOR<TV> iterator(grid,axis,face_index);
    const typename T_FACE_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Face(axis,face_index);lookup.Set_Reference_Point(grid.Face(FACE_INDEX<TV::m>(axis,face_index)));
    return AVERAGING_UNIFORM<TV,T_FACE_LOOKUP>::Average_Face_To_Face_Vector_Helper(grid,iterator,lookup);}

//#####################################################################
};
}
#endif
