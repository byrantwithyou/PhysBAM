//#####################################################################
// Copyright 2005-2006, Avi Robinson-Mosher, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __AVERAGING_UNIFORM__
#define __AVERAGING_UNIFORM__

#include <Grid_PDE/Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <Grid_PDE/Interpolation/INTERPOLATION_UNIFORM_FORWARD.h>
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
    {return Average_Face_To_Face_Vector(grid,FACE_INDEX<TV::m>(axis,face_index),u_face);}

    TV Face_To_Face_Vector(const GRID<TV>& grid,const FACE_INDEX<TV::m>& face_index,const T_FACE_LOOKUP& u_face) const
    {return Face_To_Face_Vector(grid,face_index.axis,face_index.index,u_face);}

    template<class T_FACE_LOOKUP_LOOKUP>
    static TV Average_Face_To_Face_Vector(const GRID<TV>& grid,const FACE_INDEX<TV::m>& face,const T_FACE_LOOKUP_LOOKUP& u_face)
    {const typename T_FACE_LOOKUP_LOOKUP::LOOKUP& lookup=u_face.Starting_Point_Face(face.axis,face.index);
    return Average_Face_To_Face_Vector_Helper(grid,face,lookup);}

    template<class T_FACE_LOOKUP_LOOKUP>
    static VECTOR<T,1> Average_Face_To_Face_Vector_Helper(const GRID<VECTOR<T,1> >& grid,const FACE_INDEX<1>& face,const T_FACE_LOOKUP_LOOKUP& u_face)
    {return VECTOR<T,1>(u_face(face));}

    template<class T_FACE_LOOKUP_LOOKUP>
    static VECTOR<T,2> Average_Face_To_Face_Vector_Helper(const GRID<VECTOR<T,2> >& grid,const FACE_INDEX<2>& face,const T_FACE_LOOKUP_LOOKUP& u_face)
    {int axis=face.axis;VECTOR<int,2> cell=face.index,cell1(cell),cell2(cell);cell1(axis)--;
    VECTOR<T,2> value;value[axis]=u_face(face);
    int other_axis=1-axis;
    value[other_axis]=(T).25*(u_face(other_axis,grid.First_Face_Index_In_Cell(other_axis,cell1))+u_face(other_axis,grid.Second_Face_Index_In_Cell(other_axis,cell1))+
                              u_face(other_axis,grid.First_Face_Index_In_Cell(other_axis,cell2))+u_face(other_axis,grid.Second_Face_Index_In_Cell(other_axis,cell2)));
    return value;}

    template<class T_FACE_LOOKUP_LOOKUP>
    static VECTOR<T,3> Average_Face_To_Face_Vector_Helper(const GRID<VECTOR<T,3> >& grid,const FACE_INDEX<3>& face,const T_FACE_LOOKUP_LOOKUP& u_face)
    {static const int axis_to_other_axis[3][2]={{1,2},{0,2},{0,1}};
    int axis=face.axis;VECTOR<int,3> cell=face.index,cell1(cell),cell2(cell);cell1(axis)--;
    VECTOR<T,3> value;value[axis]=u_face(face);
    for(int i=0;i<TV::m-1;i++){
        int other_axis=axis_to_other_axis[axis][i];
        value[other_axis]=(T).25*(u_face(other_axis,grid.First_Face_Index_In_Cell(other_axis,cell1))+u_face(other_axis,grid.Second_Face_Index_In_Cell(other_axis,cell1))+
                            u_face(other_axis,grid.First_Face_Index_In_Cell(other_axis,cell2))+u_face(other_axis,grid.Second_Face_Index_In_Cell(other_axis,cell2)));}
    return value;}

    template<class T_FACE_LOOKUP_LOOKUP>
    static MATRIX<T,TV::m> Cell_Centered_Gradient(const GRID<TV>& grid,const T_FACE_LOOKUP_LOOKUP& u,const TV_INT& index)
    {
        MATRIX<T,TV::m> du;
        for(int a=0;a<TV::m;a++){
            FACE_INDEX<TV::m> f0(a,index),f1(f0);
            f1.index(a)++;
            du(a,a)=u(f1)-u(f0);
            for(int b=0;b<TV::m;b++)
                if(a!=b){
                    FACE_INDEX<TV::m> f00(f0),f01(f0),f10(f1),f11(f1);
                    f00.index(b)--;
                    f10.index(b)--;
                    f01.index(b)++;
                    f11.index(b)++;
                    du(a,b)=(u(f01)+u(f11)-u(f00)-u(f10))/4;}}
        return du*DIAGONAL_MATRIX<T,TV::m>(grid.one_over_dX);
    }
//#####################################################################
};
}
#endif
