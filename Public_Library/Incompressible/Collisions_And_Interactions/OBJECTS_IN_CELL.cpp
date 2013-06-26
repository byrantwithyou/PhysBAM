//#####################################################################
// Copyright 2005, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/BLOCK_UNIFORM.h>
#include <Tools/Grids_Uniform/GRID.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_ID.h>
#include <Incompressible/Collisions_And_Interactions/OBJECTS_IN_CELL.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template <class T_GRID,class ID> OBJECTS_IN_CELL<T_GRID,ID>::
OBJECTS_IN_CELL()
{
}
//#####################################################################
// Destructor
//#####################################################################
template <class T_GRID,class ID> OBJECTS_IN_CELL<T_GRID,ID>::
~OBJECTS_IN_CELL() 
{
}
//#####################################################################
// Function Reset
//#####################################################################
template <class T_GRID,class ID> void OBJECTS_IN_CELL<T_GRID,ID>::
Reset(const T_GRID& grid,const int number_of_ghost_cells)
{
    object_in_cell.Resize(grid.Cell_Indices(number_of_ghost_cells),false,false);
    object_in_cell.Fill(ID(INT_MAX));
    object_list.Clean_Memory();
}
//#####################################################################
// Function Add_Object_To_Cell
//#####################################################################
template <class T_GRID,class ID> void OBJECTS_IN_CELL<T_GRID,ID>::
Add_Object_To_Cell(const T_INDEX& cell_index,const ID object_id)
{
    assert(object_id>=ID());
    if(object_in_cell(cell_index)==ID(INT_MAX)) object_in_cell(cell_index)=object_id;
    else if(object_in_cell(cell_index)<ID()) object_list(ID(~Value(object_in_cell(cell_index)))).Append_Unique(object_id);
    else if(object_in_cell(cell_index)!=object_id){
        object_list.Resize(object_list.m+1);
        object_list.Last().Resize(2);
        object_list.Last()(0)=object_in_cell(cell_index);
        object_list.Last()(1)=object_id;
        object_in_cell(cell_index)=ID(~Value(object_list.m-1));}
}
//#####################################################################
// Function Get_Objects_For_Cells_Cell
//#####################################################################
template <class T_GRID,class ID> void OBJECTS_IN_CELL<T_GRID,ID>::
Get_Objects_For_Cells_Cell(const T_INDEX& cell) const
{
    if(object_in_cell(cell)<ID()){
        const ARRAY<ID>& list=object_list(ID(~Value(object_in_cell(cell))));
        for(int item=0;item<list.m;item++){
            ID object=list(item);
            if(!operation_hash.Is_Marked_Current(object)){operation_hash.Mark(object);merge.Append(object);}}}
    else if(object_in_cell(cell)!=ID(INT_MAX)){
        ID object=object_in_cell(cell);
        if(!operation_hash.Is_Marked_Current(object)){operation_hash.Mark(object);merge.Append(object);}}
}
//#####################################################################
// Function Get_Objects_For_Cells
//#####################################################################
template <class T_GRID,class ID> void OBJECTS_IN_CELL<T_GRID,ID>::
Get_Objects_For_Cells(const T_INDEX* cells,const int number_of_cells,const ID number_of_collision_bodies,ARRAY<ID>& objects) const
{
    Get_Objects_For_Cells_Start(number_of_collision_bodies);
    for(int i=0;i<number_of_cells;i++) Get_Objects_For_Cells_Cell(cells[i]);
    Get_Objects_For_Cells_End(objects);
}
//#####################################################################
// Function Get_Objects_For_Cells
//#####################################################################
template <class T_GRID,class ID> void OBJECTS_IN_CELL<T_GRID,ID>::
Get_Objects_For_Cells(const T_BLOCK& block,const ID number_of_collision_bodies,ARRAY<ID>& objects) const
{
    Get_Objects_For_Cells_Start(number_of_collision_bodies);
    for(int i=0;i<T_GRID::number_of_cells_per_block-1;i++) Get_Objects_For_Cells_Cell(block.Cell(i));
    Get_Objects_For_Cells_End(objects);
}
//#####################################################################
// Function Get_Objects_For_Cell
//#####################################################################
template <class T_GRID,class ID> void OBJECTS_IN_CELL<T_GRID,ID>::
Get_Objects_For_Cell(const T_INDEX& cell_index,ARRAY<ID>& objects) const
{
    assert(!objects.m);
    if(object_in_cell(cell_index)<ID()) objects=object_list(ID(~Value(object_in_cell(cell_index))));
    else if(object_in_cell(cell_index)!=ID(INT_MAX)) objects.Append(object_in_cell(cell_index));
}
//#####################################################################
// Function Get_Objects_For_Cells
//#####################################################################
template <class T_GRID,class ID> void OBJECTS_IN_CELL<T_GRID,ID>::
Get_Objects_For_Cells(const T_INDEX& cell_index1,const T_INDEX& cell_index2,const ID number_of_collision_bodies,ARRAY<ID>& objects) const
{
    Get_Objects_For_Cells_Start(number_of_collision_bodies);
    Get_Objects_For_Cells_Cell(cell_index1);
    Get_Objects_For_Cells_Cell(cell_index2);
    Get_Objects_For_Cells_End(objects);
}
//#####################################################################
// Function Get_Objects_For_Cells_Start
//#####################################################################
template <class T_GRID,class ID> void OBJECTS_IN_CELL<T_GRID,ID>::
Get_Objects_For_Cells_Start(const ID number_of_collision_bodies) const
{
    if(operation_hash.operations.m!=number_of_collision_bodies) operation_hash.Initialize(number_of_collision_bodies);
    merge.Remove_All();
    operation_hash.Next_Operation();
}
namespace PhysBAM{
template class OBJECTS_IN_CELL<GRID<VECTOR<float,1> >,COLLISION_GEOMETRY_ID>;
template class OBJECTS_IN_CELL<GRID<VECTOR<float,2> >,COLLISION_GEOMETRY_ID>;
template class OBJECTS_IN_CELL<GRID<VECTOR<float,3> >,COLLISION_GEOMETRY_ID>;
template class OBJECTS_IN_CELL<GRID<VECTOR<double,1> >,COLLISION_GEOMETRY_ID>;
template class OBJECTS_IN_CELL<GRID<VECTOR<double,2> >,COLLISION_GEOMETRY_ID>;
template class OBJECTS_IN_CELL<GRID<VECTOR<double,3> >,COLLISION_GEOMETRY_ID>;
}
