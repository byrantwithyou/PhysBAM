//#####################################################################
// Copyright 2009, Nipun Kwatra, Avi Robinson-Mosher, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COLLISION_AWARE_INDEX_MAP
//##################################################################### 
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY.h>
#include <Rigids/Collisions/COLLISION_GEOMETRY_ID.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <Incompressible/Collisions_And_Interactions/DEFORMABLE_OBJECT_FLUID_COLLISIONS.h>
#include <Dynamics/Coupled_Evolution/COLLISION_AWARE_INDEX_MAP.h>
#include <Dynamics/Coupled_Evolution/IMPLICIT_BOUNDARY_CONDITION_COLLECTION.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED.h>
#include <Dynamics/Coupled_Evolution/UNIFORM_INTERIOR_ITERATOR_CELL.h>

using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COLLISION_AWARE_INDEX_MAP<TV>::
COLLISION_AWARE_INDEX_MAP(UNIFORM_COLLISION_AWARE_ITERATOR_FACE_INFO<TV>& info,IMPLICIT_BOUNDARY_CONDITION_COLLECTION<TV>& boundary_condition_collection_input)
    :last_coupling_cell(0),iterator_info(info),grid(iterator_info.grid),number_extra_cells(0),two_phase(false),
    boundary_condition_collection(boundary_condition_collection_input)
{
}
//#####################################################################
// Function Construct_Indices
//#####################################################################
template<class TV> void COLLISION_AWARE_INDEX_MAP<TV>::
Construct_Indices(const int ghost_cells)
{
    Clear(ghost_cells);

    for(UNIFORM_COLLISION_AWARE_ITERATOR_FACE_COUPLED<TV> iterator(iterator_info);iterator.Valid();iterator.Next()){
        if(Register_Constrained_Face_Index(SIDED_FACE_INDEX<d>(iterator.side,iterator.Full_Index())))
            Register_Cell_Index(iterator.Real_Cell_Index(),ghost_cells);}

    last_coupling_cell=indexed_cells.m;

    for(UNIFORM_COLLISION_AWARE_ITERATOR_FACE_UNCOUPLED<TV> iterator(iterator_info,0,!two_phase);iterator.Valid();iterator.Next_Fluid()){
        if(Register_Face_Index(iterator.Full_Index())){
            Register_Cell_Index(iterator.First_Cell_Index(), ghost_cells);
            Register_Cell_Index(iterator.Second_Cell_Index(),ghost_cells);}}

    real_cell_indices_reverse_map.Resize(indexed_cells.m);
    real_cell_indices_reverse_map.Fill(-1);
    for(int i=0;i<real_cell_indices.m;i++) real_cell_indices_reverse_map(real_cell_indices(i))=i;
    for(int axis=0;axis<d;axis++){
        if(boundary_condition_collection.periodic_boundary[axis]){
            for(CELL_ITERATOR<TV> iterator(grid,1,GRID<TV>::GHOST_REGION,2*axis);iterator.Valid();iterator.Next()){
                TV_INT cell_index=iterator.Cell_Index();
                TV_INT matching_cell_index(cell_index);matching_cell_index(axis)=1;
                cell_indices(cell_index)=cell_indices(matching_cell_index);
                cell_indices(matching_cell_index-TV_INT::Axis_Vector(axis))=cell_indices(cell_index-TV_INT::Axis_Vector(axis));
            }
            for(FACE_ITERATOR<TV> iterator(grid,0,GRID<TV>::BOUNDARY_REGION,2*axis,axis);iterator.Valid();iterator.Next()){
                FACE_INDEX<d> face_index=iterator.Full_Index();
                FACE_INDEX<d> matching_face_index(face_index);matching_face_index.index(axis)=1;
                if(face_indices(face_index)) face_indices(matching_face_index)=face_indices(face_index);
                // TODO(jontg): else what, for a constrained face?
            }
        }
    }
}
//#####################################################################
// Function Register_Face_Index
//#####################################################################
template<class TV> bool COLLISION_AWARE_INDEX_MAP<TV>::
Register_Face_Index(const FACE_INDEX<d>& face_index)
{
    int axis=face_index.axis;
    if(boundary_condition_collection.periodic_boundary[axis] && face_index.index(axis)==0)
        return false;
    assert(face_indices(face_index)<0);
    face_indices(face_index)=indexed_faces.Append(face_index);
    return true;
}
//#####################################################################
// Function Register_Constrained_Face_Index
//#####################################################################
template<class TV> bool COLLISION_AWARE_INDEX_MAP<TV>::
Register_Constrained_Face_Index(const SIDED_FACE_INDEX<d>& face_index)
{
    assert(!constraint_indices.Contains(face_index));
    constraint_indices.Insert(face_index,indexed_constraints.Append(face_index));
    return true;
}
//#####################################################################
// Function Register_Cell_Index
//#####################################################################
template<class TV> void COLLISION_AWARE_INDEX_MAP<TV>::
Register_Cell_Index(const TV_INT& index,const int ghost_cells)
{
    int& cell_index=cell_indices(index);
    if(cell_index<0){
        bool outside=!grid.Inside_Domain(index) || (!two_phase && (*iterator_info.outside_fluid)(index));

        if(outside){
            // TODO: handle ghost cells for periodic boundaries
            for(int axis=0;axis<d;axis++)
                if(boundary_condition_collection.periodic_boundary[axis] && (index(axis)<1 || index(axis)>grid.counts(axis)))
                    return;}
        else{cell_index=indexed_cells.Append(index);
            if(grid.Inside_Domain(index))
                real_cell_indices.Append(cell_index);}}
}
//#####################################################################
// Function Collect
//#####################################################################
template<class TV> void COLLISION_AWARE_INDEX_MAP<TV>::
Collect(const ARRAY<T,FACE_INDEX<d> >& faces,const ARRAY<T>& constrained_faces,ARRAY<T>& flattened_faces) const
{
    flattened_faces.Resize(Number_Faces());
    for(int i=0;i<indexed_faces.m;i++)
        flattened_faces(i)=faces(indexed_faces(i));
    for(int i=0;i<indexed_constraints.m;i++)
        flattened_faces(indexed_faces.m+i)=constrained_faces(i);
}
//#####################################################################
// Function Collect
//#####################################################################
template<class TV> void COLLISION_AWARE_INDEX_MAP<TV>::
Collect(const ARRAY<T,TV_INT>& cells,ARRAY<T>& flattened_cells) const
{
    flattened_cells.Resize(Number_Cells());
    for(int i=0;i<real_cell_indices.m;i++){int index=real_cell_indices(i);
        flattened_cells(index)=cells(indexed_cells(index));}
}
//#####################################################################
// Function Collect_Indexed_Cells
//#####################################################################
// Collect's pressures outside the domain as well as real ones
template<class TV> void COLLISION_AWARE_INDEX_MAP<TV>::
Collect_Indexed_Cells(const ARRAY<T,TV_INT>& cells,ARRAY<T>& flattened_cells) const
{
    flattened_cells.Resize(Number_Cells());
    for(int i=0;i<indexed_cells.m;i++)
        flattened_cells(i)=cells(indexed_cells(i));
}
//#####################################################################
// Function Distribute
//#####################################################################
template<class TV> void COLLISION_AWARE_INDEX_MAP<TV>::
Distribute(const ARRAY<T>& flattened_faces,ARRAY<T,FACE_INDEX<d> >& faces,ARRAY<T>& constrained_faces) const
{
    for(int i=0;i<indexed_faces.m;i++)
        faces(indexed_faces(i))=flattened_faces(i);
    for(int i=0;i<indexed_constraints.m;i++){
        faces(indexed_constraints(i).Face_Index())=flattened_faces(indexed_faces.m+i); // HACK: Send these back to the incompressible guys.
        constrained_faces(i)=flattened_faces(indexed_faces.m+i);}
}
//#####################################################################
// Function Distribute
//#####################################################################
template<class TV> void COLLISION_AWARE_INDEX_MAP<TV>::
Distribute(const ARRAY<T>& flattened_cells,ARRAY<T,TV_INT>& cells) const
{
    for(int i=0;i<real_cell_indices.m;i++){int index=real_cell_indices(i);
        cells(indexed_cells(index))=flattened_cells(index);}
}
//#####################################################################
// Function Distribute_Indexed_Cells
//#####################################################################
// Distributes pressures outside the domain as well as real ones
template<class TV> void COLLISION_AWARE_INDEX_MAP<TV>::
Distribute_Indexed_Cells(const ARRAY<T>& flattened_cells,ARRAY<T,TV_INT>& cells) const
{
    for(int i=0;i<indexed_cells.m;i++)
        cells(indexed_cells(i))=flattened_cells(i);
}
//#####################################################################
// Function Clear
//#####################################################################
template<class TV> void COLLISION_AWARE_INDEX_MAP<TV>::
Clear(const int ghost_cells)
{
    for(int i=0;i<indexed_faces.m;i++){
        const FACE_INDEX<d>& face_index=indexed_faces(i);
        face_indices(face_index)=0;}
    indexed_faces.Remove_All();

    last_coupling_cell=0;
    for(int i=0;i<indexed_cells.m;i++)
        cell_indices(indexed_cells(i))=0;
    indexed_cells.Remove_All();
    real_cell_indices.Remove_All();
    cell_indices.Resize(grid.Domain_Indices(ghost_cells+1),true,false,-1);
    face_indices.Resize(grid,ghost_cells+1,true,false,-1);
}
//#####################################################################
// Function Print
//#####################################################################
template<class TV> void COLLISION_AWARE_INDEX_MAP<TV>::
Print(int id) const
{
    char buff[1000];
    sprintf(buff, "index_map-%i.txt", id);
    OCTAVE_OUTPUT<T> oo(buff);
    oo.Write("last_coupling_cell",(T)last_coupling_cell);
    ARRAY<VECTOR<int,TV::m+1> > flat_faces;
    for(int i=0;i<indexed_faces.m;i++) flat_faces.Append(indexed_faces(i).index.Insert(indexed_faces(i).axis,0));
    oo.Write("indexed_faces",flat_faces);
    oo.Write("indexed_cells",indexed_cells);
    oo.Write("real_cell_indices",real_cell_indices);
}
//#####################################################################
// Function Register_Dirichlet_Cell
//#####################################################################
template<class TV> void COLLISION_AWARE_INDEX_MAP<TV>::
Register_Dirichlet_Cell(const TV_INT& index)
{
/*    if(dirichlet_cell_indices.Set(index,indexed_dirichlet_cells.m+1))
        indexed_dirichlet_cells.Append(index);
*/
}
//#####################################################################
namespace PhysBAM{
template class COLLISION_AWARE_INDEX_MAP<VECTOR<float,1> >;
template class COLLISION_AWARE_INDEX_MAP<VECTOR<float,2> >;
template class COLLISION_AWARE_INDEX_MAP<VECTOR<float,3> >;
template class COLLISION_AWARE_INDEX_MAP<VECTOR<double,1> >;
template class COLLISION_AWARE_INDEX_MAP<VECTOR<double,2> >;
template class COLLISION_AWARE_INDEX_MAP<VECTOR<double,3> >;
}
