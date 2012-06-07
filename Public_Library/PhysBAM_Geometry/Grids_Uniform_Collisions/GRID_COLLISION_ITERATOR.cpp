//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/SORT.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_COLLISION_ITERATOR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GRID_COLLISION_ITERATOR<TV>::
GRID_COLLISION_ITERATOR(const GRID<TV>& grid_input)
    :grid(grid_input),thickness((T)1e-3*grid.dX.Min()),first_z_face(0)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GRID_COLLISION_ITERATOR<TV>::
GRID_COLLISION_ITERATOR(const GRID<TV>& grid_input,const TRIANGULATED_SURFACE<T>& surface)
    :grid(grid_input),thickness((T)1e-3*grid.dX.Min()),first_z_face(0)
{
    Initialize(surface);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GRID_COLLISION_ITERATOR<TV>::
~GRID_COLLISION_ITERATOR()
{
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void GRID_COLLISION_ITERATOR<TV>::
Initialize(const TRIANGULATED_SURFACE<T>& ts)
{
    ARRAY<TRIANGLE_3D<T> > elements;
    for(int i=0;i<ts.mesh.elements.m;i++)
        elements.Append(ts.Get_Element(i));
    Initialize(elements);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void GRID_COLLISION_ITERATOR<TV>::
Initialize(const COLLISION_GEOMETRY<TV>& cg)
{
    ARRAY<TRIANGLE_3D<T> > elements;
    for(int i=0,n=cg.Number_Of_Simplices();i<n;i++)
        elements.Append(cg.World_Space_Simplex(i));
    Initialize(elements);
}
//#####################################################################
// Function Initialize
//#####################################################################
template<class TV> void GRID_COLLISION_ITERATOR<TV>::
Initialize(const ARRAY<TRIANGLE_3D<T> >& elements)
{
    faces.Remove_All();
    HASHTABLE<TV_INT,ARRAY<int> > elements_in_cell;
    for(int i=0;i<elements.m;i++){
        RANGE<TV> box=elements(i).Bounding_Box().Thickened(thickness);
        RANGE<TV_INT> range(grid.Index(box.min_corner),grid.Index(box.max_corner)+1);
        for(UNIFORM_GRID_ITERATOR_CELL<TV> it(grid,range);it.Valid();it.Next())
            elements_in_cell.Get_Or_Insert(it.Cell_Index()).Append(i);}

    HASHTABLE<FACE_INDEX<TV::m> > boundary_faces;

    for(typename HASHTABLE<TV_INT,ARRAY<int> >::ITERATOR it(elements_in_cell);it.Valid();it.Next()){
        VECTOR<FACE_INDEX<TV::m>,TV::m*2> faces;
        GRID<TV>::Neighboring_Faces(faces,it.Key());
        boundary_faces.Set_All(faces);}

    ARRAY<char,TV_INT> cell_status(grid.Domain_Indices(3));

    for(typename HASHTABLE<FACE_INDEX<TV::m> >::ITERATOR it(boundary_faces);it.Valid();it.Next()){
        TV X=grid.X(it.Key().First_Cell_Index()),Y=grid.X(it.Key().Second_Cell_Index());
        RAY<TV> ray1(X,TV::Axis_Vector(it.Key().axis),false);ray1.semi_infinite=false;ray1.t_max=grid.dX(it.Key().axis);
        RAY<TV> ray2(Y,-TV::Axis_Vector(it.Key().axis),false);ray2.semi_infinite=false;ray2.t_max=grid.dX(it.Key().axis);
        VECTOR<int,2> closest_element;
        for(int s=0;s<2;s++){
            if(const ARRAY<int>* list=elements_in_cell.Get_Pointer(it.Key().Cell_Index(s))){
                for(int i=0;i<list->m;i++){
                    if(INTERSECTION::Intersects(ray1,elements((*list)(i)),thickness)) closest_element(0)=(*list)(i);
                    if(INTERSECTION::Intersects(ray2,elements((*list)(i)),thickness)) closest_element(1)=(*list)(i);}}}
            
        if(closest_element(0) && !closest_element(1)) closest_element(1)=closest_element(0);
        else if(closest_element(1) && !closest_element(0)) closest_element(0)=closest_element(1);
        if(!closest_element(0)) continue;
        
        VECTOR<bool,2> inside(elements(closest_element(0)).Signed_Distance(X)<0,elements(closest_element(1)).Signed_Distance(Y)<0);
        ENTRY e={it.Key(),inside,closest_element};
        faces.Append(e);}

    Sort(faces);

    for(first_z_face=1;faces(first_z_face).face.axis<TV::m;first_z_face++){}
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GRID_COLLISION_ITERATOR<TV>::FACE_ITERATOR::
FACE_ITERATOR(GRID_COLLISION_ITERATOR<TV>& gci,int ghost_input)
    :faces(gci.faces),grid(gci.grid),cur(0),face_domain(grid.Face_Indices(ghost_input))
{
    Next();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> GRID_COLLISION_ITERATOR<TV>::FACE_ITERATOR::
~FACE_ITERATOR()
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GRID_COLLISION_ITERATOR<TV>::INTERIOR_FACE_ITERATOR::
INTERIOR_FACE_ITERATOR(GRID_COLLISION_ITERATOR<TV>& gci,int ghost_input)
    :faces(gci.faces),grid(gci.grid),cur(0),last(0),ghost(ghost_input),index(1,TV_INT(1,1,1)),face_domain(grid.Face_Indices(ghost_input))
{
    LOG::cout<<"gci.first_z_face "<<gci.first_z_face<<std::endl;
    Next_Helper();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> GRID_COLLISION_ITERATOR<TV>::INTERIOR_FACE_ITERATOR::
~INTERIOR_FACE_ITERATOR()
{
}
//#####################################################################
// Function Next_Helper
//#####################################################################
template<class TV> void GRID_COLLISION_ITERATOR<TV>::INTERIOR_FACE_ITERATOR::
Next_Helper()
{
    while(1){
        cur++;
        if(cur>=faces.m){last=index.index(index.axis)-1;return;}
        if(faces(cur).inside(1)) break;}
    index=faces(cur).face;
    while(1){
        cur++;
        PHYSBAM_ASSERT(cur<faces.m && faces(cur).face.axis==index.axis && faces(cur).inside(0));
        if(!faces(cur).inside(1)) break;}
    PHYSBAM_ASSERT(index.index.Remove_Index(index.axis)==faces(cur).face.index.Remove_Index(index.axis));
    last=faces(cur).face.index(index.axis);
    
    if(index.index(index.axis)<-ghost) index.index(index.axis)=-ghost;
    if(last>grid.counts(index.axis)+ghost) last=grid.counts(index.axis)+ghost;
    if(!face_domain(index.axis).Lazy_Inside_Half_Open(index.index)) Next_Helper();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> GRID_COLLISION_ITERATOR<TV>::CELL_ITERATOR::
CELL_ITERATOR(GRID_COLLISION_ITERATOR<TV>& gci,int ghost_input)
    :faces(gci.faces),grid(gci.grid),cur(gci.first_z_face-1),last(0),index(0,0,0),ghost(ghost_input)
{
    LOG::cout<<"gci.first_z_face "<<gci.first_z_face<<std::endl;
    Next_Helper();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> GRID_COLLISION_ITERATOR<TV>::CELL_ITERATOR::
~CELL_ITERATOR()
{
}
//#####################################################################
// Function Next_Helper
//#####################################################################
template<class TV> void GRID_COLLISION_ITERATOR<TV>::CELL_ITERATOR::
Next_Helper()
{
    while(1){
        cur++;
        if(cur>=faces.m) return;
        if(faces(cur).inside(1)) break;}
    index=faces(cur).face.index;
    while(1){
        cur++;
        PHYSBAM_ASSERT(cur<faces.m && faces(cur).inside(0));
        if(!faces(cur).inside(1)) break;}
    PHYSBAM_ASSERT(index.Remove_Index(TV::m)==faces(cur).face.index.Remove_Index(TV::m));
    last=faces(cur).face.index(TV::m)-1;
    if(index(TV::m)<-ghost) index(TV::m)=-ghost;
    if(last>grid.counts(TV::m)+ghost) last=grid.counts(TV::m)+ghost;
    if(!grid.Inside_Domain(index,ghost) || last<index(TV::m)) Next_Helper();
}
template class GRID_COLLISION_ITERATOR<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class GRID_COLLISION_ITERATOR<VECTOR<double,3> >;
#endif
