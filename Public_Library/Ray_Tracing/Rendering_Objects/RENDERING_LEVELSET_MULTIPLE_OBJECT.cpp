//#####################################################################
// Copyright 2005, Jiayi Chong, Jeong-Mo Hong, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Grids_Uniform/GRID.h>
#include <Tools/Matrices/MATRIX_3X3.h>
#include <Tools/Nonlinear_Equations/ITERATIVE_SOLVER.h>
#include <Geometry/Intersections/RAY_BOX_INTERSECTION.h>
#include <Geometry/Level_Sets/LEVELSET.h>
#include <Incompressible/Level_Sets/LEVELSET_MULTIPLE.h>
#include <Incompressible/Level_Sets/LEVELSET_MULTIPLE_ON_A_RAY.h>
#include <Ray_Tracing/Rendering_Objects/RENDERING_LEVELSET_MULTIPLE_OBJECT.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class T_LEVELSET_MULTIPLE> RENDERING_LEVELSET_MULTIPLE_OBJECT<T_LEVELSET_MULTIPLE>::
RENDERING_LEVELSET_MULTIPLE_OBJECT(GRID<TV>& grid_input,ARRAY<ARRAY<T,VECTOR<int,3> > >& phis_input)        
    :levelset_multiple(grid_input,phis_input),number_of_regions(phis_input.m)
{
    rendering_levelset_multiple_region_objects.Resize(phis_input.m);
    for(int i=0;i<number_of_regions;i++) rendering_levelset_multiple_region_objects(i)=new RENDERING_LEVELSET_MULTIPLE_REGION_OBJECT<T,T_LEVELSET_MULTIPLE>(levelset_multiple,i);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_LEVELSET_MULTIPLE> RENDERING_LEVELSET_MULTIPLE_OBJECT<T_LEVELSET_MULTIPLE>::
~RENDERING_LEVELSET_MULTIPLE_OBJECT()
{
    for(int i=0;i<rendering_levelset_multiple_region_objects.m;i++) delete rendering_levelset_multiple_region_objects(i);
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T_LEVELSET_MULTIPLE> bool RENDERING_LEVELSET_MULTIPLE_OBJECT<T_LEVELSET_MULTIPLE>::
Intersection(RAY<VECTOR<T,3> >& ray,const int lowest_priority,const RENDERING_OBJECT<T>** intersected_object) const
{
    if(priority<lowest_priority) return false;
    int region_start,region_end;
    if(Intersection(ray,region_start,region_end,small_number)){
        if(region_end==-1) *intersected_object=rendering_levelset_multiple_region_objects(region_start);
        else if(region_start==-1) *intersected_object=rendering_levelset_multiple_region_objects(region_end);
        else{
            if(region_start>0&&rendering_levelset_multiple_region_objects(region_start)->priority>rendering_levelset_multiple_region_objects(region_end)->priority)
                *intersected_object=rendering_levelset_multiple_region_objects(region_start);
            else *intersected_object=rendering_levelset_multiple_region_objects(region_end);}
        return true;}
    return false;
}
//#####################################################################
// Function Intersected_Region
//#####################################################################
template<class T_LEVELSET_MULTIPLE> int RENDERING_LEVELSET_MULTIPLE_OBJECT<T_LEVELSET_MULTIPLE>::
Intersected_Region(RAY<VECTOR<T,3> >& ray) const
{
    int region_start,region_end;        
    if(Intersection(ray,region_start,region_end,small_number)){
        if(region_end==-1) return region_start;
        else if(region_start==-1) return region_end;
        else{
            if(region_start>0&&rendering_levelset_multiple_region_objects(region_start)->priority>rendering_levelset_multiple_region_objects(region_end)->priority)
                return region_start;
            else return region_end;}
        return -1;}
    return -1;
}
//#####################################################################
// Function Inside_Region_Only
//#####################################################################
template<class T_LEVELSET_MULTIPLE> bool RENDERING_LEVELSET_MULTIPLE_OBJECT<T_LEVELSET_MULTIPLE>::
Inside_Region_Only(const VECTOR<T,3>& location,int region_check) const
{
    for(int i=0;i<number_of_regions;i++){
        bool inside_region=rendering_levelset_multiple_region_objects(i)->Inside(location);
        if(i!=region_check && inside_region) return false;}
    return true;
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T_LEVELSET_MULTIPLE> bool RENDERING_LEVELSET_MULTIPLE_OBJECT<T_LEVELSET_MULTIPLE>::
Inside(const VECTOR<T,3>& location,RENDERING_OBJECT<T>** intersected_object) const
{
    for(int i=0;i<number_of_regions;i++)
        if(rendering_levelset_multiple_region_objects(i)->Inside(location)){
            *intersected_object=(RENDERING_OBJECT<T>*)rendering_levelset_multiple_region_objects(i);return true;}
    return false;
}
//#####################################################################
// Function Object_Space_Bounding_Box
//#####################################################################
template<class T_LEVELSET_MULTIPLE> auto RENDERING_LEVELSET_MULTIPLE_OBJECT<T_LEVELSET_MULTIPLE>::
Object_Space_Bounding_Box() const -> RANGE<VECTOR<T,3> >
{
    return levelset_multiple.grid.domain;
}
//#####################################################################
// Function Integration_Step
//#####################################################################
template<class T_LEVELSET_MULTIPLE> auto RENDERING_LEVELSET_MULTIPLE_OBJECT<T_LEVELSET_MULTIPLE>::
Integration_Step(const T phi) const -> T
{
    T distance=abs(phi);
    if(distance > 3*levelset_multiple.grid.dX.Min()) return (T).5*distance;    
    else if(distance > levelset_multiple.grid.dX.Min()) return (T).25*distance;
    return (T).1*levelset_multiple.grid.dX.Min();
}
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T_LEVELSET_MULTIPLE> auto RENDERING_LEVELSET_MULTIPLE_OBJECT<T_LEVELSET_MULTIPLE>::
Generate_Triangles() const -> TRIANGULATED_SURFACE<T>*
{
    return TRIANGULATED_SURFACE<T>::Create();
}
namespace{
template<class T> T Iterative_Solver_Tolerance(){STATIC_ASSERT((T)false);}
template<> float Iterative_Solver_Tolerance(){return (float).01;}
template<> double Iterative_Solver_Tolerance(){return .001;}
}
//#####################################################################
// Function Intersection
//#####################################################################
template<class T_LEVELSET_MULTIPLE> bool RENDERING_LEVELSET_MULTIPLE_OBJECT<T_LEVELSET_MULTIPLE>::
Intersection(RAY<VECTOR<T,3> >& ray,int &region_start,int &region_end,const T thickness) const
{
    bool save_semi_infinite=ray.semi_infinite;T save_t_max=ray.t_max;int save_aggregate=ray.aggregate_id;

    T t_start,t_end; 
    RANGE<VECTOR<T,3> > box=levelset_multiple.grid.domain;
    bool exit_intersection=false;int exit_aggregate=0;T exit_t_max=0;
    int intersect_box=INTERSECTION::Intersects(ray,box,thickness);
    int outside_box=box.Outside(ray.endpoint,thickness);

    if(outside_box && !intersect_box) return false; // missed the box
    else if(outside_box){ // intersected the box from the outside
        VECTOR<T,3> point=ray.Point(ray.t_max); // intersection point with the box
        point=box.Thickened(-4*thickness).Clamp(point); // moves the point inside the box
        T phi_temp;
        region_end=levelset_multiple.Inside_Region(point,phi_temp);
        region_start=-1;    // -1 means outside of levelset_multiple_object
        return true;}
    else if(!intersect_box){t_start=0;t_end=ray.t_max;} // intersected some object inside the box
    else{ // intersects the box from inside
        t_start=0;t_end=ray.t_max;
        exit_intersection=true;exit_t_max=ray.t_max;exit_aggregate=ray.aggregate_id; // save for exiting rays
        ray.semi_infinite=save_semi_infinite;ray.t_max=save_t_max;ray.aggregate_id=save_aggregate;} // box intersection doesn't count    

    // start marching
    T t1=t_start+thickness;
    T phi1;region_start=levelset_multiple.Inside_Region(ray.Point(t1),phi1);
    T t2=t1+Integration_Step(phi1);
    // march through the line segment
    while(t2 <= t_end){
        T phi2;region_end=levelset_multiple.Inside_Region(ray.Point(t2),phi2);
        if(region_start != region_end){
            LEVELSET_MULTIPLE_ON_A_RAY<T_LEVELSET_MULTIPLE> levelset_multiple_on_a_ray(levelset_multiple,ray,region_start);
            ITERATIVE_SOLVER<T> iterative_solver;
            iterative_solver.tolerance=Iterative_Solver_Tolerance<T>()*thickness;
            ray.semi_infinite=false;
            ray.t_max=iterative_solver.Bisection_Secant_Root(levelset_multiple_on_a_ray,t1,t2);
            ray.aggregate_id=-1;
            return true;}
        else{
            t1=t2;phi1=phi2;t2=t1+Integration_Step(phi1);}}
    // check the last piece of the line segment
    t2=t_end;
    T phi2;region_end=levelset_multiple.Inside_Region(ray.Point(t_end),phi2);
    if(region_start != region_end){
        LEVELSET_MULTIPLE_ON_A_RAY<T_LEVELSET_MULTIPLE> levelset_multiple_on_a_ray(levelset_multiple,ray,region_start);
        ITERATIVE_SOLVER<T> iterative_solver;
        iterative_solver.tolerance=Iterative_Solver_Tolerance<T>()*thickness;
        ray.semi_infinite=false;
        ray.t_max=iterative_solver.Bisection_Secant_Root(levelset_multiple_on_a_ray,t1,t2);
        ray.aggregate_id=-1;
        return true;}
    else if(exit_intersection){
        region_end=-1;ray.semi_infinite=false;ray.t_max=exit_t_max;ray.aggregate_id=exit_aggregate;
        return true;}
    else  return false;// exiting ray
}
//#####################################################################
template class RENDERING_LEVELSET_MULTIPLE_OBJECT<LEVELSET_MULTIPLE<VECTOR<float,3> > >;
template class RENDERING_LEVELSET_MULTIPLE_OBJECT<LEVELSET_MULTIPLE<VECTOR<double,3> > >;
}
