#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Fluids/PhysBAM_Compressible/Conservation_Law_Solvers/BOUNDARY_OBJECT.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,class TV_DIMENSION> BOUNDARY_OBJECT<TV,TV_DIMENSION>::
BOUNDARY_OBJECT()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,class TV_DIMENSION> BOUNDARY_OBJECT<TV,TV_DIMENSION>::
~BOUNDARY_OBJECT()
{
}
//#####################################################################
// Function Get_State_At_Location
//#####################################################################
template<class TV,class TV_DIMENSION> void BOUNDARY_OBJECT<TV,TV_DIMENSION>::
Get_State_At_Location(const GRID<VECTOR<T,1> >& grid_1d,const ARRAY<TV_DIMENSION,VECTOR<int,1> >& U_1d,const T location,const VECTOR<int,2>& region_boundaries,TV_DIMENSION& u_1d)
{
    T dx=grid_1d.dX.x,one_over_dx=grid_1d.one_over_dX.x;
    int lower_cell=grid_1d.Cell(VECTOR<T,1>(location),4).x; //using 4 ghost cells as thats the maximum we can have.
    int higher_cell=lower_cell+1;
    T epsilon=location-grid_1d.Center(VECTOR<int,1>(lower_cell)).x;
    if(lower_cell>=region_boundaries.x){ //use lower_cell,lower_cell+1
        if(lower_cell>region_boundaries.y) lower_cell=higher_cell=region_boundaries.y;
        else if(higher_cell>region_boundaries.y) higher_cell=region_boundaries.y;

        u_1d=U_1d(lower_cell)+epsilon*(U_1d(higher_cell)-U_1d(lower_cell))*one_over_dx;}
    else{ //use lower_cell+1,lower_cell+2
        lower_cell++;higher_cell++;
        if(lower_cell>region_boundaries.y) lower_cell=higher_cell=region_boundaries.y;
        else if(higher_cell>region_boundaries.y) higher_cell=region_boundaries.y;

        u_1d=U_1d(lower_cell)+(dx-epsilon)*(U_1d(higher_cell)-U_1d(lower_cell))*one_over_dx;}
}
//#####################################################################
// Function Fill_Ghost_Cells_Neumann
//#####################################################################
template<class TV,class TV_DIMENSION> void BOUNDARY_OBJECT<TV,TV_DIMENSION>::
Fill_Ghost_Cells_Neumann(const GRID<VECTOR<T,1> >& grid_1d,ARRAY<TV_DIMENSION,VECTOR<int,1> >& U_1d,const T_FACE_ARRAYS_SCALAR& face_velocities,const VECTOR<int,TV::m-1>& node_lower_dimension,const int axis,
    const int ghost_cells,const bool use_exact_neumann_face_location,const VECTOR<int,2>& domain,const VECTOR<int,2>& region_boundaries,const VECTOR<bool,2>& psi_N,CONSERVATION_CALLBACKS<T>* callbacks)
{
    TV_DIMENSION u_1d;
    ARRAY<INTERVAL<T> > regions;regions.Resize(2);
    regions(0)=INTERVAL<T>(region_boundaries.x-ghost_cells,region_boundaries.x-1);
    regions(1)=INTERVAL<T>(region_boundaries.y+1,region_boundaries.y+ghost_cells);
    for(int side=0;side<2;side++) if(psi_N(side))
        for(int i=regions(side).min_corner;i<regions(side).max_corner;i++) if(i>=domain.x && i<domain.y){
            int reflection_face=side&1?regions(side).min_corner:regions(side).max_corner+1;
            if(use_exact_neumann_face_location){
                T location=grid_1d.Center(VECTOR<int,1>(i)).x;
                T boundary_face_location;callbacks->Get_Neumann_Face_Location(grid_1d,reflection_face,boundary_face_location);
                T reflected_location=2*boundary_face_location-location;
                Get_State_At_Location(grid_1d,U_1d,reflected_location,region_boundaries,u_1d);}
            else{
                int reflected_node=2*reflection_face-i-1;
                int extreme_node=side&1?region_boundaries.x:region_boundaries.y;
                bool is_outside_region=side&1?reflected_node<extreme_node:reflected_node>extreme_node;
                if(is_outside_region) reflected_node=extreme_node;
                u_1d=U_1d(reflected_node);}
            T neumann_face_velocity=face_velocities.Component(axis)(node_lower_dimension.Insert(reflection_face,axis));
            Apply_Neumann_Boundary_Condition(u_1d,neumann_face_velocity,axis);
            U_1d(i)=u_1d;}
}

namespace PhysBAM{
template class BOUNDARY_OBJECT<VECTOR<float,1>,VECTOR<float,1> >;
template class BOUNDARY_OBJECT<VECTOR<float,1>,VECTOR<float,2> >;
template class BOUNDARY_OBJECT<VECTOR<float,1>,VECTOR<float,3> >;
template class BOUNDARY_OBJECT<VECTOR<float,1>,VECTOR<float,4> >;
template class BOUNDARY_OBJECT<VECTOR<float,1>,VECTOR<float,5> >;
template class BOUNDARY_OBJECT<VECTOR<float,1>,VECTOR<float,6> >;
template class BOUNDARY_OBJECT<VECTOR<float,1>,float>;
template class BOUNDARY_OBJECT<VECTOR<float,2>,VECTOR<float,1> >;
template class BOUNDARY_OBJECT<VECTOR<float,2>,VECTOR<float,2> >;
template class BOUNDARY_OBJECT<VECTOR<float,2>,VECTOR<float,3> >;
template class BOUNDARY_OBJECT<VECTOR<float,2>,VECTOR<float,4> >;
template class BOUNDARY_OBJECT<VECTOR<float,2>,VECTOR<float,5> >;
template class BOUNDARY_OBJECT<VECTOR<float,2>,VECTOR<float,6> >;
template class BOUNDARY_OBJECT<VECTOR<float,2>,float>;
template class BOUNDARY_OBJECT<VECTOR<float,3>,VECTOR<float,1> >;
template class BOUNDARY_OBJECT<VECTOR<float,3>,VECTOR<float,2> >;
template class BOUNDARY_OBJECT<VECTOR<float,3>,VECTOR<float,3> >;
template class BOUNDARY_OBJECT<VECTOR<float,3>,VECTOR<float,4> >;
template class BOUNDARY_OBJECT<VECTOR<float,3>,VECTOR<float,5> >;
template class BOUNDARY_OBJECT<VECTOR<float,3>,VECTOR<float,6> >;
template class BOUNDARY_OBJECT<VECTOR<float,3>,float>;
template class BOUNDARY_OBJECT<VECTOR<double,1>,VECTOR<double,1> >;
template class BOUNDARY_OBJECT<VECTOR<double,1>,VECTOR<double,2> >;
template class BOUNDARY_OBJECT<VECTOR<double,1>,VECTOR<double,3> >;
template class BOUNDARY_OBJECT<VECTOR<double,1>,VECTOR<double,4> >;
template class BOUNDARY_OBJECT<VECTOR<double,1>,VECTOR<double,5> >;
template class BOUNDARY_OBJECT<VECTOR<double,1>,VECTOR<double,6> >;
template class BOUNDARY_OBJECT<VECTOR<double,1>,double>;
template class BOUNDARY_OBJECT<VECTOR<double,2>,VECTOR<double,1> >;
template class BOUNDARY_OBJECT<VECTOR<double,2>,VECTOR<double,2> >;
template class BOUNDARY_OBJECT<VECTOR<double,2>,VECTOR<double,3> >;
template class BOUNDARY_OBJECT<VECTOR<double,2>,VECTOR<double,4> >;
template class BOUNDARY_OBJECT<VECTOR<double,2>,VECTOR<double,5> >;
template class BOUNDARY_OBJECT<VECTOR<double,2>,VECTOR<double,6> >;
template class BOUNDARY_OBJECT<VECTOR<double,2>,double>;
template class BOUNDARY_OBJECT<VECTOR<double,3>,VECTOR<double,1> >;
template class BOUNDARY_OBJECT<VECTOR<double,3>,VECTOR<double,2> >;
template class BOUNDARY_OBJECT<VECTOR<double,3>,VECTOR<double,3> >;
template class BOUNDARY_OBJECT<VECTOR<double,3>,VECTOR<double,4> >;
template class BOUNDARY_OBJECT<VECTOR<double,3>,VECTOR<double,5> >;
template class BOUNDARY_OBJECT<VECTOR<double,3>,VECTOR<double,6> >;
template class BOUNDARY_OBJECT<VECTOR<double,3>,double>;
}
